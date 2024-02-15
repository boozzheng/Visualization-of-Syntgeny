import plotly.graph_objects as go
import pandas as pd
import dash
from dash import dcc
from dash import html
import math
import plotly.express as px
import random
import plotly.colors as colors
import numpy as np
import colorsys
import plotly.offline as pyo

file_path = 'test.tsv'  # 替换为你的文件路径
df = pd.read_csv(file_path, sep='\t')

unique_values_col2 = df.iloc[:, 1].unique()
ortholog_number = df.iloc[:, 1].nunique()
unique_values_col4 = df.iloc[:, 3].unique()

prefix_values_col4 = tuple(value.split('@')[0] for value in unique_values_col4)
species = prefix_values_col4
genes = unique_values_col2

a = len(unique_values_col2)
b = len(unique_values_col4)

species_orthologs = []
ortholog_start = []
ortholog_end = []
ortholog_direction =[]
ortholog_relpos = []

#generate tuples of information，包含ortholog关系
for i in range(b):
    orthologs = []
    starts = []
    ends = []
    directions = []
    relpostions = []

    for j in range(a):
        col2_value = unique_values_col2[j]
        col4_value = unique_values_col4[i]

        #find the row according to the ortholog relationship
        matching_rows = df[(df.iloc[:, 1] == col2_value) & (df.iloc[:, 3] == col4_value)]

        ortholog_values = matching_rows.iloc[:, 4].tolist()
        orthologs.append(ortholog_values if len(ortholog_values) > 1 else ortholog_values[0])
        
        start_position = matching_rows.iloc[:,7].unique().tolist()
        starts.append(start_position if len(start_position) > 1 else start_position[0])

        end_position = matching_rows.iloc[:,8].unique().tolist()
        ends.append(end_position if len(end_position) > 1 else end_position[0])

        direction = matching_rows.iloc[:,9].unique().tolist()
        directions.append(direction if len(direction) > 1 else direction[0])

        relpos = matching_rows.iloc[:,11].unique().tolist()
        relpostions.append(relpos if len(relpos) > 1 else relpos[0])

    species_orthologs.append(tuple(orthologs))
    ortholog_start.append(tuple(starts))
    ortholog_end.append(tuple(ends))
    ortholog_direction.append(tuple(directions))
    ortholog_relpos.append(tuple(relpostions))

# 递归函数，将列表中的元素转换为整数
def recursive_int(value):
    if isinstance(value, list):
        return [recursive_int(v) for v in value]
    elif isinstance(value, float) and math.isnan(value):
        return math.nan
    else:
        return int(value)

ortholog_start_int = [
    tuple(recursive_int(value) for value in tup)
    for tup in ortholog_start
    ]

ortholog_end_int = [
    tuple(recursive_int(value) for value in tup)
    for tup in ortholog_end
    ]

ortholog_relpos_int = [
    tuple(recursive_int(value) for value in tup)
    for tup in ortholog_relpos
    ]

# species_orthologs = [
#     ['a', 'b', 'c', 'd', 'e'], 
#     ['nan', 'g', 'h', 'i', 'j'],   
#     ['k', 'l', 'm', 'n', 'nan']
# ]

# ortholog_relpos_int = [
#     (-2, -1, 0, 1, 2), 
#     (float('nan'), float('nan'), 0, float('nan'), -9), 
#     (float('nan'), 2, 0, 1, -5)
#     ]

flat_data = [value for sublist in ortholog_relpos_int for value in sublist if not np.isnan(value)]
min_value = min(flat_data)
max_value = max(flat_data)

#根据relative position生成一个用于layout的positions列表
positions = [[] for _ in range(len(species_orthologs))]
for value in range(min_value, max_value+1):
    
    for row_index, row in enumerate(ortholog_relpos_int):
        has_value = False
        for col_index, item in enumerate(row):
            if item == value:
                positions[row_index].append((row_index, col_index))
                has_value = True
        if not has_value:
            positions[row_index].append((row_index, 'gap_1'))

            
#缩短positions中的长gap为小gap
shortened_positions = positions
find_first = True
col_index = 0
for i in range(len(positions[0])):
    if all(row[col_index][1] == 'gap_1' for row in shortened_positions):
        if find_first == True:
            for row in shortened_positions:
                # 将满足条件的位置的第二个成分改为 "gap2"
                row[col_index] = (row[col_index][0], 'gap_2')
            find_first = False
            col_index += 1
            
        else:
            for row in shortened_positions:
                del row[col_index]
    else:
        find_first = True
        col_index += 1
        continue

color_start_end_positions = [((sub.index(next(filter(lambda x: isinstance(x[1], int), sub))), 
                              sub.index(next(filter(lambda x: isinstance(x[1], int), reversed(sub)))))) 
                              for sub in shortened_positions]


viridis_colors = px.colors.sequential.Viridis
bar_colors = random.sample(viridis_colors, ortholog_number)
text_colors = []
for color in bar_colors:
    # 将 RGB 颜色转换为 HSL 颜色
    r, g, b = px.colors.hex_to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)
    
    # 在色相上进行旋转 180 度（取模 1）
    h = (h + 0.5) % 1
    
    # 将 HSL 颜色转换为 RGB 颜色
    rgb_color = colorsys.hls_to_rgb(h, l, s)
    
    # 将 RGB 颜色转换为十六进制颜色码并添加到 text_colors
    text_colors.append('#%02x%02x%02x' % tuple(int(255 * c) for c in rgb_color))


fig = go.Figure()
count_gap2 = 0
gap2_before = [(0,) * 2] * len(shortened_positions)

for col_index in range(len(shortened_positions[0])): 

    if shortened_positions[0][col_index][1] == 'gap_2':
        count_gap2 += 1

    for row, sp_name, (color_s, color_e) in zip(shortened_positions, species, color_start_end_positions):

        bar_length = 20

        if row[col_index][1] == 'gap_1':
           
            hover_text =''

        if row[col_index][1] == 'gap_2':
            
            hover_text =''
            bar_length = 10
        
        if isinstance(row[col_index][1], int):
            opacity = 1
            
            label_text = str(species_orthologs[row[col_index][0]][row[col_index][1]][0]) if isinstance(species_orthologs[row[col_index][0]][row[col_index][1]], list) else str(species_orthologs[row[col_index][0]][row[col_index][1]])

            
            hover_text_start = ortholog_start_int[row[col_index][0]][row[col_index][1]]

           
            hover_text_end = ortholog_end_int[row[col_index][0]][row[col_index][1]]

            hover_text = f'Start point: {hover_text_start} End point:{hover_text_end}'
           
       
        if isinstance(row[col_index][1], str):
            fig.add_trace(go.Bar(
        
            x=(bar_length,),
            y=[sp_name],
           
            hoverinfo='text', #隐藏bar信息, like(20,dipro)
          
            opacity =  0,
            orientation='h',
          
            ))
        #添加bars of orthologs
        else:
            fig.add_trace(go.Bar(
            
                x=(bar_length,),
                y=[sp_name],
               
                hovertext=[hover_text],
                hoverinfo='text', #隐藏bar信息, like(20,dipro)
            
                opacity =  opacity,
                orientation='h',
                marker=dict(color=bar_colors[row[col_index][1]], line=dict(color='wheat', width=2)),
    
            ))

            # 添加文本标注
            fig.add_annotation(
                x=bar_length*(col_index+1/2-count_gap2)+10*count_gap2,#-bar_length/2,  
                
                y=sp_name,
                xref='x',
                yref='y',
                
                text= f'<b>{label_text}</b>',
                showarrow=False,
                xanchor='center',  # 将文本放置在条形图的中间
                yanchor='top',  # 将文本放置在条形图的底部
              
                font=dict(color=text_colors[row[col_index][1]],size=10,family='Balto,sans-serif'),
            )

            if ortholog_direction[row[col_index][0]][row[col_index][1]] == '+':
                fig.add_annotation(
                    x = bar_length*(col_index+1/2-count_gap2)+10*count_gap2, #1.single arrow with the line
                
                    y = sp_name,
                    xref='x',
                    yref='y',
                    text='',
                    showarrow=True,
                    arrowhead=4,
                    arrowsize=2,
                   
                    arrowcolor='rgba(0, 0, 0, 1)',
                    axref='x',
                    ayref='y',
                    ax= bar_length*(col_index-count_gap2)+10*count_gap2,
                    ay=sp_name,
                )

            if ortholog_direction[row[col_index][0]][row[col_index][1]] == '-':
                fig.add_annotation(
                    x = bar_length*(col_index+1/2-count_gap2)+10*count_gap2,
                    y = sp_name,
                    xref='x',
                    yref='y',
                    text='',
                    showarrow=True,
                    arrowhead=4,
                    arrowsize=2,
                    arrowcolor='rgba(0, 0, 0, 1)',
                    axref='x',
                    ayref='y',
                    ax= bar_length*(col_index+1-count_gap2)+10*count_gap2,
                    ay=sp_name,
                )

            if col_index == color_s:
                gap2_before[row[col_index][0]]=(row[col_index][0],count_gap2)
                
            if col_index == color_e:
                fig.add_shape(
                    type='line',
                    x0=(color_s-gap2_before[row[col_index][0]][1])*bar_length+gap2_before[row[col_index][0]][1]*10,
                    y0=sp_name,
                    x1=(color_e+1-count_gap2)*bar_length+count_gap2*10,
                    y1=sp_name,
                    line=dict(color='rgba(0, 0, 0, 1)', width=2)
                )

     

fig.update_layout(
    barmode='stack',
    yaxis=dict(title='', autorange='reversed'),  # 修改此处为yaxis
    xaxis=dict(title='', showticklabels=False),  # 修改此处为xaxis
   
    showlegend=False,
)

#fig.show()
fig.show()
pyo.plot(fig, filename='demo_2.html')

# app = dash.Dash(__name__)


# app.layout = html.Div(children=[
#     html.H1(children='Ortholog Dashboard'),

#     dcc.Graph(
#         id='Ortholog-chart',
#         figure=fig
#     )
# ])

# if __name__ == '__main__':
#     app.run_server(debug=True)