import plotly.graph_objects as go
import pandas as pd
import plotly.offline as pyo

file_path = 'test.tsv'  
df = pd.read_csv(file_path, sep='\t')

unique_values_col2 = df.iloc[:, 1].unique()
unique_values_col4 = df.iloc[:, 3].unique()

prefix_values_col4 = tuple(value.split('@')[0] for value in unique_values_col4)
species = prefix_values_col4
genes = unique_values_col4

a = len(unique_values_col2)
b = len(unique_values_col4)

species_orthologs = []

for i in range(b):
    species_name = f"species_{i+1}"
    orthologs = []

    for j in range(a):
        col2_value = unique_values_col2[j]
        col4_value = unique_values_col4[i]

        # 根据条件赋值
        matching_rows = df[(df.iloc[:, 1] == col2_value) & (df.iloc[:, 3] == col4_value)]

        ortholog_values = matching_rows.iloc[:, 4].tolist()
        orthologs.append(ortholog_values if len(ortholog_values) > 1 else ortholog_values[0])
        

    species_orthologs.append(tuple(orthologs))


colors = ['rgb(188, 149, 117)', 
          'rgb(128, 137, 122)',
          'rgb(198, 175, 117)', 
          'rgb(118, 134, 146)',
          'rgb(127, 116, 132)']


fig = go.Figure()

for i in range(0, len(species_orthologs[0])): 

    for xd, yd in zip(species_orthologs, species):
       
        if str(xd[i]) == 'nan':
            opacity = 0
            hover_text = "No Ortholog Found"
        else:  
            opacity = 1.0
            hover_text = '<br>'.join([str(element) for element in xd[i]]) if isinstance(xd[i], list) else str(xd[i])
        
        fig.add_trace(go.Bar(
         
            x=(20,)*b,
            y=[yd],
          
            hovertext=[hover_text], 

            hoverinfo='text', 
            
            opacity =  opacity,
            orientation='h',
            marker=dict(color=colors[i], line=dict(color=colors[i], width=2))
        ))
    
        
        if i == 0:
            arrow_annotation = dict(
                x=20*len(species_orthologs[0])+2,  # 设置箭头的x位置，超出最后一个bar
                y=yd,  # 设置箭头的y位置
                xref='x',
                yref='y',
                text='',
                showarrow=True,
                arrowhead=4, 
                arrowcolor='rgba(0, 0, 0, 1)', 
                axref='x',
                ayref='y',
                ax=20*len(species_orthologs[0]), 
                ay=yd
            )

            arrow_length = 110  # 箭头长度
            fig.add_shape(dict(
                type='line',
                x0=0,
                y0=yd,
                x1=20*len(species_orthologs[0]),
                y1=yd,
                line=dict(color='rgba(0, 0, 0, 0.2)', width=2)
                
            ))

            fig.add_annotation(arrow_annotation)

# 设置布局
fig.update_layout(
    barmode='stack',
    yaxis=dict(title='', autorange='reversed'), 
    xaxis=dict(title='', showticklabels=False),           
    legend=dict(title='Species'),
    showlegend=False,
)

# 在网页上显示图表
fig.show()
pyo.plot(fig, filename='demo.html')
