#####################################################
# Copyright (C) 2024 Bo Zheng

# This script is for the visualization of syntgeny, 
# showing the results of homology searches across species.

# Contact: bo.zheng1996@qq.com
#####################################################


import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import colorsys
import webbrowser
import argparse
import os
import taxopy

def read_tsv(tsv_file):
    df = pd.read_csv(tsv_file, delimiter="\t")
    return df

def data_to_diagram(df,tsv_file,file_format):
    nb_ps = df.iloc[:, 1].unique()
    sp_list = df.iloc[:, 3].unique()

    contig_list = []
    for sp in sp_list:
        matching_rows = df[df.iloc[:, 3] == sp]
        contigs = matching_rows.iloc[:, 5].unique().tolist()

        for index, row in matching_rows.iterrows():
        
            if pd.notna(row.iloc[11]): 
                contigs.remove(row.iloc[5]) 
                contigs.insert(0, (row.iloc[5]))
                break  

        contig_list.extend([value for value in contigs if pd.notna(value)])

    sp_y = []
    taxdb = taxopy.TaxDb(nodes_dmp='nodes.dmp',
                         names_dmp="names.dmp")
    for contig in contig_list:
        matching_rows = df[df.iloc[:, 5] == contig]
        for index, row in matching_rows.iterrows():
            if pd.notna(row.iloc[11]):
                taxid = int(row.iloc[3].split('@')[1])
                sp_y.append(taxdb.taxid2name[taxid])
                break
            else:
                sp_y.append('')
                break

    gene_list = []
    for contig in contig_list:
        genes_info = []
        for ps in nb_ps:
            matching_rows = df[(df.iloc[:, 5] == contig)&(df.iloc[:, 1] == ps)]
            if not matching_rows.empty:
                if pd.notna(matching_rows.iloc[0,11]):
                    row_values = matching_rows.iloc[0, [1, 4, 7, 8, 9, 11, 6, 10]].tolist()
                else:
                    row_values = matching_rows.iloc[0, [1, 4, 7, 8, 9, 1, 6, 10]].tolist()
                genes_info.append(row_values)
        gene_list.append(genes_info)

    def reverse_contig(contig):
        r_contig = []
        for gene in contig:
            reverse_gene = gene.copy()

            rel_pos = reverse_gene[5]
            rel_pos *= -1
            reverse_gene[5] = rel_pos

            dir = reverse_gene[4]
            if dir == '+':
                dir = ('-','r')
            elif dir == '-':
                dir = ('+','r')
            reverse_gene[4]=dir
            
            r_contig.append(reverse_gene)
        return r_contig

    def compare_contigs(ref, contig):
        r_contig = reverse_contig(contig)

        score_contig = 0
        score_r_contig = 0
        
        for ref_gene in ref:
            ref_pos = ref_gene[5]

            for gene in contig:
                gene_ps = gene[5]
                gene_or = gene[0]
                if ref_pos == gene_ps:
                    if ref_pos != gene_or: 
                        score_contig -= 2
            
            for gene in r_contig:
                gene_ps = gene[5]
                gene_or = gene[0]
                if ref_pos == gene_ps:
                    if ref_pos != gene_or: 
                        score_r_contig -= 2
        
        if score_contig > score_r_contig:
            return contig
        elif score_r_contig > score_contig:
            return r_contig
        else:
            for gene in contig:
                if gene[0] == 0 and gene[4] == ref[len(ref) // 2][4]:
                    return contig
                else:
                    return r_contig
          
    for contig in gene_list[1:]:
        ref = gene_list[0]
        for gene in contig:
            if gene[0] == 0:
                result =  compare_contigs(ref,contig)
                contig[:] = result[:]
                break

    column_count = sorted(list(set(int(gene[5]) for contig in gene_list for gene in contig)))
    insert_count = 0
    for i in range(len(column_count) - 1):
        if column_count[i+insert_count] != column_count[i + 1+insert_count] - 1:
            column_count.insert(i + 1+insert_count, '\u2022\t\u2022\t\u2022')
            insert_count += 1

    gene_diagram = []
    for sublist in gene_list:
        complete_row = []
        for value in column_count:
            found = False
            for item in sublist:
                if item[5] == value:
                    complete_row.append(item)
                    found = True
                    break
            if not found:
                complete_row.append('gap')
        gene_diagram.append(complete_row)

    for sublist in gene_diagram:
        first_non_gap_index = None
        last_non_gap_index = None

        for i, item in enumerate(sublist):
            if item != 'gap':
                if first_non_gap_index is None:
                    first_non_gap_index = i
                last_non_gap_index = i

        if first_non_gap_index is not None and last_non_gap_index is not None:
            gap_count = 0
            for i in range(first_non_gap_index, last_non_gap_index + 1):
                if sublist[i] == 'gap':
                    gap_count += 1
                else:
                    if gap_count > 0:
                        for j in range(i - gap_count, i):
                            sublist[j] = ('inter')
                        sublist[i-gap_count]=('inter', gap_count)
                        gap_count = 0
           
            if len(sublist[first_non_gap_index][4]) == 1:
                contig_front = int(sublist[first_non_gap_index][7] - 1)
                sublist.insert(0, ('contig_front', contig_front))
                contig_back = int(sublist[last_non_gap_index+1][6] - sublist[last_non_gap_index+1][7])
                sublist.append(('contig_back', contig_back))
            
            elif len(sublist[first_non_gap_index][4]) == 2:
            
                contig_front = int(sublist[first_non_gap_index][6] - sublist[first_non_gap_index][7])
                sublist.insert(0, ('contig_front', contig_front))
                contig_back = int(sublist[last_non_gap_index+1][7] - 1)
                sublist.append(('contig_back', contig_back))

            gap_count = 0
            for item in sublist[1:]:
                if item == 'gap':
                    gap_count += 1
                else:
                    break
            sublist[0] = (sublist[0][0], sublist[0][1], gap_count)

            gap_count = 0
            for item in sublist[len(sublist)-2::-1]:
                if item == 'gap':
                    gap_count += 1
                else:
                    break
            sublist[len(sublist)-1] = (sublist[len(sublist)-1][0], sublist[len(sublist)-1][1], gap_count)

    def generate_color_list(values, color_palette):
        min_value = min(values)*3/2
        max_value = max(values)*3/2
        color_index = [(value - min_value) / (max_value - min_value) for value in values]

        colors = []
        for index in color_index:
            raw_color = color_palette[int(index * (len(color_palette) - 1))]
            
            r, g, b = px.colors.hex_to_rgb(raw_color)
            h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)
            
            s *= 0.8
            l *= 1.2
            
            rgb_color = colorsys.hls_to_rgb(h, l, s)
            
            colors.append('#%02x%02x%02x' % tuple(int(255 * c) for c in rgb_color))

        return colors

    ########## Select your favorite color palette. ##########
    # color_palette = px.colors.sequential.Viridis
    # color_palette = px.colors.sequential.Turbo
    # color_palette = px.colors.cyclical.HSV
    # color_palette = px.colors.cyclical.Twilight
    color_palette = px.colors.cyclical.mrybm
    # color_palette = px.colors.cyclical.mygbm
    color_list = generate_color_list(nb_ps, color_palette)
    max_nb_ps = max(nb_ps)

    bar_len = 20

    fig = go.Figure()

    color_bg = ['rgba(250,250,250,1)','rgba(220,220,220,1)']

    for col, column in enumerate(zip(*(gene_diagram))):
        bg_color_index = 0
        prev_bg_color = color_bg[0]  

        for row, (gene, sp,contig) in enumerate(zip(column, sp_y,contig_list)):

            label =  f"<b>{sp}</b><br>{contig}"

            if sp != '':
                bg_color = color_bg[bg_color_index]
                prev_bg_color = bg_color
                bg_color_index = (bg_color_index + 1) % len(color_bg) 
            else:
                bg_color = prev_bg_color  

            if gene == 'gap':
                fig.add_trace(go.Bar(
                    x=(bar_len,),
                    y=[label],
                    orientation='h',
                    opacity=0,
                    hoverinfo='none',
                ))
            
            if gene[0] == 'contig_front' or gene[0] == 'contig_back':
                if gene[0] == 'contig_front':
                    x_line = col*bar_len
                    x_text = (gene[2]+1)*bar_len/2
                    x0_rect = 2
                    x1_rect = (gene[2]+1)*bar_len-2
                else:
                    x_line = (col+1)*bar_len
                    x_text = (col+1)*bar_len-(gene[2]+1)*bar_len/2
                    x0_rect = (col+1)*bar_len-2
                    x1_rect = (col+1)*bar_len - (gene[2]+1)*bar_len+2

                fig.add_trace(go.Bar(
                    x=(bar_len,),
                    y=[label],
                    orientation='h',
                    opacity=0,
                    hoverinfo='none',
                ))
                fig.add_shape(
                    dict(
                        type="line",
                        x0=x_line,  
                        y0=row - 0.4,  
                        x1=x_line,  
                        y1=row + 0.4,  
                        line=dict(
                            color='rgba(200,60,0,1)',  
                            width=3  
                        )
                    )
                )

                fig.add_annotation(
                    x=x_text,  
                    y=row,  
                    text=f"{gene[1]}",  
                    showarrow=False,  
                    font=dict(color='black'),  
                )

                fig.add_shape(
                    type="rect",
                    x0=x0_rect, 
                    y0=row - 0.4, 
                    x1=x1_rect, 
                    y1=row + 0.4,
                    fillcolor='rgba(0, 0, 0, 0)',  
                    line=dict(
                        color='black',  
                        width=1 
                    ),
                )

            if gene[0] == 'inter':
                fig.add_trace(go.Bar(
                    x=(bar_len,),
                    y=[label],
                    orientation='h',
                    opacity=0,
                    hoverinfo='text',
                ))
                fig.add_annotation(
                    x=(col+0.5)*bar_len,  
                    y=row,  
                    text=f"{gene[1]}",  
                    showarrow=False,  
                    font=dict(color='black'),  
                )
                if isinstance(gene, tuple):
                    fig.add_shape(
                    type="rect",
                    x0=col*bar_len+1, 
                    y0=row - 0.4, 
                    x1=(col+gene[1])*bar_len-1, 
                    y1=row + 0.4,
                    fillcolor='rgba(0, 0, 0, 0)',  
                    line=dict(
                        color='black',  
                        width=1  
                    ),
                )

            if len(gene) > 5 :
                bar_color = color_list[gene[0]+max_nb_ps]
                hover_text = str(gene[1])
                
                fig.add_trace(go.Bar(
                    x=(bar_len,),
                    y=[label],
                    orientation='h',
                    opacity=0,
                    hovertext=[hover_text],
                    hoverinfo='text',
                    marker=dict(color=bar_color, line=dict(color=bg_color, width=10)),#bar_color'rgba(255,255,255,1)'
                ))

                ratio = 0.8

                if gene[4] == '+':
                    arrow_path = f"M {col*bar_len} {row+0.4} L {(col+ratio)*bar_len} {row+0.4} L {(col+1)*bar_len} {row} L {(col+ratio)*bar_len} {row-0.4} L {col*bar_len} {row-0.4} Z"
                else:
                    arrow_path = f"M {(col+1-ratio)*bar_len} {row+0.4} L {(col+1)*bar_len} {row+0.4} L {(col+1)*bar_len} {row-0.4} L {(col+1-ratio)*bar_len} {row-0.4} L {col*bar_len} {row} Z"
                
                fig.add_shape(
                    type="path",
                    path=arrow_path,  
                    line=dict(color=bar_color, width=1),
                    fillcolor=bar_color,
                )

            fig.add_shape(
                type="rect",
                x0=col*bar_len-1, 
                y0=row - 0.5, 
                x1=(col+1)*bar_len+1, 
                y1=row + 0.5,
                fillcolor=bg_color,
                layer="below",
                line=dict(width=0),
            )

    fig.update_layout(
        yaxis=dict(
            autorange='reversed',
            showgrid=False,
        ),
        xaxis=dict(
            title='',
            showgrid=False,
            showticklabels=False, 
            side = 'top',
            zeroline=False,
        ),
        plot_bgcolor='rgba(0, 0, 0, 0)',
        barmode='stack',
        bargap=0.3, 
        showlegend=False,
        margin=dict(pad = 20),
        width=len(column_count)*120,  
        height=len(sp_y)*80,  
        autosize=False,
    )

    file_name = os.path.splitext(os.path.basename(tsv_file))[0]
    output = f'{file_name}.{file_format}'
    
    if file_format == 'html':
        fig.write_html(output)
        print(f"diagramm has beed saved as {output}")
    else:
        fig.write_image(output, engine='kaleido')
        print(f"diagramm has beed saved as {output}")
    
    webbrowser.open(output)

    return fig

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process TSV file and display plot.')
    parser.add_argument('tsv_file', type=str, help='Path to the TSV file')
    parser.add_argument('file_format',nargs='?', type=str, choices=['html', 'png', 'jpg', 'jpeg', 'webp', 'svg','pdf'], default='html',
                        help='Output file format (html, png, jpg, jpeg, webp, svg, pdf)')# 
    args = parser.parse_args()
    if args.file_format not in ['html', 'png', 'jpg', 'jpeg', 'webp', 'svg','pdf']:
        print("Error：format must be 'html', 'png', 'jpg', 'jpeg', 'webp', 'svg','pdf'.")
        exit()


    return args

if __name__ == '__main__':
    args = parse_arguments()
    tsv_file = args.tsv_file
    file_format = args.file_format
    
    df = read_tsv(tsv_file)
    data_to_diagram(df,tsv_file,file_format)

