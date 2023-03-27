import plotly.express as px
import pandas as pd


def plot_w_pca(w_pca_df, pc, color_taxon, chrom, start, stop, w_size, w_step):
    
    '''
    Plot one PC for all included sampled along the chromosome
    '''

    fig = px.line(w_pca_df, x='window_mid', y='pc_' + str(pc), line_group='id', color=color_taxon, hover_name='id', 
                    hover_data=[x for x in list(w_pca_df.columns) if x not in ['window_mid', 'pc_' + str(pc)]],
                    width=(stop-start)/20000, height=500,
                    title=str('<b>Windowed PC ' + str(pc) + ' of ' + chrom + ':' + str(start) + '-' + str(stop) + '</b><br> (window size: ' + str(w_size) + ' bp, window step: ' + str(w_step) + ' bp)'), 
                    labels = dict(pc_1 = '<b>PC 1<b>', pc_2 = '<b>PC 2<b>', window_mid = '<b>Genomic position<b>'))

    fig.update_layout(template='simple_white', font_family='Arial', font_color='black',
                    xaxis=dict(ticks='outside', mirror=True, showline=True),
                    yaxis=dict(ticks='outside', mirror=True, showline=True),
                    legend={'traceorder':'normal'}, 
                    title={'xanchor': 'center', 'y': 0.9, 'x': 0.45})

    fig.update_traces(line=dict(width=0.5))

    return fig


def plot_w_stats(w_stats_df, chrom, start, stop, w_size, w_step):
    
    '''
    Plot supplementary information per window: % explained by PC1 and PC2 + % of sites that were included per window
    '''

    #w_stats_df['window_mid'] = w_stats_df.index
    
    fig = px.line(w_stats_df, x=w_stats_df.index, y=['pct_explained_pc_1', 'pct_explained_pc_2', 'pct_included_sites'], 
                    width=(stop-start)/20000, height=500,
                    title=str('<b>Per window-stats of ' + chrom + ':' + str(start) + '-' + str(stop) + '</b><br> (window size: ' + str(w_size) + ' bp, window step: ' + str(w_step) + ' bp)'),
                    labels = dict(Genomic_Position = '<b>Genomic Position<b>', value = '<b>variable [%]<b>', window_mid = '<b>Genomic position<b>', pct_explained_pc_1 = '% variance explained PC 1'))
    
    fig.update_layout(template='simple_white', font_family='Arial', font_color='black',
                    xaxis=dict(ticks='outside', mirror=True, showline=True),
                    yaxis=dict(ticks='outside', mirror=True, showline=True),
                    legend={'traceorder':'normal'}, 
                    title={'xanchor': 'center', 'y': 0.9, 'x': 0.45},
                    hovermode='x unified')
    
    fig.update_traces(line=dict(width=1.0))

    return fig

