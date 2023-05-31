## Dependencies
import pandas as pd
import numpy as np
import gzip
import plotly
import plotly.express as px
import plotly.subplots as ps


## Config
from modules import config


## Windowed pca scripts

def fetch_variant_file_samples(variant_file_path):
    '''
    Fetch sample ids from variant file (used by read_metadata())
    '''

    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    
    # vcf
    if variant_file_path.endswith('.vcf') or variant_file_path.endswith('.vcf.gz'):
        with read_func(variant_file_path, 'rt') as vcf:
            for line in vcf:
                if line.startswith('#CHROM'):
                    variant_file_sample_lst = line.strip().split('\t')[9:]
                    break

    # GT/PL/GL file
    elif variant_file_path.endswith('.tsv') or variant_file_path.endswith('.tsv.gz'):
        with read_func(variant_file_path, 'rt') as gt_file:
            variant_file_sample_lst = gt_file.readline().strip().split('\t')[2:]

    # BEAGLE file
    elif variant_file_path.endswith('.beagle') or variant_file_path.endswith('.beagle.gz'):
        with read_func(variant_file_path, 'rt') as gt_file:
            variant_file_sample_lst = gt_file.readline().strip().split('\t')[5:]

    # remove potential duplicates (in case of e.g. GL file) while preserving order
    variant_file_sample_lst = list(dict.fromkeys(variant_file_sample_lst))

    return variant_file_sample_lst


def read_metadata(metadata_path, variant_file_sample_lst, taxon=None, group=None):
    '''
    Read in metadata, optionally filter by taxon and sort by gt_file sample order?
    '''

    # read in metadata
    metadata_df = pd.read_csv(
        metadata_path,
        sep='\t',
    )

    # re-name first column to 'id' (this is the only required column and must have unique ids)
    metadata_df.columns.values[0] = 'id'

    # subset input samples to match taxon group specification if specified
    if taxon and group:
        metadata_df = metadata_df.loc[metadata_df[taxon].isin(group.split(','))]
    
    # remove individuals that are not in the genotype file
    exclude_lst = [x for x in list(metadata_df['id']) if x not in variant_file_sample_lst]
    for i in exclude_lst:
        metadata_df.drop(metadata_df[metadata_df['id'] == i].index, inplace=True)

    return metadata_df


def polarize(w_pca_df, mean_threshold, guide_samples):
    '''
    Polarize windowed PCA output: if no guide_samples specified polarize PC orientation using a 
    subset of samples with large absolute values and small variability (var_threshold  and 
    mean_threshold variables)
    Furthermore, switch the entire plot if the largest absolute value is negative (to make plots
    more comparable between chromosomes)
    '''



    # if guide_samples not manually specified, select the var_threshold samples with the least 
    # variance, and from those the mean_threshold samples with the highest absolute value accross 
    # all windows as guide samples to calibrate the orientation of all windows
    if guide_samples:
        guide_samples = guide_samples.split(',')
    else:
        var_threshold = round(len(w_pca_df)*0.4) # remove the 60% most variable samples from the 
                                                 # dataset (keep 40% to not only keep the ones 
                                                 # around 0) --> keep 40%
        guide_samples_df = w_pca_df.loc[w_pca_df.dropna(axis=1).abs().var(axis=1).sort_values()\
                                        .index[0:var_threshold]]
        guide_samples = guide_samples_df.dropna(axis=1).abs().sum(axis=1).sort_values()\
                                        .index[-mean_threshold:]

    # subset w_pca_df to guide_samples
    guide_samples_df = w_pca_df.loc[guide_samples]

    # considering only guide samples, if the negative absolute value of a window is closer to the
    # absolute value 
    # that in, switch orientation of that window
    # (1 --> switch, 0 --> keep)
    
    sample_switch_lst = []    
    for row in guide_samples_df.iterrows():

        # need to treat first window special because it has no reference/previous window
        switch = [0]
        row = list(row[1])
        prev_window = row[0] if not row[0] == None else 0 # only if the second window is None, the 
                                                          # first window can be None, in that case 
                                                          # set it to 0 for numerical comparisons
        
        # parse by row (=by sample), each window separately and compare to previous
        for window in row[1:]:

            # current window has no value --> encode it as 0 (=no effect on balance)
            if window == None:
                switch.append(0)
                continue
            
            # current window closer to switched previous window --> encode it as 1 (=switch)
            elif abs(window - prev_window) > abs(window - (prev_window*-1)):
                switch.append(1)
                prev_window = (window*-1) # also switch window which becomes next previous
            
            # current window closer to unswitched previous window --> encode it as 0 (=keep)
            else:
                switch.append(-1)
                prev_window = window # don't switch previous window (next reference)
        
        # add switch for this sample to sample_switch_lst
        sample_switch_lst.append(switch)

    # combine list of per sample per window switch weightings into an array
    sample_switch_arr = np.array(sample_switch_lst, dtype=int).transpose()

    # get the row sums from sample_switch_arr (if sum/balance > 0 this means switch)
    switch_balance_lst = list(sample_switch_arr.sum(axis=1))

    # switch individual windows according to switch_balance_lst (switch if value is positive)
    for idx, val in zip(list(w_pca_df.columns), switch_balance_lst):
        if val < 0:
            w_pca_df[idx] = w_pca_df[idx]*-1

    # switch Y axis if largest absolute value is negative
    if abs(w_pca_df.to_numpy(na_value=0).min()) > abs(w_pca_df.to_numpy(na_value=0).max()):
        w_pca_df = w_pca_df * -1

    return w_pca_df


def annotate(w_df, metadata_df, value_col_name):
    '''
    Pivot windowed pca output and annotate with metadata
    '''

    # annotate with metadata
    for column_name in metadata_df.columns:
        w_df[column_name] = list(metadata_df[column_name])

    # replace numpy NaN with 'NA' for plotting (hover_data display)
    w_df = w_df.replace(np.nan, 'NA')

    # convert to long format for plotting
    w_anno_df = pd.melt(
        w_df,
        id_vars=metadata_df.columns,
        var_name='window_mid',
        value_name=value_col_name,
    )

    return w_anno_df


def plot_per_sample_values(w_anno_df, plot_col, color_taxon, chrom, start, stop, w_size, w_step):
    '''
    Plot one PC for all included samples along the chromosome
    '''

    # ensure plot_col is str
    plot_col = str(plot_col)

    # modify plot_col for output plot display
    plot_col_str = plot_col.replace('_', ' ').upper() if plot_col.startswith('pc') \
        else plot_col.replace('_', ' ').capitalize()

    # compile dct that sets hover data for html display
    hover_data_dct = {x: True for x in w_anno_df.columns}
    hover_data_dct[plot_col] = False

    # plot
    fig = px.line(
        w_anno_df,
        x='window_mid',
        y=plot_col,
        width=(stop-start)/config.plot_scaling_factor, height=500,
        line_group='id',
        color=color_taxon,
        hover_name='id', 
        hover_data=hover_data_dct,
        title=str('<b>Windowed ' + plot_col_str + ' of ' + chrom + ':' + str(start) + '-' + 
                  str(stop) + '</b><br> (window size: ' + str(w_size) + ' bp, window step: ' + 
                  str(w_step) + ' bp)'), 
        labels = dict(window_mid = 'position'),
    )

    # adjust layout
    fig.update_layout(
        template='simple_white',
        font_family='Arial',
        font_color='black',
        xaxis=dict(ticks='outside', mirror=True, showline=True, title='<b>Genomic position<b>'),
        yaxis=dict(ticks='outside', mirror=True, showline=True, 
                   title='<b>' + plot_col_str + '<b>'),
        legend={'traceorder':'normal'}, 
        title={'xanchor': 'center', 'y': 0.9, 'x': 0.45})

    # set x axis range
    fig.update_xaxes(range=[start, stop])

    # set line width
    fig.update_traces(line=dict(width=0.5))

    return fig


def plot_w_stats(w_stats_df, chrom, start, stop, w_size, w_step, min_var_per_w):
    '''
    Plot per windowstats: % explained by PC1 and PC2 + # of variants per window
    '''

    global missing_stretches # delete
    
    # for simplicity
    go = plotly.graph_objects
    
    # initialize figure
    fig = ps.make_subplots(
        specs=[[{'secondary_y': True}]],
        x_title='<b>Genomic position<b>',
        subplot_titles=[
            '<b>Per window stats of ' + chrom + ':' + str(start) + '-' + str(stop) + 
            '</b><br> (window size: ' + str(w_size) + ' bp, window step: ' + str(w_step) + ' bp)'
        ],
    )

    # pc_1 variance explained
    fig.add_trace(
        go.Scatter(
            x=w_stats_df.index,
            y=w_stats_df['pct_explained_pc_1'],
            name='PC 1',
            mode='lines',
            line=dict(color='#4d61b0', width=1),
            fill='tozeroy',
            connectgaps=True,
        ),
        secondary_y=False,
    )

    # pc_2 variance explained
    fig.add_trace(
        go.Scatter(
            x=w_stats_df.index,
            y=w_stats_df['pct_explained_pc_2'],
            name='PC 2',
            mode='lines',
            line=dict(color='#458255', width=1),
            fill='tozeroy',
            connectgaps=True,
        ),
        secondary_y=False
    )

    # check if there are gaps:
    gaps = w_stats_df.isnull().values.any()

    # plotly has a bug: if filling the area under the curve, the 'fill' doesn't break at missing 
    # data even when specifying connectgaps=True --> therefore, plot white rectangles on top to 
    # cover the missing data stretches
    if gaps:
        missing_stretches, stretch = [], []
        pc_1_max = max(w_stats_df['pct_explained_pc_1'])
        for w_mid, n_variants in zip(w_stats_df.index, w_stats_df['n_variants']):
            if n_variants >= min_var_per_w:
                if stretch:
                    missing_stretches.append(stretch)
                    stretch = []
            else: stretch.append(w_mid)
        if stretch: missing_stretches.append(stretch)
        for stretch in missing_stretches:
            fig.add_trace(
                go.Scatter(
                    x=[stretch[0], stretch[-1], stretch[-1], stretch[0]], 
                    y=[0, 0, pc_1_max, pc_1_max], 
                    fill='toself',
                    mode='none',
                    fillcolor='white',
                    hoverinfo='skip',
                    showlegend=False,
                )
            )
        
        # fill only regions between min_var_per_w and n_variants if n_variants < min_var_per_w this 
        # requires some hacks, such as adding a dummy datapoint at ± 0.0001 around missing 
        # stretches to delimit grey filled areas
        w_stats_gaps_df = w_stats_df.loc[w_stats_df['n_variants'] < min_var_per_w][['n_variants']]
        gap_edges = [x[0]-0.0001 for x in missing_stretches] + \
                    [x[-1]+0.0001 for x in missing_stretches]
        gap_edges_df = pd.DataFrame(
            [min_var_per_w] * len(gap_edges),
            gap_edges, 
            columns=['n_variants'],
        )
        w_stats_gaps_df = pd.concat([w_stats_gaps_df, gap_edges_df]).sort_index()
        fig.add_trace(
            go.Scatter(
                x=w_stats_gaps_df.index,
                y=w_stats_gaps_df['n_variants'],
                mode='lines',
                line=dict(color='rgba(0, 0, 0, 0)'),
                hoverinfo='skip',
                showlegend=False,
            ),
            secondary_y=True,
        )
        # (invisible) horizontal line just to fill the grey area in gaps
        fig.add_trace(
            go.Scatter(
                x=[start, stop],
                y=[min_var_per_w, min_var_per_w],
                mode='lines',
                line=dict(color='rgba(0, 0, 0, 0)'),
                fill='tonexty',
                hoverinfo='skip',
                showlegend=False
            ),
            secondary_y=True,
        )

    # horizontal line to show min_var_per_w threshold
    fig.add_trace(
        go.Scatter(
            x=[start, stop],
            y=[min_var_per_w, min_var_per_w],
            mode='lines',
            line=dict(color='#595959', dash='dot', width=1),
            hoverinfo='skip',
            showlegend=False
        ),
        secondary_y=True,
    )

    # add annotation for min_var_per_w line
    fig.add_trace(go.Scatter(
        x=[stop],
        y=[min_var_per_w-0.05*min_var_per_w],
        mode='lines+text',
        text=['min # of variants threshold '],
        textposition='top left',
        textfont=dict(color=['#595959']),
        showlegend=False,
        ),
        secondary_y=True,
    )

    # number of variants per window
    fig.add_trace(
        go.Scatter(
            x=w_stats_df.index,
            y=w_stats_df['n_variants'],
            name='# variants',
            mode='lines',
            line=dict(color='#595959', dash='dot', width=1)
        ),
        secondary_y=True,
    )

    # set x axis range
    fig.update_xaxes(
        range=[start, stop],
    )
    
    # set y axes ranges and titles
    fig.update_yaxes(
        rangemode='tozero',
        title_text='<b>% variance explained<b>',
        secondary_y=False
    )
    fig.update_yaxes(
        rangemode='tozero',
        title_text='<b># variants per window</b>',
        secondary_y=True
    )

    # adjust layout
    fig.update_layout(
        template='simple_white',
        font_family='Arial', font_color='black',
        autosize=False,
        width=(stop-start)/config.plot_scaling_factor, height=500,
        xaxis=dict(ticks='outside', mirror=True, showline=True),
        yaxis=dict(ticks='outside', mirror=True, showline=True),
        legend={'traceorder':'normal'},
        title={'xanchor': 'center', 'y': 0.9, 'x': 0.45},
        hovermode='x unified',
    )

    return fig
