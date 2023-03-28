import gzip
import sys
import pandas as pd
import numpy as np
import plotly.express as px


def parse_arguments():

    '''
    parse command line arguments and print help message if the number of arguments is different from what is expected
    '''

    # declare all variables global

    global variant_file_path, metadata_path, output_prefix, chrom, start, stop, w_size, w_step, pc, taxon, group, color_taxon, guide_samples

    # fetch arguments
    
    _, variant_file_path, metadata_path, output_prefix, region, w_size, w_step, pc, taxon, group, color_taxon, guide_samples = sys.argv


    # print help message if incorrect number of arguments was specified

    if len(sys.argv) != 13:
        print('\nUsage:', file=sys.stderr)
        print('\tpython windowed_pca.py <variant file> <metadata> <output prefix> <region> <window size> <window step size> <filter column name> <filter column value> \ \n                         <color column name> <guide samples>\n', file=sys.stderr)
        print('\t\t<variant file>\tstr\tpath to uncompressed or gzipped variant file (VCF or genotype file; specifications --> README)', file=sys.stderr)
        print('\t\t<metadata>\t\tstr\tpath to the metadata file (specifications --> README)', file=sys.stderr)
        print('\t\t<output prefix>\t\tstr\tprefix that will be used for all output files, can also be a directory to be created', file=sys.stderr)
        print('\t\t<region>\tint\tchromosome and target coordinates in bp; use format "chr:start-stop" (i.e. chr1:1-chrom_length to analyze the entire chr1)', file=sys.stderr)
        print('\t\t<window size>\t\tint\tsize of the sliding window in bp, e.g. "1000000"', file=sys.stderr)
        print('\t\t<window step>\t\tint\tstep size of the sliding window in bp, e.g. "10000"', file=sys.stderr)
        print('\t\t<pc\t\tint\tprincipal component to use (either "1" or "2")', file=sys.stderr)
        print('\t\t<filter column name>\tstr\tset a metadata column name to be used to select individuals to be included in the analysis e.g. "genus" (see filter column value)', file=sys.stderr)
        print('\t\t<filter column value>\t\tstr\tselect a value to be filtered for in the defined filter column. Setting <filter column name> to "genus" and <filter column value> to "Homo" would include all individuals of the genus Homo in the output, and ignore all others. A comma-separated list of include values can be provided, to include for example a specific subset of genera ("Homo,Pan") ', file=sys.stderr)
        print('\t\t<color column name>\tstr\tselect a metadata column that will serve to partition included individuals into color groups in the output plots. If selecting e.g. "genus", all individuals from the same genus will have the same color in the output plots. If specifying a comma-separated list of column names (e.g. "genus,species"), two versions of each output plot will be produced, that differ only in the color scheme', file=sys.stderr)
        print('\t\t<guide samples>\tstr\t[optional]list of samples to use for polarization, e.g. "ind1,ind2,ind3" (details --> README)', file=sys.stderr)
        sys.exit()

    
    chrom = region.split(':')[0]
    start = region.split(':')[1].split('-')[0]
    stop = region.split(':')[1].split('-')[1]

    start, stop, w_size, w_step = int(start), int(stop), int(w_size), int(w_step)
    
    output_prefix = output_prefix.lower()


## GET RID OF ORDER BY GT FILE, OR NOT? --> TEST IF CURRENT VERSION CAUSES PROBLEMS
def read_metadata(variant_file_path, metadata_path, taxon=None, group=None):
    '''
    Read in metadata, optionally filter by taxon ?and sort by gt_file sample order?

    '''

    # fetch sample names from genotype file header
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as gt_file:
        samples_lst = gt_file.readline().strip().split('\t')[2:]

    # read in metadata
    metadata_df = pd.read_csv(metadata_path, sep='\t')

    # re-name first column to 'id' (this is the only required column and must have unique ids)
    metadata_df.columns.values[0] = 'id'

    # subset input samples to match taxon group specification if specified
    if taxon and group:
        metadata_df = metadata_df.loc[metadata_df[taxon].isin(group.split(','))]
    
    # remove individuals that are not in the genotype file
    exclude_lst = [x for x in list(metadata_df['id']) if x not in samples_lst]
    for i in exclude_lst:
        metadata_df.drop(metadata_df[metadata_df['id'] == i].index, inplace=True)
    
    # # get index of samples kept after filtering
    # sample_idx_lst = sorted([samples_lst.index(x) for x in list(metadata_df['id'])])
    # keep_id_lst = [samples_lst[x] for x in sample_idx_lst]

    # # sort metadata by VCF sample order 
    # metadata_df['id'] = pd.Categorical(metadata_df['id'], categories = keep_id_lst, ordered = True)
    # metadata_df.sort_values('id', inplace=True)

    return metadata_df


def polarize(w_pca_df, var_threshold, mean_threshold, guide_samples): ## IMPROVE GUIDESAMPLE SETTINGS

    '''
    Polarize windowed PCA output: if no guide_samples specified polarize PC orientation using a subset of samples 
    with large absolute values and small variability
    '''

    # if $guide_samples not manually specified, select the $var_threshold samples with the least variance, and 
    # from those the $mean_threshold with the highest absolute value accross all windows as guide samples to 
    # calibrate the orientation of all windows
    if guide_samples: # check if this makes sense from #######
        guide_samples = mean_threshold.split(',')
        guide_samples_df = w_pca_df.loc[guide_samples]
    else:
        guide_samples = list(w_pca_df.dropna(axis=1).abs().var(axis=1).sort_values(ascending=True).index[0:var_threshold])
        guide_samples_df = w_pca_df.loc[guide_samples]
        guide_samples = list(guide_samples_df.dropna(axis=1).abs().sum(axis=1).sort_values(ascending=False).index[0:mean_threshold])

    # to #######
    guide_samples_df = guide_samples_df.loc[guide_samples]

    # considering all guide samples, if the negative absolute value of each window is closer 
    # that in, switch orientation of that window
    # (1 --> switch, 0 --> keep)
    
    rows_lst = []    
    for row in guide_samples_df.iterrows():
        row = list(row[1])
        prev_window = row[0] if not row[0] == None else 0 # only if current window is None, prev_window can be None, in that case set it to 0 to enable below numerical comparisons
        out = [0]
    
        for window in row[1:]:
            if window == None:
                out.append(0)
                continue

            elif abs(window - prev_window) > abs(window - (prev_window*-1)):
                out.append(1)
                prev_window = (window*-1)
    
            else:
                out.append(-1)
                prev_window = window
    
        rows_lst.append(out)

    # sum up values from each row and save to switch_lst
    rows_arr = np.array(rows_lst, dtype=int).transpose()
    switch_lst = list(rows_arr.sum(axis=1))

    # switch individual windows according to switch_lst (switch if value is negative)
    for idx, val in zip(list(w_pca_df.columns), switch_lst):
        if val < 0:
            w_pca_df[idx] = w_pca_df[idx]*-1

    # switch Y axis if largest absolute value is negative
    if abs(w_pca_df.to_numpy(na_value=0).min()) > abs(w_pca_df.to_numpy(na_value=0).max()):
        w_pca_df = w_pca_df * -1

    return w_pca_df


def annotate(w_pca_df, metadata_df, pc):
    '''
    Pivot windowed pca output and annotate with metadata
    '''

    # annotate with metadata
    for column_name in metadata_df.columns:
        w_pca_df[column_name] = list(metadata_df[column_name])

    # replace numpy NaN with 'NA' for plotting (hover_data display)
    w_pca_df = w_pca_df.replace(np.nan, 'NA')

    # convert to long format for plotting
    w_pca_anno_df = pd.melt(w_pca_df, id_vars=metadata_df.columns, var_name='window_mid', value_name=pc)

    return w_pca_anno_df


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

    fig.update_xaxes(range=[start, stop])

    fig.update_traces(line=dict(width=0.5))

    return fig


def plot_w_stats(w_stats_df, chrom, start, stop, w_size, w_step):
    
    '''
    Plot supplementary information per window: % explained by PC1 and PC2 + % of sites that were included per window
    '''
    
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
    
    fig.update_xaxes(range=[start, stop])
    
    fig.update_traces(line=dict(width=1.0))

    return fig

