#!/usr/bin/env python
#
# Moritz Blumer | 2021-12-07
#
# Conduct sliding window PCA using scikit-allel 
#
# To add PC2 output, but it uncomment lines tagged with '# uncomment for PC 2'
# The minimum number of variants per window is set to 300 and can be changed below (otherwise NA will be returned)
min_var_per_window = 300



# import packages

import allel
import gzip
import os
import sys
import numpy as np
import pandas as pd
import plotly.express as px



def parse_arguments():

    '''
    parse command line arguments and print help message if the number of arguments is different from what is expected
    '''

    # declare all variables global

    global gt_matrix_path, metadata_path, output_prefix, chrom, chrom_len, window_size, window_step, taxon, group, color_taxon, var_threshold, mean_threshold


    # fetch arguments
    
    _, gt_matrix_path, metadata_path, output_prefix, chrom, chrom_len, window_size, window_step, taxon, group, color_taxon, var_threshold, mean_threshold = sys.argv


    # print help message if incorrect number of arguments was specified

    if len(sys.argv) != 13:
        print('\nUsage:', file=sys.stderr)
        print('\tpython windowed_pca.py <genotype matrix> <metadata> <chromosome name> <chromosome length> <window size> <window step size> <filter column name> <filter column value> \ \n                         <color column name> <variance threshold> <mean threshold> <output prefix>\n', file=sys.stderr)
        print('\t\t<genotype matrix>\tstr\tpath to the genotype matrix file produced as described in the README', file=sys.stderr)
        print('\t\t<metadata>\t\tstr\tpath to the metadata file produced as described in the README', file=sys.stderr)
        print('\t\t<output prefix>\t\tstr\tprefix that will be used for all output files, can also be a directory to be created', file=sys.stderr)
        print('\t\t<chromosome name>\tstr\tname of the chromosome, e.g. "chr1"', file=sys.stderr)
        print('\t\t<chromosome length>\tint\tlength of the chromosome in bp, e.g. "32123123"', file=sys.stderr)
        print('\t\t<window size>\t\tint\tsize of the sliding window in bp, e.g. "1000000"', file=sys.stderr)
        print('\t\t<window step>\t\tint\tstep size of the sliding window in bp, e.g. "10000"', file=sys.stderr)
        print('\t\t<filter column name>\tstr\tset a metadata column name to be used to select individuals to be included in the analysis e.g. "genus" (see filter column value)', file=sys.stderr)
        print('\t\t<filter column value>\t\tstr\tselect a value to be filtered for in the defined filter column. Setting <filter column name> to "genus" and <filter column value> to "Homo" would include all individuals of the genus Homo in the output, and ignore all others. A comma-separated list of include values can be provided, to include for example a specific subset of genera ("Homo,Pan") ', file=sys.stderr)
        print('\t\t<color column name>\tstr\tselect a metadata column that will serve to partition included individuals into color groups in the output plots. If selecting e.g. "genus", all individuals from the same genus will have the same color in the output plots. If specifying a comma-separated list of column names (e.g. "genus,species"), two versions of each output plot will be produced, that differ only in the color scheme', file=sys.stderr)
        print('\t\t<variance threshold>\tint\trelevant to correct random switching along PC axes, see code for details, if unsure, use "9"', file=sys.stderr)
        print('\t\t<mean threshold>\tint\trelevant to correct random switching along PC axes, see code for details, if unsure, use "3"\n', file=sys.stderr)
        sys.exit()

    # convert str to int where necessary

    chrom_len, window_size, window_step, var_threshold, mean_threshold = int(chrom_len), int(window_size), int(window_step), int(var_threshold), int(mean_threshold)

    output_prefix = output_prefix.lower()



def prepare_data(gt_matrix_path, metadata_path, taxon=None, group=None):

    # obtain genotype matrix rownames (positions) and column names (samples) as lists

    pos_lst = list(pd.read_csv(gt_matrix_path, header=None, usecols=[1], skiprows=1, sep='\t')[1])

    with gzip.open(gt_matrix_path, 'rt') as f:
        samples_lst = f.readline().strip().split('\t')[4:]


    # read in metadata and genus to species column for easier filtering

    metadata_df = pd.read_csv(metadata_path, sep='\t', header=None, names=['primary_id', 'simple_id', 'supplier_id', 'seq_depth', 'clade', 'genus', 'species', 'sex', 'location', 'sublocation'])

    metadata_df['species'] = metadata_df['genus'] + '_' + metadata_df['species']


    # subset input samples to match taxon group specification if specified

    if taxon and group:
        metadata_df = metadata_df.loc[metadata_df[taxon].isin(group.split(','))]


    # remove individuals that are not in the callset

    exclude_lst = [x for x in list(metadata_df['primary_id']) if x not in samples_lst]

    for i in exclude_lst:
        metadata_df.drop(metadata_df[metadata_df['primary_id'] == i].index, inplace=True)


    # get index of samples kept after filtering
    
    keep_idx_lst = sorted([samples_lst.index(x) for x in list(metadata_df['primary_id'])])
    
    keep_id_lst = [samples_lst[x] for x in keep_idx_lst]


    # sort metadata by VCF sample order 
    
    metadata_df['primary_id'] = pd.Categorical(metadata_df['primary_id'], categories = keep_id_lst,ordered = True)
    
    metadata_df.sort_values('primary_id', inplace=True)


    # read in genotype matrix file line by line, skip header line (if idx > 0) and keep only infomative rows

    print('\n[INFO] Reading genotype matrix \n', file=sys.stderr)

    rows = []
    informative_idx_lst = []

    with gzip.open(gt_matrix_path, 'rt') as gt_matrix_file:

        for idx, line in enumerate(gt_matrix_file):

            if idx > 0:

                line = line.strip().split('\t')[4:]

                gts = [line[idx] for idx in keep_idx_lst]

                if not len(set(gts)) == 1:

                    informative_idx_lst.append(idx-1) # because the header line does not count
                    rows.append(gts)

                if idx % 100000 == 0:
                    print('[INFO] Read ' + str(idx) + ' lines', file=sys.stderr)


    # remove uninformative sites from pos_lst (otherwise PCA will fail)

    pos_lst = [pos_lst[x] for x in informative_idx_lst]
    
    pos_arr = np.array(pos_lst, dtype=int)
    
    pos_idx_arr = np.array(informative_idx_lst, dtype=int)


    # convert to array

    print('\n[INFO] Converting to numpy array (this can take minutes)', file=sys.stderr)

    gt_arr = np.array(rows, dtype=np.int8) # change to np.int8


    # free up memory
    del rows

    # return genotype matrix and metadata
    return [gt_arr, pos_arr, pos_idx_arr, metadata_df]



def compile_windo_arrays(chrom_len, window_size, window_step):
    '''
    returns three arrays based on the specified chromosome length and window parameters: 1. window starts, 2. window ends, 3. window mids
    '''

    window_start_arr = np.array(range(1, chrom_len-window_size, window_step), dtype=int)

    windows_stop_arr = np.array(window_start_arr + window_size, dtype=int)

    windows_mid_arr = np.array(window_start_arr + (0.5 * window_size), dtype=int)

    return window_start_arr, windows_stop_arr, windows_mid_arr



def window_pca(gt_arr, pos_arr, window_start, window_stop, min_var_per_window=min_var_per_window):
    '''
    Conduct a PCA based on the provided window coordinates
    - if window_idx_arr is empty return array of None, because no PCA can be conducted without data
    - if window_idx_arr contains only one line (== only one pos), also return array of None, be because 
      at least two informative positions are required to get two PCs. To retain the possibility to 
      generate PC2 windowed PCAs, ignore windows aith only one position (this is only the case once 
      on chr19 so far)
    '''
    window_idx_arr = np.where((pos_arr >= window_start) & (pos_arr < window_stop))[0]

    window_gt_arr = np.take(gt_arr, window_idx_arr, 0)

    if len(window_idx_arr) <= min_var_per_window:

        empty_array = [None] * window_gt_arr.shape[1]
        
        print('[INFO] Skipped window ' + str(window_start) + '-' + str(window_stop) + ' with ' + str(window_gt_arr.shape[0]) + ' variants (threshold is ' + str(min_var_per_window) + ' variants per window)', file=sys.stderr)
        
        return empty_array, empty_array, None, None, window_gt_arr.shape[0]
    
    else:
        
        pca = allel.pca(window_gt_arr, n_components=2, copy=True, scaler='patterson', ploidy=2)
        return pca[0][: , 0], pca[0][: , 1], pca[1].explained_variance_ratio_[0]*100, pca[1].explained_variance_ratio_[1]*100, window_gt_arr.shape[0]



def do_pca(gt_arr, pos_arr, window_start_arr, windows_mid_arr, windows_stop_arr, metadata_df):
    '''
    do windowed PCAs and return separate data frames for PC1 and PC2
    '''
    # initiate empty data frames for PC1 and PC2
    pc_1_df = pd.DataFrame(columns=windows_mid_arr, index = list(metadata_df['primary_id']))
    pc_2_df = pd.DataFrame(columns=windows_mid_arr, index = list(metadata_df['primary_id']))
    pc_1_pct_explained_lst = []
    pc_2_pct_explained_lst = []
    n_variants_lst = []

    # iterrate and conduct PCAs
    num = 1
    print('[INFO] Processed 0 of ' + str(len(window_start_arr)) + ' windows', file=sys.stderr)
    for window_start, window_mid, window_stop in zip(window_start_arr, windows_mid_arr, windows_stop_arr):
        pc_1_df[window_mid], pc_2_df[window_mid], pc_1_pct_explained, pc_2_pct_explained, n_variants = window_pca(gt_arr, pos_arr, window_start, window_stop)
        pc_1_pct_explained_lst.append(pc_1_pct_explained)
        pc_2_pct_explained_lst.append(pc_2_pct_explained)
        n_variants_lst.append(n_variants)
        if num % 500 ==0:
            print('[INFO] Processed ' + str(num) + ' of ' + str(len(window_start_arr)) + ' windows', file=sys.stderr)
        num += 1

    print('[INFO] Processed all windows', file=sys.stderr)

    pc_1_pct_explained_arr = np.array(pc_1_pct_explained_lst, dtype=float)
    pc_2_pct_explained_arr = np.array(pc_2_pct_explained_lst, dtype=float)

    # compile a data frame of additional info (% variance explained for PC_1 and PC_1, the % of sites per window)
    additional_info_df = pd.DataFrame(np.array([windows_mid_arr, pc_1_pct_explained_arr, pc_2_pct_explained_arr, n_variants_lst/np.array(window_size/100)]).transpose(), columns=['Genomic_Position', '% explained PC 1', '% explained PC 2', '% included sites'], dtype=float)

    return pc_1_df, pc_2_df, additional_info_df



def calibrate_annotate(pc_df, metadata_df, pc, var_threshold=9, mean_threshold=3):

    '''
    - take a pc_df and adjust window orientation using a selection of a few samples with high absolute values and small variability
    - then annotate the df with metadata
    '''

    # select the 9 samples with the least variance, and from those the 3 with the highest absolute value accross 
    # all windows as guide samples to calibrate the orientation of all windows
    guide_samples = list(pc_df.dropna(axis=1).abs().var(axis=1).sort_values(ascending=True).index[0:var_threshold])
    guide_samples_df = pc_df.loc[guide_samples]
    guide_samples = list(guide_samples_df.dropna(axis=1).abs().sum(axis=1).sort_values(ascending=False).index[0:mean_threshold])
    guide_samples_df = guide_samples_df.loc[guide_samples]

    # for each guide sample, determine whether the positive or negative absolute value of each window is closer 
    # to the value in the previous window. If the negative value is closer, switch that windows orientation
    # (1 --> switch, 0 --> keep)
    rows_lst = []
    for row in guide_samples_df.iterrows():
        row = list(row[1])
        last_window = row[0] if not row[0] == None else 0 # only if the first window is None, last_window can be None, in that case set it to 0 to enable below numerical comparisons
        out = [0]
        for window in row[1:]:
            if window == None:
                out.append(0)
                continue
            elif abs(window - last_window) > abs(window - (last_window*-1)):
                out.append(1)
                last_window = (window*-1)
            else:
                out.append(-1)
                last_window = window
        rows_lst.append(out)

    # sum up values from each row and save to switch_lst
    rows_arr = np.array(rows_lst, dtype=int).transpose()
    switch_lst = list(rows_arr.sum(axis=1))

    # switch individual windows according to switch_lst (switch if value is negative)
    for idx, val in zip(list(pc_df.columns), switch_lst):
        if val < 0:
            pc_df[idx] = pc_df[idx]*-1

    # switch Y axis if largest absolute value is negative
    if abs(pc_df.to_numpy(na_value=0).min()) > abs(pc_df.to_numpy(na_value=0).max()):
        pc_df = pc_df * -1

    # annotate with metadata
    pc_df['primary_id'] = list(metadata_df['primary_id'])
    pc_df['species'] = list(metadata_df['species'])
    pc_df['genus'] = list(metadata_df['genus'])
    pc_df['clade'] = list(metadata_df['clade'])
    pc_df['location'] = list(metadata_df['location'])
    pc_df['sublocation'] = list(metadata_df['sublocation'])
    pc_df['supplier_id'] = list(metadata_df['supplier_id'])
    pc_df['sex'] = list(metadata_df['sex'])
    pc_df['seq_depth'] = list(metadata_df['seq_depth'])

    # replace numpy NaN with 'NA' for plotting (hover_data display)
    pc_df = pc_df.replace(np.nan, 'NA')
    
    # convert to long format for plotting
    pc_df = pd.melt(pc_df, id_vars=['primary_id', 'species', 'genus', 'clade', 'location', 'sublocation', 'sex', 'seq_depth', 'supplier_id'], var_name='window_mid', value_name=pc)

    return pc_df



def plot_pc(pc_df, pc, color_taxon, chrom, chrom_len, window_size, window_step):
    
    '''
    Plot one PC for all included sampled along the chromosome
    '''

    fig = px.line(pc_df, x='window_mid', y=pc, line_group='primary_id', color=color_taxon, hover_name='primary_id', 
                    hover_data={'window_mid': False, 'species': True, 'genus': True, 'clade': True, 'location': True, 'sublocation': True, 'sex': True, 'seq_depth': True, 'supplier_id': True, 'primary_id': False, pc: False}, 
                    width=chrom_len/20000, height=500,
                    title=str('<b>Windowed PCA of ' + chrom + '</b><br> (chromosome length: ' + str(chrom_len) + ' bp, window size: ' + str(window_size) + ' bp, window step: ' + str(window_step) + ' bp)'), 
                    labels = dict(pc_1 = '<b>PC 1<b>', pc_2 = '<b>PC 2<b>', window_mid = '<b>Genomic position<b>'))

    fig.update_layout(template='simple_white', font_family='Arial', font_color='black',
                    xaxis=dict(ticks='outside', mirror=True, showline=True),
                    yaxis=dict(ticks='outside', mirror=True, showline=True),
                    legend={'traceorder':'normal'}, 
                    title={'xanchor': 'center', 'y': 0.9, 'x': 0.45})

    fig.update_traces(line=dict(width=0.5))

    #fig.show()

    return fig



def plot_additional_info(additional_info_df, chrom, chrom_len, window_size, window_step):
    
    '''
    Plot supplementary information per window: % explained by PC1 and PC2 + number of sites that were included per window
    '''
    
    fig = px.line(additional_info_df, x='Genomic_Position', y=['% explained PC 1', '% explained PC 2', '% included sites'], 
                    width=chrom_len/20000, height=500,
                    title=str('<b>Explained variance and of proportion of variants for windowed PCAs of ' + chrom + '</b><br> (chromosome length: ' + str(chrom_len) + ' bp, window size: ' + str(window_size) + ' bp, window step: ' + str(window_step) + ' bp)'),
                    labels = dict(Genomic_Position = '<b>Genomic Position<b>', value = '<b>variable [%]<b>', window_mid = '<b>Genomic position<b>'))
    
    fig.update_layout(template='simple_white', font_family='Arial', font_color='black',
                    xaxis=dict(ticks='outside', mirror=True, showline=True),
                    yaxis=dict(ticks='outside', mirror=True, showline=True),
                    legend={'traceorder':'normal'}, 
                    title={'xanchor': 'center', 'y': 0.9, 'x': 0.45},
                    hovermode='x unified')
    
    fig.update_traces(line=dict(width=1.0))
    
    #fig.show()

    return fig



def save_results(additional_info_df, pc_1_df, pc_2_plot=None):
    '''
    plotting and saving results (HTMLs, PDFs, TSVs)
    '''
    # pc_df(s)

    for c_taxon in color_taxon.split(','):

        pc_1_plot = plot_pc(pc_1_df, 'pc_1', c_taxon, chrom, chrom_len, window_size, window_step)
        pc_1_plot.write_html(str(output_prefix + chrom + '.pc_1.' + str(c_taxon) + '.html').lower())
        pc_1_plot.write_image(str(output_prefix + chrom + '.pc_1.' + str(c_taxon) + '.pdf').lower(), engine='kaleido', scale=2.4)
        pc_1_df.to_csv(str(output_prefix + chrom + '.pc_1.tsv').lower(), sep='\t', index=False)

        # pc_2_plot = plot_pc(pc_2_df, 'pc_1', c_taxon, chrom, chrom_len, window_size, window_step)                                  # uncomment for PC 2
        # pc_2_plot.write_html(str(output_prefix + chrom + '.pc_2.' + str(c_taxon) + '.html').lower())                               # uncomment for PC 2
        # pc_2_plot.write_image(str(output_prefix + chrom + '.pc_2.' + str(c_taxon) + '.pdf).lower()', engine='kaleido', scale=2.4)  # uncomment for PC 2
        # pc_2_df.to_csv(str(output_prefix + chrom + '.pc_2.tsv').lower(), sep='\t', index=False)                                    # uncomment for PC 2


    # supplementary df

    supplementary_plot = plot_additional_info(additional_info_df, chrom, chrom_len, window_size, window_step)
    supplementary_plot.write_html(str(output_prefix + chrom + '.supplemenary_info.html').lower())
    supplementary_plot.write_image(str(output_prefix + chrom + '.supplemenary_info.pdf').lower(), engine='kaleido', scale=2.4)
    additional_info_df = additional_info_df.fillna(value='NA')
    additional_info_df.to_csv(str(output_prefix + chrom + '.supplementary_info.tsv').lower(), sep='\t', index=False)



def main():

    '''
    run main analysis
    '''

    # parse command line arguments

    parse_arguments()
    
    
    # make output directory
    
    if output_prefix.endswith('/'):
        if not os.path.exists(output_prefix):
            os.makedirs(output_prefix)


    # conduct PCA if output TSVs do not exist, else skip and instead read them in

    if not os.path.exists(output_prefix + chrom + '.pc_1.tsv') and not os.path.exists(output_prefix + chrom + '.supplementary_info.tsv'):
    # if not os.path.exists(output_prefix + chrom + '.pc_1.tsv') and not os.path.exists(output_prefix + chrom + '.supplementary_info.tsv') and not os.path.exists(output_prefix + chrom + '.pc_2.tsv'): # uncomment for PC 2


        # compile window position arrays

        window_start_arr, windows_stop_arr, windows_mid_arr = compile_windo_arrays(chrom_len, window_size, window_step)


        # load genotype array for specified samples and chromosome into memory (this usually takes some time)

        gt_arr, pos_arr, pos_idx_arr, metadata_df = prepare_data(gt_matrix_path, metadata_path, taxon, group)


        # run PCA in windows

        pc_1_df, pc_2_df, additional_info_df = do_pca(gt_arr, pos_arr, window_start_arr, windows_mid_arr, windows_stop_arr, metadata_df)


        # free up space

        del gt_arr, pos_arr, pos_idx_arr

        # calibrate and annotate

        pc_1_df = calibrate_annotate(pc_1_df, metadata_df, 'pc_1', var_threshold=var_threshold, mean_threshold=mean_threshold)
        #pc_2_plot = pc_2_df = calibrate_annotate(pc_2_df, metadata_df, 'pc_2', var_threshold=5, mean_threshold=1) # uncomment for PC 2


    else:
        
        print('\n[INFO] Reading in existing data\n', file=sys.stderr)

        pc_1_df = pd.read_csv(output_prefix + chrom + '.pc_1.tsv',  sep='\t', index_col=None)
        pc_1_df.fillna('NA', inplace=True)
    
        # pc_2_df = pd.read_csv(output_prefix + chrom + '.pc_2.tsv', sep='\t', index_col=None)  # uncomment for PC 2
        # pc_2_df.fillna('NA', inplace=True)                                                    # uncomment for PC 2

        additional_info_df = pd.read_csv(output_prefix + chrom + '.supplementary_info.tsv',  sep='\t', index_col=None)
        additional_info_df.fillna(np.nan, inplace=True)

    # print INFO

    print('\n[INFO] Processing data and plotting\n', file=sys.stderr)


    # plot and save HTML and PDF plots and TSVs

    save_results(additional_info_df, pc_1_df)
    # save_results(additional_info_df, pc_1_df, pc_2_df) # uncomment for PC 2


    # print INFO

    print('[INFO] Done\n', file=sys.stderr)

main()

