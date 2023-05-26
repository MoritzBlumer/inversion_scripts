#!/usr/bin/env python
#
# Moritz Blumer | 2021-12-07 (2023-03-15)
#
# Conduct sliding window PCA with scikit-allel 



## File info
__author__ = 'Moritz Blumer, 2023'
__version__ = '2.0'
__email__ = 'lmb215@cam.ac.uk'


## Dependencies
import sys, os
import numpy as np
import pandas as pd
import allel


## Config
from modules import config


## Private functions

def parse_arguments():
    '''
    Parse command line arguments & print help message if # of arguments is incorrect
    '''

    global variant_file_path, metadata_path, output_prefix, chrom, start, stop, w_size, w_step, \
        min_maf, pc, taxon, group, color_taxon, guide_samples

    # print help message if incorrect number of arguments was specified
    if len(sys.argv)!= 13:
        print(
            '\n   python windowed_pca.py <variant file> <metadata> <output prefix> <region>\n\
                                <window size> <window step size> <pc> <filter column name>\n\
                                <filter column value> <color column name>\n\
                                <guide samples>\n\n\
            <variant file>           str    path to uncompressed or gzipped variant file\n\
                                            (VCF or genotype file; details -> README)\n\
            <metadata>               str    path to the metadata file (details -> README)\n\
            <output prefix>          str    prefix for output files\n\
            <region>                 int    target region in format "chr:start-stop"\n\
                                            (i.e. chr1:1-chrom_length to analyze the\n\
                                            entire chr1)\n\
            <window size>            int    sliding window size in bp, e.g. "1000000"\n\
            <window step>            int    sliding window step size in bp, e.g. "10000"\n\
            <minor allel frequency>  float  minor allel frequency threshold; specify\n\
                                            "None" to disable filter [default is 0.01]\n\
            <pc>                     int    principal component to use ("1" or "2")\n\
            <filter column name>     str    metadata column name to filter for\n\
                                            individuals to includede in the analysis,\n\
                                            e.g. "genus" (see <filter column value>)\n\
            <filter column value>    str    value to be filtered for in filter column;\n\
                                            Setting <filter column name> to "genus" and\n\
                                            <filter column value> to "Homo" would\n\
                                            include all individuals of the genus Homo\n\
                                            in the output, while ignoring all others.\n\
                                            (a comma-separated list of include values\n\
                                            can be provided, e.g. "Homo,Pan")\n\
            <color column name>      str    metadata column to assign colors by in the\n\
                                            output plot; if selecting "genus", all\n\
                                            individuals from the same genus will have\n\
                                            the same color in the output plots; if\n\
                                            specifying a comma-separated list like \n\
                                            "genus,species", one output plot is \n\
                                            generated for each color scheme\n\
            <guide samples>          str    list of samples to use forpolarization,\n\
                                            e.g. "ind1,ind2,ind3"; specify "None" for\n\
                                            automatic guide sample selection (details\n\
                                            --> README)\n',
        file=sys.stderr,
        )
        sys.exit()

    # fetch arguments    
    _, variant_file_path, metadata_path, output_prefix, region, w_size, w_step, min_maf, pc, \
        taxon, group, color_taxon, guide_samples = sys.argv
    
    # fetch chrom, start, stop from regions string
    chrom = region.split(':')[0]
    start = region.split(':')[1].split('-')[0]
    stop = region.split(':')[1].split('-')[1]

    # change str to int where appropriate
    start, stop, w_size, w_step, pc = int(start), int(stop), int(w_size), int(w_step), int(pc)

    # change min_maf to float if specified and update config
    config.min_maf = float(min_maf) if not min_maf == 'None' else None

    # set PC in config
    config.pc = pc

    # change output_prefix to lower case
    output_prefix = output_prefix.lower()

    # set guide_samples to NoneType if 'None' specified
    if guide_samples == 'None': guide_samples = None 


def pca(w_gt_arr, w_start, w_size):
    '''
    Conduct PCA, but if (n_variants < min_var_per_w) generate empty/dummy output instead
    '''

    # get window mid for X value
    w_mid = int(w_start + w_size/2-1)

    # count variants
    n_variants = w_gt_arr.shape[0]

    # if # variants passes specified threshold  
    if n_variants >= config.min_var_per_w:
        pca = allel.pca(
            w_gt_arr,
            n_components=2,
            copy=True,
            scaler='patterson',
            ploidy=2,
        )
        out = [
            pca[0][:, 0],
            pca[0][:, 1],
            pca[1].explained_variance_ratio_[0]*100,
            pca[1].explained_variance_ratio_[1]*100,
            n_variants,
        ]

    # else create empty output
    else:
        print(
            '[INFO] Skipped window ' + str(w_start) + '-' + str(w_start + w_size-1) + ' with ' +
            str(n_variants) + ' variants (threshold is ' + str(config.min_var_per_w) +
            ' variants per window)',
            file=sys.stderr, flush=True,
        )
        empty_lst = [None] * len(w_gt_arr[0])
        out = [
            empty_lst,
            empty_lst,
            None,
            None,
            n_variants,
        ]

    # append output
    w_mid_lst.append(w_mid)
    w_pca_lst.append(out[0]) if config.pc==1 else out[1] # depending on whether PC1 or PC2 was specified
    w_stats_lst.append([out[2], out[3], out[4]])


def windowed_pca(variant_file_path, chrom, start, stop, metadata_df, w_size, w_step, min_maf, pca):
    '''
    Window-by-window analysis
    '''
    
    global w_mid_lst, w_pca_lst, w_stats_lst

    # initialize results containers
    w_mid_lst = []
    w_pca_lst = []
    w_stats_lst = []

    # conduct windowed PCA using window_parser() function
    if variant_file_path.endswith('.vcf') or variant_file_path.endswith('.vcf.gz'):
        from modules.window_parser import win_vcf_gt
        win_vcf_gt(
            variant_file_path,
            chrom, start, stop,
            metadata_df['id'],
            w_size, w_step,
            pca,
            skip_monomorphic=True,
            min_maf = min_maf,
        )

    elif variant_file_path.endswith('.tsv') or variant_file_path.endswith('.tsv.gz'):
        from modules.window_parser import win_gt_file
        win_gt_file(
            variant_file_path,
            chrom, start, stop,
            metadata_df['id'],
            w_size, w_step,
            pca,
            skip_monomorphic=True,
            min_maf = min_maf,
        )

    # exit if no variants found
    if len(w_pca_lst) == 0:
        print(
        '\n[ERROR] No variants found, check if region was specifified correctly\n',
        file=sys.stderr, flush=True,
        )
        sys.exit()

    # compile output dataframe for windowed PCA
    w_pca_df = pd.DataFrame(
        np.transpose(w_pca_lst),
        index=list(metadata_df['id']),
        columns=w_mid_lst,
    )
    w_pca_df.index.names = ['id']

    # compile output dataframe for supplementary info (% variance explained, # sites per window)
    w_stats_df = pd.DataFrame(
        w_stats_lst,
        index=w_mid_lst,
        columns=['pct_explained_pc_1', 'pct_explained_pc_2', 'n_variants'],
    )
    w_stats_df.index.names = ['window_mid']

    return w_pca_df, w_stats_df



## Main

def main():

    # parse command line arguments
    parse_arguments()

    # make output directory if output_prefix contains '/'
    if '/' in output_prefix:
        if not os.path.exists('/'.join(output_prefix.split('/')[0:-1]) + '/'):
            os.makedirs('/'.join(output_prefix.split('/')[0:-1]) + '/')

    # compile text and stats figure output files (pc figure depends on color taxon --> see below)
    w_pca_tsv_path =        output_prefix + '.w_pc_' + str(config.pc) + '.tsv.gz'
    w_stats_tsv_path =      output_prefix + '.w_stats'         + '.tsv.gz'
    w_stats_fig_html_path = output_prefix + '.w_stats'         + '.html'
    w_stats_fig_pdf_path =  output_prefix + '.w_stats'         + '.pdf'

    # fetch sample ids from variant file
    from modules.utils import fetch_variant_file_samples
    variant_file_sample_lst = fetch_variant_file_samples(variant_file_path)

    # read metadata
    from modules.utils import read_metadata
    metadata_df = read_metadata(
        metadata_path,
        variant_file_sample_lst,
        taxon=taxon,
        group=group,
    )

    # if guide sample was specified, check if it is in the input samples
    if guide_samples:
        gs = set(guide_samples.split(','))
        missing_gs = list((gs^set(metadata_df['id']))&gs)
        if len(missing_gs) > 0:
            print(
                '\n[ERROR] Specified guide sample(s) ' + ', '.join(missing_gs) + ' are missing ' +
                'from the input samples\n',
                file=sys.stderr, flush=True,
            )
            sys.exit()

    # check if IDs are unique
    if not len(metadata_df['id']) == len(set(metadata_df['id'])):
        print(
            '\n[ERROR] Duplicate sample IDs (first metadata column)\n',
            file=sys.stderr, flush=True,
        )
        sys.exit()

    # if there is output from a previous run, use it
    if os.path.exists(w_pca_tsv_path) and os.path.exists(w_stats_tsv_path):
        print(
            '\n[INFO] Plotting data from previous run',
            file=sys.stderr, flush=True,
        )
        w_pca_df = pd.read_csv(
            w_pca_tsv_path,
            sep='\t',
            index_col=[0],
            na_values='NA',
        )
        w_pca_df.columns = [float(x) for x in w_pca_df.columns] # change column name dtype to int
        w_stats_df = pd.read_csv(
            w_stats_tsv_path,
            sep='\t',
            index_col=[0],
            na_values='NA',
        )
    
    else:
        print(
            '\n[INFO] Conducting windowed PCA\n',
            file=sys.stderr, flush=True,
        )
                
        # run windowed PCA
        w_pca_df, w_stats_df = windowed_pca(
            variant_file_path,
            chrom, start, stop,
            metadata_df,
            w_size, w_step,
            config.min_maf,
            pca,
        )

        # polarize windowed PCA output
        from modules.utils import polarize
        w_pca_df = polarize(
            w_pca_df,
            mean_threshold=config.mean_threshold,
            guide_samples=guide_samples,
        )

        # save output data before annotation if not already present
        print(
            '\n[INFO] Writing output TSVs',
            file=sys.stderr, flush=True,
        )
        w_pca_df.to_csv(
            w_pca_tsv_path,
            sep='\t',
            na_rep='NA',
            float_format='%.' + str(config.float_precision) + 'f',
            compression='gzip',
        )
        w_stats_df.to_csv(
            w_stats_tsv_path,
            sep='\t',
            na_rep='NA',
            float_format='%.' + str(config.float_precision) + 'f',
            compression='gzip',
        )

    # pivot windowed pca output and annotate with metadata
    from modules.utils import annotate
    w_pca_anno_df = annotate(
        w_pca_df,
        metadata_df,
        'pc_' + str(config.pc),
    )

    # free up memory
    del w_pca_df

    # plot windowed PCA output & save
    print(
            '\n[INFO] Generating output HTMLs & PDFs',
            file=sys.stderr, flush=True,
    )
    from modules.utils import plot_w_pca
    for c_taxon in color_taxon.split(','): 

        # compile output paths
        w_pca_fig_html_path = output_prefix + '.w_pc_' + str(config.pc) + '.' + c_taxon + '.html'
        w_pca_fig_pdf_path =  output_prefix + '.w_pc_' + str(config.pc) + '.' + c_taxon +  '.pdf'

        # plot & save
        w_pca_fig = plot_w_pca(
            w_pca_anno_df,
            config.pc,
            c_taxon,
            chrom, start, stop,
            w_size, w_step,
        )
        w_pca_fig.write_html(
            w_pca_fig_html_path
        )
        w_pca_fig.write_image(
            w_pca_fig_pdf_path,
            engine='kaleido', scale=2.4,
        )

    # free up memory
    del w_pca_fig

    # plot window stats & save
    from modules.utils import plot_w_stats
    w_stats_fig = plot_w_stats(
        w_stats_df,
        chrom, start, stop,
        w_size, w_step,
        config.min_var_per_w,
    )
    w_stats_fig.write_html(
        w_stats_fig_html_path
    )
    w_stats_fig.write_image(
        w_stats_fig_pdf_path,
        engine='kaleido', scale=2.4
    )

    print(
        '\n[INFO] Done\n',
        file=sys.stderr, flush=True,
    )

if __name__ == '__main__':
    main()
