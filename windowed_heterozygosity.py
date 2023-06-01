#!/usr/bin/env python
#
# Moritz Blumer | 2023-05-31
#
# Calculate per sample-heterozygosity rate in windows



## File info
__author__ = 'Moritz Blumer, 2023'
__version__ = '1.0'
__email__ = 'lmb215@cam.ac.uk'


## Dependencies
import sys, os
import numpy as np
import pandas as pd


## Config
from modules import config


## Private functions

def parse_arguments():
    '''
    Parse command line arguments & print help message if # of arguments is incorrect
    '''

    global variant_file_path, metadata_path, output_prefix, chrom, start, stop, w_size, w_step, \
        het_rate_scale, taxon, group, color_taxon

    # print help message if incorrect number of arguments was specified
    if len(sys.argv)!= 11:
        print(
            '\n   python windowed_het.py <variant file> <metadata> <output prefix> <region>\n\
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
            <het rate scale>         int    scale heterozygosity to the specified number of\n\
                                            bases in each window, e.g. "100000" would scale\n\
                                            the number of heterozygous sites in each 1000000 bp\n\
                                            window to 100000 (divide by 10)\n\
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
                                            generated for each color scheme\n',
        file=sys.stderr,
        )
        sys.exit()

    # fetch arguments    
    _, variant_file_path, metadata_path, output_prefix, region, w_size, w_step, het_rate_scale, \
        taxon, group, color_taxon = sys.argv
    
    # fetch chrom, start, stop from regions string
    chrom = region.split(':')[0]
    start = region.split(':')[1].split('-')[0]
    stop = region.split(':')[1].split('-')[1]

    # change str to int where appropriate
    start, stop, w_size, w_step, het_rate_scale = int(start), int(stop), int(w_size), int(w_step), \
        int(het_rate_scale)

    # change output_prefix to lower case
    output_prefix = output_prefix.lower()


def count_hets(w_gt_arr, w_start, w_size):
    '''
    Count the number of heterozygous sites, but if (n_variants < min_var_per_w) generate 
    empty/dummy output instead
    '''

    # get window mid for X value
    w_mid = int(w_start + w_size/2-1)

    # count variants
    n_variants = w_gt_arr.shape[0]

    # if # variants passes specified threshold count hets, else set to None 
    if n_variants >= config.min_var_per_w:
        n_hets_lst = list(np.sum(w_gt_arr == 1, axis=0))
    else: 
        n_hets_lst = [None] * w_gt_arr.shape[1]

    # count # missing sites
    n_miss_str = ','.join([str(x) for x in list(np.sum(w_gt_arr == -1, axis=0))])
    
    # append output
    w_mid_lst.append(w_mid)
    w_het_lst.append(n_hets_lst)
    w_stats_lst.append([n_variants, n_miss_str])



def windowed_het(variant_file_path, chrom, start, stop, metadata_df, w_size, w_step, count_hets):
    '''
    Window-by-window analysis
    '''
    
    global w_mid_lst, w_het_lst, w_stats_lst

    # initialize results containers
    w_mid_lst = []
    w_het_lst = []
    w_stats_lst = []

    # calculate windowed heterozygosity rate using window_parser() function
    if variant_file_path.endswith('.vcf') or variant_file_path.endswith('.vcf.gz'):
        from modules.window_parser import win_vcf_gt
        win_vcf_gt(
            variant_file_path,
            chrom, start, stop,
            metadata_df['id'],
            w_size, w_step,
            count_hets,
            skip_monomorphic=False,
            min_maf = 0,
        )

    elif variant_file_path.endswith('.tsv') or variant_file_path.endswith('.tsv.gz'):
        from modules.window_parser import win_gt_file
        win_gt_file(
            variant_file_path,
            chrom, start, stop,
            metadata_df['id'],
            w_size, w_step,
            count_hets,
            skip_monomorphic=False,
            min_maf = 0,
        )

    # exit if no variants found
    if len(w_het_lst) == 0:
        print(
        '\n[ERROR] No variants found, check if region was specifified correctly\n',
        file=sys.stderr, flush=True,
        )
        sys.exit()

    # compile output dataframe for windowed het rate
    w_het_df = pd.DataFrame(
        np.transpose(w_het_lst),
        index=list(metadata_df['id']),
        columns=w_mid_lst,
    )
    w_het_df.index.names = ['id']

    # scale by specified number of bp
    scale_factor = het_rate_scale / w_size
    w_het_df = w_het_df * scale_factor

    # compile output dataframe for supplementary info (# sites per window, # missing sites per 
    # sample per window)
    w_stats_df = pd.DataFrame(
        w_stats_lst,
        index=w_mid_lst,
        columns=['n_variants', 'n_miss_per_sample'],
    )
    w_stats_df.index.names = ['window_mid']

    return w_het_df, w_stats_df



## Main

def main():

    # parse command line arguments
    parse_arguments()

    # make output directory if output_prefix contains '/'
    if '/' in output_prefix:
        if not os.path.exists('/'.join(output_prefix.split('/')[0:-1]) + '/'):
            os.makedirs('/'.join(output_prefix.split('/')[0:-1]) + '/')

    # compile text and stats figure output files (pc figure depends on color taxon --> see below)
    w_het_tsv_path =        output_prefix + '.w_het'  + '.tsv.gz'
    w_stats_tsv_path =      output_prefix + '.w_stats' + '.tsv.gz'
    w_stats_fig_html_path = output_prefix + '.w_stats' + '.html'
    w_stats_fig_pdf_path =  output_prefix + '.w_stats' + '.pdf'

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

    # check if IDs are unique
    if not len(metadata_df['id']) == len(set(metadata_df['id'])):
        print(
            '\n[ERROR] Duplicate sample IDs (first metadata column)\n',
            file=sys.stderr, flush=True,
        )
        sys.exit()

    # if there is output from a previous run, use it
    if os.path.exists(w_het_tsv_path) and os.path.exists(w_stats_tsv_path):
        print(
            '\n[INFO] Plotting data from previous run',
            file=sys.stderr, flush=True,
        )
        w_het_df = pd.read_csv(
            w_het_tsv_path,
            sep='\t',
            index_col=[0],
            na_values='NA',
        )
        w_het_df.columns = [float(x) for x in w_het_df.columns] # change column name dtype to int
        w_stats_df = pd.read_csv(
            w_stats_tsv_path,
            sep='\t',
            index_col=[0],
            na_values='NA',
        )

    else:
        print(
            '\n[INFO] Calculating windowed heterozygosity rate\n',
            file=sys.stderr, flush=True,
        )

        # calculate windowed heterozygosity rate
        w_het_df, w_stats_df = windowed_het(
            variant_file_path,
            chrom, start, stop,
            metadata_df,
            w_size, w_step,
            count_hets,
        )

        # save output data before annotation if not already present
        print(
            '\n[INFO] Writing output TSVs',
            file=sys.stderr, flush=True,
        )
        w_het_df.to_csv(
            w_het_tsv_path,
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

    # pivot windowed heterozygosity rate output and annotate with metadata
    from modules.utils import annotate
    w_het_anno_df = annotate(
        w_het_df,
        metadata_df,
        'heterozygosity_rate',
    )

    # free up memory
    del w_het_df

    # plot windowed heterozygosity rate output & save
    print(
            '\n[INFO] Generating output HTMLs & PDFs',
            file=sys.stderr, flush=True,
    )
    from modules.utils import plot_per_sample_values
    for c_taxon in color_taxon.split(','): 

        # compile output paths
        w_het_fig_html_path = output_prefix + '.w_het.' + c_taxon + '.html'
        w_het_fig_pdf_path =  output_prefix + '.w_het.' + c_taxon +  '.pdf'

        # plot & save
        w_het_fig = plot_per_sample_values(
            w_het_anno_df,
            'heterozygosity_rate',
            c_taxon,
            chrom, start, stop,
            w_size, w_step,
        )
        w_het_fig.write_html(
            w_het_fig_html_path
        )
        w_het_fig.write_image(
            w_het_fig_pdf_path,
            engine='kaleido', scale=2.4,
        )

    # free up memory
    del w_het_fig

    # plot window stats & save
    from modules.utils import plot_het_w_stats
    w_stats_fig = plot_het_w_stats(
        w_stats_df,
        chrom, start, stop,
        w_size, w_step,
        config.min_var_per_w,
        variant_file_sample_lst,
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