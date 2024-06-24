#!/usr/bin/env python
#
# Moritz Blumer | 2024-06-23
#
# Plot BAM contig or whole genome alignments to a reference sequence (using pysam)


## File info
__author__ = 'Moritz Blumer, 2023'
__email__ = 'lmb215@cam.ac.uk'


## Dependencies
import sys, os
import pysam
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import math

def parse_arguments():
    '''
    Parse command line arguments & print help message if # of arguments is incorrect
    '''

    global input_path, output_path, faidx_path, seq_prefix

    # print help message if incorrect number of arguments was specified
    if len(sys.argv) < 5:
        print(
            '   python plot_bam_alignments.py <input path> <output path> <reference seq name> \
                                    <query seq name>\n\n\
            <input path>          str  path to input SAM/BAM\n\
            <output path>         str  path to output PDF\n\
            <faidx path>          str  FAIDX of reference (used to fetch chromosomes and sizes)\n\
            <seq prefix>          str  prefix of sequences in FAIDX to use (e.g. "chr")',
        file=sys.stderr,
        )
        sys.exit()

    # fetch arguments    
    _, input_path, output_path, faidx_path, seq_prefix = sys.argv

    # check if input file was specified and if it exists:
    if not os.path.isfile(input_path):
        print('\n[ERROR] Input file does not exist.', file=sys.stderr)
        sys.exit()


def init():
    '''
    Open input BAM and check whether in_bam was simplified (simplify_bam.py) or inversion-refined 
    (bam_refine_inv_wga.py) by checking the PG lines in the header. 
    If one of these programs has been was run on BAM, read-length info is not retained in the 
    CIGAR and is instead written to the 9th (TLEN) field.
    '''

    # open input file and open output file:
    in_bam = pysam.AlignmentFile(input_path, "rb")

    # check if in_bam was simplified (bam_simplify.py) or inversion-refined (bam_refine_inv_wga.py)
    # by checking the PG lines in the header. If one of these programs was run on BAM, read-length
    # info is not retained in the CIGAR and is instead written to the 9th (TLEN) field.
    q_seq_len_from_tlen = False
    header = in_bam.header.as_dict()
    for line in header['PG']:
        if line['ID'] == 'simplify_bam.py' or line['ID'] == 'bam_refine_inv_wga.py' \
            or line['ID'] == 'bam_simplify.py':
            q_seq_len_from_tlen = True
    
    return in_bam, q_seq_len_from_tlen


def parse_alignments(in_bam, q_name, q_seq_len_from_tlen):
    '''
    Fetch all query alignments to specified reference sequence
    '''

    # initiate objects:
    alignments_lst = []
    q_seq_len_lst = []

    # iterate input BAM:
    for alignment in in_bam.fetch():

        if alignment.query_name == q_name and not alignment.flag & 0x100:

            # consider hard-clipped portions:
            soft_clip_start = 0
            soft_clip_end = 0
            hard_clip_start = 0
            hard_clip_end = 0

            # no hard clipping:
            if alignment.cigartuples[0][0] == 4:
                soft_clip_start = alignment.cigartuples[0][1]
            if alignment.cigartuples[-1][0] == 4:
                soft_clip_end = alignment.cigartuples[-1][1]

            # hard clipping
            if alignment.cigartuples[0][0] == 5:
                hard_clip_start = alignment.cigartuples[0][1]
                if alignment.cigartuples[1][0] == 4:
                    soft_clip_start = alignment.cigartuples[1][1]
            if alignment.cigartuples[-1][0] == 5:
                hard_clip_end = alignment.cigartuples[-1][1]
                if alignment.cigartuples[-2][0] == 4:
                    soft_clip_end = alignment.cigartuples[-2][1]
            
            # fetch alignment strand
            strand = '-' if alignment.is_reverse else '+'

            # reverse CIGAR (hard- and soft-clipped bases) for (-) alignments:
            if strand == '-':
                soft_clip_start, soft_clip_end = soft_clip_end, soft_clip_start
                hard_clip_start, hard_clip_end = hard_clip_end, hard_clip_start

            # add alignment dict to alignments
            alignments_lst.append(
                {
                    'q_start': hard_clip_start + soft_clip_start,
                    'q_end': hard_clip_start + soft_clip_start + alignment.query_alignment_length,
                    'q_midpoint': hard_clip_start + soft_clip_start \
                                    + (0.5 * alignment.query_alignment_length),
                    'q_strand': strand,
                    'r_start': alignment.reference_start,
                    'r_end': alignment.reference_end,
                    'r_midpoint': alignment.reference_start + 0.5 * alignment.reference_length,
                    'r_name': alignment.reference_name
                }
            )
            if q_seq_len_from_tlen:
                q_seq_len = alignment.template_length
            else: 
                q_seq_len = hard_clip_start + soft_clip_start + alignment.query_alignment_length \
                    + soft_clip_end + hard_clip_end
            q_seq_len_lst.append(q_seq_len)

    # check if the same (unaligned) query sequence length can be reconstructed from all query 
    # alignments (otherwise sth. is wrong)
    if len(set(q_seq_len_lst)) != 1:
        print('\n[ERROR] CIGAR string corrupted.\n', file=sys.stderr)
        sys.exit()
    else:
        q_seq_len = q_seq_len_lst[0]

    return alignments_lst, q_seq_len


def plot_alignments(ax, alignments_lst, r_name, r_start, r_stop, q_seq_len, fwd_aln_col, \
                    rev_aln_col, diff_chrom_col):
    '''
    Plot all query alignments to the reference either forward or reverse
    '''

    # calculate query offset (so that query and reference midpoint line up)
    r_len = r_stop - r_start
    q_offset = (r_start + (r_len/2))-(q_seq_len/2)

    # set axis limits
    min_x = min([r_start, q_offset])
    max_x = max([r_stop, q_seq_len-abs(q_offset)])
    ax_offset = 0#0.05 * max_x
    ax.set_xlim(min_x - ax_offset, max_x + ax_offset)
    ax.set_ylim(-10, 215)

    # set axis labels
    #ax.set_xlabel('Genomic Position', fontsize=10)
    ax.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    ax.axes.get_yaxis().set_visible(False)

    # fill entire reference & query sequence
    ax.fill(
        (
            q_offset ,
            q_offset + q_seq_len,
            q_offset + q_seq_len,
            q_offset,
        ),
        (200, 200, 150, 150),
        color=fwd_aln_col, alpha=1,
        linewidth=0,
    )
    ax.fill(
        (
            0,
            r_len,
            r_len,
            0,
        ),
        (50, 50, 0, 0),
        color=fwd_aln_col, alpha=1,
        linewidth=0,
    )

    for alignment in alignments_lst:

        # aligned to target chromosome
        if alignment['r_name'] == r_name:
            color = fwd_aln_col if alignment['q_strand'] == '+' else rev_aln_col
            if alignment['q_strand'] == '-': alignment['q_start'], alignment['q_end'] = alignment['q_end'], alignment['q_start']

            # plot query rectangles
            ax.fill(
                (
                    q_offset + alignment['q_start'],
                    q_offset + alignment['q_end'],
                    q_offset + alignment['q_end'],
                    q_offset + alignment['q_start'],
                ),
                (200, 200, 150, 150),
                color=color, alpha=1,
                linewidth=0,
            )
            
            # plot ref rectangles
            ax.fill(
                (
                    alignment['r_start'],
                    alignment['r_end'],
                    alignment['r_end'],
                    alignment['r_start'],
                ),
                (50, 50, 0, 0),
                color=color, alpha=1,
                linewidth=0,
            )

            # plot connection polygons:
            ax.fill(
                (
                    q_offset + alignment['q_start'],
                    q_offset + alignment['q_end'],
                    alignment['r_end'],
                    alignment['r_start'],
                ),
                (150, 150, 50, 50),
                color=color, alpha=0.3,
                linewidth=0,
            )
            
            # # plot connection lines:
            # ax.plot(
            #     (
            #         q_offset + alignment['q_midpoint'],
            #         alignment['r_midpoint'],
            #     ),
            #     (150, 50),
            #     color='k', alpha=0.3,
            #     linewidth=0.005,
            #     solid_capstyle='round',
            # )
            
        # aligned to different chromosome
        if alignment['r_name'] != r_name:
            color = diff_chrom_col
            ax.fill(
                (
                    q_offset + alignment['q_start'],
                    q_offset + alignment['q_end'],
                    q_offset + alignment['q_end'],
                    q_offset + alignment['q_start'],
                ),
                (200, 200, 150, 150),
                color=color, alpha=1,
                linewidth=0,
            )
            
    # plot start and end of target and query as tick lines
    ax.fill(
        (
            r_start - 100000,
            r_start,
            r_start,
            r_start - 100000,
        ),
        (51, 51, -1, -1),
        color = 'k', alpha=1,
        linewidth=0,
    )
    ax.fill(
        (
            r_stop + 100000,
            r_stop,
            r_stop,
            r_stop + 100000,
        ),
        (51, 51, -1, -1),
        color = 'k', alpha=1,
        linewidth=0
    )
    ax.fill(
        (
            q_offset - 100000,
            q_offset,
            q_offset,
            q_offset - 100000
        ),
        (201, 201, 149, 149),
        color = 'k', alpha=1,
        linewidth=0,
    )
    ax.fill(
        (
            q_offset + q_seq_len + 100000,
            q_offset + q_seq_len,
            q_offset + q_seq_len,
            q_offset + q_seq_len + 100000
        ),
        (201, 201, 149, 149),
        color = 'k', alpha=1,
        linewidth=0,
    )

    # remove frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    # reformat x labels to Mbp
    #ax.set_xticks([x for x in ax.get_xticks()])
    ax.set_xticklabels([int(x/1000000) for x in ax.get_xticks()])

    return ax


## Main

# parse command line arguments
parse_arguments()

# read in chromosome sizes from reference FAIDX file and calculate offsets in bp (FAIDX offset is 
# in bytes)
chrom_sizes_dct = pd.read_csv(
    faidx_path, sep='\t', 
    usecols=[0, 1], 
    names=['chrom', 'size'], 
    index_col=0
).to_dict(orient='index')

# fetch chromosomes
chroms=[x for x in chrom_sizes_dct.keys() if x.startswith('chr')]

# plot settings
n_rows = int(math.ceil(len(chroms)/2))
n_cols = 2
n_plots = n_rows * n_cols
fwd_aln_col = '#757575'
rev_aln_col = '#ed5f5f'
diff_chrom_col = '#757575'
plt.rcParams["figure.figsize"] = (20, 30)
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"

# plot
fig, axs = plt.subplots(n_rows, n_cols, sharey=True)

for i, chrom in zip(range(n_plots), chroms):

    # print info
    print('\n[INFO] Working on ' + chrom, file=sys.stdout, flush=True)

    # fetch ax
    col = i // n_rows
    row = i % n_rows
    ax = axs[row, col]

    # set variables
    r_name = chrom
    r_start = 1
    r_stop = chrom_sizes_dct[chrom]['size']

    # plot ax
    in_bam, q_seq_len_from_tlen = init()
    alignments_lst, q_seq_len = parse_alignments(in_bam, chrom, q_seq_len_from_tlen)
    ax = plot_alignments(ax, alignments_lst, chrom, r_start, r_stop, q_seq_len, fwd_aln_col, \
                         rev_aln_col, diff_chrom_col)
    ax.set_title(chrom, loc='left', fontweight='bold', pad=-10)

# change to tight layout
fig.tight_layout()

# safe
plt.savefig(output_path, transparent=True, dpi=600)

# print info
print('\n[INFO] Done.\n', file=sys.stdout, flush=True)
