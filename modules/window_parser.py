## Dependencies
import sys, os
import gzip
import numpy as np
from numba import njit, prange, set_num_threads


## Config
from modules import config


## Shared functions

def init_win(w_start, w_stop, w_idx, win, w_size, w_step):
    '''
    Initialize new window by shifting one w_step and dropping obsolete variants from previous 
    window
    '''

    w_start = w_start + w_step
    w_stop = w_start + w_size-1
    w_idx += 1
    win = [x for x in win if x[0] >= w_start]

    print(
        '[INFO] Processed ' + '' + str(w_idx) + ' of ' + str(config.n_windows) + ' windows',
            file=sys.stderr, flush=True,
        )

    return w_start, w_stop, w_idx, win


@njit()
def gt_min_maf_filter(w_gt_arr, min_maf):
    '''
    Drop SNPs with minor allele frequency below specified value
    '''

    # allele count
    n_alleles = 2 * w_gt_arr.shape[1]

    # calculate allel frequencies and multiple with -1 if AF > 0.05 (because input data may not be 
    # polarized by major/minor allel)
    afs = np.sum(w_gt_arr, axis=1) / n_alleles
    afs[afs > 0.5] = 1 - afs[afs > 0.5]
    
    # keep only sites where AF >= min_maf
    w_gt_arr = w_gt_arr[afs >= min_maf]

    return w_gt_arr


def gl_min_maf_filter(w_gl_arr, min_maf, n_threads):
    '''
    Drop SNPs with estimated minor allele frequency below specified value, using PCAngsd code
    '''

    # import minor allele frequency estimation function from PCAngsd
    from pcangsd.shared import emMAF
    from pcangsd.reader_cy import filterArrays


    # PCAngsd
    min_maf_arr = emMAF(w_gl_arr, 200, 1e-4, n_threads)

    # filter by minor allel frequency
    if min_maf > 0.0:
        min_maf_mask = (min_maf_arr >= min_maf) & (min_maf_arr <= (1 - min_maf))

        # Update arrays
        m = np.sum(min_maf_mask)
        tmp_mask = min_maf_mask.astype(np.uint8)
        filterArrays(w_gl_arr, min_maf_arr, tmp_mask)
        w_gl_arr = w_gl_arr[:m,:]
    
    return w_gl_arr, min_maf_arr


def gt_process_win(win, w_start, w_size, min_maf, func):
    '''
    Remove POS info and convert to numpy array, return empty array if there are no variants
    Call target function: func(w_gt_arr, w_start, w_size)
    '''

    # non-empty: trim off pos info, convert to numpy arr, apply min_maf filter
    if win:
        w_gt_arr = np.array([x[1:] for x in win], dtype=np.int8)
        w_gt_arr = gt_min_maf_filter(w_gt_arr, min_maf) if min_maf else w_gt_arr

    # empty: convert to empty numpy arr
    else:
        w_gt_arr = np.empty((0,0))
    
    # call target function
    func(w_gt_arr, w_start, w_size)


def pl_fetch_target_pls(gt_fields, pl_idx, sample_idx_lst):
    '''
    Extract single list of PLs for the target individuals, only return 'biallellic' (~ 3 PL values)
    '''
    pls = [gt_fields[idx].split(':')[pl_idx].split(',') for idx in sample_idx_lst]
    pls = [x for pl in pls if len(pl) == 3 for x in pl]

    return pls


@njit(parallel=True)
def pl_convert_to_gl(w_pl_arr):
    '''
    Convert a 2D array of phred scaled genotype likelihoods (PL) to normalized genotype likelihoods
    (GL) and return a 2D array, omitting the third GL, which is the expected input for PCAngsd
    '''
    
    # fetch dimensions
    n_rows, n_cols = w_pl_arr.shape

    # reshape pl_arr to have separation by sample as first dimension (to vectorize normalization)
    # --> pl_arr_3d dimensions: (samples, variants, 3 pl_values)
    w_pl_arr_3d = w_pl_arr.reshape(n_rows, -1, 3).transpose(1, 0, 2)

    # unphred
    w_pl_arr_3d = np.power(10, -w_pl_arr_3d/10)

    # preallocate output w_gl_arr (PCAngsd expects only first two GL values for memory efficiency)
    w_gl_arr = np.zeros((n_rows, n_cols*2//3), dtype=np.float32)

    # normalize all GLs partitioned by sample in parallel using numba prange
    for idx in prange(0, w_pl_arr_3d.shape[0]): # jit (prange)
        ind_arr = w_pl_arr_3d[idx]
        ind_arr = 1-(ind_arr / ind_arr.sum(axis=1).reshape(n_rows, 1))
        w_gl_arr[:,idx*2], w_gl_arr[:,idx*2+1] = ind_arr[:,0], ind_arr[:,1]

    return w_gl_arr


def pl_process_win(win, w_start, w_size, min_maf, func, n_threads=config.n_threads):
    '''
    remove POS info and convert to numpy array, return empty array if there are no variants
    '''

    # mute STDOUT by redirecting STDOUT tp /dev/null
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w') 

    # non-empty: trim off pos info, convert to numpy arr, apply min_maf filter, convert PL to GL
    if win:
        w_pl_arr = np.array([x[1:] for x in win], dtype=np.int16)
        w_gl_arr = pl_convert_to_gl(w_pl_arr)
        w_gl_arr, min_maf_arr = gl_min_maf_filter(w_gl_arr, min_maf, n_threads)

    # empty: convert to empty numpy arr & apply target function
    else:
        w_gl_arr = np.empty((0,0))
    
    # call target function
    func(w_gl_arr, min_maf_arr, w_start, w_size, n_threads)
    
    # restore sys.stdout
    sys.stdout = old_stdout


def gl_process_win(win, w_start, w_size, min_maf, func, n_threads=config.n_threads):
    '''
    remove POS info and convert to numpy array, return empty array if there are no variants
    '''

    # mute STDOUT by redirecting STDOUT tp /dev/null
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, 'w') 

    # non-empty: trim off pos info, convert to numpy arr, apply min_maf filter, convert PL to GL
    if win:
        w_gl_arr = np.array([x[1:] for x in win], dtype=np.float32)
        w_gl_arr, min_maf_arr = gl_min_maf_filter(w_gl_arr, min_maf, n_threads)

    # empty: convert to empty numpy arr & apply target function
    else:
        w_gl_arr = np.empty((0,0))
    
    # call target function
    func(w_gl_arr, min_maf_arr, w_start, w_size, n_threads)
    
    # restore sys.stdout
    sys.stdout = old_stdout



## Parsing called genotypes (GT) from genotype file

def win_gt_file(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
                skip_monomorphic=False, min_maf=config.min_maf):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) genotype file
    '''

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open uncompressed or gzip-compressed input file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        # fetch sample ids from header and derive sample index positions
        variant_file_sample_lst = variant_file.readline().strip().split('\t')[2:]
        sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]

        # initiate first window
        w_start = start
        w_stop = w_start + w_size-1
        w_idx = 0
        win = []

        # traverse input file
        for line in variant_file:
            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            
            # skip other than the specified chromosome
            if q_chrom != chrom: continue
            
            # fetch genotypes
            gts = [line[2:][idx] for idx in sample_idx_lst]
            
            # skip monomorphic sites if specified
            if skip_monomorphic and len(set(gts)) == 1: continue

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: gt_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                

            # append pos (and genotypes) to current window if larger than window start
            if pos > w_start: win.append([pos] + gts)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                gt_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )


## Parsing called genotypes (GT) from  VCF

def win_vcf_gt(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
            skip_monomorphic=False, min_maf=config.min_maf):
    '''
    Apply a target function to windows of called variants (GT field) in an (optionally gzipped) VCF
    '''

    # define genotype encoding
    gt_code_dct = {
        '0/0': 0,
        '0|0': 0,
        '0/1': 1,
        '0|1': 1,
        '1/0': 1,
        '1|0': 1,
        '1/1': 2,
        '1|1': 2,
        './.': -1,
        '.|.': -1,
        '0/.': -1,
        '0|.': -1,
        './0': -1,
        '.|0': -1,
        '1/.': -1,
        '1|.': -1,
        './1': -1,
        '.|1': -1
    }

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open iput file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        for line in variant_file:
            if line.startswith('#CHROM'):

                # fetch sample ids from header and derive sample index positions
                variant_file_sample_lst = line.strip().split('\t')[9:]
                sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]
                break

    # initiate first window
    w_start = start
    w_stop = w_start + w_size-1
    w_idx = 0
    win = []

    with read_func(variant_file_path, 'rt') as variant_file:

        # traverse input file
        for line in variant_file:
            if line.startswith('#'): continue

            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            filter = line[6]

            # skip other than the specified chromosome
            if q_chrom != chrom: continue

            # keep only PASS sites
            if filter != 'PASS': continue

            # fetch genotypes
            gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]
            gts = [gt_code_dct[x] for x in gts]
            
            # skip monomorphic sites if specified
            if skip_monomorphic and len(set(gts)) == 1: continue            

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: gt_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)

            # append pos (and genotypes) to current window
            if pos > w_start: win.append([pos] + gts)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                gt_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )


## Parsing phred-scaled genotype likelihoods (PL) from PL file

def win_pl_file(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
                min_maf=config.min_maf):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) PL file
    '''

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open uncompressed or gzip-compressed input file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        # fetch sample ids from header
        variant_file_sample_lst = variant_file.readline().strip().split('\t')[2:]

        # derive sample index positions for first of 3 columns per sample, then add the other 2
        sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]
        sample_idx_lst = [[i, i+1, i+2] for i in sample_idx_lst]
        sample_idx_lst = [x for i in sample_idx_lst for x in i] # flatten

        # remove duplicates (GL/PL file has 3 columns/sample) while preserving order
        variant_file_sample_lst = list(dict.fromkeys(variant_file_sample_lst))

        # initiate first window
        w_start = start
        w_stop = w_start + w_size-1
        w_idx = 0
        win = []

        # traverse input file
        for line in variant_file:
            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            
            # skip other than the specified chromosome
            if q_chrom != chrom: continue
            
            # fetch genotypes
            pls = [line[2:][idx] for idx in sample_idx_lst]

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: pl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
            # append pos (and genotypes) to current window if larger than window start
            if pos > w_start: win.append([pos] + pls)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                pl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )


## Parsing phred-scaled genotype likelihoods (PL) from VCF file

def win_vcf_pl(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
            min_maf=config.min_maf):
    '''
    Apply a target function to windows of phred-scaled genotype likelohoods (PL) in an (optionally 
    gzipped) VCF
    '''

    # set number of threads for numba parallel=True
    set_num_threads(config.n_threads)

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open iput file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        for line in variant_file:
            if line.startswith('#CHROM'):

                # fetch sample ids from header and derive sample index positions
                variant_file_sample_lst = line.strip().split('\t')[9:]
                sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]
                break

    # initiate first window
    w_start = start
    w_stop = w_start + w_size-1
    w_idx = 0
    win = []

    with read_func(variant_file_path, 'rt') as variant_file:

        # traverse input file
        for line in variant_file:
            if line.startswith('#'): continue

            line = line.strip().split('\t')
            q_chrom = line[0]
            
            # skip other than the specified chromosome
            if q_chrom != chrom: continue

            # fetch more fields
            pos = int(line[1])
            format = line[8].split(':')

            # keep only sites with PL field
            if not 'PL' in format: continue

            # extract PLs for target samples (using index of 'PL' in format field)
            pls = pl_fetch_target_pls(line[9:], format.index('PL'), sample_idx_lst)
            
            # skip non-biallellic sites
            if pls == []: continue     

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: pl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)

            # append pos (and genotypes) to current window
            if pos > w_start: win.append([pos] + pls)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                pl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )



## Parsing genotype likelihoods (GL) from GL file

def win_gl_file(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
                min_maf=config.min_maf):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) GL file
    '''

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open uncompressed or gzip-compressed input file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        # fetch sample ids from header
        variant_file_sample_lst = variant_file.readline().strip().split('\t')[2:]

        # derive sample index positions for first of 3 columns per sample, then add the other 2
        sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]
        sample_idx_lst = [[i, i+1] for i in sample_idx_lst]
        sample_idx_lst = [x for i in sample_idx_lst for x in i] # flatten

        # remove duplicates (GL/PL file has 3 columns/sample) while preserving order
        variant_file_sample_lst = list(dict.fromkeys(variant_file_sample_lst))

        # initiate first window
        w_start = start
        w_stop = w_start + w_size-1
        w_idx = 0
        win = []

        # traverse input file
        for line in variant_file:
            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            
            # skip other than the specified chromosome
            if q_chrom != chrom: continue
            
            # fetch genotypes
            gls = [line[2:][idx] for idx in sample_idx_lst]

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: gl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
            # append pos (and genotypes) to current window if larger than window start
            if pos > w_start: win.append([pos] + gls)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                gl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )


## Parsing genotype likelihoods (GL) from a BEAGLE file

def win_beagle(variant_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
                min_maf=config.min_maf):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) BEAGLE file as
    output by ANGSD -GL 1 -out genolike -doGlf 2
    '''

    # calculate total number of windows
    config.n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open uncompressed or gzip-compressed input file
    read_func = gzip.open if variant_file_path.endswith('.gz') else open
    with read_func(variant_file_path, 'rt') as variant_file:

        # fetch sample ids from header
        variant_file_sample_lst = variant_file.readline().strip().split('\t')[3:]

        # derive sample index positions for first of 3 columns per sample, then add the second
        sample_idx_lst = [variant_file_sample_lst.index(x) for x in target_sample_lst]
        sample_idx_lst = [[i, i+1] for i in sample_idx_lst]
        sample_idx_lst = [x for i in sample_idx_lst for x in i] # flatten

        # remove duplicates (GL/PL file has 3 columns/sample) while preserving order
        variant_file_sample_lst = list(dict.fromkeys(variant_file_sample_lst))

        # initiate first window
        w_start = start
        w_stop = w_start + w_size-1
        w_idx = 0
        win = []

        # traverse input file
        for line in variant_file:
            line = line.strip().split('\t')
            q_chrom = line[0].rsplit('_', 1)[0]
            pos = int(line[0].rsplit('_', 1)[1])
            
            # skip other than the specified chromosome
            if q_chrom != chrom: continue
            
            # fetch genotypes
            gls = [line[3:][idx] for idx in sample_idx_lst]

            # case: pos exceeds current window
            while w_stop < pos:
                
                # apply min_maf filter if specified and if window contains variants: apply function
                if win: gl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
            # append pos (and genotypes) to current window if larger than window start
            if pos > w_start: win.append([pos] + gls)

            # if end of window is reached: apply min_maf filter, function & initialize
            if w_stop <= pos:
                gl_process_win(win, w_start, w_size, min_maf, func)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                
    # print exit message
    print(
        '\n[INFO] Processed all windows',
        file=sys.stderr, flush=True,
    )