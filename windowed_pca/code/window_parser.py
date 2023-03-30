## Dependencies
import sys
import gzip



## Parsing called genotypes from genotype file or from VCF

def win_gt_file(gt_file_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
                skip_monomorphic=False):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) genotype file
    (limited to one chromosome)
    '''

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
            '[INFO] Processed ' + '' + str(w_idx) + ' of ' + str(n_windows) + ' windows',
              file=sys.stderr, flush=True,
            )

        return w_start, w_stop, w_idx, win


    # calculate total number of windows
    n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open uncompressed or gzip-compressed input file
    read_func = gzip.open if gt_file_path.endswith('.gz') else open
    with read_func(gt_file_path, 'rt') as gt_file:

        # fetch sample ids from header and derive sample index positions
        gt_file_sample_lst = gt_file.readline().strip().split('\t')[2:]
        sample_idx_lst = [gt_file_sample_lst.index(x) for x in target_sample_lst]

        # initiate first window
        w_start = start
        w_stop = w_start + w_size-1
        w_idx = 0
        win = []

        # traverse input file
        for line in gt_file:
            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            gts = [line[2:][idx] for idx in sample_idx_lst]

            # skip other than the specified chromosome
            if q_chrom != chrom: continue

            # skip monomorphic sites if specified
            if skip_monomorphic and len(set(gts)) == 1: continue

            # case: pos exceeds current window
            while w_stop < pos:
                
                # if window contains variants: apply function
                if win: func(win, w_start, w_size)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)
                

            # append pos (and genotypes) to current window if larger than window start
            if pos > w_start: win.append([pos] + gts)

            # if end of window is reached: apply function & initialize
            if w_stop <= pos:
                func(win, w_start, w_size)
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


def win_vcf(vcf_path, chrom, start, stop, target_sample_lst, w_size, w_step, func, \
            skip_monomorphic=False):
    '''
    Apply a target function to windows of variants in an (optionally gzipped) genotype file 
    (limited to one chromosome)
    '''

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
            '[INFO] Processed ' + '' + str(w_idx) + ' of ' + str(n_windows) + ' windows',
              file=sys.stderr, flush=True
        ) 

        return w_start, w_stop, w_idx, win

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
    n_windows = len(list(range(start, stop-w_size+2, w_step)))

    # open iput file
    read_func = gzip.open if vcf_path.endswith('.gz') else open
    with read_func(vcf_path, 'rt') as vcf:

        for line in vcf:
            if line.startswith('#CHROM'):

                # fetch sample ids from header and derive sample index positions
                vcf_sample_lst = line.strip().split('\t')[9:]
                sample_idx_lst = [vcf_sample_lst.index(x) for x in target_sample_lst]
                break

    # initiate first window
    w_start = start
    w_stop = w_start + w_size-1
    w_idx = 0
    win = []

    with read_func(vcf_path, 'rt') as vcf:

        # traverse input file
        for line in vcf:
            if line.startswith('#'): continue

            line = line.strip().split('\t')
            q_chrom = line[0]
            pos = int(line[1])
            filter = line[6]
            gts = [line[9:][idx].split(':')[0] for idx in sample_idx_lst]

            # skip other than the specified chromosome
            if q_chrom != chrom: continue

            # keep only PASS sites
            if filter != 'PASS': continue

            # skip monomorphic sites if specified
            if skip_monomorphic and len(set(gts)) == 1: continue            

            gts = [gt_code_dct[x] for x in gts]

            # case: pos exceeds current window
            while w_stop < pos:
                
                # if window contains variants: apply function
                if win: func(win, w_start, w_size)
                if stop < w_stop: break
                w_start, w_stop, w_idx, win = init_win(
                    w_start, w_stop,
                    w_idx, win,
                    w_size, w_step)

            # append pos (and genotypes) to current window
            if pos > w_start: win.append([pos] + gts)

            # if end of window is reached: apply function & initialize
            if w_stop <= pos:
                func(win, w_start, w_size)
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


## Parsing genotype likelihoods from BEAGLE file

# ...
