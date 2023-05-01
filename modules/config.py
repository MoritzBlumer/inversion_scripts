## GLOBAL SETINGS (ok to change)
min_var_per_w = 25  # minimum number of polymorphic variants per window
mean_threshold = 3  # if no guide samples specified, use $mean_threshold samples with largest mean 
                    # to select guide samples
float_precision = 7 # float precision for text output


## SHARED VARIABLES (don't change anything below, min_maf can be specified in the command line)
n_windows = None
min_maf = 0.01
pc = None
n_threads = 2
plot_scaling_factor = 20000
