# Conduct PCAs in in genomic windows (using scikit-allel)

![example_plot](.media/windowed_pca.png)

# PLEASE NOTE: 
This is an improved version of earlier versions. It has been tested extensively and should produce identical results, but please let me know if you notice any issues
The main improvements include:
- instead of loading the entire genotype file into memory, only the current window is read, which should greatly reduce memory requirements (and potentially speed)
- Support for local windowed PCA (regions can now be specified, e.g. chr1:1000000-8000000
- VCF support
- improved per window stats plots
- guide sample specification elevated to proper command line option
- format for the genotype file slightly changed (REF + ALT columns dropped, as they are not used in the analysis) --> keep that in mind when working with existing input data
A copy of the previous version of the script is still available in ```legacy/windowed_pca_v1.py``` 

## Overview
- useful to explore the variation/divergence landscape, and particularly to identify inversion polymorphisms in biallelic variant callsets
- generates PDF and interactive HTML plots (using plotly)
- input files: a biallelic VCF and a metadata file (details below) (instead of a VCF, a genotype file can be used as input; see below) 

## Python dependencies 
- scikit-allel (https://anaconda.org/conda-forge/scikit-allel)
- plotly (https://anaconda.org/plotly/plotly)
- gzip (https://anaconda.org/conda-forge/gzip)
- numpy (https://anaconda.org/anaconda/numpy)
- pandas (https://anaconda.org/anaconda/pandas)
--> easy to obtain using conda

## Usage
The below instructions function as a tutorial (after cloning the repo and getting the dependencies).
```test_data/``` contains a sample VCF file (```test_data/input/sample_vcf.gz```) and a corresponding metadata file (```test_data/input/metadata.tsv```). The instructions below 
guide the user through all necessary steps from preparing a genotype file to running the ```windowed_pca.py``` script.

### Preparing a metadata file
A metadata file is required to provide annotation for the HTML plots, and can also be used to control which samples to include in the windowed PCA analysis, and to assign groups that will have the same color in the plots.
The minimum requirement for the metadata file is that the first column contains unique identifiers for each sample, and is called 'id'. All additional fields are optional.
An example metadata file, which in addition to the id contains info on sequencing coverage, species and inversion state of each sample, is provided:
```
cat test_dataset/input/metadata.tsv
id      coverage        species inversion_state
ind_1   20X     species_1       inverted
ind_2   21X     species_1       inverted
ind_3   20X     species_1       inverted
ind_4   21X     species_1       heterozygous
ind_5   19X     species_1       heterozygous
ind_6   19X     species_1       heterozygous
ind_7   18X     species_2       uninverted
ind_8   26X     species_2       uninverted
ind_9   18X     species_2       uninverted
```

### Running windowed_pca.py
The python script requires 12 positional arguments, which are explained in more detail below:

```
python3 code/windowed_pca.py <variant file> <metadata> <output prefix> <region> <winow size>, <window step>, <pc>, <filter column name>, <filter column value> <color column name> <minor allel frequency> <guide samples>
```

| Argument | Type | Description |
| ----------------------- | --- | -------------------------------- | 
| **variant file**   | str | path to uncompressed or gzipped variant file (VCF or genotype file) |
| **metadata**          | str | path to the metadata file |
| **output prefix**     | str | prefix for output files |
| **region**   | str | target region in format "chr:start-stop" (i.e. chr1:1-$chrom_length to analyze the entire chr1) |
| **window size**       | int | sliding window size in bp, e.g. '1000000' |
| **window step** | int | sliding window step size in bp, e.g. '100000' |
| **minor allel frequency** | float | minor allel frequency threshold; specify "None" to disable filter |
| **pc** | int | principal component to use ('1' or '2') |
| **filter column name** | str | metadata column name to filter for individuals to includede in the analysis, e.g. 'genus' (see <filter column value>) |
| **filter column value** | str | value to be filtered for in filter column; Setting *filter column name* to 'genus' and *filter column value* to 'Homo' would include all individuals of the genus _Homo_ in the output, while ignoring all others. (a comma-separated list of include values is accepted, e.g. 'Homo,Pan') |
| **color column name** | str | metadata column to assign colors by in the output plot; if selecting 'genus', all individuals from the same genus will have the same color in the output plots; if specifying a comma-separated list like 'genus,species', one output plot is generated for each color scheme |
| **guide samples** | str | list of samples to use for polarization, e.g. 'ind1,ind2,ind3'; specify 'None' for automatic guide sample selection |

#### Sample prompt 
```
python3 code/windowed_pca.py test_dataset/input/sample_vcf.gz test_dataset/input/metadata.tsv test_dataset/output/testrun chr1:1-35000000 1000000 10000 None 1 id ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8,ind_9 inversion_state ind_7,ind_8,ind_9
```
- the variant file (```test_dataset/input/sample_vcf.gz``` this could also be ```test_dataset/input/genotype_file.tsv.gz```) and metadata file (```test_dataset/input/metadata.tsv```) are used as input files <genotype 
file> and <metadata>
- ```test_dataset/output/``` is set as the <output prefix>, which tells the script to create a new output directory ```test_dataset/output/``` (if it doesn't exist), and to create all output files therein
- <region> to analyze here is the entire chromosome 1 ('chr1:1-35000000'), but could also be a sub set of the chromosome (e.g. 'chr1:1000000-8000000')
- <window size> is set to 1 Mbp ('1000000'), because the sample dataset is downsampled to 10% of the original SNPs, and a relatively large window size is required to have enough (>100) variants per window
- <step size> is set to 100,000 bp ('100000'), resulting in overlapping windows
- <pc> is set to '1' to analyze principal component 1 (the other option would be '2' for PC 2)
- <filter column name> is set to the metadata column 'id' and <filter column value> to the a list of sample ids to include; to include e.g. only samples from species_1, set <filter column value> to 'species' and <filter column name> to 'species_1')
- <color column name> is set to 'inversion_state', resulting in all individuals with the same inversion states having the same color in the output w_pc_1 plot
- <guide samples> is set to 'ind_7,ind_8,ind_9' to use these specific samples for polarization of the PC axis for plotting


### Preparing a genotype file from a VCF file (biallelic SNPs)
A genotype file like this can be much more space efficient than a full scale VCF. It contains only the chromosome name, position and per sample genotypes encoded as (0=hom ref, 1=het, 2=hom alt, -1=mising data).
Based on the provided minimal information sample VCF, generate a genotype file in the required format.
Below is one possible way to obtain a genotype file based on a VCF, but any other approach that results in the same format works.
(A previous version included a REF and ALT column – these were not used by the script and have been dropped)
```
# set $sample_vcf and $/genotype_file variables
sample_vcf=test_dataset/input/sample.vcf.gz
genotype_file=test_dataset/input/genotype_file.tsv.gz

# define a list of samples to be included in the genotype file (samples must be subset of input VCF samples)
sample_ids='ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8,ind_9'

# derive header for $genotype_file from $sample_vcf
bcftools view -h $sample_vcf | awk '$1=="#CHROM"' | cut -f 1,2,10- | tr -d '#' | gzip -c > $genotype_file

# convert VCF rows to $genotype_file rows 
# - keep only lines without missing genotype calls
# - keep only biallelic snps that passed all filters
# - drop unnecessary VCF info (FORMAT, INFO columns)
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS $sample_vcf | bcftools query -f '%CHROM\t%POS\t[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 
's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $genotype_file
```
Have a look at the first columns of the gentype file:
```
zcat ${genotype_file}.gz | head -n 15
CHROM   POS     ind_1   ind_2   ind_3   ind_4   ind_5   ind_6   ind_7   ind_8   ind_9
chr1    10156   0       0       0       0       0       0       0       0       0
chr1    12814   0       0       0       0       0       0       0       0       0
chr1    12895   0       0       0       0       0       0       0       0       0
chr1    12957   0       0       0       0       0       0       0       0       0
chr1    55607   0       0       0       0       0       0       0       0       0
chr1    55728   0       0       0       0       0       0       0       0       0
chr1    55963   0       0       0       0       0       0       0       0       0
chr1    56234   0       0       0       0       0       0       0       0       0
chr1    56469   0       0       0       0       0       0       0       0       0
chr1    56724   2       1       1       1       0       0       2       2       1
chr1    56796   0       0       0       0       0       0       0       0       0
chr1    57369   0       0       0       0       0       0       0       0       0
chr1    57650   0       0       0       0       0       0       0       0       0
chr1    57865   0       0       0       0       0       0       0       0       0
```

### Output files
Six output files are created by default, which can be grouped as follows:
1. Windowed PCA:
    1. ```${output_prefix}.w_pc_${pc}.tsv.gz``` contains all information relevant to the windowed PCA plots: it provides the value per principal component per window per individual, as well as all metadata and is used for plotting. It is also the file to be used for any custom plotting.
    2. ```${output_prefix}.pc_${pc}.${color_column_name}.html```: interactive HTML plot of the windowed PCA results, based on ```${output_prefix}.${chromosome_name}.pc_1.tsv```. if more than one color_column_name was specified, additional versions of this plot will be produced.
    3. ```${output_prefix}.pc_${pc}.${color_column_name}.pdf```: like the HTML file(s), just as PDF(s)
2. Supplementary info:
    1. ```${output_prefix}.w_stats.tsv.gz```: contains per window stats on the PCA results per window: % explained PC 1, % explained PC 2, # of included sites per window.
    2. ```${output_prefix}.w_stats.html```: interactive HTML plot of this per windows stats
    3. ```${output_prefix}.w_stats.pdf```: same as interactive HTML, just PDF




### Notes:
- Any biallelic variants can be used as lng as the are encoded as 0 (hom ref), 1 (het), 2 (hom alt). I have used InDels smaller 20 bp (instead of SNPs) before and got nice results
- All columns in metadata will be included in hover display in HTML plots
- If output files (TSVs) from a previous run are detected (same output prefix), they will be reused for plotting instead of calculating new ones. This is useful to adjust the color scheme of the plots. To rerun everything from scratch, delete existing output files.
- The threshold for the minimum number of SNPs per window is 50 and can be adjusted at the top of windowed_pca.py.
- A major challenge for windowed PCAs is that the PCA in each window is conducted independently, and rotation is random. This results in random switching of window polarization along the component axis in the raw data, and very noisy plots. The script implements a correction step that aims at polarizing the orientation based on variance and mean of the previous window, using a specified number of the individuals with the most extreme values across the entire dataset (parameters <variance threshold> and <mean threshold>) as guide samples. This gives acceptable results for many tested datasets, but is not always ideal. After the first round of plotting, specific samples can often be visually determined which would make good guide samples (i.e. samples that are always >0 or always <0. A single guide sample can be sufficient for correct polarization.
- please contact me if you have any questions or run into problems
