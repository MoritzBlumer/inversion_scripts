# Conduct PCAs in in genomic windows (using scikit-allel)

![example_plot](windowed_pca.png)

- useful to explore the divergence landscape, and particularly to identify inversion polymorphisms in biallelic variant callsets
- generates PDF and interactive HTML plots (using plotly)
- input files: genotype matrix (easy to obtain from VCF) and a metadata file (details below)

## Python dependencies
- scikit-allel (https://anaconda.org/conda-forge/scikit-allel)
- plotly (https://anaconda.org/plotly/plotly)
- gzip (https://anaconda.org/conda-forge/gzip)
- numpy (https://anaconda.org/anaconda/numpy)
- pandas (https://anaconda.org/anaconda/pandas)

## Usage

The below instructions function as a tutorial (after cloning the repo and installing the dependencies).
'test_data' contains a sample VCF file ('test_data/input/sample_vcf.gz') and a corresponding metadata file ('test_data/input/metadata.tsv'). The instructions below guide the user through all necessary steps from preparing a genotype matrix to running the windowed_pca.py script.

### Preparing a genotype matrix from a VCF file (biallelic SNPs)
Based on the provided minimal information sample VCF, generate a genotype matrix in the required format.
Below is one possible way to obtain a genotype matrix based on a VCF, but any other approach that results in the same format works.
Columns REF and ALT are not used by windowed_pca.py and can be replaced with dummy data.
```
# set $sample_vcf and $genotype_matrix variables
sample_vcf=test_dataset/input/sample.vcf.gz
genotype_matrix=test_dataset/input/genotype_matrix.tsv.gz

# define a list of samples to be included in the genotype matrix (samples must be subset of input VCF samples)
sample_ids='ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8,ind_9'

# derive header for $genotype_matrix from $sample_vcf
bcftools view -h $sample_vcf | awk '$1=="#CHROM"' | cut -f 1,2,4,5,10- | tr -d '#' | gzip -c > $genotype_matrix

# convert VCF rows to $genotype_matrix rows 
# - keep only lines without missing genotype calls
# - keep only biallelic snps that passed all filters
# - drop unnecessary VCF info (FORMAT, INFO columns)
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS $sample_vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $genotype_matrix
```
Have a look at the first columns of the matrix:
```
zcat ${genotype_matrix}.gz | head -n 15
CHROM   POS     REF     ALT     ind_1   ind_2   ind_3   ind_4   ind_5   ind_6   ind_7   ind_8   ind_9
chr1    10156   A       T       0       0       0       0       0       0       0       0       0
chr1    12814   A       G       0       0       0       0       0       0       0       0       0
chr1    12895   G       A       0       0       0       0       0       0       0       0       0
chr1    12957   G       A       0       0       0       0       0       0       0       0       0
chr1    55607   A       T       0       0       0       0       0       0       0       0       0
chr1    55728   A       T       0       0       0       0       0       0       0       0       0
chr1    55963   G       A       0       0       0       0       0       0       0       0       0
chr1    56234   A       G       0       0       0       0       0       0       0       0       0
chr1    56469   C       T       0       0       0       0       0       0       0       0       0
chr1    56724   T       C       2       1       1       1       0       0       2       2       1
chr1    56796   G       A       0       0       0       0       0       0       0       0       0
chr1    57369   C       T       0       0       0       0       0       0       0       0       0
chr1    57650   G       A       0       0       0       0       0       0       0       0       0
chr1    57865   A       C       0       0       0       0       0       0       0       0       0
```


### Preparing a metadata file
A metadata file is required to provide annotation for the HTML plots, and can also be used to control which samples to include in the windowed PCA analysis, and to assign groups that will have the same color in the plots.
The minimum requirement for the metadata file is that the first column contains unique identifiers for each sample, and is called 'primary_id'. All additional fields are optional.
An example metadata file, which in addition to the primary_id contains info on sequencing coverage, species and inversion state of each sample, is provided:
```
cat test_dataset/input/metadata.tsv
primary_id      coverage        species inversion_state
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
The python script requires 13 positional arguments, which are explained in more detail below:

```
python3 sw_pca.py <genotype matrix> <metadata> <output prefix> <chromosome name> <chromosome length> <window size> <window step size> <filter column name> <filter column value> <color column name> <variance threshold> <mean threshold>
```

| Argument | Type | Description |
| ----------------------- | --- | -------------------------------- | 
| **genotype matrix**   | str | path to the genotype matrix file |
| **metadata**          | str | path to the metadata file |
| **output prefix**     | str | prefix that will be used for all output files ('test/' would create a new directory as the prefix and all output files would be located therein) |
| **chromosome name**   | str | name of the chromosome, e.g. 'chr1' |
| **chromosome length** | int | length of the chromosome in bp, e.g. '35000000' |
| **window size**       | int | size of the sliding window in bp, e.g. '1000000' |
| **window step** | int | step size of the sliding window in bp, e.g. '100000' |
| **filter column name** | str | set a metadata column name to be used to select individuals to be included in the analysis e.g. 'genus' (see filter column value) |
| **filter column value** | str | select a value to be filtered for in the defined filter column. Setting **filter column name** to 'genus' and **filter column value** to 'Homo' would include all individuals of the genus Homo in the output, and ignore all others |
| **color column name** | str | select a metadata column that will serve to partition included individuals into color groups in the output plots. If selecting e.g. 'genus', all individuals from the same genus will have the same color in the output plots. Specifying a comma-separated list of column names (e.g. 'genus,species'), two versions of each output plot will be produced, that only differ in the color scheme |
| **variance threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "9" |
| **mean threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "3" |


#### Sample prompt 
In the below example, the described genotype matrix (test_dataset/input/genotype_matrix.tsv.gz) and metadata file (test_dataset/input/metadata.tsv) are used as input. 'test_dataset/output/' is set as the output prefix, which tells the script to make a new output directory 'test_dataset/output/' (if it doesn't exist), and to create all output files therein. 'chr1' and '35000000' are set for chromosome name and chromosome length. window size is set to 1 Mbp ('1000000'), because the sample dataset is downsampled to 10% of the original data, and a relatively large window size is required to have enough (>100) variants per window. Step size is set to 100,000 bp ('100000'), resulting in overlapping windows. <filter column name> is set to 'primary_id' and <filter column value> to the previously defined $sample_ids variable to provide a list of samples to include. Since $sample_ids contains all samples in the metadata, all samples are included. To include e.g. only samples from species_1, set <filter column value> to 'species' and <filter column name> to 'species_1'.
```
python3 windowed_pca.py $genotype_matrix test_dataset/input/metadata.tsv test_dataset/output/ chr1 35000000 1000000 10000 primary_id $sample_ids primary_id 9 3
```

#### Output files
[...]


#### Notes:
- genotype matrix: REF/ALT fields are not used, they can be populated with dummy data
- Any biallelic variants can be used as lng as the are encoded as 0 (hom ref), 1 (het), 2 (hom alt). I have used InDels smaller 20 bp before and got nice results
- All columns in metadata will be included in hover display in HTML plots
- If output files (TSVs) from a previous run are detected (same output prefix), they will be reused for plitting instead of calculating new ones. This is useful to adjust the color scheme of the plots. To rerun everything from scratch, delete any existing output files.
- The threshold for the minimum number of SNPs per window is 100 and can be adjusted in the script. The lower the threshold, the noisier the plots.
- By default, plots are only produced for PC1, but it is easy to enable the creation of PC2 plots as well (takes 1 minute, see script for instructions)
- please contact me if you have any questions or run into problems