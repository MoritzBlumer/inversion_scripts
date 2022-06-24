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
```
# set $sample_vcf and $genotype_matrix variables
sample_vcf=test_dataset/input/sample.vcf.gz
genotype_matrix=test_dataset/input/genotype_matrix.tsv.gz

# define a list of samples to be included in the genotype matrix (samples must be subset of input VCF samples)
sample_ids="ind_1,ind_2,ind_3,ind_4,ind_5,ind_6,ind_7,ind_8,ind_9"

# derive header for $genotype_matrix from $sample_vcf
bcftools view -h $sample_vcf | awk '$1=="#CHROM"' | cut -f 1,2,4,5,10- | tr -d '#' | gzip -c > $genotype_matrix

# convert VCF rows to $genotype_matrix rows 
# - keep only lines without missing genotype calls
# - keep only biallelic snps that passed all filters
# - drop unnecessary VCF info (FORMAT, INFO columns)
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS $sample_vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' | gzip -c >> $genotype_matrix

# have a look at the first columns of the matrix
zcat ${genotype_matrix}.gz | head -n 15
```

```
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
```
cat test_dataset/input/metadata.tsv
```

```
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
python sw_pca.py <genotype matrix> <metadata> <chromosome name> <chromosome length> <window size> <window step size> <filter column name> <filter column value> <color column name> <variance threshold> <mean threshold> <output prefix>
```

| Argument | Type | Description |
| ----------------------- | --- | -------------------------------- | 
| **genotype matrix**   | str | path to the genotype matrix file |
| **metadata**          | str | path to the metadata file |
| **output prefix**     | str | prefix that will be used for all output files ('test/' would create a new directory as the prefix and all output files would be located therein) |
| **chromosome name**   | str | name of the chromosome, e.g. 'chr1' |
| **chromosome length** | int | length of the chromosome in bp, e.g. '10000000' |
| **window size**       | int | size of the sliding window in bp, e.g. '1000000' |
| **window step** | int | step size of the sliding window in bp, e.g. '100000' |
| **filter column name** | str | set a metadata column name to be used to select individuals to be included in the analysis e.g. 'genus' (see filter column value) |
| **filter column value** | str | select a value to be filtered for in the defined filter column. Setting **filter column name** to 'genus' and **filter column value** to 'Homo' would include all individuals of the genus Homo in the output, and ignore all others |
| **color column name** | str | select a metadata column that will serve to partition included individuals into color groups in the output plots. If selecting e.g. 'genus', all individuals from the same genus will have the same color in the output plots. Specifying a comma-separated list of column names (e.g. 'genus,species'), two versions of each output plot will be produced, that only differ in the color scheme |
| **variance threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "9" |
| **mean threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "3" |


#### Example prompt 
```
python3 windowed_pca.py $genotype_matrix test_dataset/input/metadata.tsv test_dataset/output/ chr1 35000000 5000000 10000 primary_id $sample_ids primary_id 9 3
```


#### Notes:
- genotype matrix: REF/ALT fields are not used, they can be filled with dummy data
- all columns in metadata will be included in hover display in HTML plots
- add description of output files