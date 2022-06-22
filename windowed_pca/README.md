# Conduct PCAs in in genomic windows (using scikit-allel)
- useful to explore the divergence landscape, and particularly to identify inversion polymorphisms in biallelic variant callsets
- generates PDF and interactive HTML plots (using plotly)
- requires a genotype matrix (0,1,2,-1) and a metadata file (details below)

## Python dependencies
- scikit-allel (https://anaconda.org/conda-forge/scikit-allel)
- plotly (https://anaconda.org/plotly/plotly)
- gzip (https://anaconda.org/conda-forge/gzip)
- numpy (https://anaconda.org/anaconda/numpy)
- pandas (https://anaconda.org/anaconda/pandas)

## Usage

### Preparing a genotype matrix from a VCF file

### Preparing a metadata file

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
python3 sw_pca.py test/input/test_gt_matrix.tsv.gz test/input/test_metadata.tsv test/results/output_prefix chr1 1000000 100000 10000 $taxon $group genus,species,location,primary_id 9 3
```

