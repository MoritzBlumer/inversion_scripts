# Conduct PCA in in genomic windows (using scikit-allel)
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
The python script requires 13 positional arguments, which are explained in more detail below:

```
python sw_pca.py <genotype matrix> <metadata> <chromosome name> <chromosome length> <window size> <window step size> <taxonomic unit> <target taxon> <color taxonomic unit> <variance threshold> <mean threshold> <output prefix>
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
| **taxonomic unit** | str | taxonomic unit (must be a column in the metadata) to be used to filter the genotype matrix, e.g. 'clade' |
| **target taxon** | str | target taxon to use, e.g. "Utaka", but can also be a comma separated list of taxa, e.g. "Utaka,Benthic" |
| **color taxonomic unit** | str | taxonomic unit used to color the plot with distinct colors, e.g. "genus". If specifying a comma separated list (e.g. "genus,species"), separate plots that only differ in coloration will be produced for each unit |
| **variance threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "9" |
| **mean threshold** | int | relevant to correct random switching along PC axes, see code for details, if unsure, use "3" |


#### Example prompt 
```
python3 sw_pca.py test/input/test_gt_matrix.tsv.gz test/input/test_metadata.tsv test/results/output_prefix chr1 1000000 100000 10000 $taxon $group genus,species,location,primary_id 9 3
```

