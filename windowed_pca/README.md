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

### Preparing a genotype matrix from a VCF file (biallelic SNPs)
```
sample_vcf=sample.vcf.gz
genotype_matrix=genotype_matrix.tsv
bcftools view -v snps -i 'F_MISSING=0' -m2 -M2 -f PASS $sample_vcf | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' | sed 's|\./\.|-1|g' | sed 's|0/0|0|g' | sed 's|1/1|2|g' | sed 's|0/1|1|g' | sed 's|1/0|1|g' > $genotype_matrix
bgzip -@ 10 $genotype_matrix
```

```
zcat ${genotype_matrix}.gz
chrom   pos     ref     alt     ind_1   ind_2   ind_3   ind_4   ind_5   ind_6
chr1    478     G       A       0       0       0       0       1       0
chr1    484     C       T       0       0       0       0       2       0
chr1    485     G       A       1       0       0       0       0       0
chr1    1221    C       T       0       0       0       0       0       0
chr1    1222    G       A       0       0       1       0       0       0
chr1    1223    C       T       0       0       0       0       0       0
chr1    1224    T       G       0      -1       0       0       2       0
chr1    1225    T       C       0       0       0       0       0       0
chr1    1234    G       A       0       0       0       0       0       0
```


### Preparing a metadata file
```
primary_id  simple_id   supplier_id seq_depth   clade   genus           species     sex location        sublocation 
ind_1       CalMal1     D19-E05     16.4        AstCal  Astatotilapia   calliptera  M   Lake_Malombe    Chimwala
ind_2       CalMal2     D19-E06     17.9        AstCal  Astatotilapia   calliptera  M   Lake_Malombe    Chimwala
ind_3       CalMal3     D19-E07     18.9        AstCal  Astatotilapia   calliptera  M   Lake_Malombe    Chimwala
ind_4       CalMal4     D19-E08     18.1        AstCal  Astatotilapia   calliptera  M   Lake_Malombe    Chimwala
ind_5       CalMal5     D19-E09     16.6        AstCal  Astatotilapia   calliptera  M   Lake_Malombe    Chimwala
ind_6       CalMsk1     D22-G01     17.5        AstCal  Astatotilapia   calliptera  M   Southwest_arm   Msaka
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
python3 sw_pca.py test/input/test_gt_matrix.tsv.gz test/input/test_metadata.tsv test/results/output_prefix chr1 1000000 100000 10000 $taxon $group genus,species,location,primary_id 9 3
```

