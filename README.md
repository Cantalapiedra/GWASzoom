# GWASzoom
An R script to plot association values of markers from 2 different datasets

The script accepts a list of markers as input. For each of those markers, retrieves all other markers and genes in its region, producing a plot showing the association values, a figure with "graphical genotypes and phenotype", and a table with markers and genes.

GWASzoom requires 10 input parameters to work:
1) The list of main markers. The region of each of these markers will be inspected, and output files will be generated for each of these markers.
2) A file in Hapmap format, with information of markers from dataset 1, including map position, genotyping, etc. The typical Hapmap format has the next columns: rs	alleles	c	pos	strand	assembly	center	protLSID	assayLSID	panelLSID	QCCode	{Genotypes,...}. Check the Hapmap format elsewhere for more information. Note that to be used with GWASzoom the header must no contain "#", which is typical of many Hapmap files. Just remove the "#" from its header before running GWASzoom.
3) A tabular file with the values to plot for each marker of dataset 1. It must have 2 columns, with a header including "SNP" and "P.value" column names. Note that this file can be obtained very easily from output produced by GWAS tools (e.g. GAPIT).
Note that the information of markers from parameter 1) should come from data in files 2) and 3).
4) A file in Hapmap format, with information of markers from dataset 2. Check parameter 2) for more info.
5) A tabular file with the values to plot for each marker of dataset 2. Check parameter 3) for more info.
6) A BED file with position information of other features (e.g. genes).
7) A file with an additional field which describes each of the features in the file of parameter 6). For example, if the file in 6) is a BED of positions of genes, the file in parameter 7) could have a description for each of those genes.
8) A number, which indicates the width of the region to inspect, in basepairs. For example, if its value is 10000, the region will be in the interval between -10000 and +10000 around the current marker.
9) Output directory, where GWASzoom will put the output files.
10) A tab-separated file with phenotypic data. Note that the genotypes in this file should match the genotypes in the previous files 3) and 5).

An example of run would be (in bash):

Rscript ./GWASzoom.R GWASmarkers.list ds1_markers.hmp ds1_markers.GWAS ds2_markers_hmp ds2_markers.GWAS genes.bed genes.txt 100000 GWASzoom_output phenotypic_data.tsv

GWASzoom generates 3 main output types, in a few files:
- A plot showing the association values of both marker datasets. This plot is shown in 2 files. The difference between them is that one file shows a smooth line of values of dataset2, whereas the other shows no trend lines of any kind.
- A file showing graphycal genotypes of all the markers, along with the phenotypic profile. Samples are clustered based on genotyping.
- A tabular file, in tab-separated plain format, showing marker and genes information.
