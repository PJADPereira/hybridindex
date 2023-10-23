# Hybrid Index
A method to determine the contribution of two genetic sources of one individual by random sampling alleles, based on a panel of diagnostic markers. The main input files are: (1) a PoPoolation sync file; (2) a marker file. 


The sync file can be produced by creating a mpileup file from a series of genome alignments with [SAMtools](https://www.htslib.org); the mpileup can then be converted to sync using [PoPoolation2](https://sourceforge.net/p/popoolation2).

A sync file has the following organization **(without the header)**:

|chrm|pos|ref|sample1|...|...|...|sampleN|
|----|---|---|-------|---|----|---|----|
|1|108645303|N	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0|
|1  |146878884|	N	|0:0:0:1:0:0	|0:0:0:2:0:0	|0:0:0:1:0:0	|0:0:0:2:0:0	|0:0:0:2:0:0|
|1  |149391231|	N	|0:0:1:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0|
|1  |160845505|	N	|0:0:0:0:0:0	|1:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0|
|1  |178761540|	N	|0:0:0:0:0:0	|1:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0|
|1  |180562339|	N	|2:0:0:0:0:0	|2:0:0:0:0:0	|1:0:0:0:0:0	|1:0:0:0:0:0	|2:0:0:0:0:0|
|1  |181011897|	N	|0:0:0:0:0:0	|0:0:0:0:0:0	|1:0:0:0:0:0	|1:0:0:0:0:0	|2:0:0:0:0:0|
|1  |190696618|	N	|1:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:0:0:0	|0:0:0:1:0:0	|0:0:0:0:0:0|


The marker file must have the following format **(without the header)**:
|chrm|pos|ref|alt|allelefreqA|allelefreqB|
|-|-|-|-|-|-|
|1| 108645303|	G|	A|	0.000004|	0.999997|
|1|	146878884|	G|	C|	0.032059|	0.999998|
|1|	149391231|	C|	G|	0.022431|	0.958408|
|1|	160845505|	A|	C|	0.000004|	0.999996|
|1|	178761540|	A|	G|	0.041247|	0.967842|
|1|	180562339|	A|	G|	0.030440|	0.963393|
|1|	181011897|	A|	T|	0.033485|	0.953533|
|1|	190696618|	G|	A|	0.970488|	0.039636|


The first two columns of the marker file have the coordinates of the diagnostic marker. Columns 3 and 4 show the reference and the alternative alleles. Columns 5 and 6 show the allele frequency of the diagnostic position in reference populations A and B. In this example, allele frequencies for in each position were calculated on low-coverage whole-genome data using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).






To use download the hybridindex.py file or clone the repository 

```
git clone https://github.com/PJADPereira/hybridindex
```

and run the contained hybridindex.py with the --help command to retrieve the possible options

```
python3 hybridindex.py --help
```
```
  usage: hybridindex.py [-h] -r replicates -d distance -m mode -ss set_initial_seed -sf sync_file -mf marker_file -of output_file

  Hybrid Index

  Required arguments:
    -sf sync file with the counts
    -mf allele frequency files of the control group
    -of output file

  Optional arguments:
    -r number of replicates (default = 100)
    -d mininum distance, in basepairs, between two markers, if mode = "all" then this argument is ignored (default = 200000)
    -m mode to select the markers, can be dynamic, static, or all. In all, all markers are considered, in static the first marker after d is considered, in dynamic, a window with half distance is used to select the next marker

```
