# hybridindex
A method to determine the contribution of two genetic sources of one individual by random sampling

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
