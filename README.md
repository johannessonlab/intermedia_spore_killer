# intermedia_spore_killer
Scripts released with the paper Svedberg et al. "Convergent evolution of complex genomic rearrangements in two fungal meiotic drive elements", (Nature Communications) 2018. https://www.nature.com/articles/s41467-018-06562-x

The paper is under review at the moment and I will update with a link on publication.

We are here depositing the scripts that were central to the data analysis in the paper: **monophylyPlot.py** and **repClusters.py**.

## monophylyPlot.py
* Author: Jesper Svedberg (jesper.svedberg@ebc.uu.se)
* Version: 0.1

The purpose of this script is to plot signals of monophyly for a subset of individuals in a population genomic dataset. It uses phylogenetic trees inferred from small windows over a chromosome. These phylogenetic trees can most easily be generated using a pipeline developed by Simon H. Martin, which can be found here: https://github.com/simonhmartin/genomics_general/

These following scripts are used:
https://github.com/simonhmartin/genomics_general/blob/master/VCF_processing/parseVCF.py
This parses SNP data in VCF format and converts them to the "geno" format.

https://github.com/simonhmartin/genomics_general/blob/master/phylo/phyml_sliding_windows.py
This splits up a chromosome into short windows and uses PhyML to generate phylogenetic trees. monophylyPlot.py can then parse the output from this script.

**monophylyPlot.py** takes a subset of individuals in a population genomic dataset and generates plots showing whether these strains form a monophyletic cluster in the phylogenetic tree, or whether just a subset of the strains form a monophyletic cluster. Based on this, recombination events can be inferred. You can also add further individuals in a sort of outgroup and plot where this outgroup clusters with the first group. Given the proper input several plots will be generated in PNG or SVG format. The script used the Python library ETE Toolkit (http://etetoolkit.org/) to parse the phylogenetic trees.

Dependencies:
* Python 2.7
* ete3
* pandas
* seaborn
* matplotlib
* numpy

At the moment this is not a script that is fairly limited in scope and is not very friendly to other datasets than what the one it was developped for. The main limitation is that the tested subgroup has to be four strains. The other limitation is that an outgroup must be specified.

Parameters:

```
  -h, --help            show help
  -t TREEFILE, --treefile TREEFILE
                        Gzipped Newick tree file from phyml_sliding_window.py. (required)
  -s STATSFILE, --statsfile STATSFILE 
                        Stats file from phyml_sliding_window.py. (required)
  -m MONO, --mono MONO  Group to check for monophyly. Names separated by
                        commas, without spaces. (required)
  -e EXTRA, --extra EXTRA
                        Extra strains to check if grouping with monophyletic
                        group. (required)
  -o OUTFORMAT, --outformat OUTFORMAT
                        Output image format. PNG or SVG. Default: PNG
  -x MARKERS, --markers MARKERS
                        List of sites to be marked by vertical lines.
  -w WINDOW, --window WINDOW
                        Number of trees in sliding window plots. Default=20
  -v STEP, --step STEP  Number of trees in sliding window plots. Default=20
  -p PREFIX, --prefix PREFIX
                        Output prefix. Default: same as statsfile
  -d DATAFILE, --datafile DATAFILE
                        CSV exported from earlier run. (Useful when rerunning the script.)
  -c, --varcolor        Allow colors for topology plots to change between
                        script runs.
```

## repClusters.py
* Author: Jesper Svedberg (jesper.svedberg@ebc.uu.se)
* Version: 0.1

This script will divide a genome or genomic region into repetitive or non-repetitive windows and outputs average GC content, Repetitive DNA fraction, cytosine methylation levels and histone methylation levels for these windows.

It requires the following input files:
* A genome in fasta format.
* RepeatMasker output.
* Cytosine methylation data in Bismark CX format.
* Histone methylation data in bedCoverage format.

Dependencies:
* Python 2.7
* BioPython
* numpy
* scipy

Parameters:
```
positional arguments:
  filename              RepeatMasker file name

optional arguments:
  -h, --help            show this help message and exit
  -s SIZE, --size SIZE  Chromosome size. Default: Last line in file.
  -w WINDOW, --window WINDOW
                        Sliding window size. Default: 2000
  -p STEP, --step STEP  Sliding window step length. Default: 2000
  -g GENOME, --genome GENOME
                        Genome fasta file.
  -m METHYLATION, --methylation METHYLATION
                        Bismark CX file
  -c COVERAGE, --coverage COVERAGE
                        Coverage file. Modified bedtools coverage format.
  -o CONTIG, --contig CONTIG
                        Specify one contig to limit the analysis to.
  -b START, --start START
                        Start coordinate.
  -e END, --end END     Stop coordinate.
  -t OUTFILE, --outfile OUTFILE
                        Output file name. Default:
                        [RepeatMasker_file].repClusters.csv
```
