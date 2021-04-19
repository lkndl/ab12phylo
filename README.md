# AB12PHYLO

![PyPI license](https://img.shields.io/pypi/l/ab12phylo?color=green)
![github version](https://img.shields.io/static/v1?label=version&message=0.5.1-beta&color=brightgreen&style=flat)
![PyPI Python version](https://img.shields.io/pypi/pyversions/ab12phylo)

[AB12PHYLO](https://github.com/lkndl/ab12phylo) is an integrated, easy-to-use pipeline for Maximum Likelihood (ML) phylogenetic tree inference from ABI trace and `FASTA` data. 
At its core, AB12PHYLO runs parallelized instances of [RAxML-NG](https://github.com/amkozlov/raxml-ng) (Kozlov et al. 2019) or [IQ-Tree](https://github.com/iqtree/iqtree2) (Nguyen et al. 2015) as well as a BLAST search in a reference database. 
It enables visual, effortless sample identification based on phylogenetic position and sequence similarity, as well as population subset selection aided by metrics like Tajima's D for estimations of ongoing evolution.  

![demo screencap of AB12PHYLO GUI](https://github.com/lkndl/ab12phylo/wiki/images/demo.gif)

AB12PHYLO was developed to identify plant pathogen populations possibly under balancing selection with *Solanum chilense*, especially in the genus *Alternaria*. With multi-gene phylogenies still a widely-used method in spite of the rise of whole-genome sequencing, future application on fungal phytopathogens or host plants might be possible.


## Documentation
Detailed installation and usage instructions can be found in the [github wiki](https://github.com/lkndl/ab12phylo/wiki). Also check the in-line help via `ab12phylo -h` and `ab12phylo-cmd -h`, or write an email to [ab12phylo@gmail.com](mailto:ab12phylo@gmail.com) for further support.

## Installation
AB12PHYLO can be installed using pip or conda:
```shell script
pip install ab12phylo
# or 
conda install -c lkndl ab12phylo
```
Check the [wiki](https://github.com/lkndl/ab12phylo/wiki/Installation) for more details, troubleshooting, installing from source or updating the package.  
As implied above, start the graphical version via `ab12phylo` from the terminal, and invoke the commandline version via `ab12phylo-cmd`.

## Quick start
ABI trace files are the main input for AB12PHYLO. Wellsplate tables can translate to original sample IDs, provided the mapping is identical for all sequenced genes. Reference data may be included in `FASTA` format, and the graphical AB12PHYLO accepts `FASTA` sequences as main input as well.  

Sequence data is derived from ABI trace files using a customisable quality control. Sequence ends are trimmed with a sliding window until a certain number (*8 out of 10 by default*) of bases reach the minimal accepted phred quality score (*between 0 and 60, 30 by default*). Bases with low phred quality are replaced by `N` only if they form a consecutive stretch that is longer than a certain threshold (*5 by default*).  
![quality control, alignment and MSA trimming](https://github.com/lkndl/ab12phylo/wiki/images/trace-trimming.png)  
Only samples available for all genes in sufficient quality are included in the analysis. Trimmed sequences are aligned into single-gene Multiple Sequence Alignments (MSAs), which are then trimmed using Gblocks.  
![quality control, alignment and MSA trimming](https://github.com/lkndl/ab12phylo/wiki/images/trim-align-trim.png)  
Single-gene MSAs are concatenated into a multi-gene MSA, which is used for ML Tree Inference.
![concatenation of single-gene MSAs](https://github.com/lkndl/ab12phylo/wiki/images/concat-single-gene-MSAs.png)  

A BLAST search for species annotation can be run on a local database, or [web BLAST](blast.ncbi.nlm.nih.gov/) results can be included.  
Phylogenetic tree inference includes finding a maximum-likelihood tree from several searches as well as bootstrapping to generate confidence values. AB12PHYLO allows some editing of the resulting tree and can calculate basic population genetics statistics, with the graphical version easier and more capable for these applications.

### A simple `ab12phylo-cmd` example
A simple real-world invocation of commandline AB12PHYLO might look like this:
```shell script
ab12phylo-cmd -abi <seq_dir> \
    -csv <wellsplates_dir> \
    -g <barcode_gene> \
    -rf <ref.fasta> \
    -bst 1000 \
    -dir <results>
```
where:
* `<seq_dir>` contains all input ABI trace files, ending in `.ab1`
* `<wellsplates_dir>` contains the `.csv` mappings of user-defined IDs to sequencer's isolate coordinates
* `<barcode_gene>` was sequenced, see [here](https://github.com/lkndl/ab12phylo/wiki/Commandline-version#genes-and-references) for more info
* `<ref.fasta>` contains full GenBank reference records [like this](https://www.ncbi.nlm.nih.gov/nuccore/AF347033.1?report=fasta&log$=seqview&format=text)
* 1000 `-bst` = `--bootstrap` trees will be generated
* `<results>` is where results will be  

## Dependencies
[Biopython](https://biopython.org/wiki/Download), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [Toytree](https://toytree.readthedocs.io/en/latest/), [Toyplot](https://toyplot.readthedocs.io/en/stable/), [matplotlib](https://matplotlib.org/), [PyYAML](https://pyyaml.org/wiki/PyYAML), [lxml](https://lxml.de/), [xmltramp2](https://pypi.org/project/xmltramp2/), [svgutils](https://github.com/btel/svg_utils), [Pillow](https://pillow.readthedocs.io/en/stable/installation.html), [Requests](https://3.python-requests.org/) and [Jinja2](https://jinja.palletsprojects.com/en/2.11.x/intro/#installation)

## External Tools
The pipeline will use existing installations of the following programs if they are found on the system `$PATH` and not considered outdated. Otherwise, AB12PHYLO will download the latest version on its initial run.
* [RAxML-NG](https://github.com/amkozlov/raxml-ng/) version >=1.0.2
* [IQ-Tree](https://github.com/iqtree/iqtree2)
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) version >=2.9
* an MSA tool: [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/), [MUSCLE](https://www.drive5.com/muscle/downloads.htm) or [T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD) *(clients for an [EMBL service](https://www.ebi.ac.uk/Tools/msa/) included)*
* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) for MSA trimming *(included)*

## References
* Alexey M. Kozlov, Diego Darriba, Tom&aacute;&scaron; Flouri, Benoit Morel, and Alexandros Stamatakis (2019)
**RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** 
*Bioinformatics, btz305* 
doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)  
  
* Nguyen,L. T., Schmidt,H. A., Von Haeseler,A., and Minh,B. Q. (2015)
**IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies.** 
*Molecular Biology and Evolution, 32, 268–274.* 
doi:[10.1093/molbev/msu300](https://doi.org/10.1093/molbev/msu300)
  


  

