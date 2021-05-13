# AB12PHYLO

![PyPI license](https://img.shields.io/pypi/l/ab12phylo?color=green)
![github version](https://img.shields.io/static/v1?label=version&message=0.5.5-beta&color=brightgreen&style=flat)
![PyPI Python version](https://img.shields.io/pypi/pyversions/ab12phylo)

[AB12PHYLO](https://github.com/lkndl/ab12phylo) is an integrated, easy-to-use pipeline for Maximum Likelihood (ML) phylogenetic tree inference from ABI trace and `FASTA` data. 
At its core, AB12PHYLO runs parallelized instances of [RAxML-NG](https://github.com/amkozlov/raxml-ng) (Kozlov et al. 2019) or [IQ-Tree](https://github.com/iqtree/iqtree2) (Nguyen et al. 2015) as well as a BLAST search in a reference database. 
It enables visual, effortless sample identification based on phylogenetic position and sequence similarity, as well as population subset selection aided by metrics like Tajima's D for estimations of ongoing evolution, or definition of haplotypes.

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
When AB12PHYLO is first run, it searches for three important non-python tools on the system: RAxML-NG, IQ-Tree and BLAST+. If they are not installed or outdated, AB12PHYLO can download the latest static binaries from GitHub or the NCBI respectively. Check the [wiki](https://github.com/lkndl/ab12phylo/wiki/Installation) for more details, troubleshooting, installing from source or updating the package.  
As implied above, start the graphical version via `ab12phylo` from the terminal, and invoke the commandline version via `ab12phylo-cmd`.

## Quick start
ABI trace files are the main input for AB12PHYLO. Wellsplate tables can translate to original sample IDs, provided the mapping is identical for all sequenced genes. Reference data may be included in `FASTA` format, and the graphical AB12PHYLO accepts `FASTA` sequences as main input as well.  

![main stages of AB12PHYLO](https://github.com/lkndl/ab12phylo/wiki/images/pipeline.png)  

**A:**
Sequence data is extracted from ABI trace files using a customisable quality control. Sequence ends are trimmed with a sliding window until a certain number (*8 out of 10 by default*) of bases reach the minimal accepted phred quality score (*between 0 and 60, 30 by default*). Bases with low phred quality are replaced by `N` only if they form a consecutive stretch that is longer than a certain threshold (*5 by default*).  

**B:**
Samples missing for a single locus are discarded for all genes. Trimmed traces as well as reference and `FASTA` sequences are aligned into single-gene Multiple Sequence Alignments (MSAs), which are then each trimmed using Gblocks 0.91b. For multi-gene analyses, the single-gene MSAs are then concatenated into a multi-gene MSA, which is used for ML Tree Inference. Trees are inferred with either RAxML-NG or IQ-Tree 2, with only the latter available for Windows. 

A BLAST search for species annotation can be run on a local database, or via the public NCBI BLAST API. However, importing XML results of a [web BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) should be preferred to remote API calls.  

AB12PHYLO allows some editing of the resulting tree and can calculate basic population genetics statistics, with the graphical version both easier and far more capable for these applications. Please have a look at the wiki pages ([graphical](https://github.com/lkndl/ab12phylo/wiki/Graphical-interface), [commandline](https://github.com/lkndl/ab12phylo/wiki/Commandline-version#results--motif-search)) for more details.

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
[Biopython](https://biopython.org/wiki/Download), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [Toytree](https://toytree.readthedocs.io/en/latest/) <= 1.2.0, [Toyplot](https://toyplot.readthedocs.io/en/stable/), [matplotlib](https://matplotlib.org/), [PyYAML](https://pyyaml.org/wiki/PyYAML), [lxml](https://lxml.de/), [xmltramp2](https://pypi.org/project/xmltramp2/), [svgutils](https://github.com/btel/svg_utils), [Pillow](https://pillow.readthedocs.io/en/stable/installation.html), [Requests](https://3.python-requests.org/), [Beautiful Soup](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup) and [Jinja2](https://jinja.palletsprojects.com/en/2.11.x/intro/#installation)

## External Tools
The pipeline will use existing installations of the following programs if they are found on the system `$PATH` and not considered outdated. Otherwise, AB12PHYLO will download the latest version on its initial run.
* [RAxML-NG](https://github.com/amkozlov/raxml-ng/) version >=1.0.2
* [IQ-Tree 2](https://github.com/iqtree/iqtree2)
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) version >=2.9
* an MSA tool: [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/), [MUSCLE](https://www.drive5.com/muscle/downloads.htm) or [T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD) *(clients for an [EMBL service](https://www.ebi.ac.uk/Tools/msa/) included)*
* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) 0.91b for MSA trimming *(included)*

## References
* Alexey M. Kozlov, Diego Darriba, Tom&aacute;&scaron; Flouri, Benoit Morel, and Alexandros Stamatakis (2019)
**RAxML-NG: A fast, scalable, and user-friendly tool for maximum likelihood phylogenetic inference.** 
*Bioinformatics, btz305* 
doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305)  
  
* Nguyen,L. T., Schmidt,H. A., Von Haeseler,A., and Minh,B. Q. (2015)
**IQ-TREE: A fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies.** 
*Molecular Biology and Evolution, 32, 268–274.* 
doi:[10.1093/molbev/msu300](https://doi.org/10.1093/molbev/msu300)
  


  

