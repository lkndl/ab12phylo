# AB12PHYLO

![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg) 
![gitlab version](https://img.shields.io/static/v1?label=version&message=0.1b.4&color=blue&style=flat)
![Python version](https://img.shields.io/static/v1?label=python&message=3.8&color=orange&style=flat&logo=python)

[AB12PHYLO](https://gitlab.lrz.de/leokaindl/ab12phylo/) is an integrated, easy-to-use pipeline for phylogenetic tree inference based on Maximum Likelihood (ML) from ABI sequencing data for multiple genes. 
At its core, AB12PHYLO runs parallelized instances of [RAxML-NG](https://github.com/amkozlov/raxml-ng) (Kozlov et al. 2019) and a BLAST search in a reference database. 
It enables visual, effortless sample identification and subset selection based on phylogenetic position and metrics like Tajima's D to estimate ongoing non-random evolution.
 
AB12PHYLO was developed to identify populations of fungal plant pathogen isolates possibly under balancing selection with *Solanum chilense*, especially in the genus *Alternaria*. As multi-gene phylogenies remain a widely-used method in spite of the rise of whole-genome sequencing, future use for fungal phytopathogens or sequencing data from host plants appears realistic not unlikely.


## Installation
 
First, please clone the AB12PHYLO GitLab repository:
 
```bash
git clone https://gitlab.lrz.de/leokaindl/ab12phylo.git
```
All [external tools](#external-tools) are in the [Bioconda](https://anaconda.org/bioconda/repo) channel, so you could install them to a python3.8 [conda](https://docs.conda.io/) environment in one go:

```bash
conda activate <your_python3_conda_env>
conda install -c bioconda "blast>=2.9.0" raxml-ng gblocks mafft clustalo muscle
```

This will require you to activate <your_python3_conda_env> anytime you want to run AB12PHYLO. Alternatively, installing BLAST+ and your preferred MSA tool will suffice and should be easy.
 
Install AB12PHYLO and all its [dependencies](#dependencies) via `pip` or `pip3` to your regular python3:
 
 ```bash
cd ab12phylo
pip install --upgrade pip
pip install .
```


## Getting Started

As AB12PHYLO is primarily a command line tool, you might want to take a look at its interface and options by running `ab12phylo -h`.


#### Test run
This pipeline comes with its own test data set. If you pass `-test`, it will read options from an auxiliary (or backup) [config file](#config-file) at `<ab12phylo_root>/ab12phylo/config/test_config.yaml` and run on these. The test run is set to `--verbose` and will run `--no_remote` BLAST search.


#### Basic options
A simple real-world invocation might look like this:

```bash
ab12phylo -abi <seq_dir> \
    -csv <wellsplates_dir> \
    -g <barcode_gene> \
    -rf <ref.fasta> \
    -bst 1000 \
    -dir <results>
```
where:
* `<seq_dir>` contains all input ABI trace files, ending in `.ab1`
* `<wellsplates_dir>` contains the `.csv` mappings of user-defined IDs to sequencer's isolate coordinates
* `<barcode_gene>` was sequenced
* `<ref.fasta>` contains full GenBank reference records [like this](https://www.ncbi.nlm.nih.gov/nuccore/AF347033.1?report=fasta&log$=seqview&format=text)
* 1000 `-bst` = `--bootstrap` trees will be generated
* `<results>` is where results will be  


#### Detailed settings
AB12PHYLO has reasonably smart defaults while allowing fine-grained access to its settings:

```bash
ab12phylo -rf <ref.fasta> \
    -dbpath <blastdb_dir> \
    -abiset <whitelist> \
    -algo <mafft-clustalo-muscle-tcoffee> \
    -gbl relaxed \
    -bst 1000 \
    -st [32,16]  \
    -s 4 \
    -v -skip
```
* **default:** AB12PHYLO will search for `.ab1` and `.csv` files in or below the current working directory
* **default:** use the `./results` subdirectory
* **default:** samples originate from `ITS1F`
* `<blastdb_dir>` is the path to a ready-to-use BLAST+ database
* only trace files in `<whitelist>` [subset](#subset-analysis) will be read
* `-algo` will generate the MSA: `mafft`, `clustalo`, `muscle` or `t_coffee`
* `-gbl` sets `Gblocks` MSA trimming mode: `skip`, `relaxed` or `strict`
* `-st`: ML tree searches from `32` random and `16` parsimony-based starting trees
* `-s` or sets the random [`--seed`](#seed) `4` for reproducibility
* `-v` or `--verbose` will print all logged events to the console
* `-skip` online BLAST for sequences not in the local BLAST+ db (read [why](#blast-api))


#### Results + Motif Search
Once ML tree inference, computing of support values and BLAST has finished, the pipeline will display a `result.html` in your web browser. This page contains a form that allows **Motif search** across node attributes, enables [subtree or subset selection](#subset-analysis) and calculates diversity metrics for you. 

If results are moved or sent, motif search will be possible by starting a CGI server in the directory via `python3 -m http.server --cgi <port>` or using `ab12phylo-view`.


#### ab12phylo-visualize + ab12phylo-view
`ab12phylo-visualize` will re-plot phylogenies and render a new `results.html`. An end user may use this to switch [support values](#support-values) or plot an [MSA visualization](#msa-visualization). 
 `ab12phylo-view` shows results of a previous run in a browser, with motif search enabled. Both commands accept a path to the AB12PHYLO results or default to `.`, and are equivalent to appending to the original `ab12phylo` call.

```bash
ab12phylo-view <result_folder>
# or
cd <results_folder>
ab12phylo-view
# or 
ab12phylo -c <my-config.yaml> -bst 1000 (...) -view
```

## Advanced use

#### Config File
AB12PHYLO comes with a config file in [YAML](https://yaml.org/) format. Adapt or replace it as you please. Also, [tmux](https://askubuntu.com/questions/8653/how-to-keep-processes-running-after-ending-ssh-session/220880#220880) and `--headless` are highly recommended for remote runs.
 

#### Seed 
An integer that is used to initialize the python3 random number generator ([RNG](https://docs.python.org/3/library/random.html)) that will be used primarily for RAxML-NG. If you set it yourself via `-s` (otherwise it's, well,  random), runs with the same seed on the same system with identical numbers of ML tree searches and Bootstraps will generate precisely the same trees and save them in the same files. This is intentional reproducibility.


#### Genes
If data for several genes is supplied, AB12PHYLO will restrict the analysis to samples that are present for all genes. If this causes a lot of samples to be dropped, it might be worth leaving out a gene entirely by setting `-g` = `--genes`.


#### Subset analysis
Alternatively, you can pass a file that serves as a whitelist if you want to analyse only a certain subset of your data via `-abiset`. It is also possible to use sample IDs: Provide a file via `-sampleset`. Defining a subset this way does not influence exclusion via `-g`, of course.


#### References
Setting references is a bit like setting genes: If a directory of reference files is supplied via `-rd` = `--ref_dir`, the pipeline will try to match the `.FASTA` files in that directory to the genes in the analysis *by their filename*. For example, `ITS1F.fasta` will be matched to sequences from the *ITS1F* gene. Alternatively, an ordered list of reference files can be passed via `-rf` = `--ref`, and file names will be ignored. 

If you're feeling neat and precise and set both the genes and individual references, be careful: The pipeline deliberately matches references and genes *by order* in this case. Accordingly, this will mess it up:

```bash
-g ITS1F OPA10 -rf ../opa10.fasta ITS1F.phy
```


#### BLAST+ Database
Sometimes, this pipeline might run headless on a server. To keep it from running head-first in a firewall when it attempts to update its BLAST+ database via FTP, please pre-supply an unzipped, ready-to-use BLAST+ database via `-dbpath` (and `-db` name). Find databases on the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/db/).


#### BLAST API
BLAST API queries are de-prioritised after just a few attempts. Accordingly, if several runs are attempted on the same data set, the `-skip`=`--no_remote` flag can be set to use data from an earlier run or leave out species annotation for samples not found in the local database. Alternatively, BLAST can be skipped entirely with `-none`=`--no_BLAST`.


#### MSA visualization
An MSA visualization can be plotted to the rectangular tree by passing `-msa_viz`. This will take some rendering time for wider alignments.


#### Support Values
For visualization, you can pick either Felsenstein Bootstrap Proportions `FBP` or [Transfer Bootstrap Expectation](https://doi.org/10.1038/s41586-018-0043-0) `TBE` support values with `-metric`. Newick tree files for both types are generated anyway, so switching support value metric only requires re-running `ab12phylo-visualize`.


#### Log File
If you're having trouble, look at the log file! It will be in your results directory and is named `ab12phylo.log`. Alternatively, you can set the `--verbose` flag and get the same information in real-time to your commandline. Your choice.


## Dependencies
[Biopython](https://biopython.org/wiki/Download), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [Toytree](https://toytree.readthedocs.io/en/latest/), [Toyplot](https://toyplot.readthedocs.io/en/stable/), [PyYAML](https://pyyaml.org/wiki/PyYAML), [DendroPy](https://dendropy.org/#installing) [lxml](https://lxml.de/), [xmltramp2](https://pypi.org/project/xmltramp2/) and [Jinja2](https://jinja.palletsprojects.com/en/2.11.x/intro/#installation)


## External Tools
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) version >=2.9
* [RAxML-NG](https://github.com/amkozlov/raxml-ng/) *(optional, currently included)*
* an MSA tool: [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/), [MUSCLE](https://www.drive5.com/muscle/downloads.htm) or [T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD) *(optional, **slow** online clients included)*
* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) for MSA trimming *(optional, currently included)*
