# AB12PHYLO

![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg) 
![gitlab version](https://img.shields.io/static/v1?label=version&message=0.1b.3&color=blue&style=flat)
![Python version](https://img.shields.io/static/v1?label=python&message=3.8&color=orange&style=flat&logo=python)

[`AB12PHYLO`](https://gitlab.lrz.de/leokaindl/ab12phylo/) is an integrated, easy-to-use pipeline for phylogenetic tree inference based on Maximum Likelihood (ML) from ABI sequencing data for multiple genes. 
At its core, AB12PHYLO runs parallelized instances of [RAxML-NG](https://github.com/amkozlov/raxml-ng) (Kozlov et al. 2019) and a BLAST search in a reference database. 
It enables visual, effortless curation of phylogenies and provides population genetics measurements.
 
AB12PHYLO is intended for research into plant pathogens and especially *Alternaria*, but future versions might allow species annotation from BLAST databases other than `ITS_RefSeq_Fungi` and unlock its use for plants or other, [real animals](https://xkcd.com/1749/).


 ## Installation
 
 First, clone the AB12PHYLO repository from GitLab:
 
```bash
git clone https://gitlab.lrz.de/leokaindl/ab12phylo.git
```
All [external tools](#external-tools) are in the [Bioconda](https://anaconda.org/bioconda/repo) channel, so you could install them to a python3.8 [conda](https://docs.conda.io/) environment in one go:

```bash
conda activate <your_python3_conda_env>
conda install -c bioconda "blast>=2.9.0" raxml-ng gblocks mafft clustalo muscle
```

This will require you to activate <your_python3_conda_env> anytime you want to run AB12PHYLO. Alternatively, installing BLAST+ and your preferred MSA tool will suffice and should be easy.
 
 Finally, install AB12PHYLO and all its [dependencies](#dependencies) via `pip` or `pip3` to your regular python3:
 
 ```bash
cd ab12phylo
pip install --upgrade pip
pip install .
```

[tmux](https://askubuntu.com/questions/8653/how-to-keep-processes-running-after-ending-ssh-session/220880#220880) 
is highly recommended for remote runs.
 

#### Dependencies
[Biopython](https://biopython.org/wiki/Download), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [Toytree](https://toytree.readthedocs.io/en/latest/), [Toyplot](https://toyplot.readthedocs.io/en/stable/), [PyYAML](https://pyyaml.org/wiki/PyYAML), [DendroPy](https://dendropy.org/#installing) and [Jinja2](https://jinja.palletsprojects.com/en/2.11.x/intro/#installation)


#### External Tools
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) version >=2.9
* [RAxML-NG](https://github.com/amkozlov/raxml-ng/) *(optional, currently included)*
* an MSA tool: [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/), [MUSCLE](https://www.drive5.com/muscle/downloads.htm) or [T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD) *(optional, **slow** online clients included)*
* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) for MSA trimming *(optional, currently included)*


## Getting Started

As AB12PHYLO is a command line tool, you might want to take a look at its available options by running `ab12phylo -h`.


#### Test run
This pipeline comes with its own test data set. If you pass `-test`, it will read its options from an auxiliary (or backup) [config file](#config-file) at `<ab12phylo_root>/ab12phylo/config/test_config.yaml` and run on these.


#### Basic options
A simple real-world invocation might look like this:

```bash
ab12phylo -abi <abi_seqs_folder> \
    -csv <lab_tables_folder> \
    -g <gene> \
    -rf <references.fasta> \
    -bst 1000 \
    -dir <result_folder>
```
where:
* `<abi_seqs_folder>` is the folder containing all input sequences as ABI trace files, but no other files ending in `.ab1`
* `<lab_tables_folder>` contains the `.csv` files mapping user-defined sample IDs to the sequencer's isolate coordinates
* sequencing data is for `<gene>`
* `<references.fasta>` is a FASTA file containing named reference sequences 
* `-bst` or `--bootstrap` sets the number of bootstrap trees to compute to 1000
* `<result_folder>` is the folder where results will be saved 


#### Detailed settings
AB12PHYLO has reasonably smart defaults while allowing fine-grained access to its settings:

```bash
ab12phylo -rf <ref.fasta> \
    -dbpath <blastdb_dir> \
    -sampleset <whitelist> \
    -msa <mafft-clustalo-muscle-tcoffee> \
    -gbl relaxed \
    -bst 1000 \
    -st [32,16]  \
    -s 4 \
    -v -skip
```
* **default:** AB12PHYLO will read in all files that end in `.ab1` or `.csv` in or below the current working directory
* **default:** result files will be in a new subdirectory `./results`
* **default:** samples originate from `ITS1F`
* `<ref.fasta>` is the path of a reference file
* `<blastdb_dir>` is the path to a ready-to-use BLAST+ database
* only samples in `<whitelist>` will be used for this analysis; one ID per line.
* `-msa` (or `--msa_algo`) defines the algorithm to generate the multiple sequence alignment: `mafft`, `clustalo`, `muscle` or `tcoffee`
* `-gbl` (or `--gblocks`) sets the MSA trimming (by `Gblocks`) mode; choices are [`skip`, `relaxed`, `strict`]
* 1000 bootstrap iterations will be run
* `-st` or `--start_trees` sets the number of random (`32`) and parsimony-based (`16`) starting trees for RAxML-NG
* `-s` or `--seed` sets the random [seed](#seed) `4` to allow reproducibility
* `-v` or `--verbose` will cause all logged events to be printed to `stdout`
* `-skip` will skip online BLAST for sequences that were not found in the local BLAST+ database (read [why](#blast-api) this is sensible)


#### Results + Motif Search
Once AB12PHYLO has finished, the package will render a `result.html` and display it in your web browser (unless you pass `--headless`).  

In the `rect` tab showing the computed phylogeny in a rectangular layout, a **Motif search** across node attributes can be run. This is intended for creation of new sample sets for further [subset analyis](#subset-analysis). Selection markup carries over to the `circ` tab.   

If the result directory is moved or sent, motif search will be possible by starting a CGI server in the directory via `python3 -m http.server --cgi 8000` or using `ab12phylo-view`:


#### ab12phylo-visualize + ab12phylo-view
`ab12phylo-visualize` allows re-plotting of trees and rendering a `results.html` from phylogenies and MSAs computed by `ab12phylo`. This can be used to switch the [support value metric](#support-values) or plot an [MSA visualization](#msa-visualization). 
 `ab12phylo-view` shows results of a previous run in a web browser. 
 Both commands expect the path to the AB12PHYLO results or will default to `.`.

```bash
ab12phylo-view <result_folder>
# or
cd <results_folder>
ab12phylo-view
```
Alternatively, append `-visualize` or `-view` to your original `ab12phylo` command.


## Advanced use

#### Config File
AB12PHYLO comes with its own default config file `<ab12phylo_root>/ab12phylo/config/config.yaml` in [YAML](https://yaml.org/) format. It can be adapted to ensure reproducibility and is recommended for remote use. If it is replaced, the path to the adapted config must be passed via `-config`.


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


#### MSA Trimming
MSA trimming by `Gblocks` can be skipped, set to a relaxed or a strict setting via `-gbl`. It's up to you!


#### MSA visualization
An MSA visualization can be plotted to the rectangular tree by passing `-msa_viz`. This will take some rendering time for wider alignments.


#### Support Values
For visualization, you can pick either Felsenstein Bootstrap Proportions `FBP` or [Transfer Bootstrap Expectation](https://doi.org/10.1038/s41586-018-0043-0) `TBE` support values with `-metric`. Newick tree files for both types are generated anyway, so switching support value metric only requires re-running `ab12phylo-visualize`.


#### Log File
If you're having trouble, look at the log file! It will be in your results directory and is named `ab12phylo.log`. Alternatively, you can set the `--verbose` flag and get the same information in real-time to your commandline. Your choice.



