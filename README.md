# AB12PHYLO

![PyPI license](https://img.shields.io/pypi/l/ab12phylo?color=green)
![github version](https://img.shields.io/static/v1?label=version&message=0.4.24-beta&color=brightgreen&style=flat)
![PyPI Python version](https://img.shields.io/pypi/pyversions/ab12phylo)

[AB12PHYLO](https://github.com/lkndl/ab12phylo) is an integrated, easy-to-use pipeline for Maximum Likelihood (ML) phylogenetic tree inference from ABI trace and `FASTA` data. 
At its core, AB12PHYLO runs parallelized instances of [RAxML-NG](https://github.com/amkozlov/raxml-ng) (Kozlov et al. 2019) and a BLAST search in a reference database. 
It enables visual, effortless sample identification based on phylogenetic position and sequence similarity, as well as population subset selection aided by metrics like Tajima's D for estimations of ongoing evolution.
 
 ![./ab12phylo/files/screen_04.png](./ab12phylo/files/screen_04.png)
 
AB12PHYLO was developed to identify plant pathogen populations possibly under balancing selection with *Solanum chilense*, especially in the genus *Alternaria*. With multi-gene phylogenies still a widely-used method in spite of the rise of whole-genome sequencing, future application on fungal phytopathogens or host plants might be possible.


## Installation
 
The recommended way to install AB12PHYLO is via [conda](https://docs.conda.io/), as all [external tools](#external-tools) are available from the [Bioconda](https://anaconda.org/bioconda/repo) channel, but any of the following approaches should work:

#### a) install to an existing python3 conda environment
```shell script
conda activate <env>
conda install -c lkndl ab12phylo
```
where `<env>` is your environment. Please do not install to `base`! Then install a combination of the external tools you would like to use via
```shell script
conda install -c bioconda "blast>=2.9.0" raxml-ng "gblocks=0.91b" mafft clustalo muscle t-coffee
```
Mind that some of these tools are not available for Windows as of 12 March 2021. 

If starting the graphical `ab12phylo` fails with something like `ValueError: Namespace Gtk not available`, `ModuleNotFoundError: No module named 'gi'` or nothing happens at all (on Windows) you are missing PyGObject, the python bindings for GTK3:
```shell script
conda install -c conda-forge pygobject gtk3  
```
If all the icons in the GUI are missing, install some:
```shell script
conda install -c conda-forge adwaita-icon-theme hicolor-icon-theme  
```

If you get an `UnsatisfiableError` in conda because of incompatible packages, please use the next approach:

#### b) install to a newly created python3 conda environment
You could run `conda create -n <env> python=3.x` with `x==6|7|8` and proceed as in **a)** or download [this file](/recipe/gtk-env.yaml) for Linux or [that file](/recipe/win-env.yaml) for Windows, then open a terminal or Anaconda Powershell in your download folder and set up the environment specified inside via
```shell script
conda env create -f gtk-env.yaml
```
and then install AB12PHYLO:
```shell script
conda activate gtk-env
conda install -c lkndl ab12phylo
```

#### c) install via pip from PyPI to your system python
This is possibly the method with the most freedom for users on Ubuntu or any Linux with GTK installed per default, but make sure you are definitely not in a conda environment, i.e. there is no `(<env>)` to the left of your shell prompt like this:
```console
(<env>) foo@bar:~$ conda deactivate
foo@bar:~$ 
```
Install an MSA tool of your choice and possibly BLAST+ by yourself, then AB12PHYLO and all its python3 [dependencies](#dependencies) via
```shell script
pip install ab12phylo
```

If starting `ab12phylo` fails with `No module named 'gi'`, [here](https://pygobject.readthedocs.io/en/latest/getting_started.html) are instructions to install PyGObject. Please be careful, python is probably important for your system!

#### d) build AB12PHYLO from source using git and pip 
Should work on any system and will always install the latest version:
```shell script
git clone https://github.com/lkndl/ab12phylo
cd ab12phylo
pip install --upgrade pip
pip install .
```
Please see **a)** and **c)** if something goes wrong.

#### e) install via pip inside a conda environment
If you are on Linux and would like to install via `pip` inside your conda `<env>`, please check `which python` and `which pip` is active inside the environment (where you will see `(<env>)` to the left of your shell prompt, not shown above).

```console
(<env>) foo@bar:~$ which python
/home/foo/anaconda3/envs/<env>/bin/python

(<env>) foo@bar:~$ which pip
/usr/bin/pip

(<env>) foo@bar:~$ which pip3
/home/foo/anaconda3/envs/<env>/bin/pip3
```

In the case outlined here, `pip` points to a version outside your conda installation, so use `pip3`. If neither points to your conda, re-start your shell and check the environment. Sometimes there is no `pip` at all, which can be fixed on Linux using your package manager (or `sudo apt-get install python3-pip` on Ubuntu), or in conda via `conda install -c conda-forge pip`.

If you see `No module named 'gi'` make sure to fix that with conda as described in **a)**.


## Getting Started

AB12PHYLO has both a graphical and a commandline interface: `ab12phylo` and `ab12phylo-cmd`. The graphical version has its own help, which you will see if you start it with `ab12phylo` from your terminal. 

The following tutorial is intended for the commandline tool. First, take a look at the options by running `ab12phylo-cmd -h`.


#### Test run
The pipeline comes with its own test data set. If you pass `-test`, it will read options from an auxiliary (or backup) config file at `<ab12phylo_root>/ab12phylo_cmd/config/test_config.yaml` and run on these. The test run is set to `--verbose` and will run `--no_remote` BLAST search.


#### Basic options
A simple real-world invocation might look like this:

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
* `<barcode_gene>` was sequenced. more [info](#genes-and-references)
* `<ref.fasta>` contains full GenBank reference records [like this](https://www.ncbi.nlm.nih.gov/nuccore/AF347033.1?report=fasta&log$=seqview&format=text)
* 1000 `-bst` = `--bootstrap` trees will be generated
* `<results>` is where results will be  


#### Detailed settings
`ab12phylo-cmd` has a lot of defaults, but still allows fine-grained access:

```shell script
ab12phylo-cmd -rf <ref.fasta> \
    -db <your_own> \
    -dbpath <your_dir> \
    -abiset <whitelist> \
    -regex_3 <"rx1" "rx2" "rx3"> \
    -algo <mafft-clustalo-muscle-tcoffee> \
    -gbl balanced \
    -local \
    -i \
    -p1

ab12phylo-cmd -p2 \
    -bst 1000 \
    -st [32,16]  \
    -s 4 \
    -v 
```
* **default:** AB12PHYLO will search for `.ab1` and `.csv` files in or below the current working directory
* **default:** use the `./results` subdirectory
* **default:** use the `GTR+Γ` model of evolution
* use `<your_own>` [BLAST+ database](#blast-database); in `<your_dir>`
* only trace files listed in the [`<whitelist>`](#results--motif-search) will be read
* plate number, gene name and well will be parsed from the `.ab1` filename using these three [RegEx](#regex).
* `-algo` will generate the MSA: `mafft`, `clustalo`, `muscle` or `t_coffee`
* `-gbl` sets `Gblocks` MSA trimming mode: `skip`, `relaxed`, `balanced` or `strict`
* `-local` skips online BLAST for sequences not in the local BLAST+ db and [read why](#blast-api)
* `-i` or `--info` shows some more run details in the console
* `-p1` run only part one, up until BLAST; also generates an MView HTML of the MSA to find funny sequences


For the second invocation:
* `-p2` run part two of AB12PHYLO, starting with RAxML-NG
* `-st`: ML tree searches from `32` random and `16` parsimony-based starting trees
* `-s` fixes the random `--seed` to `4` for reproducibility
* `-v` or `--verbose` shows all logged events in the console


## Results + Motif Search
You can import results from `ab12phylo-cmd` into the graphical `ab12phylo`. It's faster, easier and more capable when it comes to viewing and modifying trees. Open a new window/project, press `Ctrl+i` or `Import from commandline version` and select the folder with your commandline results.

Once ML tree inference, bootstrapping and BLAST has finished from the commandline, the pipeline will display a `result.html` in your web browser. This page contains a form that allows **Motif search** across node attributes and calculates diversity metrics for the matching subset/subtree. Entering a space `' '` should match all samples, and entering multiple motifs separated by commas `,` will select all leaves that match at least one motif. 

To select a subtree, enter a motif that matches several leaves or at least two separate specific leaf motifs, with one of them as far left as possible. Entering a single, specific motif for a subtree search can be used to find the index of a node for tree modifications like rooting.

Once you hit `MATCH` or `SUBTREE`, the CGI script in the package will highlight the selection in the tree, as well as compute and display diversity statistics for this 'population'. It will also write and link a file with all selected sample IDs or the paths of the original `.ab1` files in the `/query` subdirectory, which can be used as a whitelist file for a subsequent subset analysis run. Pass it via `-sampleset` if using sample IDs, or via `-abiset` if using file paths. File paths are recommended to reliably exclude outlier versions with identical base ID.

If results are moved or sent, motif search will be possible by using `ab12phylo-view` or starting a CGI server in the directory via `python3 -m http.server --cgi <port>`. Find the port on the intro tab.


## BLAST
If none of the smaller BLAST+ databases are sufficient for your search, you are better of with a [web BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi). Collate your data by running `-p1 --no_BLAST`, then upload the `<gene>/<gene>.fasta` for the gene you wish to use for species annotation. Pass the result via `--BLAST_xml`. You can pass more than one file; AB12PHYLO will use the last annotation for a sample and ignore everything but the first hit.


#### BLAST API
AB12PHYLO is perfectly capable of running online BLAST searches without a BLAST+ installation. However, this is not a suitable main BLAST strategy as BLAST API queries are de-prioritised after just a few attempts. Accordingly, if several runs are attempted on the same data set, passing the `--no_remote` flag will use data from an earlier run by default. You can also leave out species annotation by skipping BLAST entirely with `--no_BLAST`.


#### BLAST+ Database
Sometimes, this pipeline might run headless on a server. To keep it from running head-first in a firewall when it attempts to update its BLAST+ database via FTP, please pre-supply an unzipped, ready-to-use BLAST+ database via `-dbpath` (and `-db` name). Find databases on the [NCBI website](https://ftp.ncbi.nlm.nih.gov/blast/db/), but [web BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) and `--BLAST_xml` might be both easier and faster.


## Visualization

This is easier in the graphical `ab12phylo`. If you'd like to stick to the commandline:

#### ab12phylo-visualize + ab12phylo-view
`ab12phylo-visualize` will re-plot phylogenies and render a new `results.html`. An end user may use this to switch [support values](#support-values) or plot an MSA visualization with `-msa-viz`. This will take some rendering time for wider alignments.
 `ab12phylo-view` shows results of a previous run in a browser, with motif search enabled. Both commands accept a path to the AB12PHYLO results or default to `.`, and are equivalent to appending to the original `ab12phylo-cmd` call.

```shell script
ab12phylo-view <result_folder>
# or
cd <results_folder>
ab12phylo-view
# or 
ab12phylo-cmd -c <my-config.yaml> -bst 1000 (...) -view
```


#### Support Values
For visualization, you can pick either Felsenstein Bootstrap Proportions `FBP` or [Transfer Bootstrap Expectation](https://doi.org/10.1038/s41586-018-0043-0) `TBE` support values with `-metric`. Newick tree files for both types are generated anyway, so switching the support value metric only requires re-running `ab12phylo-visualize`.


#### MSA visualization
AB12PHYLO can plot an additional rectangular tree with an MSA visualization next to it by passing `-msa_viz`. This is single-threaded and takes some extra time for larger MSAs. If you run `-p1`, you can visually inspect the MView HTML of the alignment for outliers already without having to wait for the entire pipeline.


#### Tree modifications
You can `-drop_nodes`, `-replace_nodes` or subtrees with a placeholder, or `-root` the tree with an outgroup sequence. These options accept integer node indices from [motif search](#results--motif-search).

Also see the `VISUALIZATION` section of `--help`.

## Advanced use

#### Remote runs
[tmux](https://askubuntu.com/questions/8653/how-to-keep-processes-running-after-ending-ssh-session/220880#220880) and `--headless` are highly recommended for remote runs, as well as [pre-supplying a BLAST+ db](#blast-database), adapting or replacing the [YAML](https://yaml.org/) config file and setting a fixed seed for reproducibility. 


#### Genes and References
If data for several genes is supplied, AB12PHYLO will restrict the analysis to samples that are present for all genes. If this causes a lot of samples to be dropped, it might be worth leaving out a gene entirely by setting `-g` = `--genes` manually. If you somehow have `No samples shared across all genes`, there might be some trace files that seemingly belongs to another gene.

Passing references is closely inter-linked: If a directory of reference files is supplied via `-rd` = `--ref_dir`, the package will try to match the `.FASTA` files inside to genes *by their filename*. For example, `ITS1F.fasta` will be matched to trace data from the *ITS1F* gene. Alternatively, an ordered list of reference files can be passed via `-rf` = `--ref`, and file names will be ignored. 

If you're feeling this neat and precise and set both the genes and individual references, be careful: The pipeline will then deliberately match references and genes *by order*. Therefore, this will make a mess:

```shell script
-g ITS1F OPA10 -rf ../opa10.fasta ITS1F.phy
```


#### RegEx
If you provide wellsplates mappings, AB12PHYLO will parse plate number, gene name and the sequencer's isolate coordinates from the `.ab1` filename with a RegEx and fetch the user-defined ID from the corresponding `.csv` look-up table. To use your own, please consult `--help` and try out your RegEx [here](https://regex101.com/r/Yulwlf/3) or [there](https://regex101.com/r/Yulwlf/5). From a bash shell, enclose the terms with double quotes `"`. The wellsplate number or ID is also parsed from the `.csv` filename with a RegEx.

Provide a RegEx to `--regex_rev` and AB12PHYLO will look out for reverse reads, and add only their reverse complement to the data set.


#### Log File
If you're having trouble, look at the log file! It will be in your results directory and named like `ab12phylo[|-p1|-p2][-view|-viz]?.log`. Alternatively, you can set the `--verbose` flag and get the same information in real-time to your commandline. Your choice!


## Dependencies
[Biopython](https://biopython.org/wiki/Download), [NumPy](https://numpy.org/), [pandas](https://pandas.pydata.org/docs/getting_started/install.html), [Toytree](https://toytree.readthedocs.io/en/latest/), [Toyplot](https://toyplot.readthedocs.io/en/stable/), [matplotlib](https://matplotlib.org/), [PyYAML](https://pyyaml.org/wiki/PyYAML), [lxml](https://lxml.de/), [xmltramp2](https://pypi.org/project/xmltramp2/), [svgutils](https://github.com/btel/svg_utils), [Pillow](https://pillow.readthedocs.io/en/stable/installation.html), [Requests](https://3.python-requests.org/) and [Jinja2](https://jinja.palletsprojects.com/en/2.11.x/intro/#installation)


## External Tools
* [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) version >=2.9 *(optional, recommended if small or custom [database](#blast-database))*
* [RAxML-NG](https://github.com/amkozlov/raxml-ng/) *(optional, included)*
* an MSA tool: [MAFFT](https://mafft.cbrc.jp/alignment/software/), [Clustal Omega](http://www.clustal.org/omega/), [MUSCLE](https://www.drive5.com/muscle/downloads.htm) or [T-Coffee](http://www.tcoffee.org/Projects/tcoffee/index.html#DOWNLOAD) *(optional, **slow** online clients for [EMBL interface](https://www.ebi.ac.uk/Tools/msa/) included)*
* [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks.html) for MSA trimming *(optional, included)*
