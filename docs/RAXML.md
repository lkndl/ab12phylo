# A RAxML Workflow and How to Parallelize It

### Requirements
a Multiple Sequence Alignment, in our case as `.FASTA`

## Step 0: Prep
Preparations are done in `raxml.py`:

```python
def _check_msa(self):
    """
    Checks if the MSA can be understood by raxml-ng, parses it to binary and estimates effective num ...
    """
```

#### MSA check
`--check` flag in the `RAxML-NG` commandline.
Checks for actual errors and duplicate sequences. If there are duplicate sequences, an MSA without them is written: `PREFIX.raxml.reduced.phy` 

#### MSA parsing
`raxml-ng --parse`.
Compresses alignment patterns and creates a binary `.rba` file. 
Apart from dropping duplicate sequences *again*, `--parse` computes memory requirements and estimates the number of CPUs for a run on this MSA.


As we do not want duplicate sequences to be removed, we just use the *un-parsed*, original `msa.fasta` if there were any duplicates.
The pipeline relies on the number of cores per instance recommended by `--parse`, so both steps are necessary.

*****

An actual RAxML run consists of two main parts: *Finding* a likely, *parsimonious* tree topology, and computing statistical *confidence* in this topology.

## Step 1: Tree Inference
Finding the best topology and how much resources to use for it is determined primarily by the `--tree rand{%d},pars{%d}` option of RAxML-NG, translated from the `--start_trees` option of this pipeline.

An example:
```bash
raxml-ng --msa msa.fasta --tree rand{4} --prefix infer_ -seed 12
```

This will prompt 4 Maximum-Likelihood tree searches, each starting from a random topology: A hill-climbing optimisation for each starting tree towards a more parsimonious - more likely - topology.

## Step 2: Bootstrapping
Using re-sampling with replacement based on MSA columns, variants of a starting phylogeny are generated; and support values are computed based on agreement between the inferred best and current bootstrap topology.
The pipeline is intended to be used for 1000 bootstraps:

```bash
raxml-ng --msa msa.fasta --tree tree.nwk --bs-trees autoMRE{1000} --prefix bst_ --seed 16 --threads 4
```
This command would run 1000 bootstraps on the topology supplied via `--tree tree.nwk`. 
Note that `autoMRE{}` activates *bootstopping* if confidence is very high before the 1000 iterations threshold is reached, although this might be unlikely.
Also note that the `--tree` input option is now used with *a file containing one specific topology*.

## Parallelization
Developing this pipeline centers around the availability of enough computational power to run several RAxML instances in parallel, but the benefits of this customisation might not be immediately obvious. 

RAxML-NG has fine-grained parallelization already implemented:
Parallelization across alignment patterns, which can be roughly understood to increase with MSA length, is mirrored in how many threads a `raxml-ng` instance uses. 
On the other hand, very large alignments with thousands of taxa should be run on clusters with [MPI](https://de.wikipedia.org/wiki/Message_Passing_Interface) support.

But for our case of relatively narrow and low MSAs on consumer hardware, *coarse-grained* parallelism makes sense.

*****
### current implementation

As of now, the pipeline runs parallel instances like this (simplified):

```bash
raxml-ng --all --msa msa.fasta --tree rand{200},pars{200} --bs-trees autoMRE{251} --prefix results/raxml/LT91 --seed 14 --threads 3
```

It's not really *wrong*. And it *works*.

But it could be better.

### proposed implementation

#### Distribute Tree Inference

An ML tree search can converge to a local optimum, which is why it is essential to use multiple starting trees. 
These can be distributed among several threads by using distinct random seeds, so 20 ML tree searches could be split as follows:

```bash
raxml-ng --msa msa.fasta --tree rand{5},pars{5} --prefix infer1_ --seed 1

raxml-ng --msa msa.fasta --tree rand{5},pars{5} --prefix infer2_ --seed 2
```

Then we compare the best Log-Likelihood values from each thread, and pick the most parsimonious topology overall.

#### Distribute Bootstrapping

We then distribute the 1000 bootstrap iterations:

```bash
raxml-ng --msa msa.fasta --tree our_best_ML_tree.nwk --bs-trees autoMRE{500} --prefix bootstrap1_ --seed 1

raxml-ng --msa msa.fasta --tree our_best_ML_tree.nwk --bs-trees autoMRE{500} --prefix bootstrap2_ --seed 2
```

#### Re-compute Support Values

Finally, we concatenate all bootstrapping topologies into one file and re-compute support values, both Felsenstein Bootstrap Proportions (FBP) and Transfer Bootstrap Expectation (TBE). This last step is very fast and cannot be parallelized anyway.

```bash
raxml-ng --support --tree our_best_ML_tree.nwk --bs-trees all_bs_trees.nwk --prefix support_ --bs-metric fbp,tbe
```



