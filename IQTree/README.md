# Tutorial for Phylogenomics Workshop using IQ-TREE

If you haven’t installed IQ-TREE, please download and install the binary for your platform. For the next steps, the folder containing your iqtree executable should be added to your PATH enviroment variable so that IQ-TREE can be invoked by simply entering iqtree at the command-line. Alternatively, you can also copy iqtree binary into your system search.

Run the command

`iqtree`

should display something like this to the screen:

`IQ-TREE multicore version 1.7-beta6 for Mac OS X 64-bit built Nov  6 2018
Developed by Bui Quang Minh, Nguyen Lam Tung, Olga Chernomor,
Heiko Schmidt, Dominik Schrempf, Michael Woodhams.

Command-line examples (replace 'iqtree ...' by actual path to executable):`

## Input
We will use a Turtle data set to demonstrate the use of IQ-TREE throughout the workshop. Please download the following input files:

The DNA alignment (in FASTA format), which is a subset of the original Turtle data set used to assess the phylogenetic position of Turtle relative to Crocodile and Bird (Chiari et al., 2012)

    QUESTIONS:
    View the alignment in Jalview.
    Can you identify the gene boundary from the viewer? Does it roughly match the partition file?
    Is there missing data? Which taxa seem to have most missing data?
    Do you think this missing data can be problematic?


You can now start to reconstruct a maximum-likelihood (ML) tree for the Turtle data set (assuming that you are in the same folder where the alignment is stored):

`iqtree -s turtle.fa -bb 1000 -nt AUTO`

Options explained:
   
   -s turtle.fa to specify the input alignment as turtle.fa.
   -bb 1000 to perform the ultrafast bootstrap with 1000 replicates.
   -nt AUTO to determine the best number of CPU cores to speed up the analysis.

This simple command will perform three important steps:
Model selection (with ModelFinder) to select the best-fit model to the data.
Reconstruct the ML tree with the selected model.
Assess branch supports with the ultrafast bootstrap.
Once the run is done, IQ-TREE will write several output files including:

turtle.fa.iqtree: the main report file that is self-readable. You should look at this file to see the computational results. It also contains a textual representation of the final tree.

turtle.fa.treefile: the ML tree in NEWICK format, which can be visualized in FigTree.

turtle.fa.log: log file of the entire run (also printed on the screen).

turtle.fa.ckp.gz: checkpoint file used to resume an interrupted analysis.

And a few other files.

    QUESTIONS:
    Look at the report file turtle.fa.iqtree.
    What is the best-fit model? What do you know about this model?
    Visualise the tree turtle.fa.treefile in FigTree.
    Compare the tree with the published tree (Chiari et al., 2012). Are they the same or different?
    If different, where are the difference(s)?
    Look at the boostrap supports. Which branch(es) have a low support?

## Applying partition model

We now perform a partition model analysis, where one allows each partition to have its own model:

`iqtree -s turtle.fa -spp turtle.nex -bb 1000 -nt AUTO`

Options explained:

-spp turtle.nex to specify the input partition file. It also initiates an edge-linked proportional partition model (i.e. each partition has its own evolutionary rate).
QUESTIONS:
Look at the report file turtle.nex.iqtree. What are the slowest- and highest-evolving genes?
Compare the BIC score of partition model versus un-partition model done above. Which model is better?
Visualise the tree turtle.nex.treefile in Figtree and compare it with the tree from the un-partitioned model. Are they the same or different? If different, where is the difference? Which tree agrees with Chiari et al., 2012?
Look at the boostrap supports. Which branch(es) have a low support?
Further readings:

Partition models: Chernomor et al. (2016).

4) Choosing the best partitioning scheme
We now perform the PartitionFinder algorithm implemented in IQ-TREE that merges partitions to reduce potential model overfitting:

iqtree -s turtle.fa -spp turtle.nex -bb 1000 -nt AUTO -m MFP+MERGE -rcluster 10 -pre turtle.merge
Options explained:

-m MFP+MERGE to perform PartitionFinder followed by tree reconstruction.
-rcluster 10 to reduce computations by only examining the top 10% partitioning schemes (relaxed clustering algorithm).
-pre turtle.merge to set the prefix for all output files as turtle.merge.*. This is to avoid overwriting outputs from the previous analysis.
QUESTIONS:
Look at the report file turtle.merge.iqtree. How many partitions do we have now?
What is the BIC score? Is it better or worse than those of the un-partition and partition models done previously?
How does the tree look like now? How high/low are the bootstrap supports?
Further readings:

PartitionFinder: Lanfear et al., 2012.
Relaxed clustering algorithm: Lanfear et al., 2014.


5) Applying a mixture model - GHOST
We now perform a mixture model analysis, to account for heterotachy both within and between genes:

iqtree -s turtle.fa -m GTR+FO+H4 -bb 1000 -nt AUTO -pre turtle.GHOST4

Options explained:

-m GTR+FO+H4 to specify the model of sequence evolution. GTR specifies the General Time reversible model, +FO specifies that base frequencies are to be inferred by maximum likelihood (as opposed to empirical values from the alignment), +H4 means we use the GHOST model with 4 classes, with GTR model parameters linked across the 4 classes.
QUESTIONS:
Visualise the tree turtle.GHOST4.treefile in Figtree and compare it with the trees from the earlier analyses. Are they the same or different? If different, where is the difference? Which tree agrees with Chiari et al., 2012?
Look at the boostrap supports. Which branch(es) have a low support?
Further readings:

GHOST model: Crotty et al. (2019).

6) Tree topology tests
We now want to know whether the trees inferred for the Turtle data set have significantly different log-likelihoods or not. This can be conducted with Shimodaira-Hasegawa test (Shimodaira and Hasegawa, 1999), or expected likelihood weights (Strimmer and Rambaut, 2002).

First, concatenate the trees constructed by single and partition models into one file:

For Linux/MacOS:

cat turtle.fa.treefile turtle.nex.treefile >turtle.trees
For Windows:

type turtle.fa.treefile turtle.nex.treefile >turtle.trees
Now pass this file into IQ-TREE via -z option:

iqtree -s turtle.fa -spp turtle.nex.best_scheme.nex -z turtle.trees -zb 1000 -n 0 -wpl -pre turtle.test
Options explained:

-spp turtle.nex.best_scheme.nex to provide the partition model found previously to avoid model selection again.
-z turtle.trees to provide a set of trees.
-zb 1000 to specify 1000 boostrap replicates for tree topology tests.
-n 0 to avoid tree search and just perform tree topology tests.
-wpl to print partition-wise log likelihoods for both trees. This will be used in the next section.
-pre turtle.test to set the prefix for all output files as turtle.test.*.
QUESTIONS:
Look at the report file turtle.test.iqtree. There is a new section called USER TREES.
Do the two trees have significantly different log-likelihoods?
HINTS:

The KH and SH tests return p-values, thus a tree is rejected if its p-value < 0.05 (marked with a - sign).
bp-RELL and c-ELW return posterior weights which are not p-value. The weights sum up to 1 across the trees tested.
7) Concordance factors
So far we have assumed that gene trees and species tree are equal. However, it is well known that different genomic loci might lead to different trees. Therefore, we now want to quantify the agreement between gene trees and species tree, so-called concordance factors.

This feature is only available in the beta version 1.7-beta6, please download it from https://github.com/Cibiv/IQ-TREE/releases/tag/v1.7-beta6. Command iqtree in this section needs to point to the installed beta version.
You first need to compute the gene trees, one for each partition separately:

iqtree -s turtle.fa -S turtle.nex -pre turtle.loci -nt 2
Options explained:

-S turtle.nex to specify the partition file and tell IQ-TREE to infer separate trees for every partition.
All output files are similar to a partition analysis, except that the tree turtle.loci.treefile now contains a set of gene trees.

The gene concordance factor (gCF) and site concordance factor (sCF) for a reference tree are defined as follows. For every branch of a reference tree, gCF is defined as the percentage of “decisive” gene trees containing that branch and sCF is the percentage of decisive alignment sites supporting that branch.

You can now compute gCF and sCF for the tree inferred under partition model

iqtree -t turtle.nex.treefile --gcf turtle.loci.treefile -s turtle.fa --scf 100
Options explained:

-t turtle.nex.treefile to specify a species tree.
--gcf turtle.loci.treefile to specify a gene-trees file.
--scf 100 to draw 100 random quartets when computing sCF.
Once finished this run will write several files:

turtle.nex.treefile.cf.tree: tree file where branches are annotated with bootstrap/gCF/sCF values.
turtle.nex.treefile.cf.stat: a table file with various statistics for every branch of the tree.
Similarly, you can compute gCF and sCF for the tree under unpartitioned model:

iqtree -t turtle.fa.treefile --gcf turtle.loci.treefile -s turtle.fa --scf 100
QUESTIONS:
Visualise turtle.nex.treefile.cf.treein FigTree.
How do gCF and sCF values look compared with bootstrap supports?
Visualise turtle.fa.treefile.cf.tree. How do these values look like now on the contradicting branch?
FINAL QUESTION: In your opinion, which tree do you think is the true tree, given all analyses done in this Tutorial?
8) Identifying most influential genes (optional)
Now we want to investigate the cause for such topological difference between trees inferred by single and partition model. One way is to identify genes contributing most phylogenetic signal towards one tree but not the other.

How can one do this? Well, we can look at the gene-wise log-likelihood differences between the two trees, called T1 and T2. Those genes having the largest lnL(T1)-lnL(T2) will be in favor of T1. Whereas genes showing the largest lnL(T2)-lnL(T1) are favoring T2.

With the -wpl option done above, IQ-TREE will write partition-wise log-likelihoods into turtle.test.partlh file.

QUESTIONS:
Import this file into MS Excel. Compute the partition wise log-likelihood differences between two trees.
What are the two genes that most favor the tree inferred by single model?
Have a look at the paper by (Brown and Thomson, 2016). Compare the two genes you found with those from this paper. What is special about these two genes?
