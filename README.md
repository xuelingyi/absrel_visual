# absrel_visual
This R script can be used to visualize the json file output by hyphy-aBSREL (https://pmc.ncbi.nlm.nih.gov/articles/PMC4408413/), similar to the web tool http://vision.hyphy.org/aBSREL. 

# How to run
```
# This script is written in R 4.4.1

Rscript absrel_visual.R ${my.path} ${my.gene} NULL

# my.path is the path to the folder where the json file is located (NOTE: no "/" at the end). The outputs will also be put in this directory.
# ${my.gene} is the name of the gene/transcript, and the json file should be named as ${my.gene}.absrel.json 
# NULL means that the alignment heatmap (only for hyphy version >= 2.5, see details below) will use default R colors
# alternatively, replace "NULL" with "taylor" will use taylor colors, but note this may be bright
```
# More details
This script will parse the json file output by aBSREL and generate the following files in the provided directory:
* **absrel_T1.tsv**: a table similar to the branch-site table (or Table1) from hyphy vision. Columns are:
  * branches: all branches in the tree
  * transcript: the name of the gene/transcript 
  * tested: whether the branch is the tested foreground or not
  * P_corrected: the corrected P value
  * sites_nonsynonymous: the number of codons with nonsynonymous mutations on this branch (only reported if the hyphy version is 2.5 or higher)
  * sites_ER2: the number of codons with Evidence Ratio >2 on this branch (see below for ER definition; only reported if the hyphy version is 2.5 or higher)
  * rate_class: the number of rate classes aBSREL output for this branch
  * LRT: the likelihood ratio test (LRT) aBSREL output for this branch
  * w1 - w3: the omega (w) value of rate classes 1 - 3 for this branch. NA if the rate class does not exist. An error will be reported if the branch has more than three rate classes.
  * w1_percent	- w3_percent: the percentage of codon sites assigned to each rate class. NA if the rate class does not exist,

* **absrel_T2.tsv**: This table is generated only if the hyphy version is at least 2.5. This table is similar to the Table3 from hyphy vision. It only includes the branches that have omega >1. Columns are:
  * branches: all branches in the tree
  * transcript: the name of the gene/transcript
  * site: the codon site  
  * parent: the codon at this site of the parent branch
  * parent_aa: the protein at this site of the parent branch
  * codon: the codon at this site of this branch
  * codon_aa: the protein at this site of this branch
  * ER: the Evidence Ratio of this site on this branch. Calculated as ER = exp(log (L[site | selection allowed) - log (L[site | selection not allowed])), following https://github.com/veg/hyphy/issues/989 "Higher ERs mean more signal of selection."
  * substitution: the type of substitution from the parent branch to this branch, including no substitution, nonsynonymous, synonymous, or deletion.

* **absrel_tree.pdf**: This plot is generated only if at least one branch has the corrected p <=0.2. It will plot the absrel-estimated tree with branch lengths estimated by Nucleotide GTR and branches colors by the corrected P values. Names of the internal branches having corrected P <= 0.2 are also labeled.
  * If hyphy version is 2.5 or higher, this will also include a heatmap showing the amino acid alignments of the codon sites that have ER>2 on at least one branch.
  * The heatmap cells will be labeled by "the codon site"_"ER"_"ERmaxValue"_ where the ERmaxValue is the maximum ER of this site among the branches with corrected p <=0.2.
  * Note: a site can have high ER on a non-significant branch, and a significant branch can have sites with ER=1 

* **absrel_tested_alignment.pdf**: This plot is generated only if the hyphy version is at least 2.5. This plot will show the aligned codon sites of the branches with w>=1 (same branches in absrel_T2.tsv). The sites with no substitutions among these branches are removed. Colors indicate different amino acids in the taylor color. The sites with ER>2 are highlighted.

    

