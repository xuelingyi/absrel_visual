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
  * 	w1_percent	- w3_percent: the percentage of codon sites assigned to each rate class. NA if the rate class does not exist,

* **absrel_tree.pdf**: the plot of the absrel-estimated tree with branches colored by the corrected P values. Names of the internal branches having corrected P < 0.2 are also labeled.
  * If hyphy version is 2.5 or higher, a heatmap will also be plotted to show the amino acid alignments of the codon sites that have
     * at least one branch 

* absrel_tested_alignment.pdf



   T3 summary and alignment plot only for the branches with w>=1

  calculates the Evidence Ratio (ER, https://github.com/veg/hyphy/issues/989)
   "log likelihood ratio ratios": ER = 2 log (L[site | selection allowed) - 2 log (L[site | selection not allowed])
   "They have a loose interpretation; site ER that is above 0 for the optimized null setting contributes something to the signal. Typically you want to look at ER of ~2 or higher as suggestive. However, there is no rigorous statistical interpretation here, i.e. ER = 2 means that there is some probability that site is under selection. Higher ERs mean more signal of selection."


        include the maximum ER value of the site on any significant branch

          only keep node labels of the significant or rapidly evolving internal branches





     heatmap color indicating aa: only sites with at least one branch ER>2, and only give the max ER of the significant branches 
     i.e., a site can have high ER on a non-significant branch

    

