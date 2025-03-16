# absrel_visual
This R script can be used to visualize the json file output by hyphy-aBSREL (https://pmc.ncbi.nlm.nih.gov/articles/PMC4408413/), similar to the web tool http://vision.hyphy.org/aBSREL. 

 
# How to run
```
# This script is written in R 4.4.1
# required libraries:
# rjson, ggtree, tidytree, phytools, ggplot2, coRdon, Biostrings

## print the help manual
# Rscript absrel_visual.R -h

Rscript absrel_visual.R <gene_name> <json_file> [--alignment_file=NULL] [--output_dir=NULL] [--heatmap_color=NULL]

```

Required arguments: 
* **gene_name**: The name of the gene/transcript. This will be included in the output files. Needs to be the 1st arguement.
* **json_file**: The json file output by aBSREL. Add path if not in the working directory. Needs to be the 2nd arguement.

Optional arguements:
* **--alignment_file=** the name of the fasta file used to run absrel. Add path if not in the working directory. Defult is not used.
  * Note: only used with hyphy version >= 2.5 for plotting alignments with the tree.
* **--output_dir=** the directory to save output files. Default is the current working directory.
  * Note: do NOT include "/" at the end.
* **--heatmap_color=** The color scheme for the alignment heatmap. Default is the R color scheme. Alternatve is --heatmap_color=taylor. 
  * Note: only used if hyphy version >= 2.5 and the alignment file is provided.
  * The taylor in this plot may be very bright.
* **--help**: print the help message and exit


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

If at least one branch has the corrected p <=0.2, additional PDF file(s) will be generated:
* **absrel_tree.pdf**: It will plot the absrel-estimated tree with branch lengths in Nucleotide GTR and branches colors indicating the corrected P values. Names of the internal branches having corrected P <= 0.2 are also labeled.
  * If hyphy version is 2.5 or higher and the alignment file is provided, this will also include a heatmap showing the amino acid alignments of the codon sites that have ER>2 on at least one branch.
  * The heatmap cells will be labeled by "${the codon site}_ER_${ERmaxValue}_" where the ERmaxValue is the maximum ER of this site among the branches with corrected p <=0.2.
  * Note: a site can have high ER on a non-significant branch, and a significant branch can have sites with ER=1 

* **absrel_tested_alignment.pdf**: This plot is generated only if the hyphy version is at least 2.5. It will show the aligned codon sites of the branches with w>=1 (same branches in absrel_T2.tsv). The sites with no substitutions among these branches are removed. Amino acids are indicated by the taylor color. The sites with ER>2 are highlighted.

    

