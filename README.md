# absrel_visual
visualize the json file output by hyphy-aBSREL 

path: the path to the folder where the json file is located
gene: the name of the gene/transcript, and the json file should be named as ${t}.absrel.json 
heatmap_color: the color for the alignment heatmap (only for hyphy version >= 2.5). 
heatmap.color=NULL --> default R colors
heatmap.color=taylor --> taylor colors -- may be bright
Rscript absrel_visual.R ${path} ${gene} ${heatmap_color}


# this script is written in R 4.4.1

path=args[1]
## this is the path to the dir where the json file is stored (no "/" at the end); output will be in the same folder
t=args[2]
## this is the name of the gene; 


  ## T3 summary and alignment plot only for the branches with w>=1

  ## calculates the Evidence Ratio (ER, https://github.com/veg/hyphy/issues/989)
      ## "log likelihood ratio ratios": ER = 2 log (L[site | selection allowed) - 2 log (L[site | selection not allowed])
      ## "They have a loose interpretation; site ER that is above 0 for the optimized null setting contributes something to the signal. Typically you want to look at ER of ~2 or higher as suggestive. However, there is no rigorous statistical interpretation here, i.e. ER = 2 means that there is some probability that site is under selection. Higher ERs mean more signal of selection."


       ## include the maximum ER value of the site on any significant branch

         ## only keep node labels of the significant or rapidly evolving internal branches




         path=args[1]
## this is the path to the dir where the json file is stored (no "/" at the end); output will be in the same folder
t=args[2]
## this is the name of the transcript/gene; the json file is named as ${t}.absrel.json
color=args[3] ## NULL if R default, taylor if  colors (this may strike your eyes)




    ## heatmap color indicating aa: only sites with at least one branch ER>2, and only give the max ER of the significant branches 
    ## i.e., a site can have high ER on a non-significant branch

    

