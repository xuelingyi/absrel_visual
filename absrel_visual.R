#!/usr/bin/env Rscript

#module unload R
#module load R/4.4.1

suppressPackageStartupMessages(library(rjson, quietly = T))
suppressPackageStartupMessages(library(ggtree, quietly = T))
suppressPackageStartupMessages(library(tidytree, quietly = T))
suppressPackageStartupMessages(library(phytools, quietly = T))
suppressPackageStartupMessages(library(ggplot2, quietly = T))
suppressPackageStartupMessages(library(coRdon, quietly = T))
suppressPackageStartupMessages(library(Biostrings, quietly = T))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the user asked for help
if ("--help" %in% args || "-h" %in% args) {
  cat("Description: \nThis script visualizes the json file output by aBSREL.\n")
  cat("\n")
  cat("Usage: \nRscript absrel_visual.R <gene_name> <json_file> [--alignment_file=NULL] [--output_dir=.] [--heatmap_color=NULL]\n")
  cat("\n")
  cat("Arguements:\n")
  cat("  gene_name               (required) The name of the gene/transcript. This will be included in the output files.\n")
  cat("  json_file               (required) The json file output by aBSREL.\n")
  cat("\n")
  cat("  --alignment_file        (optional) The alignment file used to run aBSREL. If provided and hyphy version is >=2.5, the tree will be plotted with an alignment heatmap. Defult is not used.\n")
  cat("  --output_dir            (optional) The output path (no / at the end). Default is the current working directory. \n")
  cat("  --prefix                (optional) The prefix added to output file names. Default is the gene name. \n")
  cat("  --heatmap_color         (optional) The color scheme for the alignment heatmap. Default (NA) is the R default color scheme. Alternatve is taylor colors. --heatmap_color=taylor \n")
  cat("  --plot_nosignificance   (optional) generate plots even if no branch has p<=0.2. Note that the heatmap will not give ER values if there is no significant branch. Default is False. Pass True to turn this on.\n")
  cat("  --help                  (optional) Show this help message and exit\n")
  cat("\n")
    
  quit(status = 0)
}

# Process arguments
if (length(args) < 2) {
  cat("Error: Argument missing. Use --help for usage information.\n")
  quit(status = 1)
} else {
  t=args[1] 
  print(paste0("gene name ", t))
  json=args[2]
  print(paste0("json file ", json)) 
}

alignment_file=NULL
output_dir=getwd()
heatmap_color=NULL
plot_nosignificance=F
prefix=t
tree_only=F
if(length(args) > 2){
## parse provided optional arguements
  for (arg in args[3:length(args)]) {
    if (startsWith(arg, "--alignment_file")) {
      alignment_file=sub("--alignment_file=", "", arg)
      print(paste0("alignment file provided: ", alignment_file))
    }
    if (startsWith(arg, "--output_dir")) {
      output_dir=sub("--output_dir=", "", arg)
    }
    if (startsWith(arg, "--heatmap_color")) {
      heatmap_color=sub("--heatmap_color=", "", arg)
    }
    if (startsWith(arg, "--plot_nosignificance")) {
      plot_nosignificance=as.logical(sub("--plot_nosignificance=", "", arg))
    }
    if (startsWith(arg, "--prefix")) {
      prefix=sub("--prefix=", "", arg)
    }
    if (startsWith(arg, "--tree_only")) {
      plot_nosignificance=as.logical(sub("--tree_only=", "", arg))
    }
  }
}

get.parent.branch = function(my.branch, all.branches, tree, ...){
  return(all.branches[tree$edge[tree$edge[,2] %in% which(all.branches == my.branch), 1]])
}
color_map_taylor <- c("A" = "#CCFF00", "V" = "#99FF00", "I" = "#66FF00", "L" = "#33FF00", "M" = "#00FF00", "F" = "#00FF66",
                      "Y" = "#00FFCC", "W" = "#00CCFF", "H" = "#0066FF", "R" = "#0000FF", "K" = "#6600FF", "N" = "#CC00FF",
                      "Q" = "#FF00CC", "E" = "#FF0066", "D" = "#FF0000", "S" = "#FF3300", "T" = "#FF6600", "G" = "#FF9900",
                      "P" = "#FFCC00", "C" = "#FFFF00")

print("get absrel json file")
data <- fromJSON(file = json)

## the absrel-generated tree
tree = as.phylo(data$input$trees)
all.branches = c(tree$tip.label, tree$node.label)

## summarize omega and p per branch (absrel website Table 1)
T1 = as.data.frame(names(data$`branch attributes`$'0'))
names(T1) = "branches"
T1$transcript = t
T1$tested = sapply(T1$branches, FUN=function(x){
  if(x %in% names(data$tested$'0')){ return("yes")} else {
    return("no")}
})
T1$P_corrected = NA
T1$rate_class = NA #The number of w rate classes inferred for each branch (up to 3)
T1$LRT = NA
T1$w1 = NA #The estimated w value of the first w rate class for each branch
T1$w1_percent = NA # The proportion of sites evolving under the first w rate class
T1$w2 = NA
T1$w2_percent = NA
T1$w3 = NA
T1$w3_percent = NA

## prep table 3 only for higher hyphy versions
print(paste0("hyphy version ", data$analysis$version))
if(as.numeric(data$analysis$version) < 2.5){
  print("Note: T2 and alignemnts are not plotted for this hyphy version!")
  
} else {
  
  T1$sites_nonsynonymous = NA # the number of sites with nonsynonymous substitutions (different from the T1 in hyphy vision; not sure what they counted)
  T1$sites_ER2 = NA # the number of sites with ER > 2 (different from the T1 in hyphy vision; not sure what they counted)
  
  ## summarize branch-site codon and evidence ratio (absrel website Table 3)
  T2 = as.data.frame(rep(names(data$`Site Log Likelihood`$tested), each=length(data$substitutions$'0')))
  names(T2) = "branches"
  T2$transcript = t
  T2$site = 1:length(data$substitutions$'0')
  T2$parent=NA
  T2$parent_aa=NA
  T2$codon=NA
  T2$codon_aa=NA
  T2$ER = NA # substitutions	Evidence ratio
  T2$substitution = "no" # synonymous substitutions
}

print("loop through branches")
for(i in 1:nrow(T1)){
  my.branch=T1[i, "branches"]
  
  T1[i, "P_corrected"] = data$`branch attributes`$'0'[[my.branch]]$`Corrected P-value`
  T1[i, "LRT"] = data$`branch attributes`$'0'[[my.branch]]$LRT
  T1[i, "rate_class"] = data$`branch attributes`$'0'[[my.branch]]$`Rate classes`
  omega=data$`branch attributes`$'0'[[my.branch]]$`Rate Distributions`
  for(j in 1:T1[i, "rate_class"]){
    if(j > 3){
      print(paste0("more than 3 rate classes in branch ", my.branch, " in transcript ", t))
    } else {
      T1[i, paste0("w", j)] = format(round(omega[[j]][1], digits = 6), scientific = FALSE)
      T1[i, paste0("w", j, "_percent")] = format(round(omega[[j]][2], digits = 6), scientific = FALSE) 
    }
  }
  
  if(as.numeric(data$analysis$version) >= 2.5){
    ## T2 summary and alignment plot only for the branches with w>=1
    if(my.branch %in% T2$branches){
       T2[T2$branches == my.branch, "ER"] = exp(data$`Site Log Likelihood`$unconstrained[[1]] - data$`Site Log Likelihood`$tested[[my.branch]][[1]])
      
      # loop for each site
      for(s in 1:length(data$substitutions$'0')){
        ## site names start with 0
        
        ## find the codon of the parent branch; root is labeled as "0"
        parent.branch = get.parent.branch(my.branch, all.branches, tree)
        T2[T2$branches == my.branch & T2$site == s, "parent"] = data$substitutions$'0'[[as.character(s-1)]][["root"]]
        while(parent.branch != "0"){
          
          ## nodes in this list have substitutions at the site; other nodes at the site share the codon with their parents
          if(parent.branch %in% names(data$substitutions$'0'[[as.character(s-1)]])){
            T2[T2$branches == my.branch & T2$site == s, "parent"] = data$substitutions$'0'[[as.character(s-1)]][[parent.branch]]
            break
          } else {
            ## get the parent branch of the current parent branch and iterate -- until it shows up in the list or goes to the root
            parent.branch = get.parent.branch(parent.branch, all.branches, tree)
          }
        }
        
        ## find the codon of this branch
        if(my.branch %in% names(data$substitutions$'0'[[as.character(s-1)]])){
          T2[T2$branches == my.branch & T2$site == s, "codon"] = data$substitutions$'0'[[as.character(s-1)]][[my.branch]]
        } else {
          T2[T2$branches == my.branch & T2$site == s, "codon"] = T2[T2$branches == my.branch & T2$site == s, "parent"]
        }
        
        ## if deletion
        if(!all(unlist(strsplit(T2[T2$branches == my.branch & T2$site == s, "codon"], split="")) %in% c("A", "a", "T", "t", "G", "g", "C", "c"))){
          T2[T2$branches == my.branch & T2$site == s, "substitution"] = "delete"
        } else {
          ## test if parent to branch codon is a non-synonymous substitution
          T2[T2$branches == my.branch & T2$site == s, "codon_aa"] = Biostrings::GENETIC_CODE[T2[T2$branches == my.branch & T2$site == s, "codon"]]
          T2[T2$branches == my.branch & T2$site == s, "parent_aa"] = Biostrings::GENETIC_CODE[T2[T2$branches == my.branch & T2$site == s, "parent"]]
          if(T2[T2$branches == my.branch & T2$site == s, "codon_aa"] != T2[T2$branches == my.branch & T2$site == s, "parent_aa"]){
            T2[T2$branches == my.branch & T2$site == s, "substitution"] = "nonsynonymous"
          } else {
            ### synonymous if the nucleotides differ (otherwise no substitution)
            if(T2[T2$branches == my.branch & T2$site == s, "codon"] != T2[T2$branches == my.branch & T2$site == s, "parent"]){
              T2[T2$branches == my.branch & T2$site == s, "substitution"] = "synonymous"
            }
          }
        }
      }
    }
    
    T1[i, "sites_nonsynonymous"] = nrow(T2[T2$branches == my.branch & T2$substitution == "nonsynonymous",])
    T1[i, "sites_ER2"] = nrow(T2[T2$branches == my.branch & T2$ER > 2,])
  }
}

print(paste0("save outputs in ", output_dir))
print("save Table 1")
write.table(T1, paste0(output_dir, "/", prefix, "_absrel_T1.tsv"), sep="\t", quote = F, row.names = F)

# only branches with w>1 at sites with substitutions
if(as.numeric(data$analysis$version) >= 2.5){
  print("save Table 2 and plot codon alignments")
  write.table(T2, paste0(output_dir, "/", prefix, "_absrel_T2.tsv"), sep="\t", quote = F, row.names = F)
  
  # only summarize and plot data if at least one branch with corrected p <=0.2; or if the user asked for it
  p.all = unlist(data$`branch attributes`$'0')[grep("Corrected P-value", names(unlist(data$`branch attributes`$'0')))]
  if(min(as.numeric(p.all)) <= 0.2 | plot_nosignificance){
    ## plot alignments of the tested (w>1) branches for each transcript that have significant branches
    T2$plot = sapply(T2$site, FUN=function(x){
      if(all(T2[T2$site ==x, "substitution"] %in% c("no", "delete"))){
        ## no branch has any substitution, do not plot this site
        return("no")
      } else {
        return("yes")
      }
    })
    T2 = T2[T2$plot =="yes",]
    T2$site_with_substitution = factor(T2$site)
    T2$ER_group = sapply(T2$ER, FUN=function(x){
      if(x < 2){
        return("2")
      } else {
        if(x<10){
          return("10")
        } else {
          if(x<50){
            return("50")
          } else {
            return("large")
          }
        }
      }
    })
    T2$empty = "      "
    per=25
    
    pdf(paste0(output_dir, "/", prefix, "_absrel_tested_alignment.pdf"), width = per*0.53, height = 0.3*length(unique(T2$branches)))
    for(seq in 1:ceiling(length(unique(T2$site))/per)){
      start = sort(unique(T2$site))[per*(seq-1) +1]
      end = sort(unique(T2$site))[min(per*seq, length(unique(T2$site)))]
      
      print(ggplot(T2[T2$site >= start & T2$site <= end, ], aes(x = site_with_substitution, y = branches)) +
              geom_label(aes(label = empty, color = ER_group), size = 4, label.size=2.2) +
              scale_color_manual(values=c("2"="white", "10"="grey", "50"="grey40", "large"="black"),
                                 breaks=c("2", "10", "50", "large"),
                                 labels=c("<2", "[2,10)", "[10, 50)", ">=50")) + 
              geom_label(aes(label = codon, fill = codon_aa), color="black", size = 3, label.size=0, alpha=0.7) +
              scale_fill_manual(values=color_map_taylor, guide="none") +
              theme_classic())
    }
    dev.off()
  }
}


# with branches lengths estimated by Nucleotide GTR and colored by corrected p values if at least one significant branch
T1$P_corrected = as.numeric(T1$P_corrected)

if(min(T1$P_corrected) <= 0.2 | plot_nosignificance){
  print("plot the absrel-generated tree")
 
  branches.length = as.data.frame(tree$edge)
  colnames(branches.length)=c("parent", "node")
  branches.length$length = NA
  branches.color = as.data.frame(tree$edge)
  colnames(branches.color)=c("parent", "node")
  branches.color$p_corrected = NA
  
  ## loop through branches
  for(my.branch in names(data$`branch attributes`$'0')){
    ## annotate this branch in the tree
    branches.length[branches.length$node == which(all.branches == my.branch), "length"] = data$`branch attributes`$'0'[[my.branch]]$`Nucleotide GTR`
    branches.color[branches.color$node == which(all.branches == my.branch), "p_corrected"] = data$`branch attributes`$'0'[[my.branch]]$`Corrected P-value`
    
    ## only keep node labels of the significant or rapidly evolving internal branches
    if(!(my.branch %in% tree$tip.label) & T1[T1$branches == my.branch, "P_corrected"] > 0.2){
      tree$node.label[tree$node.label == my.branch] = ""
    }
  }
  
  ## plot the absrel-estimated tree 
  tree$edge.length = branches.length$length
  tree$edge.length[is.na(tree$edge.length)] = 0
  branches.color$p_corrected_group = unlist(sapply(branches.color$p_corrected, FUN=function(x){
    if(!is.na(x)){
      if(x <= 0.0001){return("***")}
      if(x > 0.0001 & x <= 0.001 ){return("**")}
      if(x > 0.001 & x <= 0.05 ){return("*")}
      if(x > 0.05 & x <= 0.2 ){return("rapid_evolving")}
      if(x> 0.2){return("not significant")}
    } else {
      return(NA)
    }
  }))

  p0 = ggtree(tree) %<+% branches.color + aes(colour=p_corrected_group) + 
    geom_tiplab(size=2, align=T) + hexpand(0.5,1) +
    geom_nodelab(size=2, hjust = 1, vjust=-0.8) +
    scale_color_manual(values=c("***"="red4", 
                                "**"="#f03b20",
                                "*"="orange", 
                                "rapid_evolving"="steelblue", 
                                "not significant"="grey40"), 
                       breaks=c("***", "**", "*", "rapid_evolving", "not significant"), 
                       labels=c("<=0.0001", "(0.0001,0.001]", "(0.001,0.05]", "(0.05, 0.2]", ">0.2"))
  
  ## add alignment heatmap if hyphy version is at least 2.5 or if the user asked for this, and the alignment is provided
  if((as.numeric(data$analysis$version) >= 2.5 | plot_nosignificance) & !is.null(alignment_file) & !(tree_only)){
    print("plot alignments with the tree")
    # get hyphy input alignments 
    ali = readBStringSet(alignment_file)
    ali.ER2 = data.frame(label = names(ali))
    
    # get codon sites with ER>2
    for(s in sort(unique(T2[T2$ER>2, "site"]))){
      
      # extract dna sequences of these sites and get the aa sequence 
      ali.ER2$site = sapply(ali.ER2$label, FUN=function(x){as.character(ali[[x]][(3*(s-1)+1):(3*s)])})
      ali.ER2$aa = sapply(ali.ER2$site, FUN=function(x){
        if(x == "NA"){
          return(x)
        } else {return(Biostrings::GENETIC_CODE[x])}
      })
      
      ## include the maximum ER value of the site on any significant branch
      names(ali.ER2)[ncol(ali.ER2)] = paste0(s, "_ER_", round(max(T2[T2$site ==s & T2$branches %in% T1[T1$P_corrected<=0.2, "branches"], "ER"])), "_")
    }
    
    ## heatmap color indicating aa: only sites with at least one branch ER>2, and only give the max ER of the significant branches 
    ## i.e., a site can have high ER on a non-significant branch
    heatmap2 = ali.ER2[, grep("_ER", names(ali.ER2))]
    rownames(heatmap2) = ali.ER2$label

    pdf(paste0(output_dir, "/", prefix, "_absrel_tree.pdf"), width = (7+ncol(heatmap2)*0.3))
    if(is.null(heatmap_color)) {
      print(gheatmap(p0, heatmap2, offset=0.0018*ncol(heatmap2), 
                     width=0.05*ncol(heatmap2), font.size=1.8, color="black",
                     colnames_angle=90, colnames_position = "top", hjust = 0) + 
              vexpand(0.1, 1) + 
              theme(legend.key.size = unit(0.8, "line"),
                    legend.text = element_text(size=6), legend.title = element_text(size=8),
                    legend.position = "inside", legend.position.inside = c(0.1, 0.6),
                    legend.background = element_blank()) +
              labs(title=paste0(t, "\n", length(tree$tip.label), " tips, ", length(ali[[1]]), "bp, ", length(ali[[1]])/3, " sites")))
    } else if(heatmap_color=="taylor"){
      print(gheatmap(p0, heatmap2, offset=0.0018*ncol(heatmap2), 
                     width=0.05*ncol(heatmap2), font.size=1.8, color="black",
                     colnames_angle=90, colnames_position = "top", hjust = 0) + 
              vexpand(0.1, 1) + 
              scale_fill_manual(values=color_taylor) +
              theme(legend.key.size = unit(0.8, "line"),
                    legend.text = element_text(size=6), legend.title = element_text(size=8),
                    legend.position = "inside", legend.position.inside = c(0.1, 0.6),
                    legend.background = element_blank()) + 
              labs(title=paste0(t, "\n", length(tree$tip.label), " tips, ", length(ali[[1]]), "bp, ", length(ali[[1]])/3, " sites")))
    } 
    
    
  } else {
    print("hyphy version < 2.5 or no alignment file provided or asked to plot tree only; print the tree without alignments")
    pdf(paste0(output_dir, "/", prefix, "_absrel_tree.pdf"), width = 7)
    print(p0 + labs(title=paste0(t, "\n", length(tree$tip.label), " tips")))
  }
  
  dev.off()
}

