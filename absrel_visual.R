#!/usr/bin/env Rscript

#module unload R
#module load R/4.4.1

library(rjson)
library(ggtree)
library(tidytree)
library(phytools)
library(ggplot2)
library(coRdon)
library(Biostrings)

args <- commandArgs(TRUE)
path=args[1]
## this is the path to the dir where the json file is stored (no "/" at the end); output will be in the same folder
t=args[2]
## this is the name of the transcript/gene; the json file is named as ${t}.absrel.json
heatmap.color=args[3] ## NULL if R default, taylor if taylor colors (this may strike your eyes)


get.parent.branch = function(my.branch, all.branches, tree, ...){
  return(all.branches[tree$edge[tree$edge[,2] %in% which(all.branches == my.branch), 1]])
}
color_map_taylor <- c("A" = "#CCFF00", "V" = "#99FF00", "I" = "#66FF00", "L" = "#33FF00", "M" = "#00FF00", "F" = "#00FF66",
                      "Y" = "#00FFCC", "W" = "#00CCFF", "H" = "#0066FF", "R" = "#0000FF", "K" = "#6600FF", "N" = "#CC00FF",
                      "Q" = "#FF00CC", "E" = "#FF0066", "D" = "#FF0000", "S" = "#FF3300", "T" = "#FF6600", "G" = "#FF9900",
                      "P" = "#FFCC00", "C" = "#FFFF00")


print("get absrel json file")
data <- fromJSON(file = paste0(path, "/", t, ".absrel.json"))

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

## the absrel-generated tree
tree = as.phylo(data$input$trees)
all.branches = c(tree$tip.label, tree$node.label)

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

print("save Table1")
write.table(T1, paste0(path, "/absrel_T1.tsv"), sep="\t", quote = F, row.names = F)

print("save Table 3 and plot codon alignments")
# only branches with w>1 at sites with substitutions
if(as.numeric(data$analysis$version) >= 2.5){
  write.table(T2, paste0(path, "/absrel_T2.tsv"), sep="\t", quote = F, row.names = F)
  
  # only summarize and plot data if at least one branch with corrected p <=0.2
  if(min(as.numeric(unlist(data$`branch attributes`$'0')[grep("Corrected P-value", names(unlist(data$`branch attributes`$'0')))])) <= 0.2){
    ## plot alignments of the tested (w>1) branches for each transcript that have significant branches
    T2$plot = sapply(T2$site, FUN=function(x){
      if(all(T2[T2$site ==x, "substitution"] == "no")){
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
    
    pdf(paste0(path, "/absrel_tested_alignment.pdf"), width = per*0.53, height = 0.3*length(unique(T2$branches)))
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


print("plot the absrel-generated tree")
# with branches lengths estimated by Nucleotide GTR and colored by corrected p values if at least one significant branch
T1$P_corrected = as.numeric(T1$P_corrected)

if(min(T1$P_corrected) <= 0.2){
  ## plot the tree
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
  
  ## add alignment heatmap for higher hyphy versions
  if(as.numeric(data$analysis$version) >= 2.5){
    # get hyphy input alignments 
    ali = readBStringSet(paste0(path, "/", t, ".noref.fa"))
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

    pdf(paste0(path, "/absrel_tree.pdf"), width = (7+ncol(heatmap2)*0.3))
    if(heatmap.color=="taylor"){
      print(gheatmap(p0, heatmap2, offset=0.002*ncol(heatmap2), 
                     width=0.05*ncol(heatmap2), font.size=1.8, color="black",
                     colnames_angle=90, colnames_position = "top", hjust = 0) + 
              vexpand(0.1, 1) + 
              scale_fill_manual(values=color_taylor) +
              theme(legend.key.size = unit(0.8, "line"),
                    legend.text = element_text(size=6), legend.title = element_text(size=8),
                    legend.position = "inside", legend.position.inside = c(0.1, 0.6),
                    legend.background = element_blank()) + 
              labs(title=paste0(t, "\n", length(tree$tip.label), " tips, ", length(ali[[1]]), "bp, ", length(ali[[1]])/3, " sites")))
    } else {
      print(gheatmap(p0, heatmap2, offset=0.002*ncol(heatmap2), 
                     width=0.05*ncol(heatmap2), font.size=1.8, color="black",
                     colnames_angle=90, colnames_position = "top", hjust = 0) + 
              vexpand(0.1, 1) + 
              theme(legend.key.size = unit(0.8, "line"),
                    legend.text = element_text(size=6), legend.title = element_text(size=8),
                    legend.position = "inside", legend.position.inside = c(0.1, 0.6),
                    legend.background = element_blank()) +
              labs(title=paste0(t, "\n", length(tree$tip.label), " tips, ", length(ali[[1]]), "bp, ", length(ali[[1]])/3, " sites")))
    }
    
  } else {
    ## lower hyphy versions, only print the tree
    pdf(paste0(path, "/absrel_tree.pdf"), width = 7)
    print(p0 + labs(title=paste0(t, "\n", length(tree$tip.label), " tips")))
  }
  
  dev.off()
}


