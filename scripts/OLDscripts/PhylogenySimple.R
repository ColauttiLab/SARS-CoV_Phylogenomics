library(phangorn)
library(Biostrings)
library(ips)
library(magrittr)
library(dplyr)
library(ggtree)

## ----------------- USER-DEFINED: ----------------
# USER-Determined non-informative tips from PhylogenyFull.R
CutTips<-c(91,60,84,114,64,113,93,92,83,94,12,88,27,26,15,85,2,75,103,97,102,157,160,105,
           99,96,150,123,110,106,108,112,111,152,109,82,104,125,124,107,7,78,55,80,19,47,
           81,90,63,72,86,76,79,62,145,101,73,87,100,69,24,48,89,77,43,30,38,52,50,46,51,
           70,25,22,35,18,57,53,58,41,59,121,37,74,68,36,32,95,33,23,49,34,98,6,54,21,65,
           115,56,42,28,29,44,39,40,134,131,31,8,66,137,158,156,161,155,130,142,133,122,
           140,132,138,14,147,135,153,141,136,148,149,143,139,118,126,116,119,120,151,154,
           117,146,127,128,129,144,159,67)

# Alignment file
polymorphicloci<-read.csv("./intermediatedata/polymorphicAlignmentFull.csv", stringsAsFactors = F,row.names=1)
## -------------------------------------------------

optimize_db_parsimony<-function(phy_object) {
  dna_dist <- dist.ml(phy_object)
  treeNJ <- NJ(dna_dist)
  treeUPGMA <- upgma(dna_dist)
  if (parsimony(treeNJ,phy_object) < parsimony(treeUPGMA,phy_object)) {
    fit = pml(treeNJ,phy_object)
  } else {
    fit = pml(treeUPGMA,phy_object)
  }
  return(fit)
} #optimize distance based methods NJ and UPGMA

optimize_likelihood<-function(phy_object, pml_object) {
  model_test <- modelTest(phy_object)
  model_test <- model_test %>% arrange(logLik)
  params <- strsplit(paste(model_test[nrow(model_test),]$Model),"\\+")
  model <- params[[1]][1]
  Gamma <- FALSE
  Inv <- FALSE
  if (length(params[[1]]) == 2) {
    if (params[[1]][2] == "I") {
      Inv = TRUE
    } else {
      Gamma = TRUE
    }
  } 
  if (length(params[[1]]) == 3) {
    Gamma <- TRUE
    Inv <- TRUE
  }
  print("----------------------------------------------------------------------")
  fit <- optim.pml(pml_object, model = model, optInv = Inv, optGamma = Gamma, rearrangement = "stochastic")
  print("----------------------------------------------------------------------")
  print(model_test[nrow(model_test),])
  print(paste0("Gamma = ",Gamma))
  print(paste0("Inv = ", Inv))
  print("----------------------------------------------------------------------")
  print(paste0("Unoptimized loglikelihood: ",pml_object$logLik))
  print("----------------------------------------------------------------------")
  print(paste0("Optimized loglikelihood: ",fit$logLik))
  print("----------------------------------------------------------------------")
  return(fit)
}  #Maximum likeliehood

bootstrap_tree<-function(fitted_model,bs_iterations,scale_bar,out){
  bs <- bootstrap.pml(fitted_model, bs=bs_iterations, optNni=TRUE, multicore=F, control = pml.control(trace=0))
  pdf(out,width=13,height=10)
  bstree<-plotBS(midpoint(fitted_model$tree),bs, p = 70, type="p")
  add.scale.bar(length = scale_bar, cex = 0.9, font = 2)
  dev.off()
  return(bstree)
} #Bootstrap

# Trim uninformative tips
CutTips<-paste0("r_",CutTips)
polyReduced<-polymorphicloci[!row.names(polymorphicloci) %in% CutTips , ]

#Tree prep
polychars<-apply(format(polyReduced), 1, paste, collapse="")
dna_strings <- DNAStringSet(polychars)
dna_bin <-as.DNAbin(dna_strings)
phy_object<-phyDat(dna_bin, type = "DNA")

#Trees
dbtree <- optimize_db_parsimony(phy_object) #Distance based tree
ml <- optimize_likelihood(phy_object, dbtree) #Stochastic rearrangement > NNI, -787 -> -779 logLik #Model:GTR+I
bstree <- bootstrap_tree(ml,100,0.001,"./intermediatedata/PreTrimBSTree") #Bootstrap

#Root
rooted.tree <- root(bstree, which(bstree$tip.label == "Asia_China_Hubei (Wuhan/WH04/2020)"))

#Save Tree
write.tree(rooted.tree,"./intermediatedata/SimpleTree.tree",tree.names=T)
#point david left off

#-------------------------------------Tree Design-----------------------------
#Regions
regions<-gsub("r_.*","Reference",rooted.tree$tip.label)
regions<-gsub("S.*","Sample",regions)
regions<-gsub("Asia.*","Wuhan",regions)
regionGroups<-split(bstree$tip.label,regions)
WtDTcol<-groupOTU(rooted.tree,regionGroups)

##Tip and Node label cleanup
WtDTcol$tip.label<-gsub("_","",WtDTcol$tip.label)
WtDTcol$tip.label<-gsub(".*\\((.*)\\)","\\1",WtDTcol$tip.label)
#Manual p filtering
#WtDTcol$node.label<-gsub("98.2","",WtDTcol$node.label) #Upon converting from bstree to WtDTCol, NA->98.2
#WtDTcol$node.label<-gsub("^[^9876].*","",WtDTcol$node.label) #Remove p scores < 60

#Color pallettes
CC<-c("#fc8d62","#0096c7","#66c2a5","#00b4d8")
Pal <-c("#377eb8","#e41a1c","#e78ac3")

nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
#nodetree$data[nodetree$data$node %in% c(19,21), "x"] = mean(nodetree$data$x)
pdf(file="./outputs/SimpleTreeNodes.pdf", width=12, height=12)
  nodetree
dev.off()

#ggtree
fullplot<-ggtree(WtDTcol,layout='circular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_nodelab(hjust=-0.2,size=2) +
#  geom_treescale(x = 0.003, y = 30, fontsize = 2, linesize = 0.1) +
  
  geom_hilight(node=69, fill=CC[1], alpha=0.3, extendto=0.001) +
  geom_hilight(node=78, fill=CC[2], alpha=0.3, extendto=0.001) +
  geom_hilight(node=126, fill=CC[3], alpha=0.3, extendto=0.001) +
  geom_hilight(node=120, fill=CC[4], alpha=0.3, extendto=0.001) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  geom_cladelabel(node=69,label="A.1",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[1],align=T) +
  geom_cladelabel(node=78,label="B.1",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[2],align=T) +
  geom_cladelabel(node=126,label="B.1.5",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[3],align=T) +
  geom_cladelabel(node=120,label="B.1.1",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[4],align=T) +

  theme_tree(legend.position='right')

#Save
pdf(file="./outputs/figures/Fig2.pdf", width=12, height=12)
  fullplot + geom_tiplab2(aes(colour=group,size=I(2)), hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()


