library(Biostrings)
library(phangorn)
library(dplyr)
library(ggtree)

#------------------------- User-defined parameters --------------------------
Iters<-10#000 # Bootstrap iterations (set low (10 or 100) for quick tree; 10000 for final published tree)
Scale<-0.001 # Scale bar to use on phylogeny
#------------------------------------------------------------------------------

# Import polymorphism data
dna_strings<-readDNAStringSet("./intermediatedata/NextFilt.afa")
# Import PANGOLIN lineages
PangoLins<-read.csv("./inputdata/PangoLins.csv")

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

#Tree prep
dna_bin <-as.DNAbin(dna_strings)
phy_object<-phyDat(dna_bin, type = "DNA")

#Trees
dbtree <- optimize_db_parsimony(phy_object) #Distance based tree
ml <- optimize_likelihood(phy_object, dbtree) #Stochastic rearrangement > NNI, -787 -> -779 logLik #Model:GTR+I
bstree <- bootstrap_tree(ml,Iters,Scale,"./intermediatedata/BSTree") #Bootstrap

#Root
rooted.tree <- root(bstree, which(bstree$tip.label == "Wuhan/WH04/2020"))

#Save Tree
write.tree(rooted.tree,"./intermediatedata/BaseTree.tree",tree.names=T)

# Prep tree output
# Define regions (for color-coding)
regions<-gsub("Wuhan.*","Wuhan",bstree$tip.label)
regions<-gsub("^S[0-9]+[im]$","Sample",regions)
regions<-gsub(".*[0-9]+","Reference",regions)
# Define PANGOLIN lineages
regions[gsub(".*\\.([0-9]+)$","\\1",bstree$tip.label) %in% gsub("EPI_ISL_","",PangoLins$GISAID.ID)]<-"PANGOLIN"
regionGroups<-split(bstree$tip.label,regions)
WtDTcol<-groupOTU(rooted.tree,regionGroups)

#Color palettes
CC <- c("#06BCC1","#07004D","#A44A3F","#C1AAC0")
Pal <- c("#06BCC1","#07004D","#A44A3F","#C1AAC0")

#nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
#nodetree$data[nodetree$data$node %in% c(19,21), "x"] = mean(nodetree$data$x)
#nodetree

#ggtree
fullplot<-ggtree(WtDTcol,layout='circular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_nodelab(hjust=-0.2,size=2) +
  geom_treescale(x = 0.003, y = 30, fontsize = 2, linesize = 0.1) +
  
  #geom_hilight(node=34, fill=CC[1], extendto=0.0004) +
  #geom_hilight(node=54, fill=CC[2], extendto=0.0004) +
  #geom_hilight(node=37, fill=CC[3], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  #geom_cladelabel(node=34,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[1],align=T,linetype=NA) +
  #geom_cladelabel(node=54,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[2],align=T,linetype=NA) +
  #geom_cladelabel(node=37,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position='right') 
fullplot
#Save
pdf(file="outputs/figures/Fig2.pdf", width=12, height=12)
 fullplot + geom_tiplab(aes(colour=group),size=3, hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()

