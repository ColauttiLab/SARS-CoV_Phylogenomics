library(phangorn)
library(Biostrings)
library(ips)
library(magrittr)
library(dplyr)
library(ggtree)

polymorphicloci<-read.csv("./intermediatedata/polymorphicAlignmentFull.csv", stringsAsFactors = F,row.names=1)

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
polychars<-apply(format(polymorphicloci), 1, paste, collapse="")
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
write.tree(rooted.tree,"./intermediatedata/PreTrimBaseTree.tree",tree.names=T)
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
CC<-c("#fc8d62","#8da0cb","#66c2a5")
Pal <-c("#377eb8","#e41a1c","#e78ac3")

nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
#nodetree$data[nodetree$data$node %in% c(19,21), "x"] = mean(nodetree$data$x)
nodetree

#ggtree
fullplot<-ggtree(WtDTcol,layout='rectangular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_nodelab(hjust=-0.2,size=2) +
#  geom_treescale(x = 0.003, y = 30, fontsize = 2, linesize = 0.1) +
  
  #geom_hilight(node=34, fill=CC[1], extendto=0.0004) +
  #geom_hilight(node=54, fill=CC[2], extendto=0.0004) +
  #geom_hilight(node=37, fill=CC[3], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  #geom_cladelabel(node=34,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[1],align=T,linetype=NA) +
  #geom_cladelabel(node=54,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[2],align=T,linetype=NA) +
  #geom_cladelabel(node=37,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position='right') 

#Save
pdf(file="./outputs/PreTrimBaseTree.pdf", width=6, height=24)
  fullplot + geom_tiplab(aes(colour=group,size=I(2)), hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()


