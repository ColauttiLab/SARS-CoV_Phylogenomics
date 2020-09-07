library(ape)
library(ggtree)
#library(Biostrings)
#library(dplyr)
#library(ape)
#library(ggtree)
#library(ggplot2)
#1library(rBLAST) 
#library(genbankr) 


BaseTree<-read.tree("./intermediatedata/BaseTree.tree")

#Tree building
#------- 

#Regions
regions<-gsub("_.*","",BaseTree$tip.label)
regions<-gsub("Sample","Test Sample",regions)
regions<-gsub("r","Ancestral",regions)
regions<-gsub("Asia","Wuhan",regions)
regionGroups<-split(BaseTree$tip.label,regions)
WtDTcol<-groupOTU(BaseTree,regionGroups)

WtDTcol<-groupOTU(BaseTree,regionGroups)
WtDTcol$tip.label<-gsub("r_","",WtDTcol$tip.label)
WtDTcol$tip.label<-gsub("Sample","P",WtDTcol$tip.label)
WtDTcol$tip.label<-gsub(".*\\((.*)\\)","\\1",WtDTcol$tip.label)

CC<-c("#fc8d62","#8da0cb","#66c2a5")
Pal <-c("#377eb8","#e78ac3","#e41a1c")

nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
nodetree

fullplot<-ggtree(WtDTcol,layout='rectangular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_hilight(node=58, fill=CC[1], extendto=0.0004) +
  geom_hilight(node=45, fill=CC[2], extendto=0.0004) +
  geom_hilight(node=40, fill=CC[3], extendto=0.0004) +
  geom_hilight(node=10, fill=CC[2], extendto=0.0004) +
  geom_hilight(node=154, fill=CC[1], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
    geom_tiplab2(aes(colour=group),size=1, hjust=0, offset=0.0001, align=T, linetype=3) +#
    geom_text(aes(label=node),size=1) +
  
  geom_cladelabel(node=58,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[1],align=T,linetype=NA) +
  geom_cladelabel(node=45,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[2],align=T,linetype=NA) +
  geom_cladelabel(node=40,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.0001,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position='right') 

pdf(file="outputs/prunedtree.pdf", width=10, height=60)
  fullplot + drop.tip
    geom_tiplab2(aes(colour=group),size=3, hjust=0, offset=0.0001, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2 
dev.off()
