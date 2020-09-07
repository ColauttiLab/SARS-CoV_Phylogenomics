library(phangorn)
library(Biostrings)
library(ips)
library(magrittr)
library(dplyr)
library(ggtree)

#Functions

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

# Trim alignment
alignment <- read.dna("intermediatedata/BRalignment3bp.afa", "fasta")
#checkAlignment(test, plot = FALSE)
trimmed_alignment <- trimEnds(alignment, min.n.seq = 100) # Variable to play with (see documentation). This affects # of objs in uniq_samples by trimming
# So for min.n.seq = 100, the function will cut the sequences at the position where 100 sequences have unambigiuos bases, and it will fill sequences that fall short with N's
# So a greater min.n.seq = less filling in with N's. 
# Notice how changing min.n.seq will change max(aligned_seqs@ranges@width)
# if min.n.seq = the number of sequences, then the function will cut the at the spot where there is a single unambigious base (too strict).
#checkAlignment(test, plot = FALSE)
aligned_seqs <- trimmed_alignment %>% as.list %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

aligned_df <- data.frame(ID = gsub("Sample ", "Sample_", aligned_seqs@ranges@NAMES), 
                         # Want to get rid of substr... The trimEnds function should remove boarder differences so that substr not needed
                         # For example, without substr uniq_samples(line79) = 24objs, with uniq_samples = 23 objs. Difference due to cutting out first 11 chars
                         seqs = substr(paste(aligned_seqs),11,29879),
                         Region = gsub("\\_.*", "", aligned_seqs@ranges@NAMES), 
                         Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",aligned_seqs@ranges@NAMES),
                         stringsAsFactors = FALSE)
Samples <- aligned_df[ grep("Sample",aligned_df$ID), ] 
aligned_df <- aligned_df[!(aligned_df$ID %in% Samples$ID),]

#Sample Manipulation (Combining samples with identical sequences)
uniq_samples <- Samples[ !duplicated(Samples[ , 2] ) , ] #Remove duplicate sample sequences (Confirm with real # of duplicates)
DupSampleList <-list()
for (x in 1:nrow(uniq_samples)) { #Iterate through unique sequences and group identical sequences together
  dupSeq <- paste(uniq_samples$seqs[x])
  matches <- Samples[ which(paste(Samples$seqs)==dupSeq), ]
  if (nrow(matches)>1) {
    matches$ID[2:nrow(matches)] <- gsub("Sample_", ", ", matches$ID[2:nrow(matches)])
    name = paste(matches$ID, collapse="")
    tempDF <- data.frame(ID=name, seqs=paste(matches$seqs[1]), Region=name,Country=name)
    DupSampleList[[x]] <- tempDF
  } else {
    DupSampleList[[x]] <- matches
  }
}
DupSampleList<-DupSampleList[lengths(DupSampleList) != 0]
Samples <- dplyr::bind_rows(DupSampleList)

#-------------------------------------Sorting into reference samples-----------------------------
uniq_seqs <- aligned_df[ !duplicated(aligned_df[ , 2] ) , ] #Remove duplicate sequences
AncestorList <- list()
for (x in 1:nrow(uniq_seqs)) { 
  dupSeq <- paste(uniq_seqs$seqs[x])
  matches <- aligned_df[ which(paste(aligned_df$seqs)==dupSeq), ]
  uniq_matches <- matches[ !duplicated(matches[ , 4] ) , ]
  datalist <- list()
  for (y in 1:nrow(uniq_matches)) {
    country=paste(uniq_matches$Country[y])
    datalist[[y]] <- matches[ which(matches$Country==country),]
  }
  if (nrow(matches) >= 2) {
    AncestorList[[x]] <- datalist
  }
}
AncestorList<-AncestorList[lengths(AncestorList) != 0]

#Combining with samples
RefDF <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("ID", "Region", "Country", "N"))
for (x in 1:length(AncestorList)) {
  Data <- AncestorList[[x]]
  name <- paste0("r_",x)
  tempdf <- data.frame(ID=name, seqs=paste(Data[[1]][1,]$seqs), Region=name, Country=name)
  Samples <- rbind(Samples, tempdf)
  for (y in 1:length(Data)){
    RegionDF <- Data[[y]]
    lng <- nrow(RegionDF)
    regiondf <- data.frame(ID=name, Region=paste(RegionDF$Region[1]), Country=paste(RegionDF$Country[1]), N=lng)
    RefDF <- rbind(RefDF, regiondf)
  }
}

#Add root
rootdf <- aligned_df[grep("Wuhan.*2020", aligned_df$ID),]
Samples <- rbind(Samples, rootdf)

#Convert to dna matrix
SampleMt <- data.frame(strsplit(Samples$seqs,""), stringsAsFactors = F)
colnames(SampleMt) <- Samples$ID
SampleMt <- data.frame(t(SampleMt),stringsAsFactors = F)
for (x in 1:ncol(SampleMt)) { #Replace "-" in samples 19,21 with nt copy from Wuhan base 
  if (paste(SampleMt["Sample_19",x]) == "N" || paste(SampleMt["Sample_19",x]) == "-") {
    SampleMt["Sample_19",x] <- SampleMt[161,x] # Change this number based on position (row number) of root sequence
  }
  if (paste(SampleMt["Sample_21",x]) == "N" || paste(SampleMt["Sample_21",x]) == "-") {
    SampleMt["Sample_21",x] <- SampleMt[161,x]
  }
}

#Remove monomorphic sites
polymorphicloci <-SampleMt[,c(names(Filter(function(x) length(unique(x)) != 1, SampleMt)))]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))
write.csv(polymorphicloci,"./intermediatedata/polymorphicAlignment.csv")

#Might want to trim ends here
#View(polymorphicloci)
#example: polymorphicloci <- polymorphicloci[,6:(ncol(polymorphicloci)-15)]

#Tree prep
polychars<-apply(format(polymorphicloci), 1, paste, collapse="")
dna_strings <- DNAStringSet(polychars)
dna_bin <-as.DNAbin(dna_strings)
phy_object<-phyDat(dna_bin, type = "DNA")

#Trees
dbtree <- optimize_db_parsimony(phy_object) #Distance based tree
ml <- optimize_likelihood(phy_object, dbtree) #Stochastic rearrangement > NNI, -787 -> -779 logLik #Model:GTR+I
bstree <- bootstrap_tree(ml,1000,0.001,"./intermediatedata/BSTree") #Bootstrap

#Root
rooted.tree <- root(bstree, which(bstree$tip.label == "Wuhan/WH04/2020"))

#Save Tree
write.tree(rooted.tree,"./intermediatedata/BaseTree.tree",tree.names=T)
#point david left off

#-------------------------------------Tree Design-----------------------------
#Regions
regions<-gsub("_.*","",bstree$tip.label)
regions<-gsub("Sample","Test Sample",regions)
regions<-gsub("r","Ancestral",regions)
regions<-gsub("Asia","Wuhan",regions)
regionGroups<-split(bstree$tip.label,regions)
WtDTcol<-groupOTU(rooted.tree,regionGroups)

##Tip and Node label cleanup
WtDTcol$tip.label<-gsub("r_","",WtDTcol$tip.label)
WtDTcol$tip.label<-gsub("Sample","P",WtDTcol$tip.label)
WtDTcol$tip.label<-gsub(".*\\((.*)\\)","\\1",WtDTcol$tip.label)
#Manual p filtering
WtDTcol$node.label<-gsub("98.2","",WtDTcol$node.label) #Upon converting from bstree to WtDTCol, NA->98.2
WtDTcol$node.label<-gsub("^[^9876].*","",WtDTcol$node.label) #Remove p scores < 60

#Color pallettes
CC<-c("#fc8d62","#8da0cb","#66c2a5")
Pal <-c("#377eb8","#e41a1c","#e78ac3")

nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
nodetree$data[nodetree$data$node %in% c(19,21), "x"] = mean(nodetree$data$x)

#ggtree
fullplot<-ggtree(WtDTcol,layout='rectangular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_nodelab(hjust=-0.2,size=2) +
  geom_treescale(x = 0.003, y = 30, fontsize = 2, linesize = 0.1) +
  
  #geom_hilight(node=34, fill=CC[1], extendto=0.0004) +
  #geom_hilight(node=54, fill=CC[2], extendto=0.0004) +
  #geom_hilight(node=37, fill=CC[3], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  geom_cladelabel(node=34,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[1],align=T,linetype=NA) +
  geom_cladelabel(node=54,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[2],align=T,linetype=NA) +
  geom_cladelabel(node=37,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position='right') 

#Save
pdf(file="outputs/figures/Fig2.pdf", width=12, height=12)
  fullplot + geom_tiplab(aes(colour=group),size=3, hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()

