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

# Align sequences
aligned_seqs <- readDNAStringSet("06022020alignment.fasta")
aligned_df <- data.frame(ID = gsub("Sample ", "Sample_", aligned_seqs@ranges@NAMES), 
                         seqs = substr(paste(aligned_seqs),11,29683),
                         Region = gsub("\\_.*", "", aligned_seqs@ranges@NAMES), 
                         Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",aligned_seqs@ranges@NAMES),
                         stringsAsFactors = FALSE)
Samples <- aligned_df[ grep("Sample",aligned_df$ID), ] 
aligned_df <- aligned_df[!(aligned_df$ID %in% Samples$ID),]
Samples$ID <- gsub("2019-nCoV_MN908947\\|","",Samples$ID)

#Sample Manipulation (Combing same sequece samples)
uniq_samples <- Samples[ !duplicated(Samples[ , 2] ) , ] #Remove duplicate sequences in samples 23->20
DupSampleList <-list()
for (x in 1:nrow(uniq_samples)) { #Iterate through uniq sequences and group same sequences toegther
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


#-------------------------------------Error Correcting-----------------------------
# Might be OUTDATED.
ErrorDF <- aligned_df[grep("Denmark/ALAB-SSI-1272|Denmark/ALAB-SSI-246|Denmark/ALAB-SSI-660|Denmark/ALAB-SSI-662|England/CAMB-71BDC|England/CAMB-71DE5|England/CAMB-738E2|England/CAMB-73DC5|England/CAMB-74465|England/CAMB-789F9|England/CAMB-79AF5", aligned_df$ID),]
ErrorList <- list()
for (x in 1:nrow(ErrorDF)) {
  ErrorList[[x]] <- gregexpr("N|K",paste(ErrorDF$seqs[x])) 
}
#1:5
substr(paste(aligned_df$seqs[669]),19405,19405) #N->C Europe_46
ErrorDF$seqs[1:5] <- gsub("N","C",paste(ErrorDF$seqs[1:5]))

#6:7
substr(paste(aligned_df$seqs[669]),25123,25123) #K->T Europe_47
ErrorDF$seqs[6:7] <- gsub("K","T",paste(ErrorDF$seqs[6:7]))

#8:11
substr(paste(aligned_df$seqs[1064]),10995,10995) #N->C Europe_63
ErrorDF$seqs[8:11] <- gsub("N","C",paste(ErrorDF$seqs[8:11]))

#Remove and Readd
aligned_df<-aligned_df[!(aligned_df$ID %in% ErrorDF$ID),]
aligned_df<-rbind(aligned_df, ErrorDF)

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

#-------------------------------------Filtering-----------------------------
keepDF <- data.frame(ID=c("r_1","r_68","r_69","r_70","r_46","r_39","r_16","r_72","r_53","r_67","r_57","r_3","r_14")) #Relevant ancestral nodes
samplesIDs <- data.frame(ID=Samples$ID[1:18])
keepDF <- rbind(samplesIDs,keepDF)
Samples <- merge(keepDF, Samples, by.y="ID")

#Add root
rootdf <- aligned_df[grep("Wuhan", aligned_df$ID),]
Samples <- rbind(Samples, rootdf)

#Convert to dna matrix
SampleMt <- data.frame(strsplit(Samples$seqs,""), stringsAsFactors = F)
colnames(SampleMt) <- Samples$ID
SampleMt <- data.frame(t(SampleMt),stringsAsFactors = F)
for (x in 1:ncol(SampleMt)) { #Replace "-" in samples 19,21 with nt copy from Wuhan base 
  if (paste(SampleMt["Sample_19",x]) == "-") {
    SampleMt["Sample_19",x] <- SampleMt["Asia_Wuhan_Hubei (Wuhan-Hu-1/2019)",x]
  }
  if (paste(SampleMt["Sample_21",x]) == "-") {
    SampleMt["Sample_21",x] <- SampleMt["Asia_Wuhan_Hubei (Wuhan-Hu-1/2019)",x]
  }
}

#Remove monomorphic sites
polymorphicloci <-SampleMt[,c(names(Filter(function(x) length(unique(x)) != 1, SampleMt)))]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))
write.csv(polymorphicloci,"polymorphicAlignment.csv")

#Tree prep
polychars<-apply(format(polymorphicloci), 1, paste, collapse="")
dna_strings <- DNAStringSet(polychars)
dna_bin <-as.DNAbin(dna_strings)
phy_object<-phyDat(dna_bin, type = "DNA")

#Trees
dbtree <- optimize_db_parsimony(phy_object) #Distance based tree
ml <- optimize_likelihood(phy_object, dbtree) #Stochastic rearrangement > NNI, -787 -> -779 logLik #Model:GTR+I
bstree <- bootstrap_tree(ml,1000,0.001,"bstree0602") #Bootstrap

#Root
rooted.tree <- root(bstree, which(bstree$tip.label == "Asia_Wuhan_Hubei (Wuhan-Hu-1/2019)"))

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
  
  geom_hilight(node=34, fill=CC[1], extendto=0.0004) +
  geom_hilight(node=54, fill=CC[2], extendto=0.0004) +
  geom_hilight(node=37, fill=CC[3], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  geom_cladelabel(node=34,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[1],align=T,linetype=NA) +
  geom_cladelabel(node=54,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[2],align=T,linetype=NA) +
  geom_cladelabel(node=37,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position='right') 

#Truncate sample 19,21 branches 
fullplot$data[fullplot$data$node %in% c(19,21), "x"] = mean(fullplot$data$x)

#Save
pdf(file="testtest.pdf", width=9, height=7)
fullplot + geom_tiplab(aes(colour=group),size=3, hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()

