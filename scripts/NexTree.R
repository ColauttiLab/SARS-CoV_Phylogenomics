library(Biostrings)
library(phangorn)
library(dplyr)
library(ggtree)

#------------------------- User-defined parameters --------------------------
Iters<-1#0000 # Bootstrap iterations (set low (10 or 100) for quick tree; 10000 for final published tree)
Scale<-0.001 # Scale bar to use on phylogeny
RemoveR<-NA # Refere
Comb<-10 # Used to shorten reference sequence table; Regions with < Comb are combined into 'Other' category
Exclude<-c("B.6","B.7","B.2","B.5","B.4","B.3","A.4","A.5","A.6","A.3","A.2","B.1.6","B.1.22","B.1.19") # Node names to exclude (run with NA to see full phylogeny)
#------------------------------------------------------------------------------

# Import polymorphism data
dnaIn<-readDNAStringSet("./intermediatedata/NextFilt.afa")
# Import PANGOLIN lineages
#PangoLins<-read.csv("./inputdata/PangoLins.csv")

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

## Index PANGOLIN and patient sequences
ref_seqs <- grep("^[A-z ]+\\.[0-9]{6}$",names(dnaIn)) # Reference Sequences
pang_seqs <- grep("^[AB]$|^[AB]\\.[0-9].*$",names(dnaIn)) # Pangolin Sequences
pat_seqs <- grep("S[0-9]{2}[im]",names(dnaIn)) # Patient Sequences
Wuhan <- grep("Wuhan",names(dnaIn))

#keep_seq<-grep("[A-z]+.[0-9]{6}$",names(dna_strings),invert=T) # Patient samples, PANGOLIN + Wuhan refs
#uniq_keep<-keep_seq[ !duplicated(keep_seq) ]
#ref_seq<-dna_strings[grep("[A-z]+.[0-9]{6}$",names(dna_strings))] # All other reference genomes

Dist<-dist.dna(as.DNAbin(dnaIn),model="N",pairwise.deletion=T)
DistMat<-as.matrix(Dist)

# Find unique patient samples with unique distance matrices
uniq_pat_seqs<-intersect(grep(FALSE,duplicated(DistMat[,ref_seqs])),pat_seqs)

# Find closest sequence(s) for each patient sample
keep_seqs<-DNAStringSet(NULL)
SeqTable<-data.frame(NULL)
#SeqTable<-data.frame(ID=NA,Region=NA,N=NA,Variants=NA)
r<-1
for(i in 1:length(uniq_pat_seqs)){
  Idx <- uniq_pat_seqs[i] # Index linking patient sample to DistMat
  
  # use ref_seqs as columns to compare only reference samples
  # use Idx to cycle through rows of unique patient sequences (based on distance matrices, above)
  BestSeq <- DistMat[Idx,] == min(DistMat[Idx,ref_seqs]) # Keep most similar sequences
  BestSeq <- dnaIn[BestSeq] # Convert to sequences
  BestSeq <- BestSeq[grep("^[A-z ]+\\.[0-9]{6}$",names(BestSeq))] # Remove non-reference sequences


  Complete <- rowSums(letterFrequency(BestSeq,c("A","G","C","T")))/rowSums(alphabetFrequency(BestSeq)) # Proportion of non-ambiguous bases

  # Tally Locations
  SeqTemp<-table(gsub("^([A-z ]+)\\.[0-9]+","\\1",names(BestSeq)))
  SeqTemp<-data.frame(ID=paste0("r",r),
                      Region=names(SeqTemp),
                      N=as.numeric(SeqTemp))
  SeqTable<-rbind(SeqTable,SeqTemp)
  
  # Check if sequence already exists in keep_seqs; if not, add new reference ID
  if(!BestSeq[Complete == max(Complete)][1] %in% keep_seqs){
    # Save sequence
    BestSeq<-BestSeq[Complete == max(Complete)][1]
    keep_seqs<-DNAStringSet(c(keep_seqs,BestSeq)) # Add best representative sequence ((lowest proportion of ambiguous bases))
    keep_seqs@ranges@NAMES[r]<-paste0("r",r) # Combine/simplify reference name
    r<-r+1
  } 
  Idx <- BestSeq <- Complete <- Refs <- SeqTemp <- dupr <- NA #Cleanup
}

# Combine duplicate reference/region counts
SeqTable<-SeqTable %>%
  group_by(ID,Region) %>%
  summarize(N=sum(N))

# Simplify table and sequences (remove duplicates)
#SeqTable<-SeqTable[ !duplicated(SeqTable$Refs) , ]
#keep_seqs<-keep_seqs[keep_seqs@ranges@NAMES %in% SeqTable$ID]
#for(Row in 1:nrow(SeqTable)){
#  SeqTable$ID[Row]<-paste0("r",Row)
#  keep_seqs@ranges@NAMES[Row]<-paste0("r",Row)
#}

# Combine Regions with only a few refs (see user-defined variable 'Comb', above)
SeqTable_short<-SeqTable[SeqTable$N < Comb,] %>%
  arrange(Region) %>%
  group_by(ID) %>%
  summarize(N=sum(N)) %>%
  mutate(Region="Other")

# Replace rare regions with 'Other'
SeqTable_short<- rbind(SeqTable[SeqTable$N >= Comb,],
                      SeqTable_short) %>%
  arrange(ID)



# ---------------- Reference region summary table -----------------
# Output to table
write.csv(SeqTable_short,"./outputs/RefSeqs.csv")
write.csv(SeqTable,"./outputs/RefSeqs_long.csv")
# ---------------------------------



# Phylogeny Prep
## Sequence setup
dna_bin <- as.DNAbin(c(dnaIn[c(Wuhan,pat_seqs,pang_seqs)],keep_seqs)) # Combine sequences
if(length(Exclude) >0 ){
  dna_bin<-dna_bin[!names(dna_bin) %in% Exclude]
}
phy_object<-phyDat(dna_bin, type = "DNA")

# Trees
dbtree <- optimize_db_parsimony(phy_object) #Distance based tree
ml <- optimize_likelihood(phy_object, dbtree) #Stochastic rearrangement > NNI, -787 -> -779 logLik #Model:GTR+I
bstree <- bootstrap_tree(ml,Iters,Scale,"./intermediatedata/BSTree") #Bootstrap

#Root
rooted.tree <- root(bstree, which(bstree$tip.label == "Wuhan/WH04/2020"))

# ---------------- Tree file output -----------------
#Save Tree
write.tree(rooted.tree,"./intermediatedata/BaseTree.tree",tree.names=T)
# ---------------------------------

# Prep tree output
# Define regions (for color-coding)
#regions<-gsub("Wuhan.*","Wuhan",bstree$tip.label)
regions<-gsub("^S[0-9]+[im]$","Sample",bstree$tip.label)
regions<-gsub("r[0-9]+","Reference",regions)
# Define PANGOLIN lineages
regions<-gsub("[AB]|^[AB].[0-9].*","PANGOLIN",regions)
#regions[gsub(".*\\.([0-9]+)$","\\1",bstree$tip.label) %in% gsub("EPI_ISL_","",PangoLins$GISAID.ID)]<-"PANGOLIN"
regionGroups<-split(bstree$tip.label,regions)
WtDTcol<-groupOTU(rooted.tree,regionGroups)

#Color palettes
CC <- c("#06BCC1","#07004D","#A44A3F","#C1AAC0")
Pal <- c("#47682C","#119DA4","#D64933","#1F2041") #Group Colours: PANGOLIN,Reference,Sample,Wuhan

#nodetree<- ggtree(WtDTcol,layout="rectangular",alpha=0.3) + geom_text2(aes(label=node), hjust=-.3, size=1) + geom_tippoint(aes(colour=group),size=0.5,shape=19) +geom_tiplab(size=0.3)
#nodetree$data[nodetree$data$node %in% c(19,21), "x"] = mean(nodetree$data$x)
#nodetree

#ggtree
fullplot<-ggtree(WtDTcol,layout='circular',alpha=0.5,linetype=1,size=0.1) +
  scale_colour_manual(values=Pal) +
  
  geom_nodelab2(hjust=-0.2,size=2) +
#  geom_treescale(x = 0.01, y = 38, fontsize = 2, linesize = 0.1) +
  
  #geom_hilight(node=34, fill=CC[1], extendto=0.0004) +
  #geom_hilight(node=54, fill=CC[2], extendto=0.0004) +
  #geom_hilight(node=37, fill=CC[3], extendto=0.0004) +
  
  geom_tippoint(aes(colour=group),size=1,shape=19) +
  
  #geom_cladelabel(node=34,label="S Clade",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[1],align=T,linetype=NA) +
  #geom_cladelabel(node=54,label="G Clade 1",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[2],align=T,linetype=NA) +
  #geom_cladelabel(node=37,label="G Clade 2",hjust=0.5,offset=0.0001,offset.text=0.01,colour=CC[3],align=T,linetype=NA) +
  
  theme_tree(legend.position="none") 
fullplot
#Save
pdf(file="outputs/figures/Fig2.pdf", width=12, height=8)
   fullplot + geom_tiplab2(aes(colour=group), hjust=0, align=T, linetype=3, linesize=0.5) #dont forget ab tiplab2
dev.off()


