library(phangorn)
library(Biostrings)
library(ips)
library(magrittr)
library(dplyr)
library(ggtree)

# ------------ User-Defined Parameters ------------
#Trim5<-175 #Number of bp in the alignment to trim from 5' end 
#Trim3<-325 #Number of bp in the alignment to trim from 3' end
TrimSeq<-10000  #Variable to play with (see ?trimEnds documentation). This affects # of objs in uniq_samples by trimming
# So for min.n.seq = 100, the function will cut the sequences at the position where 100 sequences have unambigiuos bases, and it will fill sequences that fall short with N's
# So a greater min.n.seq = less filling in with N's. 
# Notice how changing min.n.seq will change max(aligned_seqs@ranges@width)
# if min.n.seq = the number of sequences, then the function will cut the at the spot where there is a single unambigious base (too strict).
# --------------------------------------------------

# Trim alignment
alignment <- read.FASTA("./intermediatedata/BRaligned.afa")
# Remove  NextStrain sequences with improper alignment length (?!?)
x<-sapply(alignment,length)
alignment<-alignment[x==max(x)]
x<-NULL

#TrimFun<-function(x){
#  return(x[Trim5:(length(alignment)-Trim3)])
#}
#trimmed_alignment<-lapply(alignment,TrimFun)

#checkAlignment(test, plot = FALSE)
trimmed_alignment <- trimEnds(as.matrix(alignment), min.n.seq = TrimSeq) ##checkAlignment(test, plot = FALSE)
aligned_seqs <- trimmed_alignment %>% as.list %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

aligned_df <- data.frame(ID = aligned_seqs@ranges@NAMES, 
                         # Want to get rid of substr... The trimEnds function should remove boarder differences so that substr not needed
                         # For example, without substr uniq_samples(line79) = 24objs, with uniq_samples = 23 objs. Difference due to cutting out first 11 chars
                         seqs = paste(aligned_seqs),
                         Region = gsub(".*\\|([A-z]+)$", "\\1", aligned_seqs@ranges@NAMES), 
                         Country = gsub(".*\\/([A-z]+)\\/.*","\\1",aligned_seqs@ranges@NAMES),
                         stringsAsFactors = FALSE)

# Remove QGLO sequences already published in GISAID
aligned_df <- aligned_df[grep("QGLO",aligned_df$ID,invert=T),]

# Convert - to N 
aligned_df$seqs<-gsub("-","N",aligned_df$seqs)

# Simplify sample names
aligned_df$ID <-gsub("hCoV-19\\/[A-z _-]+\\/[^|]+\\|([^|]+).*","\\1",aligned_df$ID)

# Modify labels for IonTorrent vs Nanopore sequences
IonSeq<-grep("^2019-nCoV_.*",aligned_df$ID)
NanoSeq<-grep("nanopolish",aligned_df$ID)

aligned_df$ID[IonSeq]<-gsub("^2019-nCoV_[^|]*\\|Sample[ _](.*)","S\\1i",aligned_df$ID[IonSeq]) 
aligned_df$ID[NanoSeq]<-gsub("^Sample_([0-9]+)_.*ARTIC\\/nanopolish.*","S\\1m",aligned_df$ID[NanoSeq])
aligned_df$ID<-gsub("^S([0-9]{1})([im])","S0\\1\\2",aligned_df$ID) # Add leading 0

# Separate sample from reference sequences
Samples <- aligned_df[c(IonSeq,NanoSeq), ] 
aligned_df <- aligned_df[-c(IonSeq,NanoSeq),]

#Combine patient samples with identical sequences
uniq_samples <- Samples[ !duplicated(Samples[ , 2] ) , ] #Combine duplicate sample sequences (Confirm with real # of duplicates)
DupSampleList <-list()
for (x in 1:nrow(uniq_samples)) { #Iterate through sequences and group unique sequences together
  dupSeq <- paste(uniq_samples$seqs[x])
  matches <- Samples[ which(paste(Samples$seqs)==dupSeq), ]
  print(paste(x,"=",matches$ID))
  if (nrow(matches)>1) {
    name <- paste(matches$ID, collapse=",")
    tempDF <- data.frame(ID=name, seqs=paste(matches$seqs[1]), Region=name,Country=name)
    DupSampleList[[x]] <- tempDF
  } else {
    DupSampleList[[x]] <- matches
  }
  #dupSeq<-matches<-tempDF<-name<-NA
}
DupSampleList<-DupSampleList[lengths(DupSampleList) != 0]
Samples <- dplyr::bind_rows(DupSampleList)

# Combine reference samples with identical sequences
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
  dupSeq<-matches<-tempDF<-name<-NA
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
rootdf <- aligned_df[grep("Wuhan", aligned_df$ID),]
Samples <- rbind(Samples, rootdf)

#Convert to dna matrix
SampleMt <- data.frame(strsplit(Samples$seqs,""), stringsAsFactors = F)
colnames(SampleMt) <- Samples$ID
SampleMt <- data.frame(t(SampleMt),stringsAsFactors = F)


#Might want to trim ends here
#View(polymorphicloci)
#example: polymorphicloci <- polymorphicloci[,6:(ncol(polymorphicloci)-15)]

#Write reference sequence table to file
write.csv(RefDF,"./outputs/ReferenceTable.csv",row.names=F)

#Remove monomorphic sites and write to file:
polymorphicloci <-SampleMt[,c(names(Filter(function(x) length(unique(x)) != 1, SampleMt)))]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))
write.csv(polymorphicloci,"./intermediatedata/PolymorphicAlignmentFull.csv")
# Write polymorphic loci including only patient samples (excluding reference sequences)
polySamp<-polymorphicloci[grep("r_",row.names(polymorphicloci),invert=T),]
write.csv(polySamp,"./intermediatedata/PolymorphicAlignmentSimple.csv")
