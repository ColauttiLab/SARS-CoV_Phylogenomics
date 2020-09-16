library(Biostrings)

#------------------------- User-defined parameters
# Genomes with more than this number of substitutions vs patient samples are removed
MisMatch<-1 
# Start position from 5' (bp before this number are deleted); NA or 0 for no trim
Trim5<-359 
# End position from 3' (bp after this number are deleted); NA 0 for no trim
Trim3<-30478 
## IMPORTANT: Trim5 and Trim3 here should match Distcalc.R
#-------------------------------------------------

#Load data (see mafft.sh for details on how these files were created)
## aligned sequences
# Alignment from Nextstrain (msa)
alignIn <- readDNAStringSet("./intermediatedata/BRaligned.afa") 
## Distance calculations
dif_mat <- read.csv("./intermediatedata/DistMat.csv",row.names=1)
names(dif_mat)<-c("Sample","Ref","Dist")
## PANGOLIN lineages
PangoLins<-read.csv("./inputdata/PangoLins.csv")

# remove misaligned
x<-lengths(alignIn)
alignment<-alignIn[x==30958]

# Set trim sites if not user-defined
if(Trim5 %in% c(NA,0)){
  Trim5<-1
}
if(Trim3 %in% c(NA,0)){
  Trim3<-max(x)
}


x<-NULL



# Simplify distance matrix data
## 1. Filter dissimilar
keep<-dif_mat$Ref[dif_mat$Dist <= MisMatch] 
## 2. Remove patient samples already uploaded to GISAID ('QGLO')
keep<-keep[grep("QGLO",keep,invert=T)]
## 3. Add Wuhan basal (root of phylogeny)
keep<-c(keep,dif_mat$Ref[grep("Wuhan.WH04.2020",dif_mat$Ref)][1])
keep[length(keep)]<-"Wuhan/WH04/2020"
## 4. Add PANGOLIN lineages
keep<-c(keep,dif_mat$Ref[gsub(".*(EPI_ISL_[0-9]+).*","\\1",dif_mat$Ref) %in% PangoLins$GISAID.ID])
## 5. Remove duplicate Ref IDs
keep<-keep[!duplicated(keep)]
## 6. Add samples
keep<-c(keep,names(alignment[grep("nanopore|Sample",names(alignment))]))

# Subset alignment to include only root, patient samples, and filtered references
## Index patient samples from this study
sub_align <- alignment[ names(alignment) %in% keep]

# Parse to data frame and simplify names
# Data frame for NextStrain Alignment msa
aligned_df <- data.frame(ID = gsub("Sample[ _]", "S", sub_align@ranges@NAMES), 
                         seqs = substr(paste(sub_align),Trim5,Trim3),
                         Region = gsub("\\_.*", "", sub_align@ranges@NAMES), 
                         Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",sub_align@ranges@NAMES),
                         stringsAsFactors = FALSE)
# Replace PANGOLIN ID with Lineage name
for(i in 1:nrow(PangoLins)){
  print(paste(aligned_df$ID[grep(PangoLins$GISAID.ID[i],aligned_df$ID)],"<-",PangoLins$lineage[i]))
  aligned_df$ID[grep(PangoLins$GISAID.ID[i],aligned_df$ID)]<-PangoLins$lineage[i]
}
# Add Wuhan if missing (e.g. renamed to PANGOLIN lineage A)
if(length(grep("Wuhan.WH04.2020",aligned_df$ID)) == 0 ){
  Wuhan<-aligned_df[aligned_df$ID == "A.1",]
  Wuhan$ID<-"Wuhan/WH04/2020"
  aligned_df<-rbind(aligned_df,Wuhan)
}

# Simplify IDs
aligned_df$ID<-gsub("^2019-nCoV.*(S[0-9]{1,2})$","\\1i",aligned_df$ID) # IonTorrent
aligned_df$ID<-gsub("^(S[0-9]{1,2})_NB.*nanopolish.*","\\1m",aligned_df$ID) # MinION
aligned_df$ID<-gsub("^hCoV-19.([A-z ]+).*EPI_ISL_([0-9]+).*","\\1.\\2",aligned_df$ID) # Reference Sequences
aligned_df$ID<-gsub("^S([0-9]){1}([im])$","S0\\1\\2",aligned_df$ID)# add leading 0s

#Convert to dna matrix
FullMt <- data.frame(strsplit(aligned_df$seqs,""), stringsAsFactors = F)
names(FullMt) <- aligned_df$ID
FullMt <- data.frame(t(FullMt),stringsAsFactors = F)

#Convert patient samples only (for variant graph)
SimpleMt <- FullMt[grep("^S[0-9]+[im]$|Wuhan.WH04.2020",row.names(FullMt)),]

# Identify Remove monomorphic sites (of patient samples relative to Wuhan reference) =
keepLoci<-c(names(Filter(function(x) sum(unique(x) %in% c("A","T","G","C")) > 1, SimpleMt)))

# Remove monomorphic sites and write to file
## Write polymorphic loci including only patient samples (excluding reference sequences)
## sequences from NextSeq msa, excluding reference samples (input for VariantMap.R)
polySamp <- SimpleMt[,keepLoci]
colnames(polySamp) <- gsub("X","",colnames(polySamp))
write.csv(polySamp,"./intermediatedata/PolymorphicAlignmentSimple.csv")

## sequences from NextSeq msa 
polymorphicloci <-FullMt[,keepLoci]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))
write.csv(polymorphicloci,"./intermediatedata/PolymorphicAlignmentFull.csv")


# Write to Fasta
# NextSeq aligned
seqinr::write.fasta(sequences=paste(polymorphicloci[1,],collapse=""),names=row.names(polymorphicloci)[1],
                    file.out="./intermediatedata/NextFilt.afa",as.string=T)
for(i in 2:nrow(polymorphicloci)){
  seqinr::write.fasta(sequences=paste(polymorphicloci[i,],collapse=""),names=row.names(polymorphicloci)[i],
                      file.out="./intermediatedata/NextFilt.afa",as.string=T,
                      open="a")
}

