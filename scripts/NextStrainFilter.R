library(Biostrings)

#------------------------- User-defined parameters --------------------------
MisMatch<-5 # Genomes with more than this number of substitutions vs patient samples are removed
Trim5<-360 # Start position from 5' (bp before this number are deleted); NA or 0 for no trim
Trim3<-30420 # End position from 3' (bp after this number are deleted); NA 0 for no trim
#------------------------------------------------------------------------------

#Load data (see mafft.sh for details on how these files were created)
alignment <- readDNAStringSet("./intermediatedata/BRaligned.afa") # Alignment from Nextstrain (msa)
dif_mat <- read.csv("./intermediatedata/DistMat.csv",row.names=1)
names(dif_mat)<-c("Sample","Ref","Dist")

# Import PANGOLIN lineages
PangoLins<-read.csv("./inputdata/PangoLins.csv")

# Convert User-defined trim to sequence index
if(Trim5 %in% c(NA,0)) {
  Trim5<-1
}

if(Trim3 %in% c(NA,0)) {
  x<-lengths(alignment)
  Trim3<-max(x)
} 

# Simplify distance matrix data
## 1. Filter dissimilar (pt 1)
## Preliminary filter; smaller number speeds up code but may lose relevant sequences see below for final filter with MisMatch variable
keep<-dif_mat$Ref[dif_mat$Dist <= 20] 
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
  aligned_df$ID[grep(PangoLins$GISAID.ID[i],aligned_df$ID)]<-PangoLins$lineage[i]
}
# Add Wuhan if missing (e.g. renamed to PANGOLIN lineage A)

# Simplify IDs
aligned_df$ID<-gsub("^2019-nCoV.*(S[0-9]{1,2})$","\\1i",aligned_df$ID) # IonTorrent
aligned_df$ID<-gsub("^(S[0-9]{1,2})_NB.*nanopolish.*","\\1m",aligned_df$ID) # MinION
aligned_df$ID<-gsub("^hCoV-19.([A-z ]+).*EPI_ISL_([0-9]+).*","\\1.\\2",aligned_df$ID) # Reference Sequences
aligned_df$ID<-gsub("^S([0-9]){1}([im])$","S0\\1\\2",aligned_df$ID)# add leading 0s
aligned_df$ID[grep("Wuhan.406801",aligned_df$ID)] <- "Wuhan/WH04/2020" # Wuhan (root)

#Convert to dna matrix
FullMt <- data.frame(strsplit(aligned_df$seqs,""), stringsAsFactors = F)
names(FullMt) <- aligned_df$ID
FullMt <- data.frame(t(FullMt),stringsAsFactors = F)

#Convert patient samples only (for variant graph)
SimpleMt <- FullMt[grep("^S[0-9]+[im]$|Wuhan.WH04.2020",row.names(FullMt)),]

#Remove monomorphic sites and write to file:
# sequences from NextSeq msa 
polymorphicloci <-FullMt[,c(names(Filter(function(x) sum(unique(x) %in% c("A","T","G","C")) > 1, FullMt)))]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))
write.csv(polymorphicloci,"./intermediatedata/PolymorphicAlignmentFull.csv")

# Write polymorphic loci including only patient samples (excluding reference sequences)
# sequences from NextSeq msa, excluding reference samples (input for VariantMap.R)
polySamp <- SimpleMt[,c(names(Filter(function(x) sum(unique(x) %in% c("A","T","G","C")) > 1, SimpleMt)))]
colnames(polySamp) <- gsub("X","",colnames(polySamp))
write.csv(polySamp,"./intermediatedata/PolymorphicAlignmentSimple.csv")

# Write to Fasta
# NextSeq aligned
seqinr::write.fasta(sequences=paste(polymorphicloci[1,],collapse=""),names=row.names(polymorphicloci)[1],
                    file.out="./intermediatedata/NextFilt.afa",as.string=T)
for(i in 2:nrow(polymorphicloci)){
  seqinr::write.fasta(sequences=paste(polymorphicloci[i,],collapse=""),names=row.names(polymorphicloci)[i],
                      file.out="./intermediatedata/NextFilt.afa",as.string=T,
                      open="a")
}

