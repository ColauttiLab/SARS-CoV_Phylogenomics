library(Biostrings)

#------------------------- User-defined parameters --------------------------
Trim5<-359 # Start position from 5' (bp before this number are deleted); NA or 0 for no trim
Trim3<-30478 # End position from 3' (bp after this number are deleted); NA 0 for no trim
## IMPORTANT: Trim5 and Trim3 here should match Distcalc.R
#------------------------------------------------------------------------------

#Load data (see mafft.sh for details on how these files were created)
## aligned sequences
alignIn <- readDNAStringSet("./intermediatedata/BRaligned.afa") # Alignment from Nextstrain (msa)

# remove misaligned
x<-lengths(alignIn)
alignment<-alignIn[x==30958] # Remove misaligned

# Set trim sites if not user-defined
if(Trim5 %in% c(NA,0)){
  Trim5<-1
}
if(Trim3 %in% c(NA,0)){
  Trim3<-max(x)
}

x<-NULL

# Subset alignment to include only root, patient samples, and filtered references
## Index patient samples from this study
sub_align <- alignment[ grep("EPI_ISL_402124|MN908947",names(alignment)) ]

# Parse to data frame and simplify names
# Data frame for NextStrain Alignment msa
aligned_df <- data.frame(ID = gsub("Sample[ _]", "S", sub_align@ranges@NAMES), 
                         seqs = substr(paste(sub_align),Trim5,Trim3),
                         Region = gsub("\\_.*", "", sub_align@ranges@NAMES), 
                         Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",sub_align@ranges@NAMES),
                         stringsAsFactors = FALSE)
# Add Wuhan if missing (e.g. renamed to PANGOLIN lineage A)

# Simplify IDs
aligned_df$ID<-gsub("^2019-nCoV.*(S[0-9]{1,2})$","\\1i",aligned_df$ID) # IonTorrent
aligned_df$ID<-gsub("^(S[0-9]{1,2})_NB.*nanopolish.*","\\1m",aligned_df$ID) # MinION
aligned_df$ID<-gsub("^hCoV-19.Wuhan.*","Wuhan",aligned_df$ID) # Reference Sequences
aligned_df$ID<-gsub("^S([0-9]){1}([im])$","S0\\1\\2",aligned_df$ID)# add leading 0s

#Convert to dna matrix
FullMt <- data.frame(strsplit(aligned_df$seqs,""), stringsAsFactors = F)
names(FullMt) <- aligned_df$ID
FullMt <- data.frame(t(FullMt),stringsAsFactors = F)

# Identify Remove monomorphic sites (of patient samples relative to Wuhan reference) =
keepLoci<-c(names(Filter(function(x) sum(unique(x) %in% c("A","T","G","C")) > 1, FullMt)))

# Remove monomorphic sites and write to file
## Write polymorphic loci including only patient samples (excluding reference sequences)
## sequences from NextSeq msa, excluding reference samples (input for VariantMap.R)
polymorph <- FullMt[,keepLoci]
colnames(polymorph) <- gsub("X","",colnames(polymorph))

# Create variant map (check against Wuhan reference to only ID variant sites)
for (x in 1:ncol(polymorph)) {
  refNt <- polymorph[grep("Wuhan",rownames(polySamp)),x] #Change index number
  polymorph[which(polymorph[,x] == refNt),x]<-" "
}

suborder <- integer()
for (y in 1:nrow(polymorph)) {
  nsub<-length(which(polymorph[y,]==" "))
  suborder<-c(suborder,nsub)
}

polymorph$sample<-rownames(polymorph)
long<- polymorph %>% 
  pivot_longer(-sample, names_to = "position", values_to = "substitution")

long<-long[grep("[ATGC]",long$substitution),]

SampOrd<-order(as.numeric(gsub("[SW]([0-9]*).*","\\1",unique(long$sample))),decreasing=T) #NA okay here
PosOrder<-order(as.numeric(unique(long$position)))

long$position <- factor(long$position,levels=unique(long$position)[PosOrder],order=T) # Use this to map variable sites only
#long$position<-as.numeric(long$position) # Use this to map sites on the full-length genome
long$sample <- factor(long$sample,levels=unique(long$sample)[SampOrd],order=T)

palette <- c("#007FFF","#99001C","#008000","#fffc00")

# Variant plot with equal spacing along y-axis
VarPlot<-ggplot(long, aes(x=position, y=sample)) + 
  geom_tile(aes(colour=I("white"), fill=substitution),size=0.7) + 
  scale_fill_manual(values=palette) +
  #  scale_colour_manual(values=palette) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

VarPlot
