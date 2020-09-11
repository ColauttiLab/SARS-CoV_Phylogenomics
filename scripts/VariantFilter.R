library(Biostrings)

# Import data
RefAlign <- readDNAStringSet("./intermediatedata/RefWIM.afa")  # Alignment file containing reference genome for IonTorrent & MinION

# change to data.frame
RefAlign_df <- data.frame(ID = gsub("Sample[ _]", "S", RefAlign@ranges@NAMES), 
                          seqs = paste(RefAlign),
                          Region = gsub("\\_.*", "", RefAlign@ranges@NAMES), 
                          Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",RefAlign@ranges@NAMES),
                          stringsAsFactors = FALSE)

# Remove Reference sequences
RefAlign_df<-RefAlign_df[grep("ref",RefAlign_df$ID,invert=T),]

# Simplify ID names
RefAlign_df$ID<-gsub("^2019-nCoV.*(S[0-9]{1,2})$","\\1i",RefAlign_df$ID) # IonTorrent
RefAlign_df$ID<-gsub("^(S[0-9]{1,2})_NB.*nanopolish.*","\\1m",RefAlign_df$ID) # MinION
RefAlign_df$ID<-gsub("^hCoV-19.([A-z ]+).*EPI_ISL_([0-9]+).*","\\1.\\2",RefAlign_df$ID) # Reference Sequences
RefAlign_df$ID<-gsub("^S([0-9]){1}([im])$","S0\\1\\2",RefAlign_df$ID)# add leading 0s
RefAlign_df$ID[grep("Wuhan.406801",RefAlign_df$ID)] <- "Wuhan/WH04/2020" # Wuhan (root)

# Convert to DNA matrix
RefMt <- data.frame(strsplit(RefAlign_df$seqs,""), stringsAsFactors = F)
names(RefMt) <- RefAlign_df$ID
RefMt <- data.frame(t(RefMt),stringsAsFactors = F)

# Remove monomoprhic sites and save as csv
polyRef <-RefMt[,c(names(Filter(function(x) sum(unique(x) %in% c("A","T","G","C")) > 1, RefMt)))]
colnames(polyRef) <- gsub("X","",colnames(polyRef))
write.csv(polyRef,"./intermediatedata/PolyAlignRef.csv")
