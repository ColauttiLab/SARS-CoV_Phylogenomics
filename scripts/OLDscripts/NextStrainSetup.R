#06/02/2020 script
library(seqinr)
library(Biostrings)
library(dplyr)
library(rBLAST) 
library(plyr)
library(seqinr)

#------------------------- User-defined parameters --------------------------
MisMatch<-3 # Database genomes with more than this number of mismatches from patient samples are excluded from analysis
Bargs<-"-strand plus -word_size 100 -gapopen 0 -gapextend 0 -max_hsps 1000 -culling_limit 10 -max_target_seqs 1000" # BLAST arguments (e.g. retain 100k hits, remove alignments< 10k)
#------------------------------------------------------------------------------

#-------------------------------------Functions-------------------------------------
##From phyloBlast package (WIP)
blast_seq<-function(refseq ,blastDB, mismatches, hits = NA){
  BR<-predict(blastDB, refseq, BLAST_args=Bargs) # Limit matches to >=15kb and retain as many as possible 
  BRhitLen<-length(BR$SubjectID)
  if( !is.na(hits) ) {
    if (BRhitLen > hits ) {
      BRhitLen=hits #Limit hit length
      BR<-as.vector(BR[1:BRhitLen,]) 
    }
  }
  BR <- BR[BR$Mismatches <= mismatches, ]
  return(BR)
} #Perform blast search and return dataframe of blast hits that meet criteria

retrieve_seqs <- function(dataframe, dbpath) {
  write.table(dataframe$SubjectID, file = "./intermediatedata/GBACC.txt", row.names = F, col.names = F, quote = F)
  system('./scripts/retrieveSeq.sh')
  Seqs <- readDNAStringSet("./intermediatedata/seqs.txt")
  seqDF <- data.frame(SubjectID = paste(Seqs@ranges@NAMES), Seqs = paste(Seqs), stringsAsFactors = FALSE)
  return(seqDF)
} #Interacts with shell script to return blast result sequences

#Load data
#samples <- readDNAStringSet("./inputdata/IonTorrent_MinION.fasta")
nextstrain_db <- blast("./intermediatedata/nextstrain_sequences/nextstrainDB")
insamp<-readDNAStringSet("./intermediatedata/BRaligned.afa")
samples<-insamp[grep("^2019-nCoV_.*|nanopolish",names(insamp))]
#Create DF with suffix i for IonTorrent and m for MinION
#seqDF <- data.frame(Seqs = paste(samples), ID = gsub('2019-nCoV_MN908947\\|Sample[ _](.*)','S\\1i',paste(samples@ranges@NAMES)), stringsAsFactors = F)
seqDF <- data.frame(Seqs = paste(samples), ID = paste(samples@ranges@NAMES), stringsAsFactors = F)
seqDF$ID<-gsub("^2019-nCoV_[^|]*\\|Sample[ _](.*)","S\\1i",seqDF$ID)  #Ion Torrent
seqDF$ID<-gsub("^Sample_([0-9]+)_.*ARTIC\\/nanopolish.*","S\\1m",seqDF$ID) #MinION
seqDF$ID<-gsub("S([0-9]{1})([im])","S0\\1\\2",seqDF$ID) # Add leading 0 so all samples are 2-digit
seqDF$Seqs<-gsub("-","N",seqDF$Seqs) # replace gaps/missing with N for BLAST

#blast each sequence to local database
datalist = list()
for (i in 1:nrow(seqDF)) {
  temp <- DNAStringSet(paste(seqDF[i,]$Seqs))
  tmp_blast_results <-blast_seq(temp, nextstrain_db, mismatches = MisMatch)
  tmp_blast_results$Sample<-seqDF[i,]$ID
  tmp_blast_results$Qlen <- nchar(paste(seqDF[i,]$Seqs)) # Length of query sequence
  datalist[[i]] <- tmp_blast_results 
  temp<-tmp_blast_results<-NA
}

#Transform list into dataframe
f_blast_results <- dplyr::bind_rows(datalist)

#Remove duplicate blast results by accession code
uniq_blast_results <- f_blast_results[ !duplicated(f_blast_results[ , 2] ) , ]

# Remove redundant QGLO sequences (Samples already uploaded to GISAID)
uniq_blast_results<-uniq_blast_results[grep("QGLO",uniq_blast_results$SubjectID,invert=T),]

# Retreive seqs
blast_results_seq <- retrieve_seqs(uniq_blast_results)
blast_results_seq$SubjectID <- gsub(" $", "", blast_results_seq$SubjectID) # Remove space at end of subjectID
blast_results <- merge(uniq_blast_results, blast_results_seq, by.y = "SubjectID")

# Manual trim, directly to size of query sequence
blast_results$Seqs <- substr(blast_results$Seqs, blast_results$S.start, (blast_results$Qlen +blast_results$S.start))

# Add Wuhan Samples
nextstrain_samples <- readDNAStringSet("./inputdata/nextstrain_sequences.fasta")
wuhan_samples <- nextstrain_samples[ grep("Wuhan.WH04.2020|Wuhan.Hu-1.2019", nextstrain_samples@ranges@NAMES) ]
wuhan_df <- data.frame(SubjectID=paste(wuhan_samples@ranges@NAMES), Seqs=paste(wuhan_samples))
blast_results <- rbind.fill(blast_results,wuhan_df)

# Add metadata
names(nextstrain_metadata)[names(nextstrain_metadata) == 'strain'] <- 'SubjectID'
nextstrain_final <- merge(blast_results, nextstrain_metadata, by.y = "SubjectID")

# Rename including region_country_district (strain)
nextstrain_final$ID <-paste0(nextstrain_final$region, "_", nextstrain_final$country, "_", 
                             nextstrain_final$division, " (", nextstrain_final$SubjectID,")")

# Add samples to br dataframe
finalDF<-rbind.fill(nextstrain_final,seqDF) #Add target DF to main DF

#Write To Fasta
seqinr::write.fasta(sequences=finalDF$Seqs[1],names=paste(finalDF$ID[1]),
                    file.out="./intermediatedata/nextstrain_br.fasta",as.string=T)
for(i in 2:nrow(finalDF)){
  seqinr::write.fasta(sequences=finalDF$Seqs[i],names=paste(finalDF$ID[i]),
                      file.out="./intermediatedata/nextstrain_br.fasta",as.string=T,
                      open="a")
}

#Next: Align using MAFFT (run scripts/mafftalign.sh:
#Then: Pre-trim phylogeny

## Moved June 4 code to new file "NextStrainAnalysis_Main.R"