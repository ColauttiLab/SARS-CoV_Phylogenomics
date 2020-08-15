#06/02/2020 script
library(seqinr)
library(Biostrings)
library(dplyr)
library(rBLAST) 
library(plyr)

#-------------------------------------Functions-------------------------------------
##From phyloBlast package (WIP)
blast_seq<-function(refseq ,blastDB, mismatches, hits = NA){
  BR<-predict(blastDB, refseq)
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
  write.table(dataframe$SubjectID, file = "outputs/GBACC.txt", row.names = F, col.names = F, quote = F)
  system('scripts/retrieveSeq.sh')
  Seqs <- readDNAStringSet("outputs/seqs.txt")
  seqDF <- data.frame(SubjectID = paste(Seqs@ranges@NAMES), Seqs = paste(Seqs), stringsAsFactors = FALSE)
  return(seqDF)
} #Interacts with shell script to return blast result sequences

#Load data
samples <- readDNAStringSet("inputdata/QGLO_sequences08.fasta")
nextstrain_metadata <- read.csv("inputdata/metadata_0815.csv")
nextstrain_db <- blast("inputdata/nextstrainDB/nextstrainDB")

#Create DF
seqDF <- data.frame(Seqs = paste(samples), ID = gsub('2019-nCoV_MN908947\\|','',paste(samples@ranges@NAMES)), stringsAsFactors = F)

#blast each sequence to local database
datalist = list()
for (i in 1:1) {
  temp <- DNAStringSet(paste(seqDF[i,]$Seqs))
  blast_results <-blast_seq(temp, nextstrain_db, mismatches = 2)
  blast_results$Sample<-seqDF[i,]$ID
  blast_results$Qlen <- nchar(paste(seqDF[i,]$Seqs)) # Length of query sequence
  datalist[[i]] <- blast_results[2:nrow(blast_results),] # Remove top hit (QGOL seqs already uploaded)
}

#Transform list into dataframe
f_blast_results <- dplyr::bind_rows(datalist)

#Remove duplicate blast results by accession code
uniq_blast_results <- f_blast_results[ !duplicated(f_blast_results[ , 2] ) , ]

# Retreive seqs
blast_results_seq <- retrieve_seqs(uniq_blast_results)
blast_results_seq[,1] <- gsub(" ", "", blast_results_seq[,1])
blast_results <- merge(uniq_blast_results, blast_results_seq, by.y = "SubjectID")

# Remove blast results with alignment length < 20000
blast_results <- blast_results[blast_results$Alignment.Length > 20000,]
# Manual trim, directly to size of query sequence
blast_results$Seqs <- substr(blast_results$Seqs, blast_results$S.start, (blast_results$Qlen +blast_results$S.start))

# Add Wuhan Samples
nextstrain_samples <- readDNAStringSet("inputdata/sequences_0815.fasta")
wuhan_samples <- nextstrain_samples[ grep("Wuhan", nextstrain_samples@ranges@NAMES) ]
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

#To Fasta
seqinr::write.fasta(paste(finalDF$Seqs), paste0(finalDF$ID), "outputs/nextstrain_br.fasta")

#Align using MAFFT:
# % mafft --auto input > output

## Moved June 4 code to new file "NextStrainAnalysis_Main.R"