#06/02/2020 script
library(phangorn)
library(seqinr)
library(Biostrings)
library(dplyr)
library(ape)
library(ggtree)
library(ggplot2)
library(rBLAST) 
library(genbankr) 
library(seqRFLP)
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

blast_trim <- function(dataframe){ # add function for trim start and end
  dataframe$Seqs <- substr(dataframe$Seqs, dataframe$S.start, dataframe$S.end)
  return(dataframe)
} #Trims sequences based on blast return query

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
}  #Maximum likelihood

bootstrap_tree<-function(fitted_model,bs_iterations,scale_bar,out){
  bs <- bootstrap.pml(fitted_model, bs=bs_iterations, optNni=TRUE, multicore=F, control = pml.control(trace=0))
  pdf(out,width=13,height=10)
  bstree<-plotBS(midpoint(fitted_model$tree),bs, p = 70, type="p")
  add.scale.bar(length = scale_bar, cex = 0.9, font = 2)
  dev.off()
  return(bstree)
} #Bootstrap

#-------------------------------------May 15-------------------------------------
##Blasting sequences on nextstrain database (excluding samples 19,21)
#Load data
samples <- readDNAStringSet("inputdata/QGLO_sequences.fasta")
nextstrain_metadata <- read.csv("inputdata/nextstrain_metadata_051520.csv")
nextstrain_db <- blast("inputdata/nextstrain_sequences_051520/nextstrain_sequences_051520.fasta")

#Remove outliers 
# samples<- samples[-c(8, 10), ] 

#Create DF
seqDF <- data.frame(Seqs = paste(samples), ID = gsub('2019-nCoV_MN908947\\|','',paste(samples@ranges@NAMES)))

#Loot blast_seq function on each sequences and store in list
datalist = list()
for (i in 1:nrow(seqDF)) {
  temp <- DNAStringSet(paste(seqDF[i,]$Seqs))
  blast_results <-blast_seq(temp, nextstrain_db, mismatches = 2)
  blast_results$Sample<-seqDF[i,]$ID
  datalist[[i]] <- blast_results
}

#Transform list into dataframe
f_blast_results <- dplyr::bind_rows(datalist)

#Remove duplicate blast results by accession code
uniq_blast_results <-  f_blast_results[ !duplicated(f_blast_results[ , 2] ) , ]

# Retreive seqs
blast_results_seq <- retrieve_seqs(uniq_blast_results)
blast_results_seq[,1] <- gsub(" ", "", blast_results_seq[,1])
blast_results <- merge(uniq_blast_results, blast_results_seq, by.y = "SubjectID")

#Trim blast results to search queries S.start and S.send
trimmed_br <- blast_trim(blast_results)

# Add metadata
names(nextstrain_metadata)[names(nextstrain_metadata) == 'strain'] <- 'SubjectID'
nextstrain_final <- merge(trimmed_br, nextstrain_metadata, by.y = "SubjectID")

# Rename including region_country_district (strain)
nextstrain_final$ID <-paste0(nextstrain_final$region, "_", nextstrain_final$country, "_", 
                             nextstrain_final$division, " (", nextstrain_final$SubjectID,")")

# Add samples to br dataframe
finalDF<-rbind.fill(nextstrain_final,seqDF) #Add target DF to main DF

# Add Wuhan Samples
nextstrain_samples <- readDNAStringSet("~/ColauttiLabScratch/COVID-19/Data/05_15_2020/nextstrain_sequences_05152020.fasta")#readDNAStringSet("inputdata/nextstrain_sequences_05152020.fasta")
wuhan_samples <- nextstrain_samples[ grep("Wuhan", nextstrain_samples@ranges@NAMES) ]
wuhan_df <- data.frame(ID=paste(wuhan_samples@ranges@NAMES), Seqs=paste(wuhan_samples))
finalDF<-rbind.fill(finalDF,wuhan_df)

#To Fasta
seqinr::write.fasta(paste(finalDF$Seqs), paste0(finalDF$ID), "intermediatedata/nextstrain_br.fasta")

#Align using MAFFT:
# % mafft --auto input > output

## Moved June 4 code to new file "NextStrainAnalysis_Main.R"