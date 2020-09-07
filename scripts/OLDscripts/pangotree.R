library(Biostrings)
library(phangorn)
library(ips)
library(magrittr)
library(dplyr)
library(ggtree)
#library(ape)

#------------------------- User-defined parameters --------------------------
Trim5<-241 # Number of bp to trim from 5' 
Trim3<-29688 # Number of bp to trim from 3'
#------------------------------------------------------------------------------

#Load data
# Pangolin lineages from: Rambaut et al. PANGOLIN paper (https://doi.org/10.1101/2020.04.17.046086)
# available at https://github.com/cov-lineages/lineages/tree/master/lineages/data
alignment<-readDNAStringSet("./intermediatedata/WuhanPIM.afa")

# Convert User-defined trim to sequence index
if(Trim5 %in% c(NA,0)) {
  Trim5<-1
}

if(Trim3 %in% c(NA,0)) {
  x<-lengths(alignment)
  Trim3<-max(x)
} 

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

# Parse to data frame and simplify names
aligned_df <- data.frame(ID = gsub("Sample[ _]", "S", alignment@ranges@NAMES), 
                         # Want to get rid of substr... The trimEnds function should remove boarder differences so that substr not needed
                         # For example, without substr uniq_samples(line79) = 24objs, with uniq_samples = 23 objs. Difference due to cutting out first 11 chars
                         seqs = substr(paste(alignment),Trim5,Trim3),
                         Region = gsub("\\_.*", "", alignment@ranges@NAMES), 
                         Country = gsub(".*\\_(.*)\\_.*\\(.*","\\1",alignment@ranges@NAMES),
                         stringsAsFactors = FALSE)

# Simplify IDs
aligned_df$ID<-gsub("^2019-nCoV.*(S[0-9]{1,2})$","\\1i",aligned_df$ID) # IonTorrent
aligned_df$ID<-gsub("^(S[0-9]{1,2})_NB.*nanopolish.*","\\1m",aligned_df$ID) # MinION
aligned_df$ID[grep("Wuhan",aligned_df$ID)] <- "Wuhan/WH04/2020" # Wuhan (root)

#Convert to dna matrix
FullMt <- data.frame(strsplit(aligned_df$seqs,""), stringsAsFactors = F)
names(FullMt) <- aligned_df$ID
FullMt <- data.frame(t(FullMt),stringsAsFactors = F)
#Convert patient samples only (vor variant graph)
SimpleMt <- FullMt[grep("^S[0-9]+[im]$|Wuhan",row.names(FullMt)),]

#Remove monomorphic sites and write to file:
polymorphicloci <-FullMt[,c(names(Filter(function(x) length(unique(x)) != 1, FullMt)))]
colnames(polymorphicloci) <- gsub("X","",colnames(polymorphicloci))






# Simplify sample names
# Add labels for IonTorrent vs Nanopore sequences
IonSeq<-grep("S[0-9]+i$",aligned_df$ID)
NanoSeq<-grep("*nanopolish.*",aligned_df$ID)
aligned_df$ID[NanoSeq]<-gsub("(S[0-9]+).*","\\1m",aligned_df$ID[NanoSeq])
# Add leading zeros for 1-digit patient ID
aligned_df$ID[c(NanoSeq,IonSeq)] <- gsub("S([0-9]{1}[im])","S0\\1",aligned_df$ID[c(NanoSeq,IonSeq)])












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