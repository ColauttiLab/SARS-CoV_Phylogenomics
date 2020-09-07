library(seqinr)
library(Biostrings)
library(dplyr)
library(rBLAST) 
library(plyr)
library(ape)

#------------------------- User-defined parameters --------------------------
MisMatch<-3 # Database genomes with more than this number of mismatches from patient samples are excluded from analysis
MinLen<-10000 # Minimum alignment length 
#------------------------------------------------------------------------------

#Load data
alignment <- readDNAStringSet("./intermediatedata/BRalignedOLD.afa")
# remove misaligned
x<-lengths(alignment)
alignment<-alignment[x==max(x)]
x<-NULL

# Index patient samples from this study
QGLO<-grep("QGLO",names(alignment))
Samples<-grep("Sample.*",names(alignment))
Refs<-grep(".*",names(alignment[-c(Samples,QGLO)]))


distcalc<-

# Patient samples
#patseq<-alignment[c(QGLO,IonSeq,NanoSeq)]

# GISAID genomes with patient samples removed
#refseq<-alignment[-c(QGLO,IonSeq,NanoSeq)]

Rows<-length(Samples)*length(Refs)

DistMat<-data.frame(Query=rep(NA,Rows),Subject=rep(NA,Rows),Dist=rep(NA,Rows))

dCalc<-function(p,r){
  return(c(names(alignment[p]),
                    names(alignment[r]),
                    stringDist(alignment[c(p,r)],method="hamming")))
}

system.time({
Row<-1
for(i in 1:length(Samples)){
  for(j in 1:100){
    DistMat[Row,]<-dCalc(i,j)
    Row<-Row+1
  }
}

})




samples <- readDNAStringSet("./inputdata/IonTorrent_MinION.fasta")
bhits <- read.delim("./intermediatedata/BlastResults.out",header=F)
names(bhits)<-c("QueryID","SubjectID","Perc.Ident","Alignment.Length","Mismatches",
                "Gap.Openings","Q.start","Q.end","S.start","S.end",
                "E","Bits")
# Remove dissimilar (non-informative) sequences
filtered_bhits <- bhits[ bhits$Mismatches <= MisMatch , ]
# Remove partial genome matches
filtered_bhits <- filtered_bhits[ filtered_bhits$Alignment.Length >= MinLen , ]

#Create data.frame of patient sequence data with suffix i for IonTorrent and m for MinION
#seqDF <- data.frame(Seqs = paste(samples), ID = gsub('2019-nCoV_MN908947\\|Sample[ _](.*)','S\\1i',paste(samples@ranges@NAMES)), stringsAsFactors = F)
#seqDF <- data.frame(Seqs = paste(samples), ID = paste(samples@ranges@NAMES), stringsAsFactors = F)
#seqDF$ID <- gsub("^2019-nCoV_[^|]*\\|Sample[ _](.*)","S\\1i",seqDF$ID)  #Ion Torrent
#seqDF$ID <- gsub("^Sample_([0-9]+)_.*ARTIC\\/nanopolish.*","S\\1m",seqDF$ID) #MinION
#seqDF$Seqs <- gsub("-","N",seqDF$Seqs) # replace gaps/missing with N for BLAST

#Remove duplicate blast results by accession code
uniq_bhits <- filtered_bhits[ !duplicated(filtered_bhits$SubjectID ) , ]

# Remove non-informative sequences
## First create simplified vector names for comparison
alnames <- gsub(".*EPI_ISL_([0-9]+).*","\\1",names(alignment))
blnames <- paste(uniq_bhits$SubjectID)

# Keep BLAST hits
aligned_keep <- which(alnames %in% blnames)

# add original patient samples
patient_seqs <- which(alnames %in% names(samples))

# add Wuhan Samples
Wuhan <- grep("Wuhan.WH04.2020|Wuhan.Hu-1.2019",names(alignment))

# Use indices to subset original aligned fasta file
blast_results_seq <- alignment[c(Wuhan,patient_seqs,aligned_keep)]

# remove QGLO sequences in alignment (already published in the GISAID database)
blast_results_out <- blast_results_seq[grep("ON_QGLO-",names(blast_results_seq),invert=T)]

# Save Output
write.dna(x=blast_results_out,file="./intermediatedata/nextstrain_br.fasta",format="fasta",nbcol=1,colw=60)
