library(Biostrings)
library(ape)
library(foreach)
library(doParallel)

#------------------------- User-defined parameters --------------------------
Csave<-4 # Number of cores to keep open, to use for other tasks (>=1 for personal computer, 0 for cluster)
#------------------------------------------------------------------------------

#Load data
alignIn<-readDNAStringSet("./intermediatedata/BRaligned.afa")

# remove misaligned
x<-lengths(alignIn)
alignIn<-alignIn[x==max(x)]
x<-NULL

DistMat<-NULL
iQGLO<-grep("QGLO",names(alignIn))
iSamples<-grep("Sample.*",names(alignIn))
iRefs<-grep(".*",names(alignIn[-c(iSamples,iQGLO)]))

# Set up clusters
Clusters<-makeCluster(detectCores()[1]-Csave)
registerDoParallel(Clusters)

# Parallel for loop for calculating substitutions
DistMat<-foreach(r=1:length(iRefs), .combine=rbind) %dopar% {
  library(ape) # Must be loaded locally in each thread in order for dist.dna to work
  #Load data
  alignment <- alignIn
  dMat<-NULL
  
  # Index patient samples from this study
  QGLO<-iQGLO
  Samples<-iSamples
  Refs<-iRefs
  
  for(p in 1:length(Samples)){
    Dist<-dist.dna(as.DNAbin(alignment[c(Samples[p],Refs[r])]),model="N",pairwise.deletion=T)[[1]]

    dMat<-rbind(dMat,c(names(alignment[Samples[p]]),
                       names(alignment[Refs[r]]),
                       Dist))
  }
  return(dMat)
}

write.csv(DistMat[order(DistMat[,1]),],"./intermediatedata/DistMat.csv")



