library(tidyr)
library(dplyr)
library(ggplot2)

polymorph<-read.csv("./intermediatedata/PolyAlignRef.csv", stringsAsFactors = F)
colnames(polymorph) <- gsub("X","",colnames(polymorph))
colnames(polymorph)[1] <- "sample"

# Create variant map (check against Wuhan reference to only ID variant sites)
for (x in 2:ncol(polymorph)) {
  refNt <- polymorph[grep("Wuhan",polymorph$sample),x] #Change index number
  polymorph[which(polymorph[,x] == refNt),x]<-" "
}

suborder <- integer()
for (y in 1:nrow(polymorph)) {
  nsub<-length(which(polymorph[y,]==" "))
  suborder<-c(suborder,nsub)
}

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

pdf(file="./outputs/figures/fig1Pos.pdf",width=12,height=10)
  VarPlot
dev.off()

clong<-long
clong$position<-as.numeric(paste(clong$position))

# Variant plot with y-axis scaled to genome
VarPlotScale<-ggplot(clong, aes(x=position, y=sample),  fill=substitution) + 
  geom_tile(aes(fill=substitution, width = 100), colour = "white") + 
  scale_fill_manual(name = "Nucleotide", breaks = c('A','G','C','T','other'),values=palette) +
  scale_x_continuous(breaks=seq(0,30000,5000)) +
  xlab("Position") +
  ylab("Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(size=8))

pdf(file="./outputs/figures/fig1Scale.pdf",width=12,height=10)
 VarPlotScale
dev.off()
