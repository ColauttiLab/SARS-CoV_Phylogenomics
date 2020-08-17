library(tidyr)
library(dplyr)
library(ggplot2)

polymorph<-read.csv("intermediatedata/polymorphicAlignment.csv", stringsAsFactors = F)
colnames(polymorph) <- gsub("X","",colnames(polymorph))
colnames(polymorph)[1] <- "sample"
polymorph$sample[160] <- "Wuhan-Hu-1/2019" #Change index number based on which file
polymorph$sample[161] <- "Wuhan/WH04/2020" #Change index number based on which file


for (x in 2:ncol(polymorph)) {
  refNt <- polymorph[161,x] #Change index number
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


#long$position <- factor(long$position,levels=c(unique(long$position)))
long$position<-as.numeric(long$position)
long$sample <- factor(long$sample,levels=polymorph$sample[order(suborder)])

palette <- c("#007FFF","#99001C","#008000","#fffc00")

VarPlot<-ggplot(long, aes(x=position, y=sample)) + 
  geom_tile(aes(colour=substitution, fill=substitution),size=1.2) + 
  scale_fill_manual(values=palette) +
  scale_colour_manual(values=palette) +
  theme_classic() 

VarPlot

pdf(file="outputs/fig2.pdf",width=12,height=10)
  VarPlot
dev.off()

