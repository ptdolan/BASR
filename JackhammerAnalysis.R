library(data.table)
library(seqinr)
library(ggplot2)

hits<-fread("~/Documents/GitHub/BASR/JackHMMER.txt")
seqs<-read.fasta("~/Documents/GitHub/BASR/JackHMMER_sequences_all_final.fa")

ggplot(hits)+
  geom_histogram(aes(`E-value`))+scale_x_log10()

ggplot(hits)+
  geom_histogram(aes(as.numeric(`Domain Bit-score`)))+scale_x_log10()

ggplot(hits)+
  geom_point(aes(as.numeric(`Domain Bit-score`),Score))+scale_x_log10()+scale_y_log10()

ggplot(hits)+
  geom_segment(aes(col=Score,x=`Query Ali. Start`, xend=`Query Ali. End`, y=`Target Ali. Start`, yend=`Target Ali. End`),alpha=0.2)+coord_fixed()+scale_color_viridis_c()

ggplot(hits)+
  geom_segment(aes(col=Score,x=`Query Ali. Start`, xend=`Query Ali. End`, y=`E-value`, yend=`E-value`),alpha=0.2)+
  scale_color_viridis_c()+scale_y_log10()

ggplot(hits)+
  geom_point(aes(`Query Ali. Start`-`Query Ali. End`,Score,color=Acc))+scale_color_viridis_c()


ggplot(hits[hits$`Target Name` %in% names(table(hits$`Target Name`)[table(hits$`Target Name`)>1]) ])+
  geom_point(aes(`Query Ali. Start`-`Query Ali. End`,Score,color=Acc))+scale_color_viridis_c()

ggplot(hits[hits$`Target Name` %in% names(table(hits$`Target Name`)[table(hits$`Target Name`)>1]) ])+
  geom_segment(aes(col=Score,x=`Query Ali. Start`, xend=`Query Ali. End`, y=`Target Name`, yend=`Target Name`),alpha=0.2)+coord_fixed()+scale_color_viridis_c()

ggplot(hits[targetLength>450&targetCoverage>.90,`Target Name`])+
  geom_segment(aes(col=Score,x=`Target Ali. Start`, xend=`Target Ali. End`, y=`Target Name`, yend=`Target Name`),alpha=0.2)+scale_color_viridis_c()

hits[,hitCover:=`Query Ali. End`-`Query Ali. Start`]
hits[,targetCover:=`Target Ali. End`-`Target Ali. Start`]
hits[,targetCoverage:=sum(targetCover)/`Target Length`,by=`Target Name`]





