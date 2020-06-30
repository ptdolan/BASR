library(data.table)
library(seqinr)
library(ggplot2)

hits<-fread("~/Documents/GitHub/BASR/JackHMMER.txt")
seqs<-read.fasta("~/Documents/GitHub/BASR/JackHMMER_sequences_all_final.fa")

ggplot(hits)+
  geom_histogram(aes(`E-value`,fill=Description),position='stack')+scale_x_log10()

ggplot(hits)+
  geom_histogram(aes(`E-value`))+scale_x_log10()

ggplot(hits)+
  geom_histogram(aes(as.numeric(`Domain Bit-score`)))+scale_x_log10()

ggplot(hits)+
  geom_point(aes(as.numeric(`Domain Bit-score`),Score))+scale_x_log10()+scale_y_log10()

ggplot(hits)+
  geom_segment(aes(col=`Domain Bit-score`,x=`Query Ali. Start`, xend=`Query Ali. End`, y=`Target Ali. Start`, yend=`Target Ali. End`),alpha=0.2)+coord_fixed()
