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

ggplot(hits[FILTER,])+
  geom_segment(aes(col=Score,x=`Target Ali. Start`, xend=`Target Ali. End`, y=`Target Name`, yend=`Target Name`),alpha=1)+scale_color_viridis_c()+ggpubr::theme_pubr()+theme(axis.text.y = element_text(size = 2))


listcompare<-function(X)
{
  return(length(unique(X[,unique(`Target Ali. Start`:`Target Ali. End`)]))/X$`Query Length`)
}

hits[,`targetCoverage`:=listcompare(.SD),by=`Target Name`]

#ggplot(hits[targetCover>450&targetCoverage>.90&`Target Length`<600])+
#  geom_segment(aes(col=Score,x=`Query Ali. Start`, xend=`Query Ali. End`, y=`Target Name`, yend=`Target Name`),alpha=0.2)+scale_color_viridis_c()
FILTER<-with(hits,targetCoverage>.8&`Target Length`)
table(FILTER)
Before<-ggplot(hits)+
  stat_bin(aes(fill=as.factor(round(targetCoverage,digits = 1)),x=Score),position = 'stack')+
  scale_fill_viridis_d()+
  theme_bw()+ggtitle("Before Filter")+
  labs(fill = "Query Coverage")

After<-ggplot(hits[FILTER])+
  stat_bin(aes(fill=as.factor(round(targetCoverage,digits = 1)),x=Score),position = 'stack')+
  scale_fill_viridis_d(begin = 0.6)+
  theme_bw()+ggtitle("After Filter")+
  labs(fill = "Query Coverage")
cowplot::plot_grid(Before,After)

seqs<-read.fasta("~/Documents/GitHub/BASR/JackHMMER_sequences_all_final.fa")
seqList={}

for(i in seqs){
  if(attr(i,"name")%in%hits$`Target Name`[FILTER]){
    print("match")
    seqList<-append(seqList,i)
    }
}

write.csv(file="~/Documents/GitHub/BASR/Data/Filtered_JackHmmer.csv",hits[FILTER],quote = F,row.names = F)


