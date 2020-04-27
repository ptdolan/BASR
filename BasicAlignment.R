#Gisaid Analysis

library(Biostrings)
library(DECIPHER)

inputSeqs<-readDNAStringSet(c("~/Research Projects/SARS-COV2/PhyloAnaylsis_Hong/gisaid_cov2020_sequences.fasta"))#Vector of FASTA files. 

inputSeqs<-inputSeqs[(width(inputSeqs)>27000)&(!duplicated(inputSeqs))]#set cutoff and remove redudant sequences
names(inputSeqs)<-lapply(X = names(inputSeqs),function(X){paste(strsplit(X," ")[[1]],collapse= "_")})
inputSeqs_Gapless<-RemoveGaps(inputSeqs)
mAlign<-AlignSeqs(inputSeqs_Gapless)
SandMAlign<-StaggerAlignment(mAlign)
writeXStringSet(SandMAlign,filepath = "~/Research Projects/SARS-COV2/PhyloAnaylsis_Hong/Aligned.aln")

system(command = "FastTree -nt -gtr < ~/Research Projects/SARS-COV2/PhyloAnaylsis_Hong/Aligned.aln > GISAID.tree")
