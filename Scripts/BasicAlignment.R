#Basic Alignment.R
#Alignment of fa files. Usage (in command line): `Rscript BasicAlignment.R <myfasta.fa> [<other optional .fa's ... >]
library(Biostrings)
library(DECIPHER)
inputSeqs<-readDNAStringSet(commandArgs(trailingOnly=TRUE))#Vector of FASTA files. 

inputSeqs<-inputSeqs[]#Filter by some criteria. i
names(inputSeqs)<-lapply(X = names(inputSeqs),function(X){paste(strsplit(X," ")[[1]],collapse= "_")}) #avoids some name parsing issues with other programs. 
inputSeqs_Gapless<-RemoveGaps(inputSeqs)
MAlign<-AlignSeqs(processors = 8,inputSeqs_Gapless,refinements=10)
SMAlign<-StaggerAlignment(MAlign)
writeXStringSet(SMAlign,filepath = "basicAln_output.aln")

#system(command = "FastTree -nt -gtr < ./basicAln_output.aln > basicAln.tree")

