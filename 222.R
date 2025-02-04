library(seqinr)
n=read.fasta(file"se.fasta",as.string=TRUE)








library(seqinr)
e=c(2,3,4,5.6)
base("A","T","G","C")
barplot(e,xlab="Nucelotide",ylab="percentage",main="base",
        names.arg=base,col= c("Green","Red","Blue","Yellow"))