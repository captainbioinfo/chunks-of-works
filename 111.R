d=c("ATGTCGCTTTCTAACAAGCTGACGCTGGACAAGCTGGACGTTAAAGGGAAGCGGGTCGTTATGAGAGTCGT")
m=strsplit(d, '')
print(m)
nchar(d)
l=substr(d,start=4,stop=9)
print(l)
v=rep(c('atgc'),times=20)
print(v)
grepl("t",v)

getwd()

library(seqinr)
f=open.fasta(file="se.fasta",as.string=TRUE)




