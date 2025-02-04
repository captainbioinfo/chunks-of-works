dna="ATGATGCGATGC"
length=nchar(dna)
cat(length)
countA=nchar(gsub("[^A]","",dna))
print(countA)

pera=(countA/length)*100
print(pera)
barplot(pera,xlab=pera,ylab=100,col=c("red","yellow"))  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
