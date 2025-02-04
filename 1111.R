
DN=ATGACGTGACGTGCAGTGACGTGACGTGA
print(DN) # print the sequence
nchar(DN) # count the nucleotide in seq
q=substring(DN
          ,1,620) # extract the sequence between 1-620 and assign to variable
s= strsplit(q,split=NULL) # it split the q sequence and stor in s variable
print(s) # it print the output (splitted vactors) 
p=s[[1]]
#--------count of each base type
t_count=sum(p=="t")
a_count=sum(p=="a")
g_count=sum(p=="g")
c_count=sum(p=="c")
Total=c(t_count, a_count, g_count, c_count)
count= paste('the total of T,A,G and C is:',total)
#--------percentage count
per_t= t_count/620*100
Per_a= a_count/620*100
Per_g= g_count/620*100
per_c= c_count/620*100
percent=c(per_t,Per_a,Per_g,per_c)
print(count)
print(percent)
#---------------barplot of ATGC %----------------
base=c('T','A','G','C')
barplot(percent,xlab="Nucelotide",ylab="percentage",main="Base percentage",
        names.arg=base,col= c("Green","Red","Blue","Yellow"))
#----------formation of table through data frame-----------
df=data.frame(base,Total,percent)
print(df)