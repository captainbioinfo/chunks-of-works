u<-c(2,4,6,67,7,74,43)
b<-seq(0,100,5)
print(b)
##getting the number using index number
print(b[4])
# addition of vector
c<-u + 2
c
f<-u+c
f
### multipliation of two vector by dot product
library(geometry)#after installing the package
## calculate the dot product of the u and b
g<-dot(c,f)#like in kinemaatics
g
d<-(2:10)
d
#matrices
j<-matrix(d,nrow=3,ncol=3,byrow=F)
j#print=ctrl +enter
##bind vectors to form a matrix
l<-c(1:10)
m<-c(91:100)
lm<-rbind(l,m)#used to bind the two vector into a matrix
lm
l 
m


####data frame into a matrix
#get
df<-data.frame(1:20)
print(df)
class(df)
j<-data.matrix(df)
print(j)
class(j)
colnames(j)##to get the column name of the matrix



###to print the specific value of the column or row
jt<-matrix(1:10,nrow=3,ncol=3,byrow=F)
print(jt[3, 2])# 3 is the row no. and 2 is the col no. 
print(jt[1,])# to print all the value
#### to transpose a matrix
jl<-t(jt)
#after transpose 
jl
###matrix addition
sum.jlt<-jl + jt
### matrix multiplication
matmulti<-jl %*%jt### to matrix matrix multiplication
print(matmulti)
matmul<-jl %%jt ##((lesslikely too be used)) for element element wise multiplication
print(matmul)
## 








