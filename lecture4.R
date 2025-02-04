sys<-c(12,3,6,7)
ages<-c(1,43,7,36)
ab<-c(234,436,567,74)
head(sys)
summary(sys)
length(sys)
min(sys)
max(sys)
median(sys)
mean(sys)
# varience
var(sys) 
sd(sys)#stranded deviation fuction -under root of varience
pp<-data.frame(ages,sys,ab)
summary(pp)
########
data.file<-read.csv("DATA_FSB_SET_1.csv",header=TRUE)
data<-read.csv("DATA_FSB_SET_1.csv",header=T)# inplace of TRUE we can also write T

summary(data.file)
##subset isused to vies or selet a single column of the table


umar<-subset(data.file,Age=="20")
View(umar)
##finding mean of the ages(umar)with the sleep_hour column)

mean.pw.umar<-mean(umar$Sleep_Hours)
## frequence histogram 
hist(umar$Sugar_Blood)
#### can choose the csv file 
f<-read.csv(file.choose())
sys## press contro then enter to view
View(sys) # view or see the data
remove("ds1") # to remove any content
rm("mat")#same as above
