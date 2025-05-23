getwd()
# DNase -- it contain default data of ELISA,, which show the different enzymetic activity of enzymes
# it is default set of data which is already stored  in the R 
### the dnase sample contain- each concentration has 2 readings and each reads has 11 runs

#experiments has many repeats these repeats are called Runs
d<-DNase
# scatter plot 
#plot (variable name x axis $column name,variable name for y axis $ column name)
plot(d$conc,d$density)
###
## more featured plot
plot(d$conc,d$density,
     xlab = "concentration (ng/ml)",
     ylab = "absorbance",
     pch = 3,# for the cross sign ,, see the R documents
     col ="blue",# 
     cex = 2,# size of symbols
     cex.axis = 1.5,#size of fonts(numerical vaue) on the axis
     cex.lab = 1.5)# labels fonts size
## calculation of mean
### the dnase sample contain- each concentration has 2 readings and each reads has 11 runs
 
d.avg<-aggregate(density ~ conc,d,mean)
d.avg# it has mean of the concentration values 
##plotting
plot(d.avg$conc,d.avg$density,
     xlab = "concentration (ng/ml)",
     ylab = "absorbance",
     pch = 1,# for the cross sign ,, see the R documents
     col ="blue",# 
     cex = 1,# size of symbols
     cex.axis = 1.5,#size of fonts(numerical vaue) on the axis
     cex.lab = 1.5,
     xlim = c(0,14),## for correcting the scale value of the plot
     ylim = c(0,2))
########################
row.mean <-rowMeans(d[ ,2:3])
row.mean
### plotting
plot(d$conc,row.mean,
     xlab = "concentration (ng/ml)",
     ylab = "absorbance",
     pch = 1,# for the cross sign ,, see the R documents
     col ="blue",# 
     cex = 1,# size of symbols
     cex.axis = 1.5,#size of fonts(numerical vaue) on the axis
     cex.lab = 1.5,
     xlim = c(0,14),## for correcting the scale value of the plot
     ylim = c(0,2))
################################
plot(d.avg$conc,d.avg$density,
     xlab = "concentration (ng/ml)",
     ylab = "absorbance",
     pch = 1,# for the cross sign ,, see the R documents
     col ="blue",# 
     cex = 1,# size of symbols
     cex.axis = 1.5,#size of fonts(numerical vaue) on the axis
     cex.lab = 1.5,
     xlim = c(0,14),## for correcting the scale value of the plot
     ylim = c(0,2),
     type = "b",
     lwd = 2)

#####################
################################################
