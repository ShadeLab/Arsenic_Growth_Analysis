library(grofit)
library(reshape2)
library(dplyr)
library(ggplot2)

#read in raw data
data=data.frame(read.table("20160519_AsIII_I16,46,48.asc.txt"), row.names=TRUE)


#extract the time information, remove the extraneous "s"
time=row.names(data)
time=as.numeric(sub("s", "",time))
row.names(data)=NULL

#define negative controls and subtract from dataset
neg0=data.frame(data[,1])
neg1=data.frame(data[,2])
neg3=data.frame(data[,3])
neg5=data.frame(data[,4])
neg7=data.frame(data[,5])
neg10=data.frame(data[,6])
neg14=data.frame(data[,7])
neg20=data.frame(data[,8])
neg25=data.frame(data[,9])
neg=data.frame(cbind(neg0,neg1,neg3,neg5,neg7,neg10,neg14,neg20,neg25))
negA=data.frame(cbind(neg, rep(neg0,3)))
negB=data.frame(cbind(neg, rep(neg1,3)))
negC=data.frame(cbind(neg, rep(neg3,3)))
negD=data.frame(cbind(neg, rep(neg5,3)))
negE=data.frame(cbind(neg, rep(neg7,3)))
negF=data.frame(cbind(neg, rep(neg10,3)))
negG=data.frame(cbind(neg,neg14,neg25,neg14))
negH=data.frame(cbind(neg,neg20,neg25,neg20))
neg=data.frame(cbind(negA,negB,negC,negD,negE,negF,negG,negH))
data= data - neg

#make a time matrix based on number of observations and variables
n=nrow(data)
m=ncol(data)
time=matrix(rep(time, m), c(n, m))
time=t(time)
time=data.frame(time)


#transpose the data
data=t(data)


#read in platemap
platemap=read.csv("20160519_Platemap_MIC_Grofit.csv")

#combine data with platemap and remove well column
data=cbind(platemap, data)

#remove well column as it is no longer necessary
data=data[,-1]


#make control file to incorporate grofit options
control=grofit.control(nboot.gc=100, parameter = 31, nboot.dr=100)

#OPTION: only run table 
results=gcFit(time, data, control)

#make a data frame of the results
results=data.frame(summary(results))

#save results table
write.csv(results, "20160520_results")

#Find EC50
EC50=drFit(results, control)
EC50=data.frame(summary(EC50))

#save EC50 results
write.csv(EC50, "20160520_EC50")

#extract maximum growth parameter from dataset
mu=data.frame(results$mu.spline, results$reliability)

#add information to mu
mu=cbind(platemap, mu)

#add replicate number to data
reps=data.frame(read.table("20160520_replicates.csv", header=TRUE))
mu=cbind(mu, reps)

#remove control values with replicate numbers
normalization=subset(mu, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
mu=inner_join(mu, normalization, by = c("Strain", "Replicate"))

#normalize data
norm=mu$results.mu.spline.x/mu$results.mu.spline.y
norm=data.frame(norm)
#add information to normalized data
mu=cbind(mu, norm)

#group data
grouped=group_by(mu, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability.x=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~Strain)

#Extract maximum OD590 and reliability data
A=data.frame(results$A.spline, results$reliability)
A=cbind(platemap, A)

#add replicate number to data
reps=data.frame(read.table("20160520_replicates.csv", header=TRUE))
A=cbind(A, reps)

#remove control values with replicate numbers
normalization=subset(A, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
A=inner_join(A, normalization, by = c("Strain", "Replicate"))

#normalize data
norm=A$results.A.spline.x/A$results.A.spline.y

#add information to normalized data
A=cbind(A, norm)

#group data
grouped=group_by(A, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability.x=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~Strain) 

#extract time to exponential growth parameter (lambda) from dataset
lambda=data.frame(results$lambda.spline, results$reliability)

#add information to mu
lambda=cbind(platemap, lambda)

#add replicate number to data
reps=data.frame(read.table("20160520_replicates.csv", header=TRUE))
lambda=cbind(lambda, reps)

#remove control values with replicate numbers
normalization=subset(lambda, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
lambda=inner_join(lambda, normalization, by = c("Strain", "Replicate"))

#normalize data
norm=lambda$results.lambda.spline.x/lambda$results.lambda.spline.y

#add information to normalized data
lambda=cbind(lambda, norm)

#group data
grouped=group_by(lambda, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability.x=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1.5) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~Strain) +
  ylim(0, 5)

#extract integral growth parameter from dataset
int=data.frame(results$integral.spline, results$reliability)

#add information to int
int=cbind(platemap, int)

#add replicate number to data
reps=data.frame(read.table("20160520_replicates.csv", header=TRUE))
int=cbind(int, reps)

#remove control values with replicate numbers
normalization=subset(int, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
int=inner_join(int, normalization, by = c("Strain", "Replicate"))

#normalize data
norm=int$results.integral.spline.x/int$results.integral.spline.y
norm=data.frame(norm)

#add information to normalized data
int=cbind(int, norm)

#group data
grouped=group_by(int, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability.x=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~Strain) + 
  ylim(0,1)

