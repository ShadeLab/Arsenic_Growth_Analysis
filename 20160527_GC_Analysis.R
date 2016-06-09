library(grofit)
library(reshape2)
library(dplyr)
library(ggplot2)

#read in raw data
data=data.frame(read.table("20160527_AsIII_I23,45,59,A05,08,16,27.asc.txt"), row.names=TRUE)


#extract the time information, remove the extraneous "s"
time=row.names(data)
time=as.numeric(sub("s", "",time))
row.names(data)=NULL

#define negative controls and subtract from dataset
neg0=data.frame(data[,58])
aneg0=data.frame(rep(neg0,12))
neg0=data.frame(rep(neg0,10))
neg1=data.frame(data[,70])
aneg1=data.frame(rep(neg1,12))
neg1=data.frame(rep(neg1,10))
neg3=data.frame(data[,82])
aneg3=data.frame(rep(neg3,12))
neg3=data.frame(rep(neg3, 10))
neg5=data.frame(data[,94])
aneg5=data.frame(rep(neg5,12))
neg5=data.frame(rep(neg5,10))
neg7=data.frame(data[,59])
neg7=data.frame(rep(neg7,2))
neg=data.frame(cbind(aneg0,aneg1,aneg3,aneg5,neg0,neg7,neg1,neg7,neg3,neg7,neg5,neg7))
data= data - neg

#remove columns (wells) with less than 5 data points >=0 because GroFit will not read these
data2=data[apply(data,2,function(x)all(with(data, ave((sum(x>=0)) >5))))]
data3=data.frame(apply(data,2,function(x)all(with(data, ave((sum(x>=0)) >5)))))

#make a time matrix based on number of observations and variables
n=nrow(data2)
m=ncol(data2)
time=matrix(rep(time, m), c(n, m))
time=t(time)
time=data.frame(time)


#transpose the data
data=t(data)


#read in platemap
platemap=read.csv("20160527_platemap.csv")

#combine data with platemap and remove well column
data=cbind(platemap, data)

#re-remove all wells with less than 5 data points >=0 because GroFit will not read these
data=cbind(data3,data)
data=data[!(data$apply.data..2..function.x..all.with.data..ave..sum.x....0.....=="FALSE"),]
data=data[,-c(1,2)]

#make control file to incorporate grofit options
control=grofit.control(nboot.gc=100, parameter = 31, nboot.dr=100)

#OPTION: only run table 
results=gcFit(time, data, control)

#make a data frame of the results
results=data.frame(summary(results))

#save results table
write.csv(results, "20160527_results")

#Find EC50
EC50=drFit(results, control)
EC50=data.frame(summary(EC50))

#save EC50 results
write.csv(EC50, "20160527_EC50")

#extract maximum growth parameter from dataset
mu=data.frame(results$mu.spline, results$reliability)

#read reduced platemap
platemap=data.frame(read.csv("20160527_platemap_grofit.csv"))
mu=cbind(platemap,mu)

#add replicate number to data
reps=data.frame(read.table("20160527_replicates.csv", header=TRUE))
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
reps=data.frame(read.table("20160527_replicates.csv", header=TRUE))
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
reps=data.frame(read.table("20160527_replicates.csv", header=TRUE))
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
  facet_wrap(~Strain)

#extract integral growth parameter from dataset
int=data.frame(results$integral.spline, results$reliability)

#add information to int
int=cbind(platemap, int)

#add replicate number to data
reps=data.frame(read.table("20160527_replicates.csv", header=TRUE))
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
  facet_wrap(~Strain)


