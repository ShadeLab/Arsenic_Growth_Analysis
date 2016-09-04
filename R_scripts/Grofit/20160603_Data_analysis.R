library(grofit)
library(reshape2)
library(dplyr)
library(ggplot2)

#read in raw data
data=data.frame(read.table("20160603_AsV_A12,31,I48,16.asc"), row.names=TRUE)

#extract the time information, remove the extraneous "s"
time=row.names(data)
time=as.numeric(sub("s", "",time))
row.names(data)=NULL

#define negative controls and subtract from dataset
neg0=data.frame(data[,1])
neg0=data.frame(rep(neg0,12))
neg10=data.frame(data[,13])
neg10=data.frame(rep(neg10,12))
neg50=data.frame(data[,25])
neg50=data.frame(rep(neg50,12))
neg100=data.frame(data[,37])
neg100=data.frame(rep(neg100,12))
neg150=data.frame(data[,49])
neg150=data.frame(rep(neg150,12))
neg200=data.frame(data[,61])
neg200=data.frame(rep(neg200,10))
neg250=data.frame(data[,73])
neg250=data.frame(rep(neg250,10))
neg300=data.frame(data[,85])
neg300=data.frame(rep(neg300,10))
neg=data.frame(cbind(neg0,neg10,neg50,neg100,neg150,neg200,data[,1],data[,37],neg250,data[,13],data[,49],neg300,data[,25],data[,1]))
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
platemap=read.csv("20160603_platemap.csv")

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
write.csv(results, "20160603_results")


#extract maximum growth parameter from dataset
mu=data.frame(results$mu.spline, results$reliability)

#read reduced platemap
platemap=data.frame(read.csv("20160603_platemap.csv"))

#reduce platemap to match reduced data then join
platemap=cbind(data3, platemap)
platemap=platemap[!(platemap$apply.data..2..function.x..all.with.data..ave..sum.x....0.....=="FALSE"),]
platemap=platemap[,-1]
mu=cbind(platemap,mu)

#add replicate number to data
reps=data.frame(read.table("20160603_replicates.csv", header=TRUE))
reps=cbind(data3, reps)
reps=reps[!(reps$apply.data..2..function.x..all.with.data..ave..sum.x....0.....=="FALSE"),]
reps=reps[,-1]
mu=cbind(mu, reps)

#remove control values with replicate numbers
normalization=subset(mu, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
mu=inner_join(mu, normalization, by = c("Strain", "reps"))

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
A=cbind(A, reps)

#remove control values with replicate numbers
normalization=subset(A, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
A=inner_join(A, normalization, by = c("Strain", "reps"))

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
lambda=cbind(lambda, reps)

#remove control values with replicate numbers
normalization=subset(lambda, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
lambda=inner_join(lambda, normalization, by = c("Strain", "reps"))

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
int=cbind(int, reps)

#remove control values with replicate numbers
normalization=subset(int, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-2]
normalization=normalization[,-2]

#add negatives back
int=inner_join(int, normalization, by = c("Strain", "reps"))

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


