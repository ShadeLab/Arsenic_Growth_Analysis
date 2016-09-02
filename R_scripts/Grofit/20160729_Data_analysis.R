library(grofit)
library(reshape2)
library(dplyr)
library(ggplot2)

#read in raw data
data=data.frame(read.table("20160729_AsV_AsIII_MIC_I2742.asc"), row.names=TRUE)

#extract the time information, remove the extraneous "s"
time=row.names(data)
time=as.numeric(sub("s", "",time))
row.names(data)=NULL

#define negative controls and subtract from dataset
neg0=data.frame(data[,5])
neg0=data.frame(rep(neg0,5))
neg1=data.frame(data[,17])
neg1=data.frame(rep(neg1,5))
neg3=data.frame(data[,29])
neg3=data.frame(rep(neg3,5))
neg5=data.frame(data[,41])
neg5=data.frame(rep(neg5,5))
neg7=data.frame(data[,53])
neg7=data.frame(rep(neg7,5))
neg10=data.frame(data[,65])
neg10=data.frame(rep(neg10,5))
neg14=data.frame(data[,77])
neg14=data.frame(rep(neg14,5))
neg20=data.frame(data[,89])
neg20=data.frame(rep(neg20,5))
vneg0=data.frame(data[,10])
vneg0=data.frame(rep(vneg0,5))
vneg10=data.frame(data[,22])
vneg10=data.frame(rep(vneg10,5))
vneg50=data.frame(data[,34])
vneg50=data.frame(rep(vneg50,5))
vneg100=data.frame(data[,46])
vneg100=data.frame(rep(vneg100,5))
vneg150=data.frame(data[,58])
vneg150=data.frame(rep(vneg150,5))
vneg200=data.frame(data[,70])
vneg200=data.frame(rep(vneg200,5))
vneg250=data.frame(data[,82])
vneg250=data.frame(rep(vneg250,5))
vneg300=data.frame(data[,94])
vneg300=data.frame(rep(vneg300,5))
x=data.frame(data[,11:12])
neg=data.frame(cbind(neg0,vneg0,x,neg1,vneg10,x,neg3,vneg50,x,neg5,vneg100,x,neg7,vneg150,x,neg10,vneg200,x,neg14,vneg250,x,neg20,vneg300,x))
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
platemap=read.csv("20160729_platemap.csv")

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
write.csv(results, "20160729_results.csv")

#extract maximum growth parameter from dataset
mu=data.frame(results$mu.spline, results$reliability)

#read reduced platemap
platemap=data.frame(read.csv("20160729_platemap.csv"))

#reduce platemap to match reduced data then join
platemap=cbind(data3, platemap)
platemap=platemap[!(platemap$apply.data..2..function.x..all.with.data..ave..sum.x....0.....=="FALSE"),]
platemap=platemap[,-1]
mu=cbind(platemap,mu)

#remove control values with replicate numbers
normalization=subset(mu, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-3]
normalization=normalization[,-4]

#summarize controls
norm=normalization %>%
  group_by(Arsenic, Strain) %>%
  summarise(N=length(results.mu.spline),
            Average=mean(results.mu.spline))

#add negatives back
mu=inner_join(mu, norm, by = c("Arsenic", "Strain"))

#normalize data
norm=mu$results.mu.spline/mu$Average
norm=data.frame(norm)

#add information to normalized data
mu=cbind(mu, norm)

#group data
grouped=group_by(mu, Arsenic, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average, group_by(Arsenic.x, Strain, Concentration))) +
  geom_line(size=1) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(Strain~Arsenic, scales="free_x")

#Extract maximum OD590 and reliability data
A=data.frame(results$A.spline, results$reliability)
A=cbind(platemap, A)

#remove control values with replicate numbers
normalization=subset(A, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-3]
normalization=normalization[,-4]

#summarize controls
norm=normalization %>%
  group_by(Arsenic, Strain) %>%
  summarise(N=length(results.A.spline),
            Average=mean(results.A.spline))

#add negatives back
A=inner_join(A, norm, by = c("Strain", "Arsenic"))

#normalize data
norm=A$results.A.spline/A$Average

#add information to normalized data
A=cbind(A, norm)

#group data
grouped=group_by(A, Arsenic, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(Strain~Arsenic, scales="free_x") 

#extract time to exponential growth parameter (lambda) from dataset
lambda=data.frame(results$lambda.spline, results$reliability)

#add information to mu
lambda=cbind(platemap, lambda)

#remove control values with replicate numbers
normalization=subset(lambda, Concentration %in% 0)

#remove unimportant columns
normalization=normalization[,-1]
normalization=normalization[,-3]
normalization=normalization[,-4]

#summarize controls
norm=normalization %>%
  group_by(Arsenic, Strain) %>%
  summarise(N=length(results.lambda.spline),
            Average=mean(results.lambda.spline))

#add negatives back
lambda=inner_join(lambda, normalization, by = c("Strain", "Arsenic"))

#normalize data
norm=lambda$results.lambda.spline.x/lambda$results.lambda.spline.y

#add information to normalized data
lambda=cbind(lambda, norm)

#group data
grouped=group_by(lambda, Arsenic, Strain, Concentration)

#remove tagged(bad) wells from data analysis
grouped2=grouped[!(grouped$results.reliability=="FALSE"),]

#Calculate average and stdev
stats=summarise(grouped2, Average=mean(norm), StDev=sd(norm))

#plot data
ggplot(data=stats, aes(x=Concentration, y=Average)) +
  geom_line(size=1.5) +
  geom_point(aes(color=Concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=stats$Average + stats$StDev, ymin=stats$Average - stats$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(Strain~Arsenic, scales="free_x")
