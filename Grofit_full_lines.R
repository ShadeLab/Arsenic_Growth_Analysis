library(ggplot2)
library(dplyr)
library(reshape2)
library(stats)

#read in data by calling on files with "results" 
#in the name, and combine all data from separate experiments
data=data.frame(do.call(rbind, lapply(list.files(pattern="*results"), read.csv)))

#remove unnecessary first column and rename a data file. 
#we will need this structure for EC50 later
data=data[-1]
orig.data=data

#fix problem spline with models. (I did this manually by saving
#the data as a .csv and adding a column "quality" that describes
#whether or not I need to use the model)
write.csv(orig.data, "20160610_data.csv")
data=read.csv("20160610_data_model.csv")
data$mu.spline <- with( data, ifelse( Quality == 0, mu.model, mu.spline ))
data$A.spline <- with( data, ifelse( Quality == 0, A.model, A.spline ))
data$lambda.spline <- with( data, ifelse( Quality == 0, lambda.model, lambda.spline ))
data$integral.spline <- with( data, ifelse( Quality == 0, integral.model, integral.spline ))

#remove all model data since we are using spline
data=data[,-grep("model", names(data))]
data=data[,-grep(".bt", names(data))]
data=data[,-grep("log.", names(data))]
data=data[,-grep("nboot", names(data))]
data=data[,-grep("Quality", names(data))]

#add binary column to describe yes/no growth
data$growth=as.integer(data$reliability=="TRUE")

#remove unnecessary first column
data=data[,-1]

#Average data for each variable before melting
stats=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu=mean(mu.spline),
            A=mean(A.spline), 
            lambda=mean(lambda.spline), 
            integral=mean(integral.spline),
            growth=mean(growth))

#Call samples with 1/3 growing no growth
stats$growth[stats$growth<0.4]=0

#Flag samples with 2/3 growth so we can remove poor replicate
stats$flag[which(stats$growth==0 | stats$growth==1)] =0
stats$flag[is.na(stats$flag)]=1

#remove unnecessary columns from stats
stats=stats[,-c(4:8)]

#combine flagged sites with original data
data=inner_join(data, stats, by=c("TestId", "AddId", "concentration"), copy=TRUE)

#remove rows that are called false AND are flagged
data=data[-which(data$reliability=="FALSE" & data$flag==1),]

#remove flag and reliability columns as they are no longer necessary 
data=data[,!colnames(data) %in% grep("reliability",colnames(data), value=TRUE)]
data=data[,!colnames(data) %in% grep("flag",colnames(data), value=TRUE)]

#remove negative control data
data=data[!(data$TestId=="NEG"),]

#adjust all negative lambda values to 1000 (just below lowest positive lambda)
#this only affect acinetobacter (A2705, A2716, I2759)
data$lambda.spline[data$lambda.spline<0]=1000

##normalize data to growth without arsenic
#extract data without arsenic
normalization=subset(data, concentration %in% 0)

#remove concentration information (all 0)
normalization=normalization[,-3]

#average data
normalization=normalization %>%
  group_by(TestId, AddId) %>%
  summarise(mu.spline=mean(mu.spline),
           A.spline=mean(A.spline),
           lambda.spline=mean(lambda.spline))

#join normalization data with full dataset
data=inner_join(data, normalization, by = c("TestId", "AddId"))

#normalize data to growth without arsenic
data$mu=data$mu.spline.x/data$mu.spline.y
data$A=data$A.spline.x/data$A.spline.y
data$lambda=data$lambda.spline.x/data$lambda.spline.y

#remove columns with spline (no longer needed)
data=data[,!colnames(data) %in% grep("spline",colnames(data), value=TRUE)]

##calculate means and stdev for all samples
data=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu.sd=sd(mu),
            A.sd=sd(A),
            lambda.sd=sd(lambda),
            mu=mean(mu), 
            A=mean(A), 
            lambda=mean(lambda), 
            growth=mean(growth))

#Call samples with 1/3 growing no growth
data$growth[data$growth<0.4]=0

#Call parameters NA for no growth
data$A[data$growth==0]=NA
data$mu[data$growth==0]=NA
data$lambda[data$growth==0]=NA

#melt data for plotting
data.long=melt(data, id=c("TestId", "AddId", "concentration"), measure=c("mu", "A", "lambda"))

#split arsenate and arsenite
as5=data[-which(data$AddId=="III"),]
as3=data[-which(data$AddId=="V"),] 

#plot data
setEPS()
postscript("arsenate.grofit.mu.eps")
ggplot(data=as5, aes(x=concentration, y=mu, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as5$mu+as5$mu.sd, ymin=as5$mu-as5$mu.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()

setEPS()
postscript("arsenate.grofit.A.eps")
ggplot(data=as5, aes(x=concentration, y=A, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as5$A+as5$A.sd, ymin=as5$A-as5$A.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()

setEPS()
postscript("arsenate.grofit.lambda.eps")
ggplot(data=as5, aes(x=concentration, y=lambda, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as5$lambda+as5$lambda.sd, ymin=as5$lambda-as5$lambda.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()

setEPS()
postscript("arsenite.grofit.mu.eps")
ggplot(data=as3, aes(x=concentration, y=mu, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as3$mu+as3$mu.sd, ymin=as3$mu-as3$mu.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()

setEPS()
postscript("arsenite.grofit.A.eps")
ggplot(data=as3, aes(x=concentration, y=A, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as3$A+as3$A.sd, ymin=as3$A-as3$A.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()

setEPS()
postscript("arsenite.grofit.lambda.eps")
ggplot(data=as3, aes(x=concentration, y=lambda, color=concentration)) +
  geom_line(size=1, color="black") +
  geom_point(aes(color=concentration), size=2) +
  geom_point(shape=1, size=2, color="black") +
  geom_errorbar(aes(ymax=as3$lambda+as3$lambda.sd, ymin=as3$lambda-as3$lambda.sd), color="black") +
  facet_wrap(~TestId) +
  labs(x="Concentration (mM)", y="Average") + 
  scale_color_gradientn(colors=rainbow(6)) +
  theme_bw()
dev.off()
