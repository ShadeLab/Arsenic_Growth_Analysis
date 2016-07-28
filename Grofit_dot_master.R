setwd("/Users/dunivint/Documents/GitHubRepos/Arsenic_Growth_Analysis/")

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

#Average data for each variable before melting
data=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu=mean(mu.spline),
            A=mean(A.spline), 
            lambda=mean(lambda.spline), 
            integral=mean(integral.spline),
            growth=mean(growth),
            mu.sd=sd(mu.spline),
            A.sd=sd(A.spline),
            lambda.sd=sd(lambda.spline),
            integral.sd=sd(integral.spline))

#Call samples with 1/3 growing no growth
data$growth[data$growth<0.4]=0

#Call parameters NA for no growth
data$A[data$growth==0]=NA
data$mu[data$growth==0]=NA
data$integral[data$growth==0]=NA
data$lambda[data$growth==0]=NA

#remove negative control data
data=data[!(data$TestId=="NEG"),]

#adjust all negative lambda values to 1000 (just below lowest positive lambda)
data$lambda[data$lambda<0]=1000

#normalize data to growth without arsenic
normalization=subset(data, concentration %in% 0)
normalization=normalization[,-3]
data=inner_join(data, normalization, by = c("TestId", "AddId"))
full.data=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu=mu.x/mu.y,
            lambda=lambda.x/lambda.y,
            mu.sd=mu.sd.x/mu.y,
            lambda.sd=lambda.sd.x/lambda.y)

#remove rows with no arsenic (controls)
data=full.data[!(full.data$concentration=="0"),]

#remove growth and integral because we do not need them for analysis
data=full.data[,!colnames(full.data) %in% grep("gr", colnames(full.data), value=TRUE)]

#read in genus information
genus=read.csv("isolate_genus.csv")

#merge genus with data
data=inner_join(genus, data, by="TestId")

#separate arsenate and arsenite
as5=data[-which(data$AddId=="III"),]
as3=data[-which(data$AddId=="V"),]

#remove concentrations 20 and 25mM from arsenite because no isolates grow
as3=as3[!as3$concentration %in% grep("2", as3$concentration, value=TRUE),]
data=data[!data$concentration %in% grep("25", data$concentration, value=TRUE),]
data=data[!data$concentration %in% grep("20 ", data$concentration, value=TRUE),]

as3.long=melt(as3, id=c("TestId", "Genus", "concentration", "AddId"), measure=c("mu", "lambda"))
as5.long=melt(as5, id=c("TestId", "Genus", "concentration", "AddId"), measure=c("mu", "lambda"))
data.long=melt(data, id=c("TestId", "Genus", "concentration", "AddId"), measure=c("mu", "lambda"))


#plot & save arsenite graph
setEPS()
postscript("arsenite.grofit.eps")
ggplot(data=as3.long, aes(x=concentration, y=value, color=Genus, group_by(c("lambda", "mu", "A")))) +
  geom_point(size=1.5) +
  facet_grid(variable~Genus, scales="free_y") +
  labs(x="Concentration (mM)", y="Average") + 
  theme_bw()
dev.off()

#plot & save arsenate graph
setEPS()
postscript("arsenate.grofit.eps")
ggplot(data=as5.long, aes(x=concentration, y=value, color=Genus, group_by(c("lambda", "mu", "A")))) +
  geom_point(size=1.5) +
  facet_grid(variable~Genus, scales="free_y") +
  labs(x="Concentration (mM)", y="Average") + 
  theme_bw()
dev.off()


setEPS()
postscript("arsenate.grofit.eps")
ggplot(data=data.long, aes(x=concentration, y=value, color=Genus, group_by(c("lambda", "mu", "A")))) +
  geom_point(size=1.5) +
  facet_grid(AddId+variable~Genus, scales="free_y") +
  labs(x="Concentration (mM)", y="Average") + 
  theme_bw()
dev.off()