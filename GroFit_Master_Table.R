library(ggplot2)
library(grofit)
library(dplyr)
library(reshape2)
library(data.table)
library(vegan)
library(gplots)

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

#remove all model data since we are using spline
data=data[,-grep("model", names(data))]
data=data[,-grep(".bt", names(data))]
data=data[,-grep("log.", names(data))]
data=data[,-grep("nboot", names(data))]
data=data[,-grep("Quality", names(data))]

#add binary column to describe yes/no growth
data$growth=as.integer(data$reliability=="TRUE")

#remove unnecessary reliability column
data=data[,-1]
data=data[,-4]

#Average data for each variable before melting
data=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu=mean(mu.spline),
            A=mean(A.spline), 
            lambda=mean(lambda.spline), 
            integral=mean(integral.spline),
            growth=mean(growth))

#Call samples with 1/3 growing no growth
data$growth[data$growth<0.4]=0

#Call parameters NA for no growth
data$A[data$growth==0]=NA
data$mu[data$growth==0]=NA
data$integral[data$growth==0]=NA
data$lambda[data$growth==0]=NA

#remove negative control data
data=data[!(data$TestId=="NEG"),]

#normalize data to growth without arsenic
normalization=subset(data, concentration %in% 0)
normalization=normalization[,-3]
data=inner_join(data, normalization, by = c("TestId", "AddId"))
full.data=data %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(mu=mu.x/mu.y,
            A=A.x/A.y,
            lambda=lambda.x/lambda.y,
            integral=integral.x/integral.y,
            growth=growth.x)

#remove rows with no arsenic (controls)
data=full.data[!(full.data$concentration=="0"),]

#widen data set so that concentrations and their 
#corresponding variables are listed horizontally
data=dcast(setDT(data), formula= TestId ~ AddId + concentration, value.var=c("mu", "A", "lambda", "integral", "growth"))

#call untested arsenic concentrations with NA growth 0 for no growth
data$growth_V_150[is.na(data$growth_V_150)]=0 
data$growth_V_200[is.na(data$growth_V_200)]=0 
data$growth_V_250[is.na(data$growth_V_250)]=0 
data$growth_V_300[is.na(data$growth_V_300)]=0 
data$growth_III_3[is.na(data$growth_III_3)]=0 
data$growth_III_5[is.na(data$growth_III_5)]=0 
data$growth_III_7[is.na(data$growth_III_7)]=0 
data$growth_III_10[is.na(data$growth_III_10)]=0 
data$growth_III_14[is.na(data$growth_III_14)]=0 
data$growth_III_20[is.na(data$growth_III_20)]=0 
data$growth_III_25[is.na(data$growth_III_25)]=0 

#make variable ID the row name
data=data.frame(data, row.names=data[,1])

#remove arsenite concentrations 20 and 25 because nothing grows there
data=data[,-grep("III_20", names(data))]
data=data[,-grep("III_25", names(data))]

#temporarily remove growth data because it produces NAs in 
#standardization
growth=data[,grep("growth", names(data))]
data=data.frame(data[,-grep("growth", names(data))])

#Standardize data
std=decostand(data, method="standardize")

#OR scaled data
scl=scale(data, center=TRUE, scale=TRUE)

#add back growth data
stdg=data.frame(cbind(std, growth))

#separate arsenate and arsenite
As3=stdg[,!colnames(stdg) %in% grep("V", colnames(stdg), value=TRUE)]
As5=stdg[,!colnames(stdg) %in% grep("III", colnames(stdg), value=TRUE)]

#Calculate euclideandistance for all datasets 
#full dataset
dist.stdg=dist(stdg)
hc.stdg=hclust(dist.stdg)
plot(hc.stdg)
str(hc.stdg)
order=hc.as3$order

matrix=as.matrix(stdg)
data.heatmap <- heatmap(matrix,  Rowv=NA, Colv=NA, scale="column", order(order))
heatmap.2(matrix, Rowv=TRUE, Colv=NA, scale="column", na.rm=TRUE, trace="none")

#arsenate

#lets try removing colums with high concentrations
As5=As5[,!colnames(As5) %in% grep("200", colnames(As5), value=TRUE)]
As5=As5[,!colnames(As5) %in% grep("250", colnames(As5), value=TRUE)]
As5=As5[,!colnames(As5) %in% grep("300", colnames(As5), value=TRUE)]
As5=As5[,!colnames(As5) %in% grep("growth", colnames(As5), value=TRUE)]
As5=As5[,!colnames(As5) %in% grep("integral", colnames(As5), value=TRUE)]

dist.as5=dist(As5)
hc.as5=hclust(dist.as5)
plot(hc.as5)
matrix=as.matrix(As5)
heatmap.2(matrix, Rowv=TRUE, Colv=NA, scale="column", na.rm=TRUE, trace="none")


#arsenite
dist.as3=dist(As3)
hc.as3=hclust(dist.as3)
plot(hc.as3)
matrix=as.matrix(As3)
heatmap.2(matrix, Rowv=TRUE, Colv=NA, scale="column", na.rm=TRUE, trace="none")


str(hc.as3)
hc.as3$order
range(hc.as3$height)

#Look at heatmap for all datasets 
#full
a=data.matrix(data)
heatmap.2(a, dendrogram = "row", trace="none", col=bluered(20), na.rm=TRUE)

#arsenate
b=as.matrix(dist.as5)
heatmap.2(b, dendrogram = "row", trace="none", col=bluered(20))

#arsenite
c=as.matrix(dist.as3)
heatmap.2(c, dendrogram = "row", trace="none", col=bluered(20))

#will need to plot arsenate and arsenic separately 
#due to different scales (concentrations), so first we will 
#split dataset
full.as5=full.data[!(full.data$AddId=="III"),]
full.as3=full.data[!(full.data$AddId=="V"),]

#plot time to exponential growth in arsenate
ggplot(data=full.as5, aes(x=concentration, y=lambda)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot time to exponential growth in arsenite
ggplot(data=full.as3, aes(x=concentration, y=lambda)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot max growth rate in arsenate
ggplot(data=full.as5, aes(x=concentration, y=mu)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot max growth rate in arsenite
ggplot(data=full.as3, aes(x=concentration, y=mu)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot max OD590 rate in arsenate
ggplot(data=full.as5, aes(x=concentration, y=A)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot max OD590 rate in arsenite
ggplot(data=full.as3, aes(x=concentration, y=A)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot growth integral in arsenate
ggplot(data=full.as5, aes(x=concentration, y=integral)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#plot growth integral in arsenite
ggplot(data=full.as3, aes(x=concentration, y=integral)) +
  geom_line(size=0.5) +
  geom_point(aes(color=concentration), size=1.5) +
  geom_point(shape=1, size=1.5, color="black") +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_wrap(~ TestId)

#set up grofit control file
control=grofit.control(parameter = 31, nboot.dr=100)

#separate arsenite and arsenate in original data format 
#because grofit's EC50 cannot handle two "treatment" 
#types at once
orig.as5=orig.data[!(orig.data$AddId=="III"),]
orig.as3=orig.data[!(orig.data$AddId=="V"),]

#Calculate integral EC50 for arsenate but first remove A2733 
#becuase its integral EC50 cannot be analyzed
orig.as5.int=orig.as5[!(orig.as5$TestId=="A2733"),]
As5.EC50=drFit(orig.as5.int, control)
As5.EC50=data.frame(summary(As5.EC50))

#Plot arsenate integral EC50s with StDev
ggplot(data=As5.EC50, aes(x=Test, y=EC50)) +
  geom_bar(stat="identity") +
  geom_errorbar(ymax=As5.EC50$EC50 + As5.EC50$sdEC50, ymin=As5.EC50$EC50 - As5.EC50$sdEC50) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Isolate", y="EC50 (integral)")

#Calculate max growth rate (mu) EC50 for arsenate
#first change grofit control to reflect mu
control.mu=grofit.control(parameter = 28, nboot.dr=100)

#remove A2716 and A2733 from dataset bc its not working
#(maybe i will fix this by using a model)
orig.as5.mu=orig.as5[!(orig.as5$TestId=="A2716"),]
orig.as5.mu=orig.as5[!(orig.as5$TestId=="A2733"),]

#test for EC50 (mu arsenate)
As5.EC50.mu=drFit(orig.as5.mu, control.mu)
As5.EC50.mu=data.frame(summary(As5.EC50.mu))

#Plot arsenate max growth rate (mu) EC50s with StDev
ggplot(data=As5.EC50.mu, aes(x=Test, y=EC50)) +
  geom_bar(stat="identity") +
  geom_errorbar(ymax=As5.EC50.mu$EC50 + As5.EC50.mu$sdEC50, ymin=As5.EC50.mu$EC50 - As5.EC50.mu$sdEC50) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Isolate", y="EC50 (max growth rate)")

#Calculate max OD590 (A) EC50 for arsenate
#first change grofit control to reflect A
control.a=grofit.control(parameter = 28, nboot.dr=100)

#remove problem A2733
orig.as5.a=orig.as5[!(orig.as5$TestId=="A2733"),]

#test for max OD590 EC50 in arsenate
As5.EC50.a=drFit(orig.as5.a, control.a)
As5.EC50.a=data.frame(summary(As5.EC50.a))

#Plot arsenate max OD590 (A) EC50s with StDev
ggplot(data=As5.EC50.a, aes(x=Test, y=EC50)) +
  geom_bar(stat="identity") +
  geom_errorbar(ymax=As5.EC50.a$EC50 + As5.EC50.a$sdEC50, ymin=As5.EC50.a$EC50 - As5.EC50.a$sdEC50) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x="Isolate", y="EC50 (max OD590)")

