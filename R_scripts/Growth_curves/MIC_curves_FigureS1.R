library(ggplot2)
library(dplyr)
library(reshape2)

#read in annotated data into data frames & keep filename
names=list.files(pattern="*annotated")
data=do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read.csv(X))}))

##remove duplicate isolates
#make list of dates with isolate duplicates
date=list("20150922_MIC_annotated.csv","20160418_MIC_annotated.csv","20150820_MIC_annotated.csv", "20160421_MIC_annotated.csv","20160425_MIC_annotated.csv","20160428_MIC_annotated.csv","20160516_MIC_annotated.csv")
strain=list("I2716", "I2746", "I2742","I2748","I2702","I2718","I2720","A2727","I2745","A2705","I2749","I2747","A2712","A2707","A2731")

#remove this list
data=data[-which(data$id %in% date & data$Strain %in% strain),]

#repeat for dates not included before (wouldve artificially removed wells)
date=list("20160501_MIC_annotated.csv","20160523_MIC_annotated.csv","20160527_MIC_annotated.csv","20160729_MIC_annotated.csv")
strain=list("I2748","I2747","A2723","I2742b")

#remove this new list
data=data[-which(data$id %in% date & data$Strain %in% strain),]

#define confidence interval of 95
conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}

#remove unnecessary column X
data=data[,-2]

#remove negative controls
data=data[-which(data$Strain=="NEG"),]
data=data[-which(data$Strain=="x"),]
data=data[-which(data$Strain=="NEG2"),]


#make OD590 a number
OD590=as.numeric(data$OD590)
data=cbind(data,OD590)
data=data[,-4]

#make time a number
time=as.numeric(data$Time)
data=cbind(data,time)
data=data[,-2]

# Group the data by the different experimental variables and calculate the
# sample size, average OD590, and 95% confidence limits around the mean
# among the replicates. Also remove all records where the Strain is NA.
stats <- data %>%
  group_by(Arsenic, Strain, Concentration, time) %>%
  summarise(N=length(OD590),
            Average=mean(OD590),
            CI95=conf_int95(OD590)) %>%
  filter(!is.na(Strain))

#Calculate upper and lower CI95s for each curve
stats$CI95_up=stats$Average+stats$CI95
stats$CI95_lo=stats$Average-stats$CI95

#remove negative CI95 because they are not biologically relevant
stats$CI95_up[stats$CI95_up<0]=0
stats$CI95_lo[stats$CI95_lo<0]=0

#separate arsenate and arsenite
as5=stats[-which(stats$Arsenic=="III"),]
as3=stats[-which(stats$Arsenic=="V"),]


# Plot the average OD590 over time for each strain in each environment (as5)
ggplot(data=as5, aes(x=time/3600, y=Average, color=Concentration, group=Concentration)) +
  geom_ribbon(aes(ymin=CI95_lo, ymax=CI95_up, color=Concentration),
              color=NA, alpha=0.3) +
  geom_line(aes(color=Concentration)) +
  scale_color_gradientn(colours=rainbow(8)) +
  facet_wrap(~Strain) +
  ylim(0,1.5) +
  labs(x="Time (Hours)", y="OD590") +
  theme_bw()

# Plot the average OD590 over time for each strain in each environment (as3)
ggplot(data=as3, aes(x=time/3600, y=Average, color=Concentration, group=Concentration)) +
  geom_ribbon(aes(ymin=CI95_lo, ymax=CI95_up),
              color=NA, alpha=0.3) +
  geom_line(aes(color=Concentration)) +
  scale_color_gradientn(colours=rainbow(8)) +
  facet_wrap(~Strain) +
  labs(x="Time (Hours)", y="OD590")+
  theme_bw()



