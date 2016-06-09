library(reshape2)
library(dplyr)
library(ggplot2)

# Read in the raw data and the platemap. You may need to first change your
# working directory with the setwd command.
rawdata <- read.table("20160531_AsV_AsIII_I23,46,47,49,59.asc")
platemap <- read.csv("20160531_platemap.csv")

#remove time temporarily from dataset
time=rawdata[1]
rawdata=rawdata[,-1]

#define negative controls and subtract from dataset
neg0=data.frame(rawdata[,7])
neg10=data.frame(rawdata[,8])
neg50=data.frame(rawdata[,9])
neg100=data.frame(rawdata[,10])
neg150=data.frame(rawdata[,11])
neg200=data.frame(rawdata[,12])
negV=data.frame(cbind(neg0,neg10,neg50,neg100,neg150,neg200))
neg=data.frame(rep(negV,9))
neg0=data.frame(rawdata[,55])
neg1=data.frame(rawdata[,56])
neg3=data.frame(rawdata[,57])
neg5=data.frame(rawdata[,58])
neg=data.frame(cbind(neg,neg0,neg1,neg3,neg5,neg5,neg5,negV))
neg0=data.frame(rep(neg0,3))
neg1=data.frame(rep(neg1,3))
neg3=data.frame(rep(neg3,3))
neg5=data.frame(rep(neg5,3))
neg=data.frame(cbind(neg,neg0,neg1,neg0,neg1,neg3,neg5,neg0,neg1,neg3,neg5))
rawdata= rawdata - neg

#add time back to data
rawdata=cbind(time, rawdata)

#add time and well infromation to data (warning message is OK in first line)
labels=read.csv("Time,wells.csv", header=FALSE)
labels=t(labels)
rawdata=t(rawdata)
rawdata=cbind(labels,rawdata)
rawdata=data.frame(rawdata, row.names=TRUE)
rawdata=data.frame(t(rawdata))

#remove s from time column
rawdata$Time=gsub('s', '', rawdata$Time)

# Reshape the data. Instead of rows containing the Time, Temperature,
# and readings for each Well, rows will contain the Time, Temperature, a
# Well ID, and the reading at that Well.
reshaped <- melt(rawdata, id="Time", variable.name="Well",
                 value.name="OD590")

# Add information about the experiment from the plate map. For each Well
# defined in both the reshaped data and the platemap, each resulting row
# will contain the absorbance measurement as well as the additional columns
# and values from the platemap.
annotated <- inner_join(reshaped, platemap, by="Well")

# Save the annotated data as a CSV for storing, sharing, etc.
write.csv(annotated, "20160531_MIC_annotated.csv")

conf_int95 <- function(data) {
  n <- length(data)
  error <- qt(0.975, df=n-1) * sd(data)/sqrt(n)
  return(error)
}

#make OD590 a number
OD590=as.numeric(annotated$OD590)
annotated=cbind(annotated,OD590)
annotated=annotated[,-3]

#make time a number
time=as.numeric(annotated$Time)
annotated2=cbind(annotated,time)
annotated=annotated2[,-1]
  
# Group the data by the different experimental variables and calculate the
# sample size, average OD600, and 95% confidence limits around the mean
# among the replicates. Also remove all records where the Strain is NA.
stats <- annotated %>%
  group_by(Arsenic, Strain, Concentration, time) %>%
  summarise(N=length(OD590),
            Average=mean(OD590),
            CI95=conf_int95(OD590)) %>%
  filter(!is.na(Strain))

# Plot the average OD600 over time for each strain in each environment
ggplot(data=stats, aes(x=time/3600, y=Average, color=Concentration, group=Concentration)) +
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, color=Concentration),
              color=NA, alpha=0.3) +
  geom_line(aes(color=Concentration)) +
  scale_color_gradientn(colours=rainbow(8)) +
  facet_wrap(Arsenic~Strain) +
  labs(x="Time (Hours)", y="Absorbance at 590 nm")
