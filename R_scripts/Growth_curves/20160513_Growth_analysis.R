library(reshape2)
library(dplyr)
library(ggplot2)

# Read in the raw data and the platemap. You may need to first change your
# working directory with the setwd command.
rawdata <- read.table("20160513_01_AsIII_AsV_A12,24,31,33,07.asc")
platemap <- read.csv("20160513_Platemap_MIC.csv")

#remove time temporarily from dataset
time=rawdata[1]
rawdata=rawdata[,-1]

#define negative controls and subtract from dataset
Ineg0=data.frame(rawdata[,1])
Ineg1=data.frame(rawdata[,2])
Ineg3=data.frame(rawdata[,3])
Ineg5=data.frame(rawdata[,4])
Ineg7=data.frame(rawdata[,5])
Ineg=data.frame(cbind(Ineg0,Ineg1,Ineg3,Ineg5,Ineg7))
Aneg0=data.frame(rawdata[,60])
Aneg10=data.frame(rawdata[,72])
Aneg50=data.frame(rawdata[,84])
Aneg100=data.frame(rawdata[,96])
Aneg150=data.frame(rawdata[,10])
Aneg=data.frame(cbind(Aneg0,Aneg10,Aneg50,Aneg100,Aneg150))
negA=data.frame(cbind(Ineg,Ineg0,Ineg1,Ineg3,Ineg5,Aneg150,Aneg0,Aneg0))
negB=data.frame(cbind(Ineg,Ineg,Aneg10,Aneg10))
negC=data.frame(cbind(Ineg,Ineg,Aneg50,Aneg50))
negD=data.frame(cbind(Ineg,Ineg,Aneg100,Aneg100))
negE=data.frame(cbind(Ineg,Ineg,Aneg0,Aneg0))
negF=data.frame(cbind(Ineg,Aneg,Aneg10,Aneg10))
negG=data.frame(cbind(Ineg,Aneg,Aneg50,Aneg50))
negH=data.frame(cbind(Ineg,Aneg,Aneg100,Aneg100))
neg=data.frame(cbind(negA,negB,negC,negD,negE,negF,negG,negH))
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
write.csv(annotated, "20160513_MIC_annotated.csv")

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

#separate arsenate and arsenite values
arsenate=subset(annotated, Arsenic %in% "V")
arsenite=subset(annotated, Arsenic %in% "III")

# Group the data by the different experimental variables and calculate the
# sample size, average OD600, and 95% confidence limits around the mean
# among the replicates. Also remove all records where the Strain is NA.Repeat this for both AsV and AsIII
stats_V <- arsenate %>%
  group_by(Strain, Arsenic, Concentration, time) %>%
  summarise(N=length(OD590),
            Average=mean(OD590),
            CI95=conf_int95(OD590)) %>%
  filter(!is.na(Strain))

stats_III <- arsenite %>%
  group_by(Strain, Arsenic, Concentration, time) %>%
  summarise(N=length(OD590),
            Average=mean(OD590),
            CI95=conf_int95(OD590)) %>%
  filter(!is.na(Strain))

# Plot the average OD600 over time for each strain in each environment. Repeat for AsIII and AsV
ggplot(data=stats_V, aes(x=time/3600, y=Average, color=Concentration, group=Concentration)) +
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, color=Concentration),
              color=NA, alpha=0.3) +
  geom_line(aes(color=Concentration)) +
  scale_color_gradientn(colours=rainbow(8)) +
  facet_wrap(~Strain, nrow=2) +
  labs(x="Time (Hours)", y="Absorbance at 590 nm")

ggplot(data=stats_III, aes(x=time/3600, y=Average, color=Concentration, group=Concentration)) +
  geom_ribbon(aes(ymin=Average-CI95, ymax=Average+CI95, color=Concentration),
              color=NA, alpha=0.3) +
  geom_line(aes(color=Concentration)) +
  scale_color_gradientn(colours=rainbow(8)) +
  facet_wrap(~Strain, nrow=2) +
  labs(x="Time (Hours)", y="Absorbance at 590 nm")
