library(reshape2)
library(dplyr)
library(ggplot2)

# Read in the raw data and the platemap. You may need to first change your
# working directory with the setwd command.
rawdata <- read.table("20160425_01_AsV_MIC_A06,07,12,24.asc")
platemap <- read.csv("20160425_Platemap_High_MIC.csv")

#remove time temporarily from dataset
time=rawdata[1]
rawdata=rawdata[,-1]

#define negative controls and subtract from dataset
neg0=data.frame(rawdata[,10])
neg0=data.frame(rep(neg0, 12))
neg10=data.frame(rawdata[,22])
neg10=data.frame(rep(neg10, 12))
neg50=data.frame(rawdata[,34])
neg50=data.frame(rep(neg50, 12))
neg100=data.frame(rawdata[,46])
neg100=data.frame(rep(neg100, 12))
neg150=data.frame(rawdata[,58])
neg150=data.frame(rep(neg150, 12))
neg200=data.frame(rawdata[,70])
neg200=data.frame(rep(neg200, 10))
neg200=data.frame(cbind(neg200,rawdata[,10],rawdata[,46]))
neg250=data.frame(rawdata[,82])
neg250=data.frame(rep(neg250,10))
neg250=data.frame(cbind(neg250,rawdata[,22],rawdata[,58]))
neg300=data.frame(rawdata[,94])
neg300=data.frame(rep(neg300,10))
neg300=data.frame(cbind(neg300,rawdata[,34],rawdata[,10]))
neg=data.frame(cbind(neg0,neg10,neg50,neg100,neg150,neg200,neg250,neg300))
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
write.csv(annotated, "20160425_MIC_annotated.csv")

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
  group_by(Strain, Concentration, time) %>%
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
  facet_wrap(~Strain, nrow=2) +
  labs(x="Time (Hours)", y="Absorbance at 590 nm")
