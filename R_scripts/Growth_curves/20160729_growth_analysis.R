library(reshape2)
library(dplyr)
library(ggplot2)

# Read in the raw data and the platemap. You may need to first change your
# working directory with the setwd command.
rawdata <- read.table("20160729_AsV_AsIII_MIC_I2742.asc")
platemap=read.csv("20160729_platemap.csv")

#remove time temporarily from dataset
time=rawdata[1]
rawdata=rawdata[,-1]

#define negative controls and subtract from dataset
neg0=data.frame(rawdata[,5])
neg0=data.frame(rep(neg0,5))
neg1=data.frame(rawdata[,17])
neg1=data.frame(rep(neg1,5))
neg3=data.frame(rawdata[,29])
neg3=data.frame(rep(neg3,5))
neg5=data.frame(rawdata[,41])
neg5=data.frame(rep(neg5,5))
neg7=data.frame(rawdata[,53])
neg7=data.frame(rep(neg7,5))
neg10=data.frame(rawdata[,65])
neg10=data.frame(rep(neg10,5))
neg14=data.frame(rawdata[,77])
neg14=data.frame(rep(neg14,5))
neg20=data.frame(rawdata[,89])
neg20=data.frame(rep(neg20,5))
vneg0=data.frame(rawdata[,10])
vneg0=data.frame(rep(vneg0,5))
vneg10=data.frame(rawdata[,22])
vneg10=data.frame(rep(vneg10,5))
vneg50=data.frame(rawdata[,34])
vneg50=data.frame(rep(vneg50,5))
vneg100=data.frame(rawdata[,46])
vneg100=data.frame(rep(vneg100,5))
vneg150=data.frame(rawdata[,58])
vneg150=data.frame(rep(vneg150,5))
vneg200=data.frame(rawdata[,70])
vneg200=data.frame(rep(vneg200,5))
vneg250=data.frame(rawdata[,82])
vneg250=data.frame(rep(vneg250,5))
vneg300=data.frame(rawdata[,94])
vneg300=data.frame(rep(vneg300,5))
x=data.frame(rawdata[,11:12])
neg=data.frame(cbind(neg0,vneg0,x,neg1,vneg10,x,neg3,vneg50,x,neg5,vneg100,x,neg7,vneg150,x,neg10,vneg200,x,neg14,vneg250,x,neg20,vneg300,x))
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
write.csv(annotated, "20160729_MIC_annotated.csv")

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
