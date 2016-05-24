library(ggplot2)
library(dplyr)

#read in data by calling on files with "results" in the name, and combine all data from separate experiments
data=data.frame(do.call(rbind, lapply(list.files(pattern="*results"), read.csv)))

#remove flagged data
data=data[!(data$reliability=="FALSE"),]

#normalize the growth parameters to growth at 0mM arsenic and average them
normalization=subset(data, concentration %in% 0)
normalization=normalization %>%
  group_by(TestId, AddId) %>%
  summarise(N=length(mu.spline), 
            mu=mean(mu.spline), 
            A=mean(A.spline),
            lambda=mean(lambda.spline),
            integral=mean(integral.spline))

#Add averaged control data to the original data table
data=inner_join(data, normalization, by=c("TestId", "AddId"))

#Calculate graphs for mu
mu=cbind(data[2:4], data.frame(data$mu.spline/data$mu.y))

#Calculate statistics for mu
mu=mu %>%
  group_by(TestId, AddId, concentration) %>%
  summarise(Average=mean(data.mu.spline.data.mu.y), 
            StDev=sd(data.mu.spline.data.mu.y))

#plot data
ggplot(data=mu, aes(x=concentration, y=Average)) +
  geom_line(size=1) +
  geom_point(aes(color=concentration), size=2.5) +
  geom_point(shape=1, size=2.5, color="black") +
  geom_errorbar(aes(ymax=mu$Average + mu$StDev, ymin=mu$Average - mu$StDev)) +
  scale_color_gradientn(colors=rainbow(6)) +
  facet_grid(~AddId) +
  facet_wrap(~TestId + AddId) +
  ylim(-0.01,1.01)




