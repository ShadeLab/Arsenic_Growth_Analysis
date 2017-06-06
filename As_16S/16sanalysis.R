#load required packages
library(tidyverse)
library(reshape2)

##########################
#Env setup (read in data)#
##########################

#label working directory
wd <- getwd()

#read in OTU table
otutable <- read_delim(paste(wd, "/data/MASTER_OTU_hdf5_filteredfailedalignments_rdp_rmCM_collapse_even321000.txt", sep = ""), delim = "\t", comment = "#") 

#read blast results
blast <- read_delim(file = paste(wd, "/data/as.cen.tabular.txt", sep = ""), delim = "\t", col_names = FALSE, 
                    col_types =   list(col_character(), col_character(), col_double(), col_integer(), col_integer(),col_integer(),col_integer(),col_integer(),col_integer(), col_integer(),col_double(),col_integer()))


#read in metadata
meta <- read.delim(paste(wd, "/data/Centralia_full_map.txt", sep = ""), sep=" ")

################################
#Absolute abundance comparisons#
################################
#add column names to blast results
colnames(blast) <- c("query", "OTU_ID", "perc.id", "align.length", "mismatches", "gaps", "q.start", "q.end", "s.start", "s.end", "e.value", "bit.score")

#separate site nubmber
cen.blast <- blast %>%
  left_join(otutable, by = "OTU_ID")

#only look at cen13
blast.96 <- cen.blast[which(cen.blast$perc.id > 96),]

#check that all isolates have hit
length(unique(blast.96$query))
#should be =25

#remove rows where C13 is zero
blast.96.present <- blast.96[-which(blast.96$C13 == 0),]

#check that all isolates still have hit
length(unique(blast.96.present$query))
#should be =25

#make query a character
blast.96.present$query <- as.character(blast.96.present$query)

#remove duplicate hits for isolate (ie keep hit with highest %identity)
blast.96.present.top <- blast.96.present[!duplicated(blast.96.present$query),]

#remove duplicate OTU ids (they mean the isolates are the same via 16S)
blast.96.present.top <- blast.96.present.top[!duplicated(blast.96.present.top$OTU_ID),]

#remove unimportant columns
blast.96.present.top <- blast.96.present.top %>%
  select(query, C03: ConsensusLineage)

#tidy data for plotting
blast.96.tidy <- melt(blast.96.present.top, id.vars = c("query", "ConsensusLineage"), variable.name = "Site", value.name = "Abundance")

#adjust site name to relate to metadata
blast.96.tidy$Site <- gsub("C", "Cen", blast.96.tidy$Site)

#join metadata with tidy blast data
blast.96.tidy <- blast.96.tidy %>%
  left_join(meta, by = "Site") 

#order Site by temperature
blast.96.tidy$Site <- factor(blast.96.tidy$Site , 
                             levels = blast.96.tidy$Site [order(blast.96.tidy$SoilTemperature_to10cm)])



#plot data
(abundance.plot <- ggplot(blast.96.tidy, aes(x = Site, y = Abundance, fill = Classification)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
  facet_wrap(~ query) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.95,vjust=0.2)))
        
ggsave(abundance.plot, filename = paste(wd, "/figures/isolate.abundance.png", sep = ""))

#get relative abundance of Centralia OTUs
colsums <- apply(otutable[,2:19], 2, sum)
#all should be = 321000

#divide abundance by 321000 to get relative abundance
blast.96.tidy.rel <- blast.96.tidy %>%
  mutate(rel.Abundance = Abundance/ 321000)
  
#plot data
(rel.abundance.plot <- ggplot(blast.96.tidy.rel, aes(x = Site, y = rel.Abundance*10000, fill = Classification)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = c("firebrick2", "yellow1", "green3")) +
    facet_wrap(~ query) +
    ylab("Relative abundance (E-4)") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, size = 12, hjust=0.95,vjust=0.2)))

ggsave(rel.abundance.plot, filename = paste(wd, "/figures/isolate.rel.abundance.png", sep = ""))

  
  
  
  
  
  
  
  
  


