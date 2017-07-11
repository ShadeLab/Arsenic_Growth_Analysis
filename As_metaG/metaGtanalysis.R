##########################
#READ IN DATA, SET UP ENV#
##########################

#read dependencies
library(phyloseq)
library(vegan)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(taxize)
library(psych)

#print working directory for future references
#note the GitHub directory for this script is as follows
#https://github.com/ShadeLab/Arsenic_Growth_Analysis/tree/master/As_metaG
wd <- print(getwd())

#change to metaG working directory
setwd("../As_metaG")
wd <- print(getwd())

#temporarily change working directory to data to bulk load files
setwd(paste(wd, "/data", sep = ""))

#read in abundance data
names <- list.files(pattern="*_45_taxonabund.txt")
data <- do.call(rbind, lapply(names, function(X) {
  data.frame(id = basename(X), read_delim(X, delim = "\t"))}))

#move back up a directory to proceed with analysis
setwd("../")

#change cen13 titles to more meaningful
data$id <- gsub("cen13_", "metagenome_", data$id)
data$id <- gsub("cen13-ct", "cultivated_", data$id)

#split columns and tidy dataset
data <- data %>%
  separate(col = id, into = c("Site", "junk"), sep = 10, remove = TRUE) %>%
  separate(col = junk, into = c("Gene", "junk"), sep = "_45_", remove = TRUE) %>%
  select(-junk)

#remove _ in front of gene name
data$Gene <- gsub("_", "", data$Gene)

#separate out rplB data for normalization
rplB <- data[which(data$Gene == "rplB"),]
data <- data[-which(data$Gene == "rplB"),]

##prep rlpB information (ie get genome estimates)
#split columns 
rplB.summarised <- rplB %>%
  group_by(Site) %>%
  summarise(Total = sum(Fraction.Abundance), rplB = sum(Abundance))

#Tidy gene data
data.tidy <- data %>%
  separate(col = Taxon, into = c("Code", "Organism"), sep = "organism=") %>%
  separate(col = Organism, into = c("Organism", "Definition"), sep = ",definition=") %>%
  select(Site, Gene, Organism:Fraction.Abundance) %>%
  group_by(Gene, Site)

#make sure abundance and fraction abundance are numbers
#R will think it's a char since it started w taxon name
data.tidy$Fraction.Abundance <- as.numeric(data.tidy$Fraction.Abundance)
data.tidy$Abundance <- as.numeric(data.tidy$Abundance)

#double check that all fraction abundances = 1
#slightly above or below is okay (Xander rounds)
summarised.total <- data.tidy %>%
  summarise(N = length(Site), Total = sum(Fraction.Abundance))

#make column for organism name and join with microbe census data and normalize to it
data.annotated <- data.tidy %>%
  left_join(rplB.summarised, by = "Site") %>%
  mutate(Normalized.Abundance.rplB = Abundance / rplB)

#summarise data to get number of genes per gene per site
data.site <- data.annotated %>%
  group_by(Gene, Site) %>%
  summarise(Count = sum(Abundance), 
            Count.rplB = sum(Normalized.Abundance.rplB))

#plot genes
ggplot(data.annotated, aes(x = Site, y = Normalized.Abundance.rplB)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Gene)

#add taxanomic information 
data.ncbi <- tax_name(query = data.annotated$Organism, 
                      get = c("genus", "class", "phylum"), db = "ncbi")


#label query "Organism" for joining purposes
data.ncbi$Organism <- data.ncbi$query

#add ncbi information to data
data.annotated <- data.annotated %>%
  left_join(data.ncbi, by = "Organism")

#replace NA in phylum with uncultured bacterium
data.annotated$phylum[is.na(data.annotated$phylum)] = "uncultured bacterium"

#make color list
color <- c("#FF7F00", "#7570B3", "#CAB2D6", "#FBB4AE", "#F0027F", "#BEBADA", "#E78AC3", "#A6D854", "#B3B3B3", "#386CB0", "#BC80BD", "#FFFFCC", "#BF5B17", "#984EA3", "#CCCCCC", "#FFFF99", "#B15928", "#F781BF", "#FDC086", "#A6CEE3", "#FDB462", "#FED9A6", "#E6AB02", "#E31A1C", "#B2DF8A", "#377EB8", "#FCCDE5", "#80B1D3", "#FFD92F", "#33A02C", "#66C2A5", "#666666", "black", "brown")

#plot genes with phylym
ggplot(data.annotated, aes(x = Site, y = Normalized.Abundance.rplB)) +
  geom_bar(stat = "identity", aes(fill = phylum)) +
  scale_fill_manual(values = color) +
  facet_wrap(~Gene) + 
  theme_classic()

#plot rplB
ggplot(rplB, aes(x = Site, y = Abundance)) +
  geom_bar(stat = "identity", aes(fill = Taxon), position = "stack") +
  scale_fill_manual(values=color)
