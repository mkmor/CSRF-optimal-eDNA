######################################
##### Date Created  ## 05 Oct 22 #####
##### Date Modified ## 24 Oct 22 #####
######################################

# This code is built to stack sample columns for raw metabarcoding output data 
# and obtain species taxonomic information.

library(dplyr)
library(tidyr)
library(tidyverse)
library(worms) 

# Confirm file pathway and file name patterns. 
## Ensure all files from the same project are in the same path with a naming pattern.
### e.g., 18s2_run7_lobster_species_table.csv

files=data_frame(filename=list.files(path = "/Users/morrisonme/Documents/Projects/1 new/",
                                     pattern = "species_table", 
                                     full.names = FALSE))

files = files %>% 
  mutate(filepath = paste0("/Users/morrisonme/Documents/Projects/1 new/", filename))

files


# Loads files as dataframes into a list.
df = lapply(files$filepath, read.csv)


# Names the dataframes in the list. 
names(df) = files$filename


# Create empty lists for the following loops.
taxon.names=list()
taxons = list()
stacked=list()


# Concatenate Genus and Species values into a new column - "TaxonName".
for (i in 1:length(df)) {
  df[[i]]$scientificName = paste(df[[i]]$Genus, df[[i]]$Species)
  taxon.names[[i]] = subset(df[[i]], Genus != "Hits",
                            select=scientificName) 
  df[[i]] = subset(df[[i]], select=-c(Group,Genus,Species,TaxonName))
}   

names(taxon.names) = files$filename


# Obtain taxonomic information from WoRMS for each dataframe in list.
for (j in 1:length(taxon.names)) {
  taxons[[j]] = subset(wormsbynames(taxon.names[[j]]$scientificName, marine_only = FALSE),
                       select=c(scientificname,kingdom,phylum,class,order,family,genus))
}


# Joins taxonomic information to metabarcoding data, trims extraneous columns,
# and transforms dataframes to long-form (i.e., one column per sample).
for (k in 1:length(df)){
  stacked[[k]] = df[[k]] %>%
    full_join(taxons[[k]], by=c('scientificName'='scientificname')) %>%
    relocate(c(scientificName,kingdom,phylum,class,order,family,genus), .before = Total) %>%
    drop_na(scientificName) 
  stacked[[k]] = cbind(stacked[[k]][1:7], stack(stacked[[k]][-(1:7)]))
}

names(stacked) = paste0("stacked_",files$filename)

# Write new CSV files to folder.
sapply(names(stacked), 
       function (x) write.csv(stacked[[x]], 
                              file=paste0("/Users/morrisonme/Documents/Projects/1 new/", x),
                              row.names = FALSE))
