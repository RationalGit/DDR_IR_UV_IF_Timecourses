library("plyr")
library("stringr")
library("data.table")
library("doMC")
library("gtools")
library("plyr")
library("ggplot2")
library("ggrepel")
library("scales")
library("RColorBrewer")
library("gridExtra")
library("cowplot")
library("ggpubr")
#Clear out existing environment
rm(list = ls())

#Set the working directory to the local input folder
setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/1,2,9,11,12,13/Intensity")
##Intensity file locations
intensity.locations <- list("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/1,2,9,11,12,13/Intensity/CompartmentalAnalysis.V4_09-20-18_04;17;42",
"/Volumes/webbert-3/HTS data/201809 DDR 27 plates/4/Intensity/CompartmentalAnalysis.V4_09-20-18_01;38;48",
"/Volumes/webbert-3/HTS data/201809 DDR 27 plates/5, 7, 10/Intensity/CompartmentalAnalysis.V4_09-20-18_01;43;11",
"/Volumes/webbert-3/HTS data/201809 DDR 27 plates/R-5-5/Intensity/CompartmentalAnalysis.V4_10-16-18_06;05;17")
##Spot file locations
spot.locations <- list("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/1,2,9,11,12,13/Spot/SpotDetector.V4_09-14-18_02;16;52",
"/Volumes/webbert-3/HTS data/201809 DDR 27 plates/5, 7, 10/Spot/SpotDetector.V4_10-18-18_03;31;52",
"/Volumes/webbert-3/HTS data/201809 DDR 27 plates/R-5-5/Spot/SpotDetector.V4_10-16-18_05;52;59")

Import.data <- function(x, input){
#Check and retrieve correct order of folders:
setwd(x) 

###For each file location - to ensure that the plates are assigned in the correct order (as some were processed in the incorrect order, and outputted with a )
##Make a list of folder names
file.list = list.files(getwd())
##Make a list of folder paths
directories <- sapply(file.list, normalizePath)
##Make a list of the header files (plate.csv) in these folders
files <- lapply(directories, list.files, pattern="Plate.csv", full.names = TRUE)
##Sort the order of this list (by name)
#files <- lapply(files, sort)

##Pull out names of columns in header files
header.names <- read.table(files[[1]], nrow = 1, stringsAsFactors = FALSE, sep = ",")
##Check that there are 11 columns
#ncol(header.names)
##List apply extraction of header file contents
headers.list <- lapply(files, function(i){
  headers.data <- read.table(i, skip = 1, stringsAsFactors = FALSE, sep = ",")
})
##Bind together header file contents  
df <- rbind.fill(headers.list)
names(df) <- header.names

##Overspill of PLateID and Plate Name are consequences of naming with commas
if(is.numeric(df$`Plate Name`)){
  if(sum(df$`Plate Name`)==54){

df[,3] <- paste(df[,3], df[,4], df[,5], df[,6], df[,7], df[,8], sep = "_")
df[,9] <- paste(df[,9], df[,10], df[,11], df[,12], df[,13], df[,14], sep = "_")
df[,4] <- df[,9]
for(i in 5:11){
  j=i+10
  df[,i] <- df[,j]
  df
}
df <- df[,1:11]
} else if(sum(df$`Plate Name`)==7*27){
  df[,3] <- paste(df[,3], df[,4], df[,5], sep = "_")
  df[,6] <- paste(df[,6], df[,7], df[,8], sep = "_")
  df[,4] <- df[,6]
  for(i in 5:11){
    j=i+4
    df[,i] <- df[,j]
    df
  }
}
}

##Order df by correctly formatted Plate Name
df <- with(df, df[order(nchar(df$`Plate Name`), df$`Plate Name`),])  
##Extract correct order of Plate descriptors (folder names)
plateOrder <- df$`Unique Plate Descriptor`
##Rearrange file list to correct order
ordered.file.list <- file.list[match(plateOrder, file.list)]
##Generate ordered filepath list
ordered.directories <- sapply(ordered.file.list, normalizePath)
##Correctly ordered list of data csv
ordered.files <- lapply(ordered.directories, list.files, pattern=input, full.names = TRUE)

##empty list 
all.data.list <- list()
##Create list of csvs for data files and add a column for the plate ID (folder name)
all.data.list <- lapply(ordered.files, read.csv)

for(i in 1:length(all.data.list)){
  if(all.data.list[[i]][1,1] == 1){
    all.data.list[[i]][1] <- i
  } else{ 
    print("script already run")}
}
all.data <- rbindlist(all.data.list, idcol = "fileName")

all.data<- all.data[,(colSums(all.data != 0)) > 0, with = FALSE]
return(all.data)

View(all.data[,1:12])
View(all.data[,109:114])
}

# Wells data: spots and intensity measurements. Step 1, create list of dataframes. Step 2: bind into single dataframes
complete.dataset.spots <- lapply(spot.locations, Import.data, input = "Well.csv")
complete.dataset.intensity <- lapply(intensity.locations, Import.data, input = "Well.csv")

all.spots <- rbindlist(complete.dataset.spots, fill = TRUE)
all.intensity <- rbindlist(complete.dataset.intensity, fill = TRUE)

# Cells data: spots and intensity measurements. Step 1, create list of dataframes. Step 2: bind into single dataframes
cells.dataset.spots <- lapply(spot.locations, Import.data, input = "Cell.csv")
cells.dataset.intensity <- lapply(intensity.locations, Import.data, input = "Cell.csv")
cells.spots <- rbindlist(cells.dataset.spots, fill = TRUE)
cells.intensity <- rbindlist(cells.dataset.intensity, fill = TRUE)

# Write files ready for next steps
setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/")
write.table(all.spots, file = "all.spots.txt")
write.table(all.intensity, file = "all.intesity.txt")
write.table(cells.spots, file = "cells.spots.txt")
write.table(cells.intensity, file = "cells.intensity.txt")
