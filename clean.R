library("plyr")
library("stringr")
library("data.table")
library("doMC")
library("gtools")
library("plyr")
library("stringr")
library("data.table")
library("doMC")
library("gtools")

# Import dataframes
setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/")
all.spots <- read.table("all.spots.txt")
all.intensity <- read.table("all.intesity.txt")

cells.spots <- read.table("cells.spots.txt")
cells.intensity <- read.table("cells.intensity.txt")

# Identify columns of interest for this analysis
colnames(cells.spots)
cells.spot.columns <- c("fileName", "PlateNumber","WellId", "Row", "Column", "CellNumber", "ObjectAreaCh1", "ObjectTotalIntenCh1", "SpotCountCh2", "SpotTotalAreaCh2", "SpotAvgAreaCh2", "SpotTotalIntenCh2", "SpotAvgIntenCh2", "TotalIntenCh2", "AvgIntenCh2", "SpotCountCh3", "SpotTotalAreaCh3", "SpotAvgAreaCh3", "SpotTotalIntenCh3", "SpotAvgIntenCh3", "TotalIntenCh3",  "AvgIntenCh3")
colnames(cells.intensity)
cells.intensity.columns <- c("fileName", "PlateNumber","WellId", "Row", "Column", "CellNumber", "ObjectAreaCh1", "ObjectTotalIntenCh1", "CircTotalIntenCh2", "CircAvgIntenCh2", "TotalIntenCh2", "AvgIntenCh2", "CircTotalIntenCh3", "CircAvgIntenCh3", "TotalIntenCh3", "AvgIntenCh3")
cells.spots <- setDT(cells.spots)[, ..cells.spot.columns]
cells.intensity <- setDT(cells.intensity)[, ..cells.intensity.columns]

# add triplicate numbers based on PlateNumber 
triplicates <- rep(1:9, each = 3)
names(triplicates) <- 1:27
cells.spots$triplicate.number <- triplicates[cells.spots$PlateNumber]
cells.intensity$triplicate.number <- triplicates[cells.intensity$PlateNumber]
# sequence along plate number and rep(1:9, each 3)
cells.spot <- setDT(cells.spot)[, triplicate.number := sequence(.N), by = rep(1:9, each = 3)]


### Check for cells data
Median.df <- function(df){
setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/")

channel.key <- read.csv("ChannelKeyCellomics.csv")
triplicate.key <- read.csv("triplicateKeyCellomics.csv")

df.sorted <- setorder(df, PlateNumber, Row, Column)
df.sorted$WellId <- as.factor(gsub("\\s", "", df.sorted$WellId))
#df.sorted <- intensity.table.all[order(PlateNumber, Row, Column),]
#df.sorted$WellId <- as.factor(gsub("\\s", "", df.sorted$WellId))

###remove spaces from 'WellId' factor to allow the key to match the channel.key table to the intensity.table data
###Testing matching of WellId factors:
#oldWellId <- df.sorted$WellId
#newWellId <- channel.key$WellId
#Reduce(intersect, list(oldWellId, newWellId))


###Create median table
df.sorted <- merge(setDT(df.sorted), triplicate.key)

if(length(grep("MEAN",colnames(df), value =TRUE))>0){
# MEANs
 df.sorted <- setDT(df.sorted)[, sapply(.SD, function(x) list(median = median(x))), .SDcols = unlist(lapply(setDT(df.sorted), is.numeric)), by = list(WellId, triplicate.number)]
# Cells
}

df.sorted <- merge(setDT(df.sorted), channel.key, by = "WellId")
# MEANs
# df.sorted$triplicate.number <- rep(1:27, each = nrow(df.sorted)/9, length.out = nrow(df.sorted))
# df.sorted$triplicate.well.ID <- paste(df.sorted$triplicate.number, df.sorted$WellId)


#View(df.sorted)
# df.sorted$Ch2 <- channel.key$Ch2[match(df.sorted$WellId, channel.key$WellId)]
# df.sorted$Ch3 <- channel.key$Ch3[match(df.sorted$WellId, channel.key$WellId)]
# df.sorted$siRNA <- siRNA.key$siRNA[match(df.sorted$WellId, channel.key$WellId)]
# df.sorted$Column <- df.sorted$Column[match(df.sorted$WellId, df.sorted$WellId)]
# df.sorted$Row <- df.sorted$Row[match(df.sorted$WellId, df.sorted$WellId)]
df.sorted$Order <- channel.key$Order[match(df.sorted$WellId, channel.key$WellId)]
df.sorted$stain.order <- channel.key$stain.order[match(df.sorted$WellId, channel.key$WellId)]

# df.sorted$treatment <- triplicate.key$treatment[match(df.sorted$triplicate.number, triplicate.key$triplicate.number)]

return(df.sorted)
}

colnames(channel.key)
colnames(cells.spots)

cells.spots$WellId
channel.key$WellId
# Remove sapces for WellIDs (prevents merge)
as.factor(gsub("\\s", "", cells.spots$WellId))

cells.spots <- Median.df(cells.spots)
cells.intensity <- Median.df(cells.intensity)

median.spots <- Median.df(all.spots)
median.intensity <- Median.df(all.intensity)

View(rev(median.spots)[,1:10])
View(median.spots[,1:20])

setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/")
write.table(median.spots, file = "median.spots.txt")
write.table(median.intensity, file = "median.intesity.txt")
write.table(cells.spots, file = "cells.spots.final.txt")
write.table(cells.intensity, file = "cells.intensity.final.txt")

### MAking G2 Cells data:
## Bimodal distribution to find minima for siNT of Untreated for each dataframe:

# df <- read.csv(header=F,"data.txt")
# colnames(df) = "X"
# 
# # bimodal
x <- setDT(Dfs.intensity.IR[[1]])[siRNA == "siNT" & treatment == "Untreated"]

View(x)
setDT(x)[siRNA == "siNT"]
colnames(ch3.Dfs.intensity.IR[[1]])
km <- kmeans(x[,"ObjectTotalIntenCh1"],centers=2)
k <- km$size
y <- x[km$cluster == 1]
min(x$ObjectTotalIntenCh1)
x$clust <- as.factor(km$cluster)
 library(ggplot2)
 ggplot(x, aes(x =ObjectTotalIntenCh1)) + 
   geom_histogram(aes(fill=clust,y=..count../sum(..count..)),
                  binwidth=500000, color="grey50")+
   stat_density(geom="line", color="red")
