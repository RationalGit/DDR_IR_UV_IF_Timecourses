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

setwd("/Volumes/webbert-3/HTS data/201809 DDR 27 plates/")
median.spots <- read.table("median.spots.txt")
median.intensity <- read.table("median.intesity.txt")


#Set the working directory to the local input folder
setwd('/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/DDR IF/DDR IF/201809 27 plate triplicate/R analysis')
a <- median.intensity
b <- median.spots
median.intensity <- a
median.spots <- b

### Function will begin here, taking 'median.intensity and median.spots' as the arguments

# to.match <- c("GENE1", "GENE2")
# to.match <- c("GENE3", "GENE2")
# to.match <- c("GENE3", "GENE1")

# median.spots <- median.spots[ -grep(paste(to.match, collapse = "|"),median.spots$siRNA),]
# median.intensity <- median.intensity[ -grep(paste(to.match, collapse = "|"),median.intensity$siRNA),]

## Subset dataframes into IR or UV
intensity.IR.df <- median.intensity[ -grep("UV",median.intensity$treatment),]
intensity.UV.df <- median.intensity[ -grep("IR",median.intensity$treatment),]

spots.IR.df <- median.spots[ -grep("UV",median.spots$treatment),]
spots.UV.df <- median.spots[ -grep("IR",median.spots$treatment),]

## Subset dataframes into relevant stains
Channel.sep <- function(df, new.name){
  new.name[[1]]<- as.character((unique(df$Ch2)))
  new.name[[2]]<- as.character((unique(df$Ch3)))
  return(new.name)
}
stains.intensity <- c(ch2.stains.intensity = list(), ch3.stains.intensity = list())
stains.spots <- c(ch2.stains.spots = list(), ch3.stains.spots = list())

all.stains.intensity <- unlist(Channel.sep(median.intensity, stains.intensity))
all.stains.intensity = all.stains.intensity[!is.na(all.stains.intensity)]

all.stains.spots <- unlist(Channel.sep(median.spots, stains.spots))
all.stains.spots = all.stains.spots[!is.na(all.stains.spots)]

#make a list of dataframes split by the stains
Dfs.intensity.IR <- split(intensity.IR.df, by = "Stain.order", drop = FALSE)
Dfs.intensity.UV <- split(intensity.UV.df, by = "Stain.order", drop = FALSE)

Dfs.spots.IR <- split(spot.IR.df, by = "Stain.order", drop = FALSE)
Dfs.spots.IR[["4"]] <- data.frame()
Dfs.spots.UV <- split(spot.UV.df, by = "Stain.order", drop = FALSE)
Dfs.spots.UV[["4"]] <- data.frame()

## Find measurments of interest (manual)
### medians of means

## Intensity details (reference of the measurement and the 'chart name')
#print(colnames(median.intensity)[grep("MEAN",colnames(median.intensity))])

graph.features.intensity = data.frame( channel = c("Ch1", "Ch2", "Ch2", "Ch3", "Ch3"), 
                                       measurement = c("ValidObjectCount.median", 
                                                       "MEAN_CircTotalIntenCh2.median", "MEAN_CircAvgIntenCh2.median",
                                                       "MEAN_CircTotalIntenCh3.median", "MEAN_CircAvgIntenCh3.median"),
                                       title =       c("Nuclei count",  "Mean Total Nuc. intensity", "Mean Nuc. intensity", 
                                                       "Mean Total Nuc. intensity", "Mean Nuc. intensity"))
## Spots details (reference of the measurement and the 'chart name')
#print(colnames(median.spots)[grep("MEAN",colnames(median.spots))])

graph.features.spots = data.frame(channel = c("Ch1", "Ch2", "Ch2", "Ch2", "Ch2", "Ch3", "Ch3", "Ch3", "Ch3"), 
                                  measurement = c("MEAN_ObjectAreaCh1.median", "MEAN_ObjectSpotTotalCountCh2.median", "MEAN_ObjectSpotAvgAreaCh2.median", "MEAN_ObjectSpotAvgIntenCh2.median", "MEAN_ObjectSpotTotalIntenCh2.median", "MEAN_ObjectSpotTotalCountCh3.median", "MEAN_ObjectSpotAvgAreaCh3.median", "MEAN_ObjectSpotAvgIntenCh3.median", "MEAN_ObjectSpotTotalIntenCh3.median"),
                                  title =       c("Mean Object area",  "Mean Spot Count", "Mean Spot area", "Mean Spot intensity", "Mean Total Spot intensity per nucleus", "Mean Spot Count", "Mean Spot area", "Mean Spot intensity", "Mean Total Spot intensity per nucleus"))

### CELLS
## Intensity details (reference of the measurement and the 'chart name')
#print(colnames(median.intensity)[grep("MEAN",colnames(median.intensity))])

graph.cells.spots = data.frame( channel = c("Ch1", "Ch1", "Ch1", 
                                            "Ch2", "Ch2", "Ch2", "Ch2", "Ch2",  
                                            "Ch3", "Ch3", "Ch3", "Ch3", "Ch3"), 
                                measurement = c("CellNumber", "ObjectAreaCh1", "ObjectTotalIntenCh1", 
                                                "SpotCountCh2", "SpotTotalAreaCh2", "SpotAvgAreaCh2", "SpotTotalIntenCh2", "SpotAvgIntenCh2",  
                                                "SpotCountCh3", "SpotTotalAreaCh3", "SpotAvgAreaCh3", "SpotTotalIntenCh3", "SpotAvgIntenCh3"),
                                title =       c("Nuclei Count", "Nuclei Area", "DAPI Intensity",
                                                "Nuclear Spot Count", "Nuclear Spot Area", "Average Spot Area", "Total Spot Intensity", "Avg. Spot Intensity", 
                                                "Nuclear Spot Count", "Nuclear Spot Area", "Average Spot Area", "Total Spot Intensity", "Avg. Spot Intensity"))
## Spots details (reference of the measurement and the 'chart name')
#print(colnames(median.spots)[grep("MEAN",colnames(median.spots))])

graph.cells.intensity = data.frame(channel = c("Ch1", "Ch1", "Ch1", 
                                               "Ch2", "Ch2",  
                                               "Ch3", "Ch3"), 
                                   measurement = c("CellNumber", "ObjectAreaCh1", "ObjectTotalIntenCh1", 
                                                   "CircTotalIntenCh2", "CircAvgIntenCh2",
                                                   "CircTotalIntenCh3", "CircAvgIntenCh3"),
                                   title =       c("Nuclei Count", "Nuclei Area", "DAPI Intensity", 
                                                   "Total Nuclear Intensity", "Average Nuclear Intensity", 
                                                   "Total Nuclear Intensity", "Average Nuclear Intensity"))

#Iterate through this list to create labelled bar plots (and pngs) for each stain
setwd('/Users/webbert/Documents/Work/Experiments/LateSubDriv Large files:folders (unsynced)/Gene workup/DDR/DDR IF/DDR IF/201809 27 plate triplicate/R analysis/Output/')  
setwd('./G2')
setwd('./All_4_genes')
save.image = FALSE

plot.title <- "(IR treated)"
plot.titles.list <- c("(IR treated)", "(UV treated)")


### Cells plots

Make.violin.plots = function(save.image = FALSE, df, plot.title, graph.features, feature, cell.cycle.phase = "All"){
  
  y_Axis = as.character(graph.features$title[graph.features$measurement == feature])
  print(y_Axis)
  ifelse(as.character(graph.features$channel[graph.features$measurement == feature]) == "Ch1",
         channel.label <- "DAPI", 
         channel.label <- df[[1,as.character(graph.features$channel[graph.features$measurement == feature])]])
  
  ### Allow separation of G1 and G2 phases (from siNT, untreated plot)
  
  if(cell.cycle.phase %in% c("G2", "G1")){
    km.df <- setDT(df)[siRNA == "siNT" & treatment == "Untreated"]
    km <- kmeans(km.df[,"ObjectTotalIntenCh1"],centers=2)
    km.N.1 <- km.df[km$cluster == 1]
    km.N.2 <- km.df[km$cluster == 2]
    G2.line <- max(min(km.N.1[,"ObjectTotalIntenCh1"]),min(km.N.2[,"ObjectTotalIntenCh1"]))
    print(nrow(df))
  }
  if(cell.cycle.phase == "G2"){
    df <- df[ObjectTotalIntenCh1 > G2.line]
    print(nrow(df))
    plot.title <- paste0(plot.title, " (~G2 only)")
  } else if(cell.cycle.phase == "G1"){
    df <- df[ObjectTotalIntenCh1 < G2.line]
    print(nrow(df))
    plot.title <- paste0(plot.title, " (~G1 only)")
  }
  
  # Colour of graph based on IR or UV
  if(plot.title=="(IR treated)"){
    col.pal = brewer.pal(n=9, "Oranges")[4:8]
  } else if(plot.title=="(UV treated)"){
    col.pal = brewer.pal(n=9, "Purples")[4:8]
  }
  ylim1 = boxplot.stats(df[[feature]])$stats[c(1, 5)]
  #df <- df[treatment %in% c("Untreated", "0.5 hr post-IR")]
  
  graph.title = paste(y_Axis, " (", channel.label, ") ", plot.title, sep = "")
  
  working.plot = ggplot(df,
                        aes(x = reorder(siRNA, Order), fill = reorder(treatment, triplicate.number), y = df[[feature]])) +
    labs(x = "siRNA", y = y_Axis, fill = "Plate Treatment", title = graph.title) + #S9.6 #pRPA2 S33
    geom_boxplot(position = "dodge") +
    
    # colour of bars
    scale_fill_manual(values = col.pal) +
    #facet_wrap(~plateTreatment +reorder(siRNA, Row), nrow=3) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    coord_cartesian(ylim = ylim1*1.05)
  
  #+
  #ylim(0, newLim) +
  working.plot <- working.plot + theme(text = element_text(size = 7), 
                                       axis.text.x = element_text(size = 7,
                                                                  angle = 45, vjust = 1, hjust=1),
                                       axis.title.x = element_blank(),
                                       axis.text.y = element_text(size = 7),
                                       axis.title.y = element_text(face = "bold"),
                                       plot.title = element_text(size = 8, hjust = 0),
                                       aspect.ratio = 3/4,
                                       legend.text = element_text(size = 6),
                                       legend.key.size = unit(0.5, "line"),
                                       legend.box.margin = margin(0,0,0,-20),
                                       plot.margin = margin(0,0,0,0))
  
  return(working.plot)
  # if(save.image){
  #   file_name = paste(gsub("/", "x", graph.title), ".png", sep="")
  #   png(file_name)
  #   print(working.plot)
  #   dev.off()
  #   print("saved")
  #   print(working.plot)
  # 
  # } else{
  #   print(working.plot)
  # }
}

Make.violin.plots(save.image = FALSE, 
                  df = Dfs.spots.IR[[6]], 
                  plot.title = plot.title, 
                  graph.features = graph.cells.spots, 
                  "SpotCountCh2"
)

limits: 
  
  
  compare_means(ch3.Dfs.intensity.IR[[6]][,"ObjectTotalIntenCh1"] ~ ch3.Dfs.intensity.IR[[6]][,"siRNA"] + ch3.Dfs.intensity.IR[[6]][,"treatment"],  data = ToothGrowth, method = "anova")
##intensity IR:
# if(plot.title == "(IR treated)"){
#   ##IR:
#   for (i in 1:length(ch3.Dfs.spots.IR)){
#     this.df <- ch3.Dfs.spots.IR[[i]]
#     #yl <- yLimitsIR[i]
#     makePlots(save.image, df=this.df, plot.title, graph.features, "MEAN_CircAvgIntenCh2.median")
#   } if(plot.title == "(UV treated)"){
#     ##UV:
#     for (i in 1:length(ch3.Dfs.spots.UV)){
#       this.df <- ch3.Dfs.spots.UV[[i]]
#       #yl <- yLimitsUV[i]
#       makePlots(save.image, df=this.df, plot.title, graph.features, "MEAN_CircAvgIntenCh2.median")
#     }
#   }
# }

### Function to make plots, 1 at a time
Make.plots = function(save.image = FALSE, df, plot.title, graph.features, feature){
  # Use measurement from dataframe of features
  y_Axis = as.character(graph.features$title[graph.features$measurement == as.character(feature)])
  # Choice of label for plot title: 
  if( as.character(graph.features$channel[graph.features$measurement == as.character(feature)]) == "Ch2"){
    channel.label <- df[[1,"Ch2"]]
  } else{
    channel.label <- df[[1,"Ch3"]]
  } 
  # Graph title 
  graph.title = paste(y_Axis, " (", channel.label, ") ", plot.title, sep = "")
  
  # Colour of graph based on IR or UV
  if(plot.title=="(IR treated)"){
    col.pal = brewer.pal(n=9, "Oranges")[4:8]
  } else if(plot.title=="(UV treated)"){
    col.pal = brewer.pal(n=9, "Purples")[4:8]
  }
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  working.plot = ggplot(df,
                        aes(x = reorder(siRNA, Order), fill = reorder(treatment, triplicate.number), y = df[[feature]])) + 
    # labels
    labs(x = "siRNA", y = y_Axis, fill = "Plate Treatment", title = graph.title) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_bar(stat="identity", position = "dodge") + 
    
    # SEM error bars
    geom_errorbar(aes(x = reorder(siRNA, Order), ymin = (df[[feature]]-df[[gsub("MEAN", "SE", feature)]]), ymax = (df[[feature]]+df[[gsub("MEAN", "SE", feature)]])), width = .2, position = position_dodge(.9)) +
    
    # colour of bars
    scale_fill_manual(values = col.pal) +
    # text orientation
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  print(df[[(gsub("MEAN", "SD", feature))]])
  print(df[[feature]])
  print(working.plot)
  ###Saving in function
  #if(save.image){
  #  ggsave(
  #    paste(gsub("/", "x", graph.title), ".png", sep=""),
  #    working.plot,
  #    width = 6,
  #    height = 6,
  #    scale = 1,
  #    dpi = 180
  #  )
  #    dev.off()
  #    print("saved")
  #    print(working.plot)
  
  #file_name = paste(gsub("/", "x", graph.title), ".png", sep="")
  #png(file_name, width = 960, height = 960, bg = "transparent", pointsize = 12)
  #print(working.plot)
  
  
}

Make.plots(save.image, Dfs.intensity.IR[[8]], plot.title, graph.features.intensity, "MEAN_CircTotalIntenCh2.median")


for(i in 1:2){
  plot.title = plot.titles.list[i]
  if(i ==1){df = Dfs.spots.IR
  } else if(i ==2){df = Dfs.spots.UV
  }
  plots <- lapply(df, make.plots.spot, save.image=FALSE, plot.title, feature= lapply(1:9, function(x) {as.character(graph.features.spot[x,2])}))
  #lapply()
  #lapply(names(plots), function(x) ggsave(filename=paste(x,".png",sep=""), plot=plots[[x]])))
}

for(i in 1:2){
  plot.title = plot.titles.list[i]
  if(i ==1){df = Dfs.spots.IR
  } else if(i ==2){df = Dfs.spots.UV
  }
  plots <- lapply(df, 
                  apply(
                    lapply(feature=as.character(graph.features.spot[,2]), make.plots.spot, save.image=FALSE, plot.title=plot.title, graph.features=graph.features.spot),
                    grid.arrange()))
  return(plots)
}


##spots IR:
if(plot.title == "(IR treated)"){
  ##IR:
  for (i in 1:length(ch3.Dfs.spots.IR)){
    this.df <- ch3.Dfs.spots.IR[[i]]
    #yl <- yLimitsIR[i]
    make.plots.spot(save.image, df=this.df, plot.title, graph.features.spot, as.character(graph.features.spot[1,2]))
  }                                                            
} else if (plot.title == "(UV treated)"){
  ##UV:
  for (i in 1:length(ch3.Dfs.spots.UV)){
    this.df <- ch3.Dfs.spots.UV[[i]]
    #yl <- yLimitsUV[i]
    make.plots.spot(save.image, df=this.df, plot.title, graph.features.spot, as.character(graph.features.spot[6,2]))
  }
}



lapply(Dfs.spots.IR, make.plots.spot, save.image=FALSE, plot.title, feature=as.character(graph.features.spot[1,2] ))

dfs.spot[1]
setwd('./Test/')  

## Putting into 3 x 3 grid (unfinished)

# Need to pass this function a single dataframe

# if( as.character(graph.features$channel[graph.features$measurement == as.character(feature)]) == "Ch2"){
#   channel.label <- df[[1,"Ch2"]]
# } else{
#   channel.label <- df[[1,"Ch3"]]
# } 
# Graph title 
#graph.title = paste( " (", channel.label, ") ", plot.title, sep = "")
# Spot:       4 in 2x2 grid 
# Intensity:  3 in 2 x 2 grid (Nuclei count(1) + ChX(2))
# Combine both to give 4 x 2 grid
# Need to extract S9.6 and Nucleolin out for separate processing

# Plot for generating plots for multi-plotting (title and axes changed)
Make.multi.plots = function(save.image, df, plot.title, graph.features, feature){
  # Use measurement from dataframe of features
  y_Axis = as.character(graph.features$title[graph.features$measurement == as.character(feature)])
  # Graph title 
  graph.title = y_Axis
  
  # Colour of graph based on IR or UV
  if(plot.title=="(IR treated)"){
    col.pal = brewer.pal(n=9, "Oranges")[4:8]
  } else if(plot.title=="(UV treated)"){
    col.pal = brewer.pal(n=9, "Purples")[4:8]
  } else{ print("no plot.title")}
  
  # Plot function, takes dataframe, plots siRNA (plate order) on x axis by treatment time, against measurement on y axis
  working.plot = ggplot(df,
                        aes(x = reorder(siRNA, Order), fill = reorder(treatment, triplicate.number), y = df[[feature]])) + 
    # labels
    labs(y = y_Axis, fill = "Plate Treatment", title = graph.title) + #S9.6 #pRPA2 S33    
    # graph type, stat+ position plot treatments alongside one another
    geom_bar(stat="identity", position = "dodge") + 
    # colour of bars
    scale_fill_manual(values = col.pal)
  
  # SEM error bars
  if(feature != "ValidObjectCount.median"){
    working.plot <- working.plot + geom_errorbar(aes(x = reorder(siRNA, Order), ymin = (df[[feature]]-df[[gsub("MEAN", "SE", feature)]]), ymax = (df[[feature]]+df[[gsub("MEAN", "SE", feature)]])), width = .2, position = position_dodge(.9)) 
    
  }
  
  # text orientation
  working.plot <- working.plot + theme(text = element_text(size = 7), 
                                       axis.text.x = element_text(size = 7,
                                                                  angle = 45, vjust = 1, hjust=1),
                                       axis.title.x = element_blank(),
                                       axis.text.y = element_text(size = 7),
                                       axis.title.y = element_text(face = "bold"),
                                       plot.title = element_text(size = 8, hjust = 0),
                                       aspect.ratio = 4/3,
                                       legend.text = element_text(size = 6),
                                       legend.key.size = unit(0.5, "line"),
                                       legend.box.margin = margin(0,0,0,-20),
                                       plot.margin = margin(0,0,0,0))
  
  return(working.plot)
  
  ###Saving in function
  #if(save.image){
  #  ggsave(
  #    paste(gsub("/", "x", graph.title), ".png", sep=""),
  #    working.plot,
  #    width = 6,
  #    height = 6,
  #    scale = 1,
  #    dpi = 180
  #  )
  #    dev.off()
  #    print("saved")
  #    print(working.plot)
  
  #file_name = paste(gsub("/", "x", graph.title), ".png", sep="")
  #png(file_name, width = 960, height = 960, bg = "transparent", pointsize = 12)
  #print(working.plot)
  
  
}

rearrange.ch2.2by4 <- function (df.intensity, df.spot, plot.title, graph.features.intensity, graph.features.spot) {
  plot.list.ch2 <- list()
  plot.list.ch3 <- list()
  #if((nrow(df.int) > 0)){
  # Intensity Ch2 measurements
  for(i in 1:3){
    j = i
    plot.list.ch2[[j]] <- Make.violin.plots(df.intensity, save.image=FALSE, plot.title, graph.features.intensity, feature= as.character(graph.features.intensity[i,2]))
  }
  # Spot Ch2 measurements
  if((nrow(df.spot) > 0) & !is.na(df.intensity$Ch2[1])){
    if(df.intensity$Ch2[1] == df.spot$Ch2[1]){
      for(i in 2:5){
        j = i+2
        plot.list.ch2[[j]] <- Make.violin.plots(df.spot, save.image=FALSE, plot.title, graph.features.spot, feature= as.character(graph.features.spot[i,2]))
      }
    }
  }
  g.ch2 <- arrangeGrob(grobs= plot.list.ch2, ncol= 4, top = as.character(df.intensity[,Ch2][1]))
  ggsave(file=paste(gsub("/", "x", as.character(df.intensity[,Ch2][1])), "(Ch2)", plot.title, ".png", sep = " "), 
         g.ch2,
         width = 20,
         height = 10,
         scale = 0.75,
         dpi = 180)
  # Intensity Ch3 measurements
  for(i in 1:3){
    j = i
    plot.list.ch3[[j]] <- Make.violin.plots(df.intensity, save.image=FALSE, plot.title, graph.features.intensity, feature= as.character(graph.features.intensity[c(1,4,5),2][i]))
  }
  # Spot Ch3 measurements
  if((nrow(df.spot) > 0) & !is.na(df.intensity$Ch3[1])){  
    if(df.intensity$Ch3[1] == df.spot$Ch3[1]){
      for(i in 6:9){
        j = i-2
        plot.list.ch3[[j]] <- Make.violin.plots(df.spot, save.image=FALSE, plot.title, graph.features.spot, feature= as.character(graph.features.spot[i,2]))
      }
    }
  }
  #grid.arrange(grobs= plot.list, ncol= 3, top = "")
  g.ch3 <- arrangeGrob(grobs= plot.list.ch3, ncol= 4, top = as.character(df.intensity[,Ch3][1]))
  ggsave(file=paste(gsub("/", "x", as.character(df.intensity[,Ch3][1])), "(Ch3)", plot.title, ".png", sep = " "), 
         g.ch3,
         width = 20,
         height = 10,
         scale = 0.75,
         dpi = 180)
}  

rearrange.ch2.2by4(Dfs.intensity.IR[[2]], Dfs.spots.IR[[2]], plot.titles.list[1], graph.cells.intensity, graph.cells.spots)

#####


for(i in 1:1){
  plot.title = plot.titles.list[i]
  if(i ==1){
    df.int = Dfs.intensity.IR
    df.spot = Dfs.spots.IR
  } else if(i ==2){
    df.int = Dfs.intensity.UV
    df.spot = Dfs.spots.UV
  }
  for(j in 1:1){
    if(j %in% names(df.int)){
      rearrange.ch2.2by4(df.int[[as.character(j)]], df.spot[[as.character(j)]], plot.title, graph.features.intensity, graph.features.spot)
    }
  }
}

