library(here)
library(readxl)
library(dplyr)
library(ggplot2)



quantstudio <- function(name, output, sheet = "Results", skip = 46,
                        housekeeper = "GAPDH", cont.sample = NULL) {
  

  # read in the data frame
  df <- read_excel(here(paste0(name,".xls")), sheet = sheet,
                   skip = 46, n_max = 96)
  
  # rename the first few wells, as these will be the important columns
  df <- df[,c("Well Position", "Sample Name", "Target Name", "CT", "Ct Mean", "Ct SD")]
  
  # subset only the wells we require: those renamed above, by getting the distinct rows
  # this assumes that the excel sheet was properly labelled prior to this step (either on QS5 software or in excel)
  df.sum <- distinct(df[,-c(1,4)])
  
  # will remove all entries where the sd of Ct is >0.9, which is a good indication
  # that the technical replicates are not accurate
  
  # loop iterates through rows to find those Ct values where Ct SD is > 0.9 and converts those rows to NA
  for (i in 1:nrow(df.sum)) {
    if (df.sum[i,4]>0.9) {
      df.sum[i,] <- NA
    } 
  }
  
  # then remove those rows from the data frame
  df.sum <- df.sum[!is.na(df.sum[,1]),]
  
  # write this file as a .csv
  write.csv(df.sum, here(paste0(output, "/", name, "tech_rep_filter.csv")))
  
  # loop iterates through rows to find the Ct value of the assigned housekeeper, and identify those with Ct > 22
  # then adds the sample name to the variable "weird"
  weird <- as.vector(NULL)
  for(i in 1:nrow(df.sum)) {
    if (df.sum[i,1] == housekeeper & df.sum[i,3] > 22) {
       weird <- c(weird, as.character(df.sum[i,2]))
    }
  }
  
  # loop iterates through rows to find all samples identified as having bad housekeeper values and eliminating them
  for (i in 1:nrow(df.sum)) {
    for(j in length(weird)) {
      if (df.sum[i,2] == weird[j]){
        df.sum[i,] <- NA
      }
    } 
    
  }
  
  df.sum <- df.sum[!is.na(df.sum[,1]),]
  
  # write this file as a .csv
  write.csv(df.sum, here(paste0(output, "/", name, "housekeeper_filter.csv")))
  
  # identifies all sample names and assigns to vector "samples"
  samples <- as.vector(df.sum %>% 
                         group_by(`Target Name`) %>% 
                         summarise(n_distinct(`Sample Name`)))$`Target Name`

  
  # iterates through the df and adds a tibble df of each sample to a list
  sum.list <- as.list(NULL)
  temp <- as.data.frame(NULL)
  for (i in 1:length(samples)) {
    
    for(j in 1:nrow(df.sum)) {
      if (df.sum[j, 2] == samples[i]) {
        temp <- rbind(temp, df.sum[j,])
      }
      sum.list[[samples[i]]] <- temp
    }
    temp <- as.data.frame(NULL)
  }
  
  # Now generate the Dct value for each gene, and then take the 2^-Dct to generate normalised expression
  sum.list.dct <- sum.list
  for (i in 1:length(sum.list)) {
    
    for (j in 1:nrow(sum.list[[i]])) {
      
      sum.list.dct[[i]][j,"DCT"] <- sum.list[[i]][j,3] - sum.list[[i]][sum.list[[i]][,1]==housekeeper,3]
      sum.list.dct[[i]][j,"NormExp"] <- 2^-(sum.list.dct[[i]][j,"DCT"])
      print(sum.list.dct[[i]][j,"NormExp"])
    }
  }
  
  # removes control (GAPDH) from dataset for further analysis
  sum.list.dct <- lapply(sum.list.dct, function(x) x[x$`Sample Name` != housekeeper,])
  
  #generates df containing all current information inc. dct
  df.dct <- as.data.frame(NULL)
  for (i in 1:length(sum.list.dct)) {
    for (j in 1:nrow(sum.list.dct[[i]])) {
      df.dct <- rbind(df.dct, sum.list.dct[[i]][j,])
    }
  }
  write.csv(df.dct, here(paste0(output,name,"Normalised.csv")))
  
  # print graphs
  g2 <- ggplot(df.ddct, aes(`Target Name`, NormExp, fill = `Sample Name`)) +
    geom_col() +
    facet_grid(`Sample Name` ~ ., scales = "free_y")
  jpeg(paste0(output, "/", name,"NormExp_graph.jpg"))
  plot(g2)
  dev.off()
}
