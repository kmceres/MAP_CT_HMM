#!/usr/bin/env Rscript
# usage example: Rscript 01_process_raw_data.R ../data/raw_data/ /Users/kristinaceres/Desktop/CT_HMM/data/

############################################################################################
#  This file processes raw data found in two csv files (fecal.csv, gen_info.csv)
#  The file takes two csv files, and combines them to produce one cleaned csv file to feed
#  into the HMM. 
#
############################################################################################

#combineData takes a path to a folder, and then reads in three dataframes within that folder
#  fecal.csv contains MAP colony forming unit (CFU) counts
#  gen_info.csv contains cow birthdates which is used to calculate age at testing
# It then combines relevant information to form a new dataframe. 
# It also calculates age at testing from birthdate and sampling date
combineData <- function(path){
  setwd(path)
  fecal_data<-read.csv("fecal.csv")
  gen_info<-read.csv("gen_info.csv")
 
   #load package
  library(dplyr) #for cleaning dataframes
  library(tidyr)
  
  #reformat dates 
  fecal_data$SampDate <- as.POSIXct(strptime(fecal_data$SampDate, format="%d/%m/%y"))
  gen_info$BDate<-as.POSIXct(strptime(gen_info$BDate, format="%d.%m.%Y"))
  
  
  #make new df of relevant informatio nfrom gen_inifo (birthdate), and milk_data (305ME)
  BDate <- gen_info$BDate
  JoinedID <- gen_info$JointID
  bday_df <- tibble(JoinedID, BDate)
  
  
  #fixing mistake in data entry, that would make the cows birthday after sampling
  #original sampling date was 1996-01-25, but birthday is 2003-03-16 and the isolate
  # was plated 2006-01-25
  fecal_data$SampDate[fecal_data$JoinedID == "A_1567_03" & fecal_data$SampDate == "1996-01-25"] <- "2006-01-25"
  
  combined_dfs <- left_join(fecal_data, bday_df, by = "JoinedID") %>%
    #remove all other types of result (nt, lost, no sample)
    filter(Result == "Positive" | Result == "Negative") %>%
    mutate(age_at_testing = SampDate - BDate) %>% 
    mutate(CombinedID = paste(FarmID, CowID, sep = "_"))#(CombinedID = paste(FarmID, CowID, format(SampDate, '%Y'), sep="_"))
  
    return(combined_dfs)
}

#remove values other than CFU values in Tube columns 
#set CFU value to lower bound in cases where a range is given 
#(if Tube = 200500, indicating a range from [200,500] CFU, set value to 200)

#bc  = Bacterial contamination /tube was over grown (LOST rather than Negative)
#b#= bacterial growth / MAP colonies also on the tube.
#b - = bacterial contamination/ tube still readable(negative)
#f-  = fungal growth / tube still readable (negative)
#f# = fungal growth present / Map colonies also on the tube
#fc = fungal  contamination/ tube was over grown (LOST rather than Negative)

# A sample is determined to be LOST when all four tubes become unreadable 
# due to contamination overgrowth before the 12 week reading.  
# It is designated as LOST rather than negative because the true Johnes status is undetermined.

fix_characters_totcfu<- function(combined_dfs){
  
  combined_rmv <- combined_dfs
  combined_rmv$Tube.1<-as.numeric(gsub("[^[:digit:]]","",combined_rmv$Tube.1))
  combined_rmv$Tube.2<-as.numeric(gsub("[^[:digit:]]","",combined_rmv$Tube.2))
  combined_rmv$Tube.3<-as.numeric(gsub("[^[:digit:]]","",combined_rmv$Tube.3))
  combined_rmv$Tube.4<-as.numeric(gsub("[^[:digit:]]","",combined_rmv$Tube.4))
  
  combined_rmv$Tube.1[combined_rmv$Tube.1 == 200250]<-200
  combined_rmv$Tube.2[combined_rmv$Tube.2 == 200250]<-200
  combined_rmv$Tube.3[combined_rmv$Tube.3 == 200250]<-200
  combined_rmv$Tube.4[combined_rmv$Tube.4 == 200250]<-200
  
  combined_rmv$Tube.1[combined_rmv$Tube.1 == 150200]<-150
  combined_rmv$Tube.2[combined_rmv$Tube.2 == 150200]<-150
  combined_rmv$Tube.3[combined_rmv$Tube.3 == 150200]<-150
  combined_rmv$Tube.4[combined_rmv$Tube.4 == 150200]<-150
  
  combined_rmv$Tube.1[combined_rmv$Tube.1 == 100150]<-100
  combined_rmv$Tube.2[combined_rmv$Tube.2 == 100150]<-100
  combined_rmv$Tube.3[combined_rmv$Tube.3 == 100150]<-100
  combined_rmv$Tube.4[combined_rmv$Tube.4 == 100150]<-100
  
  combined_rmv$Tube.1[combined_rmv$Tube.1 == 75100]<-75
  combined_rmv$Tube.2[combined_rmv$Tube.2 == 75100]<-75
  combined_rmv$Tube.3[combined_rmv$Tube.3 == 75100]<-75
  combined_rmv$Tube.4[combined_rmv$Tube.4 == 75100]<-75
  
  combined_rmv$Tube.1[combined_rmv$Tube.1 == 250300]<-250
  combined_rmv$Tube.2[combined_rmv$Tube.2 == 250300]<-250
  combined_rmv$Tube.3[combined_rmv$Tube.3 == 250300]<-250
  combined_rmv$Tube.4[combined_rmv$Tube.4 == 250300]<-250
  
  #NA values are set to 0 so they are not considered in final CFU count
  combined_rmv$Tube.1[is.na(combined_rmv$Tube.1) == T] <- 0
  combined_rmv$Tube.2[is.na(combined_rmv$Tube.2) == T] <- 0
  combined_rmv$Tube.3[is.na(combined_rmv$Tube.3) == T] <- 0
  combined_rmv$Tube.4[is.na(combined_rmv$Tube.4) == T] <- 0
  
  
  combined_rmv<-combined_rmv %>% mutate(totCFU= (Tube.1+Tube.2+Tube.3+Tube.4)*5.3)
  
  
  return(combined_rmv)
}

add_time <- function(df){
  #arrange data in order of sampling date by cow
  sorted_data <-df %>% arrange(CombinedID, SampDate) 
  
  #find cows with two observations on the same day, set their corrected totCFU count (cor_totCFU)
  #to the average of the two samples 
  corrected_cfu <- sorted_data  %>% 
    group_by(CombinedID, SampDate) %>% 
    mutate(dub = ifelse(n() > 1, "double", "single")) %>%
    group_by(CombinedID, SampDate, dub) %>% mutate(cor_totCFU= as.integer(mean(totCFU))) %>%
    group_by(CombinedID, SampDate, cor_totCFU) %>% distinct(cor_totCFU, .keep_all = T)
  
  #only keep cows that have two positive results in 3 consecutive samplings
    corrected_cfu <- corrected_cfu %>% group_by(CombinedID) %>%
      filter( 
        (Result == "Positive" & lag(Result) == "Positive")| 
          (Result == "Positive" & lead(Result) == "Positive") | 
          (Result == "Positive" & lag(lag(Result) == "Positive")) |
          (Result == "Positive" & lead(lead(Result) == "Positive")) 
      )
  
  ####### Might delete ^ #################################################################
  
  # #only keep cows that have at least 2 observations, add time 
  out <- corrected_cfu %>% group_by(CombinedID) %>% mutate(n = n()) %>% filter(n >= 2) %>%
    mutate(time = difftime((SampDate), (lag(SampDate)), units = "weeks")) %>% 
    mutate(time = round(time)) %>%
    mutate(time = ifelse(time == 0, "NA", time)) %>%
    #multiply by 5.3 to correct for dilution as in pradhan et al 2011
    #add 1.1 to correct for log errors for 0 counts
    mutate(cor_totCFU = (log(cor_totCFU + 1.1))) %>% 
    mutate(OD.Value = as.numeric(as.character(OD.Value))) %>%
    mutate(OD.Value = log(OD.Value)) %>%
    select(CombinedID, FarmID.x, CowID.x, SampDate, Result, cor_totCFU, time, OD.Value) 
  
  return(out)
  
}

get_neg_data = function(combined_rmv, pos_data){
  #arrange data in order of sampling date by cow
  sorted_data <-combined_rmv %>% arrange(CombinedID, SampDate) 
  
  #find cows with two observations on the same day, set their corrected totCFU count (cor_totCFU)
  #to the average of the two samples 
  corrected_cfu <- sorted_data  %>% 
    group_by(CombinedID, SampDate) %>% 
    mutate(dub = ifelse(n() > 1, "double", "single")) %>%
    group_by(CombinedID, SampDate, dub) %>% mutate(cor_totCFU= as.integer(mean(totCFU))) %>%
    group_by(CombinedID, SampDate, cor_totCFU) %>% distinct(cor_totCFU, .keep_all = T)
  
  neg_data <- anti_join(corrected_cfu, pos_data, by = "CombinedID")
  return(neg_data)
}

## main 
main <- function(){
  args <- commandArgs(trailingOnly = T)
  path1 <- args[1]
  path2 <- args[2]
  df1 <- combineData(path1)
  df2<-fix_characters_totcfu(df1)
  df3<-add_time(df2)
  neg <- get_neg_data(df2, df3)
  write.csv(neg, file=paste(path2, "neg.csv", sep = ""))
  write.csv(df3, file = paste(path2,"hmm_ready_data.csv", sep = ""))
  return(df3)
}

main()
