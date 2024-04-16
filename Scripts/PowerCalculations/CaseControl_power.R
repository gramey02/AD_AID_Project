# Power analysis script for case-control study design

#load libraries
library(dplyr)
library(MatchIt)
library(hash)
library(tibble)
library(epitools)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(pwr)

#load case-control study data & disease keywords
load("M:/Matched_data/matched_cc_clean_2023_08_04_m1b.Rdata")
cc<-cc_matched
keywords<-readxl::read_xlsx("~/EHR/keywords_2022_10_26.xlsx")


#construct a row of a dataframe for the power analysis
df<-data.frame(name = NA, #name of analysis
               a = NA, #autoimmune disorder cases who also got AD
               b = NA, #autoimmune disorder cases who did not get AD
               c = NA, #non-autoimmune individuals who also got AD
               d = NA, #non-autoimmune individuals who did not get AD
               OR = NA, #Odds Ratio that was calculated
               p1 = NA, #proportion of exposed individuals who also got AD
               p2 = NA, #proportion of non-exposed individuals who also got AD
               h = NA) #effect size calculated from p1 and p2
#now duplicate that row so we have an empty data frame to fill in
df<-do.call("rbind", replicate( 
  105, df, simplify = FALSE)) #105 because that's the total number of comparisons that we do with this study group


#first get the numbers for the overall cohort
{
  idx=1
  df$name[idx] = "Overall, full group"
  tbl<-table(cc$aid, cc$alz)
  df$a[idx]<-tbl[2,2] #tbl[4] #exposed and with the outcome
  df$b[idx]<-tbl[2,1] #tbl[2] #exposed and without the outcome
  df$c[idx]<-tbl[1,2] #tbl[3] #not exposed and with the outcome
  df$d[idx]<-tbl[1,1] #tbl[1] #not exposed and without the outcome
}

#now for female subgroup
{
  idx=2
  df$name[idx] = "Overall, female group"
  female <- cc %>% filter(gender=="Female")
  tbl<-table(female$aid, female$alz)
  df$a[idx]<-tbl[2,2] #tbl[4] #exposed and with the outcome
  df$b[idx]<-tbl[2,1] #tbl[2] #exposed and without the outcome
  df$c[idx]<-tbl[1,2] #tbl[3] #not exposed and with the outcome
  df$d[idx]<-tbl[1,1] #tbl[1] #not exposed and without the outcome
}

#now for male subgroup
{
  idx=3
  df$name[idx] = "Overall, male group"
  male <- cc %>% filter(gender=="Male")
  tbl<-table(male$aid, male$alz)
  df$a[idx]<-tbl[2,2] #tbl[4]
  df$b[idx]<-tbl[2,1] #tbl[2]
  df$c[idx]<-tbl[1,2] #tbl[3]
  df$d[idx]<-tbl[1,1] #tbl[1]
}


#now run for the disease subtype analyses
group_names<-unique(keywords$Group)
start<-grep(group_names[1], colnames(cc)) #get the column index where the first disease group starts
cur_idx<-idx + 1
for(i in 1:length(group_names)){
  print(i)
  print(group_names[i])
  #filter the dataframe to only include the cur full disease group
  df_filt<-cc %>% filter(cc[,start+i-1]==1)
  group_df<-cc %>% filter(subclass %in% df_filt$subclass) #also include the corresponding controls

  {#in the full group
    df$name[cur_idx]<-paste(group_names[i], " Subtype, full group", sep="")
    tbl<-table(group_df$aid, group_df$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
        print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
      } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  cur_idx<-cur_idx+1
  
  group_df_f<-group_df %>% filter(gender=="Female")
  {#in the female-only group
    df$name[cur_idx]<-paste(group_names[i], " Subtype, female group", sep="")
    tbl<-table(group_df_f$aid, group_df_f$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
      print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  cur_idx<-cur_idx+1
  
  group_df_m<-group_df %>% filter(gender=="Male")
  {#in the male-only group
    df$name[cur_idx]<-paste(group_names[i], " Subtype, male group", sep="")
    tbl<-table(group_df_m$aid, group_df_m$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
      print(paste("check ", group_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  cur_idx<-cur_idx+1
}


#now run for the individual disease analysis
aid_names<-unique(keywords$Col_names)
start<-grep(aid_names[1], colnames(cc))[1] #get the column index where the first individual disease starts
for(i in 1:length(aid_names)){
  #filter the dataframe to only include the cur disease group
  print(i)
  print(aid_names[i])
  df_filt<-cc %>% filter(cc[,start+i-1]==1)
  aid_df<-cc %>% filter(subclass %in% df_filt$subclass)
  
  {#in the full group
    df$name[cur_idx]<-paste(aid_names[i], " , full group", sep="")
    tbl<-table(aid_df$aid, aid_df$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  cur_idx<-cur_idx+1
  
  aid_df_f<- aid_df %>% filter(gender=="Female")
  {#in females
    df$name[cur_idx]<-paste(aid_names[i], " , female group", sep="")
    tbl<-table(aid_df_f$aid, aid_df_f$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  cur_idx<-cur_idx+1
  
  aid_df_m<-aid_df %>% filter(gender=="Male")
  {#in males
    df$name[cur_idx]<-paste(aid_names[i], " , male group", sep="")
    tbl<-table(aid_df_m$aid, aid_df_m$alz)
    print(tbl)
    #include some checks to make sure the table readout is conducive
    if(ncol(tbl)==1){
      df$d[cur_idx]<-tbl[2,1]
      df$b[cur_idx]<-tbl[1,1]
    } else if(nrow(tbl)<2){#skip if the table is not full
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else if(ncol(tbl)<2){
      print(paste("check ", aid_names[i], " table, idx: ", as.character(cur_idx), sep=""))
    } else {
      #fill in the dataframe
      df$a[cur_idx]<-tbl[2,2] #tbl[4]
      df$b[cur_idx]<-tbl[2,1] #tbl[2]
      df$c[cur_idx]<-tbl[1,2] #tbl[3]
      df$d[cur_idx]<-tbl[1,1] #tbl[1]
    }
  }
  
  cur_idx<-cur_idx+1
}


#recalculate the odds ratios
df$OR<-(df$a*df$d)/(df$b*df$c)
#calculate proportion 1 (i.e. the proportion of exposed cases that got the outcome)
df$p1<-df$a/(df$a + df$b)
#calculate proportion 2 (i.e. the proportion of non-exposed cases that got the outcome)
df$p2<-df$c/(df$c + df$d)
#calculate the effect size detectable with these two proportions (from two proportions comparison, see pwr package)
df$h<-(2*asin(sqrt(df$p1)))-(2*asin(sqrt(df$p2)))

#calculate the power you have to detect this difference
#the number of exposed and non-exposed cases is not equal here, so use the pwr.2p2n test
#note here that we are fixing the number of controls (n1) based on what we observe in our study group to get the expected number of cases
df<-df %>% add_column(power=NA)
for(i in 1:length(df$name)){
  #test will only run if each disease has a sample size of at least 2
  cur_h<-df$h[i]
  if(is.na(df$h[i])==FALSE){
    if((df$a[i] + df$b[i])>=2 && (df$c[i] + df$d[i])>=2){
      cur_power<-(pwr.2p2n.test(h=ES.h(df$p1[i], df$p2[i]), n2 = df$c[i]+df$d[i], n1 = df$a[i] + df$b[i], sig.level=0.05))[['power']]
      df$power[i]<-cur_power
    }
  }
}

#finally, fill remaining NA values in the a,b,c,d columns with zeroes
df<- df %>% mutate(a = ifelse(is.na(a), 0, a))
df<- df %>% mutate(b = ifelse(is.na(b), 0, b))
df<- df %>% mutate(c = ifelse(is.na(c), 0, c))
df<- df %>% mutate(d = ifelse(is.na(d), 0, d))

###change this code here to the proper file location
#save the resulting power table
write.csv(df, file="~/EHR/toBox/cc_power_stats.csv") 
