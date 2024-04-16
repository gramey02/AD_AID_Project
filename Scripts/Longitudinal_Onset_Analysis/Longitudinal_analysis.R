#### Script to perform temporal analysis on retrospective cohort

#load libraries
library(dplyr)
library(MatchIt)
library(hash)
library(tibble)
library(epitools)
library(ggplot2)
library(ggsci)

#load data for the longitudinal study group
load("M:/Matched_data/matched_rc_temporal_2023_08_04_m1.Rdata")
df<-rc_temporal #rename the resulting data frame

#load keywords dataframe with search terms and disorder groupings
keywords<-readxl::read_xlsx("~/EHR/keywords_2022_10_26.xlsx")

#summary statistics for plotting
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#create new column aggregating the sex and autoimmunity (aid) status variables
df<-df %>% add_column(sex_aid = paste(df$gender, df$aid, sep=", "))
df$sex_aid<-as.factor(df$sex_aid)

#plot the aid status against age of AD onset
p<- ggplot(df, aes(x=as.factor(aid), y=aaAD_years, fill=as.factor(aid))) + 
  geom_violin()+
  stat_summary(fun.data=data_summary, col="black")+
  theme_minimal()+
  ylab("Age at AD diagnosis (years)")+
  xlab("Autoimmune Disease Status")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size=16))+
  scale_fill_manual(values=c("lightgrey", "#009193"))+
  scale_y_log10(limits=c(NA, 115))+
  scale_x_discrete(labels=c(paste("Control\nN=", length((df %>% filter(aid==0))$person_id), sep=""),
                            paste("Autoimmune\nN=", length((df %>% filter(aid==1))$person_id), sep="")))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'))+
  geom_segment(aes(x=1,xend=2,y=95,yend=95))+
  annotate("text", x = 1.5, y = 97, label = "N.S.", size=3.5)
p
#save the plot
ggsave(filename="M:/Figures/temporal_overall_m1.png", plot=p, device='png', dpi=600, width=5, height=5)
ggsave(filename="~/EHR/toBox/temporal_overall_m1.png", plot=p, device='png', dpi=600, width=5, height=5)

#print the mean age at AD diagnosis for these stratifications
agg_tbl1 <- df %>% group_by(aid) %>% 
  summarise(mean_aaAD_years=mean(aaAD_years),
            .groups = 'drop')
print(agg_tbl1)

#run mann-whitney u test to see if difference is signficant
res<-wilcox.test(aaAD_years~aid, data = df, exact=FALSE)
res
res$p.value

#get the mean age at AD diagnosis by autoimmunity status and sex
agg_tbl2 <- df %>% group_by(sex_aid) %>% 
  summarise(mean_aaAD_years=mean(aaAD_years),
            .groups = 'drop')
print(agg_tbl2)

#run a mann-whitney u test to see if the difference between the distributions is significant
#get dataframes for intra- and inter-gender comparisons
df_female <-df %>% filter(gender=="Female")
df_male <-df %>% filter(gender=="Male")
df_femaleAID_maleAID<-df %>% filter(aid==1)
df_female_male<-df %>% filter(aid==0)

#within females
res1<-wilcox.test(aaAD_years~aid, data = df_female, exact=FALSE)
res1
res1$p.value
#within males
res2<-wilcox.test(aaAD_years~aid, data = df_male, exact=FALSE)
res2
res2$p.value
#between males with autoimmunity (AID) and females with AID
res3<-wilcox.test(aaAD_years~gender, data = df_femaleAID_maleAID, exact=FALSE)
res3
res3$p.value #interesting, women with autoimmunity have a younger onset age than men with autoimmunity
#between control males and females
res4<-wilcox.test(aaAD_years~gender, data = df_female_male, exact=FALSE)
res4
res4$p.value #interesting, women without autoimmunity have a younger onset age than men without autoimmunity

#plot the distributions
p<- ggplot(df, aes(x=as.factor(sex_aid), y=aaAD_years, fill=as.factor(sex_aid))) + 
  geom_violin()+
  stat_summary(fun.data=data_summary, col="black")+
  theme_minimal()+
  ylab("Age at AD diagnosis")+
  xlab("Autoimmune Disease Status")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size=16))+
  scale_fill_manual(values=c("#FFE999", "#FED439", "#CDDBF4", "#709AE1"))+
  scale_y_log10(limits=c(NA, 120))+
  scale_x_discrete(labels=c("Control,\nFemale", "Autoimmune,\nFemale", "Control,\nMale", "Autoimmune,\nMale"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'))+
  geom_segment(aes(x=1,xend=3,y=95,yend=95))+
  annotate("text", x = 2, y = 99, label = paste("** p = ", sprintf(as.numeric(res4$p.value),
                                                  fmt = '%#.2f'), sep=""), size=5)+
  geom_segment(aes(x=2,xend=4,y=103,yend=103))+
  annotate("text", x = 3, y = 107, label = paste("* p = ", sprintf(as.numeric(res3$p.value),
                                                                   fmt = '%#.2f'), sep=""), size=5)
p
#save the plot
ggsave(filename="M:/Figures/temporal_sex_m1.png", plot=p, device='png', dpi=600, width=7, height=5)
ggsave(filename="~/EHR/toBox/temporal_sex_m1.png", plot=p, device='png', dpi=600, width=7, height=5)


###---------break up the people who have more than one autoimmune disease (i.e. the same person will show up twice but with two different autoimmune diseases)
dups<-NULL
#create a lookup table to get the corresponding column name for each autoimmune disease name
h<-hash()
for(i in 1:length(keywords$AID_name)){
  h[[keywords$AID_name[i]]]<-keywords$Col_names[i]
}
h2<-hash()
for(i in 1:length(keywords$AID_name)){
  h2[[keywords$AID_name[i]]]<-keywords$Group[i]
}
#sort the dataframe so that the autoimmune disease case comes before the control
df<-df[order(df$subclass,desc(df$aid )),]
#split the patients with multiple autoimmune diseases into separate rows (so they will be duplicated)
#--this includes the controls with multiple autoimmune diseases listed
df_split<-cSplit(df, "all_AIDs", direction = "long")
#now remove the columns for each autoimmune disease and disease group, and reassign based on the split
df_split<-select(df_split, -keywords$Col_names)
df_split<-select(df_split, -unique(keywords$Group))
#add columns for individual diseases
n_col<-ncol(df_split)
for(i in 1:length(keywords$AID_name)){
  df_split<-df_split %>% add_column(temp=0)
  colnames(df_split)[colnames(df_split) == 'temp'] <- keywords$Col_names[i]
}
#add 1s to the aid columns depending on the aid names for the person
for(i in 1:length(df_split$person_id)){
  for(j in 1:length(keywords$AID_name)){
    if(grepl(df_split$all_AIDs[i], pattern=keywords$AID_name[j])){
      df_split[[n_col+j]][i]<-1
      }
    }
}
#add columns for disease groups
col_n<- ncol(df_split) #get current number of columns in df
#add a column for each disease group
for(i in 1:length(unique(keywords$Group))){
  df_split<-df_split %>% add_column(temp = 0)
  colnames(df_split)[colnames(df_split) == 'temp'] <- unique(keywords$Group)[i]
}
#create a hash to lookup the group for each AID name
group_hash<-hash()
for(i in 1:length(keywords$AID_name)){
  group_hash[[keywords$AID_name[i]]]<-keywords$Group[i]
}
#now assign each person's autoimmune diseases to disease groups
for(i in 1:length(df_split$all_AIDs)){
  group_list<-c() #vector that will hold disease group names for each person
  for(j in 1:length(keywords$AID_name)){
    if(grepl(df_split$all_AIDs[i], pattern=keywords$AID_name[j])){
      group_name<-group_hash[[keywords$AID_name[j]]]
      df_split[[group_name]][i] <- 1
      group_list<-c(group_list, group_name)
    }
  }
  group_list<-unique(group_list)
  df_split$disease_groups[i]<-paste(group_list, collapse=", ")
}

#awesome! Now patients with multiple autoimmune diseases are split up so we can count them each
#separately for different disease group and specific disease analyses

####Disease Group Distributional Analysis-------------------------------
#first look at the different distributions by disease group:

#remove the other columns to make sure you don't have a single patient represented more than once within a certain disease group
#(since a person can have multiple autoimmune diseases, each of those could all belong to one disease subtype group)
df_split_group<-select(df_split, -keywords$Col_names)
df_split_group<-select(df_split_group, -c(all_AIDs, earliest_aid, aid_dates, earliest_aid_date, numAIDs))
df_split_group<-distinct(df_split_group)

#now plot
df_split_group$group_aid <- paste(df_split_group$disease_groups, df$aid, sep=", ")
df_split_group$group_aid<-as.factor(df_split_group$group_aid)
p<- ggplot(df_split_group, aes(x=group_aid, y=aaAD_years, group=group_aid)) + 
  geom_violin()+
  stat_summary(fun.data=data_summary, col="black", position=position_dodge(0.9))+
  theme_minimal()+
  ylab("AD diagnosis age")+
  xlab("Autoimmune Disease Group")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size=16))+
  scale_fill_manual(values=c("lightgrey", "#009193"))+
  scale_y_log10(limits=c(NA, NA))+
  scale_x_discrete()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=9, face='bold', angle=45, vjust = 1.2, hjust=1.1),
        axis.title.y = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'))
  
p
#save the plot
ggsave(filename="M:/Figures/temporal_dg_overall_m1.png", plot=p, device='png', dpi=600, width=12, height=5)
ggsave(filename="~/EHR/toBox/temporal_dg_overall_m1.png", plot=p, device='png', dpi=600, width=12, height=5)

#get the mean age at AD diagnosis by disease group
agg_tbl5 <- df_split_group %>% group_by(group_aid) %>% 
  summarise(mean_aaAD_years=mean(aaAD_years),
            .groups = 'drop')
print(agg_tbl5)
#run a mann-whitney u test to see if the difference between the distributions is significant
for(i in 1:length(unique(keywords$Group))){
  print(unique(keywords$Group)[i])
  df_filt<-df_split_group %>% filter(disease_groups==unique(keywords$Group)[i])
  res<-wilcox.test(aaAD_years~aid, data = df_filt, exact=FALSE)
  print(res)
  print(res$p.value)
  print("----------------------------")
}



####Specific Disease Distributional Analysis---------------------------------------
df_split_aid<-df_split
df_split_aid$aid_status <- paste(df_split_aid$all_AIDs, df_split_aid$aid, sep=", ")
df_split_aid$aid_status <- as.factor(df_split_aid$aid_status)
#now plot
p<- ggplot(df_split_aid, aes(x=aid_status, y=aaAD_years, group = aid_status)) + 
  geom_violin()+
  stat_summary(fun.data=data_summary, col="black", position=position_dodge(0.9))+
  theme_minimal()+
  ylab("AD diagnosis age")+
  xlab("Autoimmune Disease")+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1, size=16))+
  scale_fill_simpsons()+
  scale_y_log10(limits=c(NA, NA))+
  scale_x_discrete()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=9, face='bold', angle=90, vjust = 1.2, hjust=0.8),
        axis.title.y = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'))
p
#get the mean age at AD diagnosis by disease group
agg_tbl6 <- df_split_aid %>% group_by(aid_status) %>% 
  summarise(mean_aaAD_years=mean(aaAD_years),
            .groups = 'drop')
print(n = 100, agg_tbl6)
#run a mann-whitney u test to see if the difference between the distributions is significant
for(i in 1:length(unique(keywords$AID_name))){
  print(unique(keywords$AID_name)[i])
  df_filt<-df_split_aid %>% filter(all_AIDs==unique(keywords$AID_name)[i])
  if(length(df_filt$person_id)>2){
    res<-wilcox.test(aaAD_years~aid, data = df_filt, exact=FALSE)
    print(res)
    print(res$p.value)
    print("----------------------------")
  }
}





####---------------------------------------hazard ratio analysis
library(survival)
library(survminer)
df$gender<-as.factor(df$gender)
df$aid<-as.factor(as.character(df$aid))
#first look at aid as the only variable
res.cox <- coxph(Surv(aaAD_years, alz) ~ aid, data = df)
summary(res.cox)
#now work in gender
res.cox2<-coxph(Surv(aaAD_years, alz) ~ gender + aid, data = df)
summary(res.cox2)
#plot the results
grid.draw.ggsurvplot <- function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}
fit <- survfit(Surv(aaAD_years, alz) ~ aid+gender, data = df)
p<-ggsurvplot(fit,
              xlim = c(50,100),
              break.x.by = 5,
              pval = TRUE, conf.int = FALSE,
              legend.labs = c("Control,\nFemale", "Control,\nMale", "Autoimmune,\nFemale", "Autoimmune,\nMale"),
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              #linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw(), # Change ggplot2 theme
              palette = c("#8F7000", "#214C97", "#FED439", "#709AE1"),
              risk.table.height = 0.2,
              font.main = c(18, "bold", "black"),
              font.x = c(18, "bold", "black"),
              font.y = c(18, "bold", "black"),
              font.tickslab = c(16, "plain", "black"),
              fontsize=6)+
  xlab("Time (days)")+ylab("AD probability")
  
p

#save the plot
ggsave(filename="M:/Figures/temporal_sex_survival_overall_m1.png", p, height = 6, width = 16)
ggsave(filename="~/EHR/toBox/temporal_sex_survival_overall_m1.png", p, height=6, width=16)


#compare within the female group
df_f<-df %>% filter(gender=="Female")
res.cox3<-coxph(Surv(aaAD_years, alz) ~ aid, data=df_f)
summary(res.cox3)
fit <- survfit(Surv(aaAD_years, alz) ~ aid, data = df_f)
p<-ggsurvplot(fit,
              xlim = c(50,92),
              break.x.by = 5,
              pval = TRUE, conf.int = FALSE,
              legend.labs = c("Control,\nFemale", "Autoimmune,\nFemale"),
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              #linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw(), # Change ggplot2 theme
              palette = c("#8F7000", "#FED439"),
              risk.table.height = 0.2,
              font.main = c(18, "bold", "black"),
              font.x = c(18, "bold", "black"),
              font.y = c(18, "bold", "black"),
              font.tickslab = c(16, "plain", "black"),
              fontsize=6)+
  xlab("Time (days)")+ylab("AD probability")
p

#compare within the male group
df_m<-df %>% filter(gender=="Male")
res.cox4<-coxph(Surv(aaAD_years, alz) ~ aid, data=df_m)
summary(res.cox4)
fit <- survfit(Surv(aaAD_years, alz) ~ aid, data = df_m)
p<-ggsurvplot(fit,
              xlim = c(50,92),
              break.x.by = 5,
              pval = TRUE, conf.int = FALSE,
              legend.labs = c("Control,\nMale", "Autoimmune,\nMale"),
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              #linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw(), # Change ggplot2 theme
              palette = c("#214C97", "#709AE1"),
              risk.table.height = 0.2,
              font.main = c(18, "bold", "black"),
              font.x = c(18, "bold", "black"),
              font.y = c(18, "bold", "black"),
              font.tickslab = c(16, "plain", "black"),
              fontsize=6)+
  xlab("Time (days)")+ylab("AD probability")
p

#finally, to verify that autoimmunity alone is enough to separate rates of AD diagnosis
res.cox5<-coxph(Surv(aaAD_years, alz) ~ aid, data=df)
summary(res.cox5)
fit <- survfit(Surv(aaAD_years, alz) ~ aid, data = df)
p<-ggsurvplot(fit,
              xlim = c(50,92),
              break.x.by = 5,
              pval = TRUE, conf.int = FALSE,
              legend.labs = c("Autoimmune", "Control"),
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              #linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw(), # Change ggplot2 theme
              palette = c("blue", "red"),
              risk.table.height = 0.2,
              font.main = c(18, "bold", "black"),
              font.x = c(18, "bold", "black"),
              font.y = c(18, "bold", "black"),
              font.tickslab = c(16, "plain", "black"),
              fontsize=6)+
  xlab("Time (days)")+ylab("AD probability")
p



####----------------Disease group hazard ratio analysis

#Within disease group comparison
#use df_split_group from before
#first look at aid as the only variable
df_split_group$aid<-as.factor(as.character(df_split_group$aid))
for(i in 1:length(unique(df_split_group$disease_groups))){
  df_filt<-df_split_group %>% filter(disease_groups==unique(df_split_group$disease_groups)[i])
  print(unique(df_split_group$disease_groups)[i])
  res.cox <- coxph(Surv(aaAD_years, alz) ~ aid, data = df_filt)
  print(summary(res.cox))
  print("------------------------------")
}
res.cox <- coxph(Surv(aaAD_years, alz) ~ aid + disease_groups, data = df_split_group)
summary(res.cox)
res.cox2<-coxph(Surv(aaAD_years, alz) ~ gender + aid, data = df_split_group)
summary(res.cox2)

##Across disease group comparison
df_filt<-df_split_group %>% filter(aid==1)
res.cox <- coxph(Surv(aaAD_years, alz) ~ disease_groups, data = df_filt)
summary(res.cox)
fit <- survfit(Surv(aaAD_years, alz) ~ disease_groups, data = df_filt)
p<-ggsurvplot(fit,
              pval = TRUE, conf.int = TRUE,
              risk.table = TRUE, # Add risk table
              risk.table.col = "strata", # Change risk table color by groups
              #linetype = "strata", # Change line type by groups
              surv.median.line = "hv", # Specify median survival
              ggtheme = theme_bw())+
  xlab("Age (years)")+ylab("AD probability")
p


####------------Specific Disease Hazard Ratio Analysis------------------------
#Within specific comparison
#use df_split_aid from before
#first look at aid as the only variable
df_split_aid$aid<-as.factor(as.character(df_split_aid$aid))
for(i in 1:length(unique(df_split_aid$all_AIDs))){
  df_filt<-df_split_aid %>% filter(all_AIDs==unique(df_split_aid$all_AIDs)[i])
  print(unique(df_split_aid$all_AIDs)[i])
  res.cox <- coxph(Surv(aaAD_years, alz) ~ aid, data = df_filt)
  print(summary(res.cox))
  print("------------------------------")
}