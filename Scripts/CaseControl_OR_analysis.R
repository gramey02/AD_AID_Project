#### Script to run odds ratio analysis on case-control study group

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

#load some helper functions
source('M:/Helper_Functions/get_cc_pwr_stats.R') #function to calculate the power for each level of the analysis
source('~/EHR/get_OR_results.R') #function that takes in the result and metrics data frame, and fills it with values
source('~/EHR/AD_prevalence_calc.R') #function for bootstrapping to obtain cofidence intervals

#Data loading------------------------------------------------------------------------------------
#see https://github.com/gramey02/AD_AID_Project for an example dataframe structure
load("M:/Matched_data/matched_cc_clean_2023_08_04_m1b.Rdata") #main text is matched_cc_clean_2023_08_04_m1b
cc<-cc_matched #rename the dataframe for easier typing below

#load keywords dataframe with search terms
keywords<-readxl::read_xlsx("~/EHR/keywords_2022_10_26.xlsx")

#Overall OR Analysis------------------------------------------------------------------------------
#dataframe to hold odds ratios for different cohort stratifications
cohort_type<-c('Full Cohort', 'Female Cohort', 'Male Cohort')
metrics<-data.frame("Cohort_Type" = cohort_type,
                    "OR" = 0,
                    "lower" = 0,
                    "upper" = 0,
                    "pval" = 0)

#full cohort analysis
result<-oddsratio.fisher(table(cc$aid, cc$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=1) #i=1 for the first overall comparison
cur_table<- table(cc$aid, cc$alz) #now fill in some sample sizes

#within female OR test
female <- cc %>% filter(gender=="Female")
result<-oddsratio.fisher(table(female$aid, female$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=2) #i=2 for the female comparison

#within male OR test
male <- cc %>% filter(gender=="Male")
result<-oddsratio.fisher(table(male$aid, male$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=3) #i=2 for the female comparison

#plot these results in a forest plot
p<-ggplot(data=metrics, aes(x=reorder(Cohort_Type, (OR)), y=(OR), 
                            ymin=(lower), ymax=(upper))) +
  geom_pointrange() + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cohort Type") + ylab("Odds Ratio (95% CI)") +
  theme_bw()+
  ggtitle('Case-Control Analysis\nAD Odds Ratios')+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))+
  scale_y_log10(limits=c(NA, NA))+
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1.1)
#view the plot
p
#calculate sig levels, Bonferroni correction factor = 3, for the 3 comparisons
metrics<- metrics %>% add_column(sig=NA)
for(i in 1:length(metrics$Cohort_Type)){
  if(metrics$pval[i]<=0.001/length(metrics$Cohort_Type)){
    metrics$sig[i]<-"***"
  }
  else if(metrics$pval[i]<=0.01/length(metrics$Cohort_Type)){
    metrics$sig[i]<-"**"
  }
  else if(metrics$pval[i]<=0.05/length(metrics$Cohort_Type)){
    metrics$sig[i]<-"*"
  }
}

#save data frames for later plotting
cc_metrics_overall<-metrics
cc_table_e_o<- table(cc$aid, cc$alz) #2 x 2 exposure x outcome table for everyone, "e", overall "o" (i.e. not stratified by disease/group)
cc_table_f_o<- table(female$aid, female$alz) #table for females, "f"
cc_table_m_o<- table(male$aid, male$alz) #table for males, "m"
save(cc_metrics_overall, cc_table_e_o, cc_table_f_o, cc_table_m_o, 
     file='M:/Plotting_Data/CC_overall_m1b.Rdata') #m1b representing the df for the analysis in the main text, not a sensitivity analysis
#file for main text is file='M:/Plotting_Data/CC_overall_m1b.Rdata'


#Sex-specific prevalence analysis----------------------------------------------
#match individuals based on gender, to ensure the other demographics of females and males are not too different
control_male<-cc %>% filter(aid==0) %>% filter(gender=="Male") #get the control males
control_female<-cc %>% filter(aid==0) %>% filter(gender=="Female") #get the control females
controls<-cc %>% filter(aid==0) #get all controls
colnames(controls)[colnames(controls) == "subclass"] <- "subclass1" #rename subclass columns, as previous name can throw errors
controls<-subset(controls, select = -c(prop.score, weights)) #remove some info related to prior matching

#uncomment and change to own column name for gender below
#controls$___UCSF-specific column name redacted___#<-as.factor(controls$gender) #rename gender column for shorter typing below

#conduct matching of male and female controls
control_m.out<-matchit(gender ~ birth_year + race + ethnicity+aaDeath_years,
                       data = controls, method = "nearest", ratio = 1, verbose=T)
control_m.data<-match.data(control_m.out, distance = "prop.score")

#match male and female autoimmune disorder (aid) patients
aid_female<-cc %>% filter(aid==1) %>% filter(gender=='Female') #get the aid females
aid_male<-cc %>% filter(aid==1) %>% filter(gender=='Male') #get the aid males
aid<-cc %>% filter(aid==1) #get all aid patients
colnames(aid)[colnames(aid) == "subclass"] <- "subclass1" #rename subclass column
aid<-subset(aid, select = -c(prop.score, weights)) #remove some prior matching columns

#uncomment and change to own column name for gender below
#aid$___UCSF-specific column name redacted___#<-as.factor(aid$gender) #rename gender column for shorter typing below

#conduct matching
aid_m.out<-matchit(gender ~ birth_year + race + ethnicity+aaDeath_years,
                       data = aid, method = "nearest", ratio = 1, verbose=T)
aid_m.data<-match.data(aid_m.out, distance = "prop.score")

aid_female<-aid_m.data %>% filter(aid==1) %>% filter(gender=="Female")
control_female<-control_m.data %>% filter(aid==0) %>% filter(gender=="Female")
aid_male<-aid_m.data %>% filter(aid==1) %>% filter(gender=="Male")
control_male<-control_m.data %>% filter(aid==0) %>% filter(gender=="Male")

#calculate AD prevalences in these matched groups
aid_female_prev<-
  length(unique((aid_female %>% filter(alz==1))$person_id))/length(unique(aid_female$person_id))
control_female_prev<-
  length(unique((control_female %>% filter(alz==1))$person_id))/length(unique(control_female$person_id))
aid_male_prev<-
  length(unique((aid_male %>% filter(alz==1))$person_id))/length(unique(aid_male$person_id))
control_male_prev<-
  length(unique((control_male %>% filter(alz==1))$person_id))/length(unique(control_male$person_id))

#bootstrap the dataframes to create a confidence interval for each prevalence calculation
result_list<-AD_prevalence_calc(aid_female, aid_male, control_female, control_male, n=1000)
prev_df<-result_list[[1]] #parse out the resulting confidence intervals for plotting
distributions<-result_list[[2]] #parse out the resulting bootstrapped distributions

#plot these prevalances and confidence intervals
overall_p<-ggplot(data=prev_df)+geom_bar(aes(x=aid_status, y=ad_prev, fill = gender, label = n), 
                                 stat='identity', position=position_dodge(),
                                 color='black')+
  geom_errorbar(aes(x = aid_status, ymin=lower_CI, ymax=upper_CI, fill=gender), width=.2,
                position=position_dodge(0.85), stat="identity")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_text(size=12, face='bold'),
        axis.title.y = element_text(size=12, face='bold'),
        axis.text.y = element_text(size=10, face='bold'),
        legend.text=element_text(face='bold'),
        legend.title=element_text(face="bold"))+
  scale_x_discrete(labels=c("AID", "Control"))+
  scale_y_continuous()+
  labs(fill="Sex")+
  ylab("AD prevalence")+
  scale_fill_simpsons()
overall_p

#now run statistical tests to determine if prevalences are significantly different

#first compare aid females to aid males
total_df<-select(rbind(aid_female %>% add_column(stratification = "aid_female"),
                aid_male %>% add_column(stratification = "aid_male")), 
                c(alz, stratification))
aid_test<-fisher.test(table(total_df$stratification, total_df$alz))
aid_test

#now compare control females to control males
total_df<-select(rbind(control_female %>% add_column(stratification = "control_female"),
                       control_male %>% add_column(stratification = "control_male")), 
                 c(alz, stratification))
control_test<-fisher.test(table(total_df$stratification, total_df$alz))
control_test

#compare aid females to control females
total_df<-select(rbind(aid_female %>% add_column(stratification = "aid_female"),
                       control_female %>% add_column(stratification = "control_female")), 
                 c(alz, stratification))
female_test<-fisher.test(table(total_df$stratification, total_df$alz))
female_test

#compare aid males to control males
total_df<-select(rbind(aid_male %>% add_column(stratification = "aid_male"),
                       control_male %>% add_column(stratification = "control_male")), 
                 c(alz, stratification))
male_test<-fisher.test(table(total_df$stratification, total_df$alz))
male_test

#great, now we have the significance values for the prevalence plots! save these and the prevalence dataframe for plotting
cc_prev_sig<-data.frame(comparison = c("aid female vs aid male", "control female vs control male", "aid female vs control female", "aid male vs control male"),
                        pval = c(aid_test$p.value, control_test$p.value, female_test$p.value, male_test$p.value),
                        pval_adj = 4 * c(aid_test$p.value, control_test$p.value, female_test$p.value, male_test$p.value))
#also save some sample sizes from the prevalence calculations
cc_prev_ss<-data.frame(aid = length(aid_female$person_id), control = length(control_female$person_id))
cc_prev_df<-prev_df #rename to be more specific
#save
save(cc_prev_ss, cc_prev_sig, cc_prev_df,
     file='M:/Plotting_Data/CC_prev_m1b.Rdata') #CC_prev_m1b representing the df for the analysis in the main text, not a sensitivity analysis


###Autoimmune Disorder Subtype/Group Odds Ratio Analysis-----------------------------------------------------###
group_names<-unique(keywords$Group) #get a list of disorder subtype names
metrics<-data.frame("group_names" = group_names,
                    "OR" = NA,
                    "upper_CI" = NA,
                    "lower_CI" = NA,
                    "pvals" = NA,
                    "N" = NA)
#now calculate odds ratios for each group
start<-grep(group_names[1], colnames(cc)) #get the column index where the first disease group starts
cur_idx<-length(pr$name) + 1
for(i in 1:length(metrics$group_name)){
  #filter the dataframe to only include the current disease group
  print(i)
  df_filt<-cc %>% filter(cc[,start+i-1]==1)
  group_df<-cc %>% filter(subclass %in% df_filt$subclass) #get both the cases and controls for this specific group
  metrics$N[i]<-length(group_df$person_id) #so this N includes autoimmune cases and non-autoimmune controls
  print(group_names[i])
  print(length(group_df$person_id))
  print(table(group_df$aid, group_df$alz))
  print(table(group_df$alz, group_df$gender))
  result<-oddsratio.fisher(table(group_df$aid, group_df$alz))
  metrics$OR[i]<-result$measure[2]
  metrics$lower_CI[i]<-result$measure[4]
  metrics$upper_CI[i]<-result$measure[6]
  metrics$pvals[i]<-result$p.value[4]
}

#compute some significance metrics, with Bonferroni correction factor = 8
metrics<-metrics %>% add_column(sig=NA)
for(i in 1:length(metrics$group_names)){
  if(metrics$pvals[i]<=0.001/length(metrics$group_names)){
    metrics$sig[i]<-"***"
  }
  else if(metrics$pvals[i]<=0.01/length(metrics$group_names)){
    metrics$sig[i]<-"**"
  }
  else if(metrics$pvals[i]<=0.05/length(metrics$group_names)){
    metrics$sig[i]<-"*"
  }
}


#run same comparison but for the female-specific subsets of each disorder subtype
metrics_female<-data.frame("group_names" = group_names,
                           "OR" = NA,
                           "upper_CI" = NA,
                           "lower_CI" = NA,
                           "pvals" = NA,
                           "N" = NA)
for(i in 1:length(metrics_female$group_name)){
  #filter the dataframe to only include the current disease group
  df_filt<-cc %>% filter(cc[,start+i-1]==1) %>% filter(gender=="Female")
  group_df<-cc %>% filter(subclass %in% df_filt$subclass) %>% filter(gender=="Female") #get both the cases and controls for this specific group
  metrics_female$N[i]<-length(group_df$person_id)
  print(group_names[i])
  print(table(group_df$aid, group_df$alz))
  print(table(group_df$alz))
  result<-oddsratio.fisher(table(group_df$aid, group_df$alz))
  metrics_female$OR[i]<-result$measure[2]
  metrics_female$lower_CI[i]<-result$measure[4]
  metrics_female$upper_CI[i]<-result$measure[6]
  metrics_female$pvals[i]<-result$p.value[4]
}

#compute some significance metrics, with Bonferroni correction factor = 8
metrics_female<-metrics_female %>% add_column(sig=NA)
for(i in 1:length(metrics_female$group_names)){
  if(metrics_female$pvals[i]<=0.001/length(metrics_female$group_names)){
    metrics_female$sig[i]<-"***"
  }
  else if(metrics_female$pvals[i]<=0.01/length(metrics_female$group_names)){
    metrics_female$sig[i]<-"**"
  }
  else if(metrics_female$pvals[i]<=0.05/length(metrics_female$group_names)){
    metrics_female$sig[i]<-"*"
  }
}

#run the same comparison but for the male-specific subset
metrics_male<-data.frame("group_names" = group_names,
                                    "OR" = NA,
                                    "upper_CI" = NA,
                                    "lower_CI" = NA,
                                    "pvals" = NA,
                                  "N" = NA)
for(i in 1:length(metrics_male$group_name)){
  #filter the dataframe to only include the current disease group
  df_filt<-cc %>% filter(cc[,start+i-1]==1) %>% filter(gender=="Male")
  group_df<-cc %>% filter(subclass %in% df_filt$subclass) %>% filter(gender=="Male") #get both the cases and controls for this specific group
  metrics_male$N[i]<-length(group_df$person_id)
  print(group_names[i])
  print(table(group_df$aid, group_df$alz))
  print(table(group_df$alz))
  result<-oddsratio.fisher(table(group_df$aid, group_df$alz))
  metrics_male$OR[i]<-result$measure[2]
  metrics_male$lower_CI[i]<-result$measure[4]
  metrics_male$upper_CI[i]<-result$measure[6]
  metrics_male$pvals[i]<-result$p.value[4]
}

#indicate significance
metrics_male<-metrics_male %>% add_column(sig=NA)
for(i in 1:length(metrics_male$group_names)){
  if(metrics_male$pvals[i]<=0.001/length(metrics_male$group_names)){
    metrics_male$sig[i]<-"***"
  }
  else if(metrics_male$pvals[i]<=0.01/length(metrics_male$group_names)){
    metrics_male$sig[i]<-"**"
  }
  else if(metrics_male$pvals[i]<=0.05/length(metrics_male$group_names)){
    metrics_male$sig[i]<-"*"
  }
}

#merge these dataframes
all_metrics<-rbind(metrics %>% add_column(cohort_type="Overall"),
                   metrics_female %>% add_column(cohort_type="Female"),
                   metrics_male %>% add_column(cohort_type = "Male"))

#add some important columns/metrics for plotting
#calculate -log10 pval column
all_metrics<-all_metrics %>% add_column(neg_log10_pval = -log10(all_metrics$pvals))
#create a column for pval coloring
all_metrics <- all_metrics %>% add_column(adj_pval = all_metrics$pvals*length(group_names))

#initial plot
group_p<-ggplot(data=all_metrics, aes(y=reorder(group_names, -log10(adj_pval)), x=cohort_type, 
                                            fill=ifelse(adj_pval <= 0.05, adj_pval, NA),
                                            size=OR))+
  geom_point(shape=21, color='grey')+
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "#fa1702",
                      high = "#faccc8",
                      name="adjusted pval",
                      trans="log")+
  scale_size(range=c(1,20),
             breaks=c(1, 3, 5, 10, 15, 20),
             labels=c("1", "3", "5", "10", "15", "20"),
             name="OR")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "grey", fill=NA))+
  xlab("Cohort Type")+ylab("AID Group")
group_p

#save information for more aesthetic plotting later
group_cc<-all_metrics
save(group_cc,
     file='M:/Plotting_Data/CC_group_m1b.Rdata') #CC_group_m1b representing the df for the analysis in the main text, not a sensitivity analysis



#Disorder-specific OR analysis---------------------------------------------------
aid_names<-unique(keywords$Col_names) #get a list of autoimmune disorder names
metrics<-data.frame("aid_names" = aid_names,
                    "OR" = NA,
                    "lower" = NA,
                    "upper" = NA,
                    "pval" = NA,
                    "N" = NA)

#now calculate odds ratios for each disease
start<-grep(aid_names[1], colnames(cc))[1] #get the column index where the first disorder starts
for(i in 1:length(metrics$aid_name)){
  #filter the dataframe to only include the current disorder
  print(i)
  print(aid_names[i])
  df_filt<-cc %>% filter(cc[,start+i-1]==1)
  aid_df<-cc %>% filter(subclass %in% df_filt$subclass) #get both the cases and controls for this specific disorder
  metrics$N[i]<-length(aid_df$person_id) #get the sample size
  print(length(aid_df$person_id))
  print(table(aid_df$aid, aid_df$alz))
  print(table(aid_df$aid, aid_df$gender))
  if(length(aid_df$person_id)!=0){
    #print(table(aid_df$aid, aid_df$alz))
    #print(table(aid_df$aid, aid_df$gender_source_value))
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics<-get_OR_results(metrics_df = metrics, result = result, i=i)
  }
}

#run the same comparison for the female-specific subset
metrics_female<-data.frame("aid_names" = aid_names,
                           "OR" = NA,
                           "lower" = NA,
                           "upper" = NA,
                           "pval" = NA,
                           "N" = NA)
for(i in 1:length(metrics_female$aid_name)){
  #filter the dataframe to only include the cur disease group
  print(i)
  df_filt<-cc %>% filter(cc[,start+i-1]==1) %>% filter(gender=="Female")
  aid_df<-cc %>% filter(subclass %in% df_filt$subclass) %>% filter(gender=="Female") #get both the cases and controls for this specific disorder
  if(length(aid_df$person_id)!=0){
    print(aid_names[i])
    print(table(aid_df$aid, aid_df$alz))
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics_female <- get_OR_results(metrics_df = metrics_female, result=result, i=i)
    metrics_female$N[i]<-length(aid_df$person_id)
  }
  else{metrics_female$N[i]<-0}
}

#run the same comparison but for the male-specific subset
metrics_male<-data.frame("aid_names" = aid_names,
                           "OR" = NA,
                           "lower" = NA,
                           "upper" = NA,
                           "pval" = NA,
                           "N" = NA)
for(i in 1:length(metrics_male$aid_name)){
  #filter the dataframe to only include the cur disease group
  print(i)
  df_filt<-cc %>% filter(cc[,start+i-1]==1) %>% filter(gender=="Male")
  aid_df<-cc %>% filter(subclass %in% df_filt$subclass) %>% filter(gender=="Male") #get both the cases and controls for this specific disorder
  if(length(aid_df$person_id)!=0 && (length(unique(aid_df$alz))!=1) && (length(unique(aid_df$aid))!=1)){
    print(aid_names[i])
    print(table(aid_df$aid, aid_df$alz))
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics_male <- get_OR_results(metrics_df = metrics_male, result=result, i=i)
    metrics_male$N[i]<-length(aid_df$person_id)
  }
  else{metrics_male$N[i]<-0}
}
#merge dataframes together
all_specific_metrics<-rbind(metrics %>% add_column(cohort_type = "Overall"),
                   metrics_female %>% add_column(cohort_type = "Female"),
                   metrics_male %>% add_column(cohort_type = "Male"))

#here, make a within-disease-group pvalue adjustment
#start by getting the counts of diseases within each disease grouping
count_df<-select(keywords, -c(Search_name, Remove_names, Notes, Col_names))
count_df<-count_df %>% group_by(Group) %>% 
  summarise(total_count=n(),
            .groups = 'drop')
#make a dictionary to look up number of AIDs based on disease subtype name
h1<-hash()
for(i in 1:length(count_df$Group)){
  h1[[count_df$Group[i]]]<-count_df$total_count[i]
}
#make a dictionary to look up disease subtype based on AID name
h2<-hash()
for(i in 1:length(keywords$AID_name)){
  h2[[keywords$Col_names[i]]]<-keywords$Group[i]
}

#now add in a column that denotes the pval correction for each AID
all_specific_metrics <- all_specific_metrics %>% add_column(correction_val = NA)
for(i in 1:length(all_specific_metrics$aid_names)){
  cur_aid<-all_specific_metrics$aid_names[i]
  all_specific_metrics$correction_val[i] <- h1[[h2[[cur_aid]]]]
}
#create a new column that represents the adjusted pvals after within-subtype multiple testing correction
all_specific_metrics <- all_specific_metrics %>% add_column(adj_pval = all_specific_metrics$pval*all_specific_metrics$correction_val)
#add column indicating significance
all_specific_metrics<-all_specific_metrics %>% add_column(sig=NA)
for(i in 1:length(all_specific_metrics$aid_names)){
  #print(i)
  if(is.na(all_specific_metrics$pval[i])==FALSE){
    if(all_specific_metrics$adj_pval[i]<=0.001){
      all_specific_metrics$sig[i]<-"***"
    }
    else if(all_specific_metrics$adj_pval[i]<=0.01){
      all_specific_metrics$sig[i]<-"**"
    }
    else if(all_specific_metrics$adj_pval[i]<=0.05){
      all_specific_metrics$sig[i]<-"*"
    }
  }
}


#initial disorder-stratified plot
p<-ggplot(data=all_specific_metrics, aes(y=reorder(aid_names, -adj_pval),
          x=cohort_type, 
          fill=ifelse(adj_pval <= 0.05, adj_pval, NA),
                                size=OR))+
  geom_point(shape=21, color='grey')+
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "#fa1702",
                      high = "#faccc8",
                      name="adjusted pval",
                      trans="log")+
  scale_size(range=c(3,6),
             breaks=c(1, 3, 5, 10),
             labels=c("1", "3", "5", "10"),
             name="OR")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "grey", fill=NA))+
  xlab("Stratification")+ylab("AID")+
  theme(axis.title.y = element_text(angle = 90))+
  theme(legend.key.size = unit(0.4, 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=10), #change legend title font size
        legend.text = element_text(size=8))+ #change legend text font size
  theme(legend.position = "right")
p

#save information for more asethetic plotting later
group_lookup<-h2 #table to look up what subtype group each individual autoimmune disease belongs to
ind_cc<-all_specific_metrics #ind_cc for "individual disorder case-control"
save(ind_cc,group_lookup,
     file='M:/Plotting_Data/CC_ind_m1b.Rdata') #CC_ind_m1b representing the df for the analysis in the main text, not a sensitivity analysis
