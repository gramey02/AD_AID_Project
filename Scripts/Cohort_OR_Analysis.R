#### Script to run odds ratio analysis on retrospective cohort

#load libraries
library(dplyr)
library(epitools)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(hash)
library(MatchIt)
library(pwr)
library(tibble)

#load helper functions
source('~/EHR/get_OR_results.R') #function that takes in the result and metrics data frame, and fills it with values
source('~/EHR/AD_prevalence_calc.R') #function to create confidence intervals by bootstrapping

#Data loading------------------------------------------------------------------
#see https://github.com/gramey02/AD_AID_Project for an example dataframe structure
load("M:/Matched_data/matched_rc_clean_2023_08_04_m1b.Rdata") #main text is matched_rc_clean_2023_08_04_m1b
rc<-rc_matched #rename the data frame for easier typing below

#load keywords dataframe with autoimmune disorder search terms
keywords<-readxl::read_xlsx("~/EHR/keywords_2022_10_26.xlsx")

#Overall OR Analysis------------------------------------------------------------
#dataframe to hold odds ratios for different cohort stratifications
cohort_type<-c('Full Cohort', 'Female Cohort', 'Male Cohort')
metrics<-data.frame("Cohort_Type" = cohort_type,
                    "OR" = 0,
                    "lower" = 0,
                    "upper" = 0,
                    "pval" = 0)
#full cohort analysis
result<-oddsratio.fisher(table(rc$aid, rc$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=1) #i=1 for the first overall comparison

#within female OR test
female <- rc %>% filter(gender=="Female")
result<-oddsratio.fisher(table(female$aid, female$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=2) #i=2 for the female comparison
#within male OR test
male <- rc %>% filter(gender=="Male")
result<-oddsratio.fisher(table(male$aid, male$alz))
metrics<-get_OR_results(metrics_df = metrics, result = result, i=3) #i=2 for the female comparison

#plot overall OR results
p<-ggplot(data=metrics, aes(x=reorder(Cohort_Type, (OR)), y=(OR), 
                            ymin=(lower), ymax=(upper))) +
  geom_pointrange() + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Cohort Type") + ylab("Odds Ratio (95% CI)") +
  theme_bw()+
  ggtitle('Retrospective Cohort Analysis\nAD Odds Ratios')+
  theme(plot.title = element_text(hjust = 0.5, vjust = 1))+
  scale_y_log10(limits=c(NA,NA), 
                breaks=c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4))+
  geom_hline(yintercept = 1, linetype="dashed", 
             color = "black", size=1.1)
#view the plot
p
#indicate significance levels, Bonferroni correction val = 3, for 3 comparisons
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

#save data frames and sample sizes for plotting later
rc_metrics_overall<-metrics
rc_table_e_o<- table(rc$aid, rc$alz) #2 x 2 exposure x outcome table for everyone, "e", overall "o" (i.e. not stratified by disease/group)
rc_table_f_o<- table(female$aid, female$alz) #table for females, "f"
rc_table_m_o<- table(male$aid, male$alz) #table for males, "m"
save(rc_metrics_overall, rc_table_e_o, rc_table_f_o, rc_table_m_o, 
     file='M:/Plotting_Data/RC_overall_m1b.Rdata') #m1b representing the df for the analysis in the main text, not a sensitivity analysis


#Sex-specific prevalence analysis----------------------------------------------
#isolate the controls
controls<-rc %>% filter(aid==0)
colnames(controls)[colnames(controls) == "subclass"] <- "subclass1" #rename the subclass column bc it can throw errors
controls<-subset(controls, select = -c(prop.score, weights)) #remove prior matching variables
controls$gender<-as.factor(controls$gender) #factorize the gender column
#get equal numbers of control females and control males with matching
control_m.out<-matchit(gender ~ birth_year + race + ethnicity+aaDeath_years,
                       data = controls, method = "nearest", ratio = 1, verbose=T)
control_m.data<-match.data(control_m.out, distance = "prop.score")
control_female<-control_m.data %>% filter(aid==0) %>% filter(gender=="Female")
control_male<-control_m.data %>% filter(aid==0) %>% filter(gender=="Male")

#isolate the cases
aid<-rc %>% filter(aid==1)
colnames(aid)[colnames(aid) == "subclass"] <- "subclass1" #rename the subclass column
aid<-subset(aid, select = -c(prop.score, weights)) #remove some prior matching columns
aid$gender<-as.factor(aid$gender) #factorize the gender column
#get equal numbers of case females and case males with matching
aid_m.out<-matchit(gender ~ birth_year + race + ethnicity+aaDeath_years,
                       data = aid, method = "nearest", ratio = 1, verbose=T)
aid_m.data<-match.data(aid_m.out, distance = "prop.score")
aid_female<-aid_m.data %>% filter(aid==1) %>% filter(gender=="Female")
aid_male<-aid_m.data %>% filter(aid==1) %>% filter(gender=="Male")

#calculate prevalences
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

#plot the prevlances
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
  labs(fill="Sex")+
  ylab("AD prevalence")+
  scale_fill_simpsons()
overall_p

#now compute significances for the prevalence in each group using fisher's exact test

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

#compare aid females to control females - need to account for differing numbers of columns
aid_female_sub<-select(aid_female, c(alz))
control_female_sub<-select(control_female, c(alz))
total_df<-select(rbind(aid_female_sub %>% add_column(stratification="aid_female"),
                       control_female_sub %>% add_column(stratification="control_female")),
                 c(alz, stratification))
female_test<-fisher.test(table(total_df$stratification, total_df$alz))
female_test

#compare aid males to control males
aid_male_sub<-select(aid_male, c(alz))
control_male_sub<-select(control_male, c(alz))
total_df<-select(rbind(aid_male_sub %>% add_column(stratification = "aid_male"),
                       control_male_sub %>% add_column(stratification = "control_male")), 
                 c(alz, stratification))
male_test<-fisher.test(table(total_df$stratification, total_df$alz))
male_test

#great, now we have the significance values for the prevalence plots! save these and the dataframe for plotting
rc_prev_sig<-data.frame(comparison = c("aid female vs aid male", "control female vs control male", "aid female vs control female", "aid male vs control male"),
                        pval = c(aid_test$p.value, control_test$p.value, female_test$p.value, male_test$p.value),
                        pval_adj = 4 * c(aid_test$p.value, control_test$p.value, female_test$p.value, male_test$p.value))
#also save some sample sizes
rc_prev_ss<-data.frame(aid = length(aid_female$person_id), control = length(control_female$person_id))
rc_prev_df<-prev_df #rename to be more specific
save(rc_prev_ss, rc_prev_sig, rc_prev_df,
     file='M:/Plotting_Data/RC_prev_m1b.Rdata') #m1b representing the df for the analysis in the main text, not a sensitivity analysis


#Autoimmune Disorder Subtype Odds Ratio Analysis------------------------------------------------------------
group_names<-unique(keywords$Group) #get a list of subtype names
metrics<-data.frame("group_names" = group_names,
                    "OR" = NA,
                    "lower" = NA,
                    "upper" = NA,
                    "pval" = NA,
                    "N" = NA)
#calculate odds ratios for each group
start<-grep(group_names[1], colnames(rc)) #get the column index where the first disease group starts
for(i in 1:length(metrics$group_names)){
  #filter the dataframe to only include the current disease group and corresponding controls
  df_filt<-rc %>% filter(rc[,start+i-1]==1)
  print(group_names[i])
  print(length(df_filt$person_id))
  metrics$N[i]<-length(df_filt$person_id)
  print(table(df_filt$aid, df_filt$alz))
  print(table(df_filt$aid, df_filt$gender))
  result<-oddsratio.fisher(table(df_filt$aid, df_filt$alz))
  metrics<-get_OR_results(metrics_df = metrics, result = result, i=i) #store the resulting Fisher's exact outputs
}

#Conduct the same comparison for the female-specific subset
metrics_female<-data.frame("group_names" = group_names,
                    "OR" = NA,
                    "lower" = NA,
                    "upper" = NA,
                    "pval" = NA,
                    "N" = NA)
for(i in 1:length(metrics_female$group_name)){
  #filter the dataframe to only include the cur disease group
  df_filt<-rc %>% filter(rc[,start+i-1]==1) %>% filter(gender=="Female")
  print(group_names[i])
  print(table(df_filt$aid, df_filt$alz))
  metrics_female$N[i]<-length(df_filt$person_id)
  result<-oddsratio.fisher(table(df_filt$aid, df_filt$alz))
  metrics_female<-get_OR_results(metrics_df = metrics_female, result = result, i=i)
}

#conduct the same comparison for the male-specific subset
metrics_male<-data.frame("group_names" = group_names,
                         "OR" = NA,
                         "lower" = NA,
                         "upper" = NA,
                         "pval" = NA,
                         "N" = NA)
for(i in 1:length(metrics_male$group_name)){
  #filter the dataframe to only include the cur disease group
  df_filt<-rc %>% filter(rc[,start+i-1]==1) %>% filter(gender=="Male")
  print(group_names[i])
  print(table(df_filt$aid, df_filt$alz))
  metrics_male$N[i]<-length(df_filt$person_id)
  result<-oddsratio.fisher(table(df_filt$aid, df_filt$alz))
  metrics_male<-get_OR_results(metrics_df = metrics_male, result = result, i=i)
}

#combine dataframes
all_group_metrics<-rbind(metrics %>% add_column(Cohort_Type = "Overall"),
                         metrics_female %>% add_column(Cohort_Type = "Female"),
                         metrics_male %>% add_column(Cohort_Type = "Male"))

#indicate significance, Bonferroni correction val = 8, for 8 subtype categories
all_group_metrics<- all_group_metrics %>% add_column(sig=NA)
for(i in 1:length(all_group_metrics$group_name)){
  if(all_group_metrics$pval[i]<=0.001/(length(group_names))){
    all_group_metrics$sig[i]<-"***"
  }
  else if(all_group_metrics$pval[i]<=0.01/(length(group_names))){
    all_group_metrics$sig[i]<-"**"
  }
  else if(all_group_metrics$pval[i]<=0.05/(length(group_names))){
    all_group_metrics$sig[i]<-"*"
  }
}
#calculation -log10 pval column
all_group_metrics<-all_group_metrics %>% add_column(neg_log10_pval = -log10(all_group_metrics$pval))
#create a column for pval coloring
all_group_metrics <- all_group_metrics %>% add_column(adj_pval = all_group_metrics$pval*length(group_names))

#initital plot
group_p<-ggplot(data=all_group_metrics, aes(y=reorder(group_names, -log10(adj_pval)), x=Cohort_Type, 
                                fill=ifelse(adj_pval <= 0.05, adj_pval, NA),
                                size=OR))+
  geom_point(shape=21, color='grey')+
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "#fa1702",
                       high = "#faccc8",
                       name="adjusted pval",
                      trans="log")+
  scale_size(range=c(5,13),
             breaks=c(1, 3, 5, 10),
             labels=c("1", "3", "5", "10"),
             name="OR")+
  theme_minimal()+
  theme(panel.border = element_rect(color = "grey", fill=NA))+
  xlab("Cohort Type")+ylab("AID Group")
group_p

#save information for plotting later
group_rc<-all_group_metrics
save(group_rc,
     file='M:/Plotting_Data/RC_group_m1b.Rdata') #m1b representing the df for the analysis in the main text, not a sensitivity analysis)



#Specific Autoimmune Disorder OR Analysis---------------------------------------------------
aid_list<-keywords$Col_names #get list of autoimmune disorder (aid) names
#create a dataframe to hold the odds ratio metrics for each disorder
metrics_stratified<-data.frame(AID=aid_list,
                               "OR" = 0,
                               "lower" = 0,
                               "upper" = 0,
                               "pval" = NA,
                               "N" = NA)
#loop through each aid in the overall dataframe, calculating ORs
start<-grep(aid_list[1], colnames(rc))[1] #get the index where the first disorder starts in the data frame
for(i in 1:length(metrics_stratified$AID)){
  print(i)
  aid_name<-metrics_stratified$AID[i]
  #filter the total df for the current aid
  aid_df<-rc %>% filter(rc[,start+i-1]==1)
  metrics_stratified$N[i]<-length(aid_df$person_id)
  if(length(aid_df$person_id)==0 || sum(aid_df$alz)==0 || sum(aid_df$aid)==0){}
  else{
    #fill in the sample size
    print(aid_name)
    print(length(aid_df$person_id))
    print(table(aid_df$aid, aid_df$alz))
    print(table(aid_df$aid, aid_df$gender))
    #calculate odds ratios
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics_stratified<-get_OR_results(metrics_df = metrics_stratified, result = result, i=i)
  }
}

#conduct the same comparison for the female-specific subset
female<-rc %>% filter(gender=="Female")
metrics_stratified_female<-data.frame(AID=aid_list,
                                      "OR" = 0,
                                      "lower" = 0,
                                      "upper" = 0,
                                      "pval" = NA,
                                      "N" = NA)
for(i in 1:length(metrics_stratified_female$AID)){
  print(i)
  aid_name<-metrics_stratified_female$AID[i]
  print(aid_name)
  #filter the total df for the current aid
  aid_df<-female %>% filter(female[,start+i-1]==1)
  #fill in the sample size
  print(table(aid_df$aid, aid_df$alz))
  N_cases<-length((aid_df %>% filter(aid==1))$person_id)
  N_controls<-length((aid_df %>% filter(aid==0))$person_id)
  metrics_stratified_female$N[i]<-length(aid_df$person_id)
  if((N_cases==0) || (N_controls==0) || sum(aid_df$alz)==0){print("pass")}
  else{
    #calculate odds ratios
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics_stratified_female<-get_OR_results(metrics_df = metrics_stratified_female, result = result, i=i)
  }
}

#conduct the same comparison for the male-specific subset
male<-rc %>% filter(gender=="Male")
metrics_stratified_male<-data.frame(AID=aid_list,
                                      "OR" = 0,
                                      "lower" = 0,
                                      "upper" = 0,
                                      "pval" = NA,
                                      "N" = NA)
for(i in 1:length(metrics_stratified_male$AID)){
  print(i)
  aid_name<-metrics_stratified_male$AID[i]
  print(aid_name)
  #filter the total df for the current aid
  aid_df<-male %>% filter(male[,start+i-1]==1)
  #fill in the sample size
  print(table(aid_df$aid, aid_df$alz))
  N_cases<-length((aid_df %>% filter(aid==1))$person_id)
  N_controls<-length((aid_df %>% filter(aid==0))$person_id)
  metrics_stratified_male$N[i]<-length(aid_df$person_id)
  if((N_cases==0) || (N_controls==0) || sum(aid_df$alz)==0){print("pass")}
  else{
    #calculate odds ratios
    result<-oddsratio.fisher(table(aid_df$aid, aid_df$alz))
    metrics_stratified_male<-get_OR_results(metrics_df = metrics_stratified_male, result = result, i=i)
  }
}
#merge data frames
all_specific_metrics<-rbind(metrics_stratified %>% add_column(cohort_type = "Overall"),
                            metrics_stratified_female %>% add_column(cohort_type="Female"),
                            metrics_stratified_male %>% add_column(cohort_type="Male"))

#here, make a within-disease-group pvalue adjustment
#start by getting the counts of diseases within each disease grouping
count_df<-select(keywords, -c(Search_name, Remove_names, Notes, Col_names))
count_df<-count_df %>% group_by(Group) %>% 
  summarise(total_count=n(),
            .groups = 'drop')
#make a dictionary to look up number of AIDs based on group name
h1<-hash()
for(i in 1:length(count_df$Group)){
  h1[[count_df$Group[i]]]<-count_df$total_count[i]
}
#make a dictionary to look up disease group based on AID name
h2<-hash()
for(i in 1:length(keywords$AID_name)){
  h2[[keywords$Col_names[i]]]<-keywords$Group[i]
}

#add in a column that denotes the pval correction for each individual disorder
all_specific_metrics <- all_specific_metrics %>% add_column(correction_val = NA)
for(i in 1:length(all_specific_metrics$AID)){
  cur_aid<-all_specific_metrics$AID[i]
  all_specific_metrics$correction_val[i] <- h1[[h2[[cur_aid]]]]
}
#create a new column that represents the adjusted pvals after multiple testing correction
all_specific_metrics <- all_specific_metrics %>% add_column(adj_pval = all_specific_metrics$pval*all_specific_metrics$correction_val)
#indicate significance
all_specific_metrics<-all_specific_metrics %>% add_column(sig=NA)
for(i in 1:length(all_specific_metrics$AID)){
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


#initial plot
specific_p<-ggplot(data=all_specific_metrics, aes(y=reorder(AID,-adj_pval),
                                                  #factor(AID, level=c('IBD','Psor','RA','PMR',
                                                  #     'PA','AIT','T1D')), 
                                 x=cohort_type, 
                                 fill=ifelse(adj_pval <= 0.05, adj_pval, NA),
                                size=OR))+
  geom_point(shape=21, color='grey')+
  scale_x_discrete(position = "top") +
  scale_fill_gradient(low = "#fa1702",
                      high = "#faccc8",
                      name="adjusted pval",
                      trans="log")+
  scale_size(range=c(3,8),
             breaks=c(1, 3, 5, 10),
             labels=c("1", "3", "5", "10"),
             name="OR")+
  #geom_text(aes(label = round(OR, digits=2)), color = "black", size=4, 
  #          fontface='bold', hjust=0.5, vjust=-1.5)+
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
  
specific_p


#save this information for plotting later
ind_rc<-all_specific_metrics #ind_cc for "individual disease case-control"
save(ind_rc,
     file='M:/Plotting_Data/RC_ind_m1b.Rdata') #m1b representing the df for the analysis in the main text, not a sensitivity analysis
