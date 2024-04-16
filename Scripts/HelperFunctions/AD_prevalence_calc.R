AD_prevalence_calc <- function(aid_female, aid_male, control_female, control_male, n=1000){
  #note: n is number of bootstrap iterations
  set.seed(103114) #save final sys.time seed #also move to top!!!
  l<-list(aid_female, aid_male, control_female, control_male)
  true_prevs<-rep(0,4)
  distributions<-list(aid_female, aid_male, control_female, control_male)
  lower_CI<-rep(0, 4)
  upper_CI<-rep(0,4)
  i<-1
  for(df in l){
    true_prevs[i]<-length((df %>% filter(alz==1))$person_id)/length(df$person_id)
    cur_df_prevs<-c()
    for(j in 1:1000){
      #print(j)
      samp <- df[sample(nrow(df),size=nrow(df),replace=TRUE),] #create a bootstrap sample
      #calculate our test statistic, the ad prevlance, in the bootstrap sample
      prev<-length((samp %>% filter(alz==1))$person_id)/length(samp$person_id)
      cur_df_prevs<-c(cur_df_prevs, prev)
    }
    lower_CI[i]<-quantile(cur_df_prevs, probs = c(0.025, 0.95))[[1]] #2.5th percentile
    upper_CI[i]<-quantile(cur_df_prevs, probs = c(0.025, 0.95))[[2]] #97.5th percentile
    distributions[[i]]<-cur_df_prevs
    i<-i+1
  }
  #combine prevalences into a data frame with confidence intervals
  prev_df<-data.frame("group" = c("aid_female", "aid_male", "control_female", "control_male"),
                      "gender" = c("Female", "Male", "Female", "Male"),
                      "aid_status" = c("aid", "aid", "control", "control"),
                      "ad_prev" = true_prevs,
                      "lower_CI" = lower_CI,
                      "upper_CI" = upper_CI,
                      "n" = c(length(aid_female$person_id), length(aid_male$person_id),
                              length(control_female$person_id), length(control_male$person_id)))
  return(list(prev_df, distributions))
}