get_OR_results <- function(metrics_df, result, i){
  metrics_df$OR[i]<-result$measure[2]
  metrics_df$lower[i]<-result$measure[4]
  metrics_df$upper[i]<-result$measure[6]
  metrics_df$pval[i]<-result$p.value[4]
  return(metrics_df)
}