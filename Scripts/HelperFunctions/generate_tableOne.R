generate_tableOne <- function(df, study_type){
  #will need to stratify by and summarize slightly different variables depending on cohort type
  if(study_type=="cc"){
    strata_var = "alz"
    sum_var = "aid"
  }
  else{
    strata_var = "aid"
    sum_var = "alz"
  }
  
  # #print variables
  # dput(names(df))
  
  #ensure that certain variables are factors
  df$gender <- as.factor(df$gender)
  df$ethnicity <- as.factor(df$ethnicity)
  df$race <- as.factor(df$race)
  df$birth_year <- as.numeric(as.character(df$birth_year))
  df$aaDeath_years <- as.numeric(as.character(df$aaDeath_years))
  df$aid <- as.factor(as.character(df$aid))
  df$alz <- as.factor(as.character(df$alz))
  
  #identify variables that you want summarized
  myVars <- c("gender", "race", "ethnicity", "birth_year", 
              "aaDeath_years", sum_var)
  #generate and print table
  tab <- CreateTableOne(vars = myVars, strata = strata_var , data = df, smd = TRUE)
  return(tab)
}