# Power analysis script for cohort study design

#load libraries
library(dplyr)
library(epiR)
library(MatchIt)
library(hash)
library(tibble)
library(epitools)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(powerSurvEpi)
library(pwr)

#load longitudinal study group, where everyone gets AD at different ages
load("M:/Matched_data/matched_rc_temporal_2023_08_04_m1.Rdata")
df<-rc_temporal
keywords<-readxl::read_xlsx("~/EHR/keywords_2022_10_26.xlsx")

###first perform a power analysis for the age distribution comparison
#power analysis for the mann whitney U test is the same as a two-sample t test, but with sample size adjustments based on distribution of two sample populations
#use pwr.t.test (since sample sizes are equal here)

res<-pwr.t.test(n=(length((df %>% filter(aid==1))$person_id)) * ((pi)/3), sig.level = 0.05, 
           type="two.sample", alternative="two.sided", d=(cohens_d(df, aaAD_years~aid, paired=FALSE))$effsize)
print(res$power)
#note how we had to adjust the sample size by pi/3, because our samples are roughly normally distributed



###now calculate the power for the survival analysis/cox proportional hazard model
#here, treat is the number of people with an autoimmune disease & control is the number of people without
library(epiR)
res2<-epi.sscomps(treat = NA, control = NA, n = length(df$person_id), power = 0.80, 
            r = 1, design = 1, sided.test = 2, nfractional = FALSE, conf.level = 0.95)
print(res2)

longitudinal_power<-data.frame("N_autoimmune" = length((df %>% filter(aid==1))$person_id),
                               "N_control" = length((df %>% filter(aid==0))$person_id),
                               "distributional_power" = res$power,
                               "cox_ph_ss_for_80_power_stanford" = res3$hazard)

###change these paths here
write.csv(longitudinal_power, file = "~/EHR/toBox/longitudinal_power_stats.csv")
save(res, res2, res3, file = "~/EHR/toBox/longitudinal_results.Rdata")
