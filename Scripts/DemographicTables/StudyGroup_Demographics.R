#load libraries
library(dplyr)
library(MatchIt)
library(hash)
library(tibble)
library(epitools)
library(ggplot2)
library(tableone)

#Table One Generation-----------------------------------------------------------

#load the case-control study group data
load('M:/Matched_data/matched_cc_clean_2023_08_04_m1b.Rdata') ####change file name here for other case-control (cc) study group
cc_df<-cc_matched
#load the function that generates a demographic table for different study designs
source("~/EHR/generate_tableOne.R")
#generate table for case-control cohort
cc_tab<-generate_tableOne(df = cc_df, study_type = "cc")

#load the cohort study group data and generate the correpsonding table
load('M:/Matched_data/matched_rc_clean_2023_08_04_m1b.Rdata') ####change file name here for other cohort (rc) study group
rc_df<-rc_matched
rc_tab <- generate_tableOne(df = rc_df, study_type = "rc")
rm(m.data5, exact_data, m.out5)

#to export, print the quote-delimited tables and copy into txt file
print(cc_tab, quote = TRUE, noSpaces = TRUE, smd = TRUE)
print(rc_tab, quote = TRUE, noSpaces = TRUE, smd = TRUE)
