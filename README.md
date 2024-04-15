# AD_AID_Project
This repository contains the means to conduct a cross-trait analysis on electronic health record study groups. Its first use is detailed in...
___insert publication here___

...where it was used to quantify risk associations between autoimmune disorders and Alzheimer's Disease (AD) in both a case-control (CC) and cohort study design.

## Data Availability
* While the original data input for the following code was curated from the University of California, San Francisco (UCSF) Observational Medical Outcomes Partnership (OMOP) Common Data Model (CDM), it can be easily applied to other electronic health record databases for other traits of interest.

## Required packages


## Scripts
* Figures.Rmd - generates all figures in the manuscript
* PowerAnalysis.Rmd - generates statistical power analysis quantifying necessary sample sizes and power for different statistical tests used througout the study
* StudyGroup_Demographics - generates tables containing information on demographics of the study groups, including sself-reported race and ethnicity, gender, age, etc.
* CC_OR_Analysis.R - quantifies odds ratios and differences in AD prevalence between autoimmune patients and non-autoimmune controls in the case-control study group
* Cohort_OR_Analysis.R - quantifies odds ratios and differences in AD prevalence between autoimmune patients and non-autoimmune controls in the cohort study group
* DataCleaning_CC.R - data cleaning pipeline for the case-control study group
* DataCleaning_Cohort.R - data cleaning pipeline for the cohort study group 

## Figures
* Fig 1 - 

## Tables
