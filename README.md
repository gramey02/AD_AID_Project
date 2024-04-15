# AD_AID_Project
This repository contains the means to conduct a cross-trait risk analysis on electronic health record study groups. It was first used to quantify risk associations between autoimmune disorders and Alzheimer's Disease (AD) in both a case-control (CC) and cohort study design.

## Data Availability
The original data input for the following code was curated from the University of California, San Francisco (UCSF) Observational Medical Outcomes Partnership (OMOP) Common Data Model (CDM), and is not publicly available due to patient data sharing privacy. Analagous data can be extracted from other electronic health record databases using the medical billing codes in the files below, or the pipeline can be adapted for other traits of interest.

## Required packages
Required R packages can be found in the _requirements.txt_ file

## Scripts
* Figures.Rmd - generates all figures in the manuscript
* PowerAnalysis.Rmd - generates statistical power analysis quantifying necessary sample sizes and power for different statistical tests used througout the study
* StudyGroup_Demographics - generates tables containing information on demographics of the study groups, including sself-reported race and ethnicity, gender, age, etc.
* CC_OR_Analysis.R - quantifies odds ratios and differences in AD prevalence between autoimmune patients and non-autoimmune controls in the case-control study group
* Cohort_OR_Analysis.R - quantifies odds ratios and differences in AD prevalence between autoimmune patients and non-autoimmune controls in the cohort study group
* DataCleaning_CC.R - data cleaning pipeline for the case-control study group
* DataCleaning_Cohort.R - data cleaning pipeline for the cohort study group

_note: due to proprietary OMOP naming conventions at different institutions, database column names and other UCSF-specific information has been redacted_

## Figures

## Tables
