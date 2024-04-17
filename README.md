# Cross-trait Risk Analysis: Alzheimer's & Autoimmunity
This repository contains the means to conduct a cross-trait risk analysis on electronic health record study groups. It was first used to quantify risk associations between autoimmune disorders and Alzheimer's Disease (AD) in both a case-control (CC) and cohort study design.


![Fig1](https://github.com/gramey02/AD_AID_Project/assets/94878687/a1fbe3e5-3c48-4d1a-8897-1596c69854b6)


## Data Availability
The original data input for the following code was curated from the 
University of California, San Francisco (UCSF) Observational Medical 
Outcomes Partnership (OMOP) Common Data Model (CDM), and is not publicly 
available due to patient data sharing privacy. Analagous data can be 
extracted from other electronic health record databases using the medical 
billing codes in the files below, disease keyword searches in the 
'keywords' file below, or 
the pipeline can be adapted for other traits/billing codes of interest.

* [Alzheimer's Disease billing codes](https://docs.google.com/spreadsheets/d/1bzQN4iUvpV92Ke8re1JsipjxPmY9MTlKL0ZoBdUmBPg/edit#gid=0)
* [Autoimmune disorder billing codes](https://docs.google.com/spreadsheets/d/1d-O7TLsyBrxEE4MqsEuf6322mG9UcLf68ndhuOm6mu4/edit#gid=0)
* [Autoimmune Disorder and Alzheimer's Disease
Keywords](https://docs.google.com/spreadsheets/d/1ImZNCqbBNpE3UKrMWOxgX_v3An8_SYkZOIEB8gsSo5Y/edit?usp=sharing)

## Required packages
Required R packages can be found in the _requirements.rtf_ file

## Scripts
### General Workflow

* __DataCleaning__ - data cleaning pipeline for case-control and cohort study groups
* __Matching__ - matching pipeline to generate matched case-control, cohort, and longitudinal study group data sets
* __OddsRatio_Analysis__ - quantifies odds ratios in the case-control and cohort study groups. Input: matched data, output: dataframes for the Figures.Rmd plotting file
* __Longitudinal_Onset_Analysis__ - quantifies distributional differences between AD onset age between different group stratifications. Also generates survival curves and hazard ratio estimates
* __Figures.Rmd__ - generates all figures in the manuscript
* __PowerCalculations__ - generates statistical power analyses that quantify the necessary sample sizes and power for different comparisons of the odds ratio and longitudinal analyses
* __HelperFunction__ - functions that streamline code in the previous files

_note: due to proprietary OMOP naming conventions at different institutions, database column names and other UCSF-specific information has been redacted from these scripts_

## Figures
Figures from the manuscript, including all main text figures and 
supplementary figures.

## Tables
Tables from the mansucript.

## Data Structure Example
This is an example data frame that represents the structure (e.g. column 
formats) that would be ideal for input into the data cleaning scripts.


## Acknowledgements
Funding sources for this research include NIA R01AG060393, NIAMS P30 AR070155, F30 Fellowship 1F30AG079504-01, and the UCSF Discovery Fellows. We would also like to acknowledge the use of the [UCSF Information Commons](https://informationcommons.ucsf.edu/) and [UCSF Research Analysis Environment](https://it.ucsf.edu/service/rae) computational research platforms. Through these platforms, the project was supported by the National Center for Advancing Translational Sciences, National Institutes of Health, through UCSF-CTSI Grant Number UL1 TR001872. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the NIH.
