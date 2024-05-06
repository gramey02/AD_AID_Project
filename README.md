# Cross-trait Risk Analysis: Alzheimer's & Autoimmunity
This repository contains the means to conduct a cross-trait risk analysis on electronic health record (EHR) study groups. It was first used to quantify risk associations between autoimmune disorders and Alzheimer's Disease (AD) in both a case-control and cohort study design.

Please find the corresponding preprint here: https://www.medrxiv.org/content/10.1101/2024.05.02.24306649v1

![Fig1](https://github.com/gramey02/AD_AID_Project/assets/94878687/ec43c2f0-989e-4d4d-8112-eae558b406b6)



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
### General Workflow:
![Code_workflow](https://github.com/gramey02/AD_AID_Project/assets/94878687/37bdb612-27b5-4a58-b741-9dea309c5688)


* __DataCleaning__ - data cleaning pipeline for case-control and cohort study groups
* __Matching__ - matching pipeline to generate matched case-control, cohort, and longitudinal study group data sets
* __OddsRatio_Analysis__ - quantifies odds ratios in the case-control and cohort study groups. Input: matched data, output: dataframes for the Figures.Rmd plotting file
* __Longitudinal_Onset_Analysis__ - quantifies distributional differences between AD onset age between different group stratifications. Also generates survival curves and hazard ratio estimates
* __Figures.Rmd__ - generates all figures in the manuscript
* __PowerCalculations__ - generates statistical power analyses that quantify the necessary sample sizes and power for different comparisons of the odds ratio and longitudinal analyses
* __HelperFunctions__ - functions that streamline code in the previous files

_note: due to proprietary OMOP naming conventions at different institutions, database column names and other UCSF-specific information has been redacted from these scripts_

## Figures
Figures from the manuscript, including all main text figures and 
supplementary figures.

## Tables
Tables from the mansucript.

## Data Structure Example
This is an example data frame that represents the structure (e.g. column 
formats) that would be ideal for input into the data cleaning scripts.

### For the autoimmunity patients and corresponding controls:
This is what the "prematched_rc_data.Rdata" input file to the data cleaning scripts should look like.

| person_id     | num_dx        | num_visits | date of death | gender | ethnicity             | race                                       | birth_year | aid  | aid_date | first_visit | last_visit | AID name                    |
|---------------|---------------|------------|---------------|--------|-----------------------|--------------------------------------------|------------|------|----------|-------------|------------|-----------------------------|
| 1             | 10            | 30         | 5/5/20        | Female | Hispanic/Latino       | Other                                      | 1940       | 1    | 1/15/55  | 1/1/10      | 7/6/18     | Type 1 Diabetes             |
| 2             | 30            | 70         | 1/30/19       | Male   | Not Hispanic/Latino   | White                                      | 1955       | 1    | 3/1/00   | 1/1/99      | 8/2/15     | Inflammatory Bowel Disease  |
| 2             | 30            | 70         | 1/30/19       | Male   | Not Hispanic/Latino   | White                                      | 1955       | 1    | 9/15/10  | 1/1/99      | 8/2/15     | Rheumatoid Arthritis        |
| 3             | 25            | 50         | 6/1/21        | Female | Not Hispanic/Latino   | Black/African American                     | 1959       | 1    | 2/1/12   | 3/4/10      | 5/6/20     | Vitiligo                    |
| ...           | ...           | ...        | ...           | ...    | ...                   | ...                                        | ...        | ...  | ...      | ...         | ...        | ...                         |
| ...           | ...           | ...        | ...           | ...    | ...                   | ...                                        | ...        | ...  | ...      | ...         | ...        | ...                         |
| 32601         | 27            | 100        | 4/3/15        | Male   | Not Hispanic/Latino   | Native Hawaiian or Other Pacific Islander  | 1935       | 0    | NA       | 5/12/04     | 6/7/14     | NA                          |
| 32602         | 15            | 40         | 7/3/19        | Female | Not Hispanic/Latino   | Unknown                                    | 1951       | 0    | NA       | 6/20/00     | 7/30/12    | NA                          |

_Notice that a person (e.g. person_id 2) can have more than one autoimmune disorder_

### For the Alzheimer's patients and corresponding controls:
This is what the "prematched_cc_data.Rdata" input file to the data cleaning scripts should look like.

| person_id     | num_dx        | num_visits | date of death | gender | ethnicity             | race                                       | birth_year | alz  | alz_date | first_visit | last_visit |
|---------------|---------------|------------|---------------|--------|-----------------------|--------------------------------------------|------------|------|----------|-------------|------------|
| 1             | 50            | 20         | 6/24/16       | Female | Hispanic/Latino       | Other                                      | 1934       | 1    | 5/3/04   | 12/4/00     | 8/29/10    |
| 2             | 20            | 16         | 10/10/20      | Male   | Not Hispanic/Latino   | Asian                                      | 1945       | 1    | 9/21/00  | 11/26/99    | 3/7/06     |
| ...           | ...           | ...        | ...           | ...    | ...                   | ...                                        | ...        | ...  | ...      | ...         | ...        |
| 3             | 11            | 30         | 8/8/2019      | Female | Not Hispanic/Latino   | White                                      | 1949       | 0    | 10/1/15  | 1/15/10     | 2/2/17     |
| ...           | ...           | ...        | ...           | ...    | ...                   | ...                                        | ...        | ...  | ...      | ...         | ...        |

AID = autoimmune disorder,
aid_date = first instance of autoimmune disorder billing code,
ALZ = Alzheimer's disease,
alz_date = first instance of Alzheimer's billing code,
dx = diagnosis


## Acknowledgements
Funding sources for this research include NIA R01AG060393, NIAMS P30 AR070155, F30 Fellowship 1F30AG079504-01, and the UCSF Discovery Fellows. We would also like to acknowledge the use of the [UCSF Information Commons](https://informationcommons.ucsf.edu/) and [UCSF Research Analysis Environment](https://it.ucsf.edu/service/rae) computational research platforms. Through these platforms, the project was supported by the National Center for Advancing Translational Sciences, National Institutes of Health, through UCSF-CTSI Grant Number UL1 TR001872. Its contents are solely the responsibility of the authors and do not necessarily represent the official views of the NIH.
