# Comparing Antimicrobial Resistance Profiles of Patients and Community Residents: Implications for Targeted Interventions

<br/>

Pei-Yu Huang *, Chloe Sze Lon Chui , Peng Wu

School of Public Health, The University of Hong Kong, Hong Kong SAR, China

<br/>

*pyh0930@hku.hk

<br/>

## Getting Started

### Folder structure
**0_Data/:** Contains the standardized MIC fold change data and susceptibility cutoffs used in the analysis.<br>
**1_Scripts/:** Contains all the scripts used for analysis and processing.<br>
**2_Plots_Tables/:** Contains the generated plots and tables.<br>
**3_PDF/: Contains** the final report in PDF format.<br>

<br/>

### Utilization
Our analytical approaches enable you to explore and compare resistance levels across various antibiotics and bacteria, irrespective of differing susceptibility cutoffs.<br/> 
To conduct this analysis, you will need the complete MIC data for all relevant antibiotics. Your data should be formatted similarly to `20240913_atlas_rif_HK_mic_sus_fold_change_meta_3Bac.rds` in the "0_Data" folder. Please feel free to reach out if you have any queries. <br/>

<br/>

## Overview

Antimicrobial resistance (AMR) poses a significant global health threat, compromising treatment effectiveness for infectious diseases. Key pathogens, such as ESBL-producing and carbapenem-resistant *E. coli*, *K. pneumoniae*, and MRSA, show alarming resistance rates, particularly in clinical settings. However, understanding their AMR profiles in the community remains limited, hindering effective treatment development.<br>

This project analyses clinical data from Pfizer's ATLAS initiative, obtained through https://amr.vivli.org, and community data from Hong Kong to compare AMR profiles of these pathogens. The insights gained aim to inform targeted interventions, enhancing antimicrobial stewardship in both clinical and community contexts.

<br/>

## Methods ##
### Antibiotic Resistance Profile

We assessed the minimal inhibitory concentration (MIC) and antibiotic resistance of 11 antibiotics for clinical isolates of *E. coli*, *K. pneumoniae*, and *S. aureus*—including clinically significant strains—collected from 2018 to 2022 in mainland China, Taiwan, and Hong Kong, as well as community isolates from Hong Kong.<br/>
The data was merged and cleaned using `1_Data_cleaning.R`. The analysis and plots was performed using `2_Figure1_Antibiogram.R` and `3_Figure2_Pathogen_prevalence.R`.


<br/>

### Susceptibility Standardization

Categorical MICs were converted to numerical values using a normal distribution approach. MIC fold changes were then calculated against CLSI/EUCAST susceptible cutoffs, which was shown in `20240911_ast_mic_cutoff_susceptible_all.xlsx` in the "0_Data" folder, to determine susceptibility.<br/>
This step was conducted using `3_MIC_fold_change_transformation.R`.

<br/>

### Statistical Analysis
Bray-Curtis dissimilarity distances and Principal Coordinates Analysis (PCoA) were used to visualize variations in AMR profiles. Statistical significance was assessed using PERMANOVA and pairwise comparisons.<br/>
This analysis was conducted using `4_Figure3_PCoA.R`.

<br/>

## Results

AMR profiles differed notably between clinical and community settings, especially for *E. coli* and *K. pneumoniae*, with resistance rates in clinical isolates up to 80 times higher. Over 98% of clinical isolates were intermediate resistant to colistin.<br/>
In community settings, higher rates of vancomycin-intermediate (VISA) and vancomycin-resistant (VRSA) strains were observed, indicating a significant public health concern. The PCoA results showed clear differences in AMR profiles, driven by higher resistance to specific antibiotics in clinical isolates.

<br/>

## Impact

This study highlights the urgent need for robust AMR surveillance and targeted interventions, as resistance profiles between clinical and community settings are converging. The findings underscore AMR as a public health priority, necessitating international collaboration to standardize detection efforts and enhance comparability across regions.

<br/>

## Acknowledgments
We thank Vivli and Pfizer for providing the ATLAS clinical AMR data essential to this research. We are grateful to our team leader, Dr. Peng Wu, for her guidance and support, and to Dr. Keiji Fukuda for initiating the AMR study in the Hong Kong community. Our appreciation also goes to the HKU-Pasteur Research Pole for providing laboratory facilities. This study was supported by the Research Impact Fund (R7033-18) and the Collaborative Research Fund (C5063-22GF) from the Research Grants Council of the Hong Kong Government.

