
# standardise the AMR using MIC fold changes
setwd("/Users/pyhuang//SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/")
source("~/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/0_vivli_rif_functions.R")

# load packages
pkgs<-c("ggplot2","dplyr","tidyr","stringr","cowplot","vegan", "phyloseq","openxlsx","dplyr")
sapply(pkgs, require, character.only = TRUE)

# input data
# China, Hong Kong , and Taiwan vivli data after 2018
library(openxlsx)
mic_raw = openxlsx::read.xlsx("WPR_atlas_rif_merged_mic_20240913_Greater_China_after2018.xlsx")
mic_df = mic_raw %>% relocate(Phenotype, .before ="Amoxicillin/Clavulanic.acid")

# ==============================================================================
# Prepare table for MIC fold change  with susceptible as baseline)
# ==============================================================================
bac_level<-c("Escherichia coli","Klebsiella pneumoniae" ,"Staphylococcus aureus")
host_level = c("Clinical", "Community")

ecoli = mic_df %>% filter(Species =="Escherichia coli")
kpneu = mic_df %>% filter(Species =="Klebsiella pneumoniae")
saurs = mic_df %>% filter(Species =="Staphylococcus aureus")

# remove NA columns
ecoli_tidy <- ecoli[,colSums(is.na(ecoli))<nrow(ecoli)]
kpneu_tidy <- kpneu[,colSums(is.na(kpneu))<nrow(kpneu)]
saurs_tidy <- saurs[,colSums(is.na(saurs))<nrow(saurs)]


#-------------------------------------------------------------
# match results with MIC data
#-------------------------------------------------------------
ecoli_keep = keep_both_antibiotic(ecoli_tidy)
kpneu_keep = keep_both_antibiotic(kpneu_tidy)
saurs_keep = keep_both_antibiotic(saurs_tidy)

#-------------------------------------------------------------
# tidy up results
#-------------------------------------------------------------
col1=colnames(ecoli_keep)
col2=colnames(kpneu_keep)
col3=colnames(saurs_keep)

ast_wide_result1 = openxlsx::read.xlsx("WPR_atlas_rif_merged_result_20240913_Greater_China_after2018.xlsx")

ecoli_r = ast_wide_result1 %>% filter(Species =="Escherichia coli")
kpneu_r = ast_wide_result1 %>% filter(Species =="Klebsiella pneumoniae")
saurs_r = ast_wide_result1 %>% filter(Species =="Staphylococcus aureus")

ecoli_r1 = ecoli_r %>% select(all_of(col1))
kpneu_r1 = kpneu_r %>% select(all_of(col2))
saurs_r1 = saurs_r %>% select(all_of(col3))

# ------------------------------------------------------------
# add a random number to the MIC value
# ------------------------------------------------------------
ecoli_num = mic_numeric (ecoli_keep, 0.2)
kpneu_num = mic_numeric (kpneu_keep, 0.2)
saurs_num = mic_numeric (saurs_keep, 0.2)


# ------------------------------------------------------------
# transform mic value to fold changes compared to mic of susceptilbe
# ------------------------------------------------------------
# change from long format to wide format and transform to MIC fold changes relative to CLSI susceptible cutoffs
ecoli_fold = convert_to_fold_changes(ecoli_num)
kpneu_fold = convert_to_fold_changes(kpneu_num)
saurs_fold = convert_to_fold_changes(saurs_num)

#clipr::write_clip(ecoli_fold)

ecoli_meta = ecoli_num%>%select(Isolate.Id:Phenotype) %>% unique %>%  tibble::column_to_rownames(var = "Isolate.Id")
kpneu_meta = kpneu_num%>%select(Isolate.Id:Phenotype) %>% unique %>%  tibble::column_to_rownames(var = "Isolate.Id")
saurs_meta = saurs_num%>%select(Isolate.Id:Phenotype) %>% unique %>%  tibble::column_to_rownames(var = "Isolate.Id")

# merge MIC fold changes and meta data
ecoli_all = merge(ecoli_meta, ecoli_fold , by="row.names")
kpneu_all = merge(kpneu_meta, kpneu_fold , by="row.names")
saurs_all = merge(saurs_meta, saurs_fold , by="row.names")


# output as a list
output_list = list(ecoli_fold, kpneu_fold, saurs_fold, ecoli_meta, kpneu_meta, saurs_meta)
names(output_list) = c("ecoli_fold", "kpneu_fold", "saurs_fold", "ecoli_meta", "kpneu_meta", "saurs_meta")
saveRDS(output_list, "20240913_atlas_rif_HK_mic_sus_fold_change_meta_3Bac.rds")

