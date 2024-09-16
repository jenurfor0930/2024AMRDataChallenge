
# load vivli data
setwd("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/")
source("~/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_intro/AMR_analysis/00_Manuscript_1/AMR_functions_20230308.R")
v_raw = read.csv("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/2024_05_28 atlas_antibiotics.csv") # 917049 rows, 135 columns

# should remove the wrong data
# 1964315
v_raw1 = v_raw %>% filter(!Isolate.Id=="1964315")

# convert country name to country code
library(devtools)
devtools::install_github("gpw13/whoville")
library(whoville)
v_df =
  v_raw1 %>%
  mutate(iso3 = names_to_code(Country),
         who_region = iso3_to_regions(iso3, region = "who_region"), .after = "Country")

# add taiwan back to WPR
v_df$who_region[v_df$Country=="Taiwan"] = "WPR"

# filtering criteria
# 1. E. coli, K. pneumoniae, S. aureus
# 2. WPR nations
# 3. known 18+ age
# extract E. coli and S. aureus and WPR (Wstern pacific region) nations
v_sub = v_df %>% filter(Species %in%c("Escherichia coli", "Staphylococcus aureus", "Klebsiella pneumoniae") & who_region=="WPR") %>% select(-State) %>%
  filter (!Age.Group %in% c("0 to 2 Years", "13 to 18 Years","3 to 12 Years", "Unknown"))# 32066 rows, 135 columns

# write out the data
openxlsx::write.xlsx(v_sub, "atlas_antibiotics_over18_WPS_ECOLI_KPNEU_SAURS.xlsx")


#table(v_sub$Species, v_sub$Phenotype)
# -----------------------------------------------------------------------
# Start from here! Read in data
# -----------------------------------------------------------------------
# atlas data
v_df = openxlsx::read.xlsx("atlas_antibiotics_over18_WPS_ECOLI_KPNEU_SAURS.xlsx") # rows 32066  cols 135

v_df1 = v_df %>% filter (Country %in% c("Hong Kong", "China", "Taiwan") ) # 35353 rows, 135 columns
unique(v_df1$Speciality)

v_df %>% filter(Country =="Hong Kong") %>% select(Clindamycin) # 1 row, 5 columns)

# -----------------------------------------------------------------------
# remove genotype data
# -----------------------------------------------------------------------
v_ast = v_df %>% select(Isolate.Id:Tebipenem_I) # phenotype data # 112
# remove all NA columns
v_ast <- v_ast[,colSums(is.na(v_ast))<nrow(v_ast)] # 105 columns; 32066 rows

# adjust for merging
v_ast1 = v_ast %>% select(-c(Family, Speciality,In...Out.Patient)) %>% mutate(Sector = 'Clinical', .after = "Study", Study  = "Atlas" )

# -----------------------------------------------------------------------
# read in RIF data
# -----------------------------------------------------------------------
#r_df = openxlsx::read.xlsx("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_intro/AMR_analysis/00_Manuscript_1/1_Raw_data/ast_result_view_20230725_rmFW_nonFit_First_Second.xlsx") # rows 32066  cols 135
r_df = readRDS("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_intro/AMR_analysis/00_Manuscript_1/1_Raw_data/AST_questionnaire_raw_MIC_data_20231020.rds")

# extract human ECOLI and SAURS isolates above 19 years old
r_ast = r_df %>% filter(a_bacteria %in%c("ECOLI", "KPNEU","SAURS") & Sector == "Human" & age_BL >=19) # 2768 rows, 85 columns

r_ast  %>% filter(a_lab_id =="923658KBS1")

# extract metadata for merging
# "Isolate.Id"                  "Study"                       "Sector"                      "Species"
#[5] "Country"                     "iso3"                        "who_region"                  "Gender"
#[9] "Age.Group"                   "Source"                      "Year"                        "Phenotype"
r_ast1 =r_ast %>% select(a_lab_id, a_bacteria, a_source, age_BL, gender, s_collection_date, a_amocla_mic:a_vancomycin_result)
r_ast1 = r_ast1 %>% rename(Isolate.Id = a_lab_id, Species = a_bacteria, Source = a_source, Gender = gender)
r_ast1 = r_ast1 %>% mutate(Study = "RIF", Sector = "Community", .after = "Isolate.Id") %>%
  mutate(Country= "Hong Kong",iso3="HKG", who_region="WPR", .after = "Species") %>% relocate (Gender, .after = who_region)
# tidy age group

r_ast2 = r_ast1 %>% mutate(Age.Group = cut(age_BL, c(0,64,84, Inf),include.lowest = TRUE,
                                           labels=c('19 to 64 Years','65 to 84 Years','85 and Over')), .after = Gender) %>% select(-age_BL)
# table(r_ast1$age_BL, r_ast1$Age.Group) # checking

# tidy date
r_ast2 = r_ast2 %>% mutate(s_collection_date = cleanDate1(s_collection_date)) %>%
  mutate(Year = substring(s_collection_date,1,4), .after = "Source") %>% select(-s_collection_date)

# tidy penotype after merge
# add a blank column for pheontype
r_ast2 = r_ast2 %>% mutate(Phenotype = NA_character_, .after=Year )

# -----------------------------------------------------------------------
# match antibiotics
# -----------------------------------------------------------------------
abx_name = openxlsx::read.xlsx("match_atlas_rif_antibiotics.xlsx")
match_abx = abx_name %>% filter(!is.na(antibiotic)) %>% filter(!is.na(atlas))

# transform to long format
v_ast_long = v_ast1 %>% pivot_longer(cols = -c(Isolate.Id: Phenotype), names_to = "atlas", values_to = "value")
v_ast_long = v_ast_long %>% mutate(MIC = ifelse(grepl("_I", atlas), "Result", "MIC"), .after =atlas) %>% mutate(atlas = gsub("_I", "", atlas))

r_ast_long = r_ast2 %>% pivot_longer(cols = -c(Isolate.Id: Phenotype), names_to = "antibiotic", values_to = "value")
r_ast_long = r_ast_long %>% mutate(MIC = ifelse(grepl("_mic", antibiotic), "MIC", "Result"), .after =antibiotic) %>%
  mutate(antibiotic = gsub("a_(.*?)_.*$", "\\1", antibiotic))


# unify abx name (vivli vs rif)
# vivli
match_abx_v = match_abx %>% select(atlas,name_full)
v_ast_long1 = merge(v_ast_long, match_abx_v) %>% select(-atlas) %>%  relocate(name_full, .after =Phenotype)

# rif
match_abx_r = match_abx %>% select(antibiotic,name_full)
r_ast_long1 = merge(r_ast_long, match_abx_r) %>% select(-antibiotic) %>%  relocate(name_full, .after =Phenotype)

# combine vivli and rif
ast_long = rbind(v_ast_long1, r_ast_long1)

ast_long$MIC = factor(ast_long$MIC, levels = c("Result","MIC"))
ast_long = ast_long %>% arrange(MIC)

# tidy up result
ast_long = ast_long %>% mutate(value1 = ifelse(value == "SUS", "Susceptible",
                                               ifelse (value == "SDD", "Susceptible Dose Dependent",
                                                       ifelse(value == "SUS/INT", "Susceptible/Intermediate",
                                                              ifelse(value == "RES", "Resistant",
                                                                     ifelse(value == "INT", "Intermediate",
                                                                            ifelse (value == "NOTSUS", "Resistant",
                                                                                    ifelse (value =="", NA, value))))))))


ast_long = ast_long %>% select(-value) %>% rename(value = value1)




# tidy up bacteria
unique(ast_long$Species)
ast_long = ast_long %>% mutate(Species = ifelse(Species == "ECOLI", "Escherichia coli",
                                                ifelse (Species == "SAURS", "Staphylococcus aureus",
                                                        ifelse (Species=="KPNEU", "Klebsiella pneumoniae", Species))))


# saveRDS file
saveRDS (ast_long, "RIF_Atlas_WPR_3Bac_ast_long_20240912.rds")

#-------------------------------------------------------------
# start from here
#-------------------------------------------------------------
# only keep Greater China and after 2018
ast_long = readRDS ("RIF_Atlas_WPR_3Bac_ast_long_20240912.rds")
ast_long = ast_long %>% filter(Country %in% c("China","Hong Kong", "Taiwan") & Year >=2018)
list_bac = c("Escherichia coli","Klebsiella pneumoniae","Staphylococcus aureus")

ast_long = ast_long %>% mutate(Sector1 = ifelse (Country == "Hong Kong" & Sector=="Community",  "Community (Hong Kong)",Country), .after=Sector)
ast_long$Sector1 = factor (ast_long$Sector1, levels =c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)"))

# Split the result and mic into two data.frames
ast_long_result = ast_long %>% filter(MIC == "Result") %>% select(-MIC)
ast_long_mic = ast_long %>% filter(MIC == "MIC") %>% select(-MIC)

# make it wide
ast_wide_result = ast_long_result %>% pivot_wider(names_from = name_full, values_from = value)
ast_wide_mic = ast_long_mic %>% pivot_wider(names_from = name_full, values_from = value)




# -----------------------------------------------------------------------
# tidy Phenotype (clincially important pathogen)
# -----------------------------------------------------------------------
ast_wide_result1 = ast_wide_result %>% mutate(Carbapenem = ifelse(Species %in% c(list_bac[1], list_bac[2]) & Imipenem == "Resistant" | Meropenem == "Resistant", "Carbapenem_Resistant", NA),
                      ESBL = ifelse(Species %in% c(list_bac[1], list_bac[2]) & Cefepime == "Resistant" | Ceftazidime == "Resistant" | Ceftriaxone == "Resistant", "ESBL", NA),
                      MRSA = ifelse(Species==list_bac[3] & Oxacillin == "Resistant", "MRSA",
                                    ifelse (Species=="SAURS" & Oxacillin == "Susceptible", "MSSA", NA)),
                      VIRSA = ifelse(Species==list_bac[3] & Vancomycin == "Resistant", "VRSA",
                                    ifelse (Species == list_bac[3] & Vancomycin == "Intermediate", "VISA", NA)), .after = Phenotype)

# merge patogen into a column
ast_wide_result2 = ast_wide_result1 %>% mutate(Pathogen = paste(ESBL, MRSA, VIRSA, Carbapenem, sep = "_")) %>% select(-c(ESBL, MRSA, VIRSA, Carbapenem))
ast_wide_result3 = ast_wide_result2 %>% mutate (Pathogen1 = ifelse (Pathogen == "NA_NA_NA_NA", NA,
                                                                    ifelse (Pathogen %in% c("ESBL_NA_NA_Carbapenem_Resistant", "NA_NA_NA_Carbapenem_Resistant"), "Carbapenem_Resistant",
                                                                            ifelse (Pathogen == "ESBL_NA_NA_NA", "ESBL",
                                                                                    ifelse (Pathogen  %in% c("NA_MRSA_NA_NA", "NA_MRSA_VISA_NA","NA_MRSA_VRSA_NA", "NA_NA_VISA_NA"), "MRSA_VIRSA",Pathogen
                                                                                            ))))) %>% select(-Pathogen) %>% rename(Pathogen = Pathogen1)

ast_wide_result3 = ast_wide_result3 %>% mutate(Pathogen = ifelse(is.na(Pathogen), "Others", Pathogen))
ast_wide_result4= ast_wide_result3 %>% mutate(Sector_Pathogen = paste0(Sector, "_", Pathogen), .before = Phenotype)


# add phenotypes back to MIC table
pheno = ast_wide_result4 %>% select(Isolate.Id, Pathogen, Sector_Pathogen)
ast_wide_mic1 = merge(ast_wide_mic, pheno) %>% relocate(Pathogen, Sector_Pathogen, .before =Phenotype)

openxlsx::write.xlsx(ast_wide_result4, "WPR_atlas_rif_merged_result_20240913_Greater_China_after2018.xlsx")
openxlsx::write.xlsx(ast_wide_mic1, "WPR_atlas_rif_merged_mic_20240913_Greater_China_after2018.xlsx")




