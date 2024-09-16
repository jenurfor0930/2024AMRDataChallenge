setwd("~/Desktop/")

pkgs<-c("viridis", "RColorBrewer", "ggplot2", "ggsci", "wesanderson","rainbow", "scales", "wesanderson",
        "readxl","dplyr","tidyr","janitor","stringr","ggpubr","ggstance","cowplot", "egg", "grid")
sapply(pkgs, require, character.only = TRUE)

#-------------------------------------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------------------------------------
ast_data_ecoli <- openxlsx::read.xlsx("~/SPH Dropbox/Proj_ACES/ACES_HK/Shared RIF/RIF/Output/Conference/2024/2024_Vivli_data_challenge/1_input_data/20240911_vivli_rif_Result_ECOLI_Greater_China_after2018.xlsx")
ast_data_kpneu <- openxlsx::read.xlsx("~/SPH Dropbox/Proj_ACES/ACES_HK/Shared RIF/RIF/Output/Conference/2024/2024_Vivli_data_challenge/1_input_data/20240911_vivli_rif_Result_KPNEU_Greater_China_after2018.xlsx")
ast_data_saurs <- openxlsx::read.xlsx("~/SPH Dropbox/Proj_ACES/ACES_HK/Shared RIF/RIF/Output/Conference/2024/2024_Vivli_data_challenge/1_input_data/20240911_vivli_rif_Result_SAURS_Greater_China_after2018.xlsx")


#-------------------------------------------------------------------------------------------------------------
# Merge data
#-------------------------------------------------------------------------------------------------------------
merged_data <- bind_rows(ast_data_ecoli, ast_data_kpneu, ast_data_saurs)



#-------------------------------------------------------------------------------------------------------------
# Identify bacteria
#-------------------------------------------------------------------------------------------------------------

# mutate(ESBL_ECOLI = ifelse(Species == "Escherichia coli"  & (Cefepime == "Resistant" | Ceftazidime == "Resistant" | Ceftriaxone == "Resistant" | Phenotype == "ESBL"), "Yes", NA)) %>%
# mutate(ESBL_KPNEU = ifelse(Species == "Klebsiella pneumoniae"  & (Cefepime == "Resistant" | Ceftazidime == "Resistant" | Ceftriaxone == "Resistant" | Phenotype == "ESBL"), "Yes", NA)) %>%
# mutate(MRSA_2 = ifelse(Species == "Staphylococcus aureus" & Oxacillin == "Resistant", "Yes", NA)) %>%
# mutate(MRSA_2 = case_when(Phenotype == "MRSA" ~ "Yes", TRUE ~ MRSA_2)) %>%


ast_result = merged_data %>%
  mutate(Carbapenems_ECOLI = ifelse(Species == "Escherichia coli" & Carbapenem == "Carbapenem_Resistant", "Yes", NA)) %>%
  mutate(Carbapenems_KPNEU = ifelse(Species == "Klebsiella pneumoniae" & Carbapenem == "Carbapenem_Resistant", "Yes", NA)) %>%
  mutate(ESBL_ECOLI = ifelse(Species == "Escherichia coli"  & ESBL == "ESBL", "Yes", NA)) %>%
  mutate(ESBL_KPNEU = ifelse(Species == "Klebsiella pneumoniae"  & ESBL == "ESBL", "Yes", NA)) %>%
  mutate(MRSA_2 = ifelse(Species == "Staphylococcus aureus"  & MRSA == "MRSA", "Yes", NA)) %>%
  mutate(VISA = ifelse(Species == "Staphylococcus aureus"  & VIRSA == "VISA", "Yes", NA)) %>%
  mutate(VRSA = ifelse(Species == "Staphylococcus aureus"  & VIRSA == "VRSA", "Yes", NA))
  

df = ast_result %>%
  select(Sector, Country, Species, Carbapenems_ECOLI, Carbapenems_KPNEU, ESBL_ECOLI, ESBL_KPNEU, MRSA_2, VRSA, VISA)

ast_result_binary <- df %>%
  mutate(across(Carbapenems_ECOLI:VISA, ~ ifelse(. == "Yes", 1, 0), .names = "binary_{.col}"))

# Summarize data to calculate percentages
total_counts <- ast_result_binary %>%
  group_by(Sector, Species) %>%
  summarise(Total = n(), .groups = 'drop')


summary_data <- ast_result_binary %>%
  group_by(Sector, Species) %>%
  summarise(
    Carbapenems_ECOLI = sum(binary_Carbapenems_ECOLI, na.rm = TRUE),
    Carbapenems_KPNEU = sum(binary_Carbapenems_KPNEU, na.rm = TRUE),
    ESBL_ECOLI = sum(binary_ESBL_ECOLI, na.rm = TRUE),
    ESBL_KPNEU = sum(binary_ESBL_KPNEU, na.rm = TRUE),
    MRSA_2 = sum(binary_MRSA_2, na.rm = TRUE),
    VRSA = sum(binary_VRSA, na.rm = TRUE),
    VISA = sum(binary_VISA, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_longer(cols = -c(Sector, Species), names_to = "Phenotype", values_to = "Count") %>%
  left_join(total_counts, by = c("Sector", "Species")) %>%  # Join total counts
  mutate(Percentage = Count / Total * 100) %>%  # Calculate percentage
  mutate(Phenotype = gsub("binary_", "", Phenotype)) 

#-------------------------------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------------------------------
custom_colors <- c("Clinical" = "#d1495b", "Community" = "#66a182")

plot <- ggplot(summary_data, aes(x = Phenotype, y = Percentage, fill = Sector)) +
  geom_bar(stat = "identity", color = "black",position = position_dodge()) +
  geom_text(aes(label = ifelse(Count > 0, paste("n =", Count), "")),
            position = position_dodge(width = 0.9), 
            vjust = -0.5,
            size = 3) +
  labs(x = "Clinically important bacteria", y = "Percentage (%)") +
  ylim(0,50) +
  scale_fill_manual(values = custom_colors, name = NULL) +
  scale_x_discrete(labels = c("Carbapenems_KPNEU" = "Carbapenem-resistant\nK. pneumoniae",
                              "Carbapenems_ECOLI" = "Carbapenem-resistant\nE. coli", 
                              "ESBL_ECOLI" = "ESBL-producing\nE. coli", 
                              "ESBL_KPNEU" = "ESBL-producing\nK. pneumoniae",
                              "MRSA_2" = "Methicillin-resistant\nS. aureus\n(MRSA)", 
                              "VRSA" = "Vancomycin-resistant\nS. aureus\n(VRSA)", 
                              "VISA" = "Vancomycin-intermediate\nS. aureus\n(VISA)")) +
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        legend.position = "bottom")

ggsave("bugdrug_prevaence.pdf", width = 27, height = 20, units = "cm")
