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

merged_data_2 <- merged_data %>%
  select(-1) %>%
  select(-c(ESBL, Phenotype, MRSA, VIRSA, Carbapenem))

ast_data_raw <- merged_data_2 %>%
  mutate(Source = "Community", Imipenem = ifelse(Imipenem == "Intermediate", "Susceptible", Imipenem))


#-------------------------------------------------------------------------------------------------------------
# Select data
#-------------------------------------------------------------------------------------------------------------
df = ast_data_raw %>%
  select(Sector, Species, 'Amoxicillin/Clavulanic.acid',
         Ampicillin, Aztreonam, Cefepime, Ceftazidime, Ciprofloxacin, Colistin, 
         Gentamicin, Imipenem, Meropenem, 'Trimethoprim/Sulfamethoxazole', Clindamycin,
         Daptomycin, Erythromycin, Levofloxacin, Linezolid, Oxacillin, Tigecycline, Vancomycin) %>%
  pivot_longer(
    cols = -c(Sector, Species),
    names_to = "Antibiotic",
    values_to = "Value")


df_filtered <- df %>%
  filter(!is.na(Value)) %>%
  group_by(Species, Sector, Antibiotic, Value) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(Species, Sector, Antibiotic) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()

# Calculate the total
total_counts <- df_filtered %>%
  group_by(Species, Sector) %>%
  summarise(total_n = sum(count), .groups = 'drop')
#-------------------------------------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------------------------------------
df_filtered$Value <- factor(df_filtered$Value, levels = rev(c("Resistant", "Intermediate", 
                                                              "Susceptible Dose Dependent", "Susceptible", "NOCLSI")))

desired_order <- c("Tigecycline", "Vancomycin", "Tetracycline", "Quinupristin/Dalfopristin",
                   "Colistin", "Linezolid", "Daptomycin", "Trimethoprim/Sulfamethoxazole", 
                   "Azithromycin", "Erythromycin", "Gentamicin", "Clindamycin",
                   "Moxifloxacin", "Ciprofloxacin", "Levofloxacin", "Aztreonam",
                   "Imipenem", "Meropenem", "Cefoxitin", "Cefepime", "Ceftazidime",
                   "Ceftriaxone", "Penicillin", "Oxacillin", "Ampicillin", "Amoxicillin/Clavulanic.acid")

df_filtered <- df_filtered %>%
  mutate(Antibiotic = factor(Antibiotic, levels = desired_order))

# colour
col = rev(c("#ff7f24","#FFD700","#e5d8bd","#DDF1DA","#c3cfb6"))

plot <- ggplot(df_filtered, aes(fill = Value, y = percentage, x = Antibiotic)) + 
  geom_bar(position = "fill", stat = "identity") +
  facet_grid(Species ~ Sector, scales = "free", space = "free") +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis labels to percentage
  labs(y = "Percentage") +
  scale_fill_manual(values = col) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text.x=element_text(angle=90, hjust=0.95,vjust=0.2, size=11),
        axis.text.y=element_text(size=11),
        axis.title.x = element_text(size=11, face = "bold"),
        legend.title=element_blank(),
        legend.text = element_text(face = "plain", size=12),
        legend.position="bottom",
        legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.spacing.y  = unit(0.007, 'cm'),
        strip.background =element_blank(),
        strip.text.x = element_text(face = "bold", size=15),
        strip.text.y = element_text(face = "bold.italic",size=15)) 


pdf(file = "AST_clinical_community.pdf", width = 10, height = 9)
plot
dev.off()
