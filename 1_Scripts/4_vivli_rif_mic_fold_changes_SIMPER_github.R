
# Use MIC folc change of all antibiotic profile of each isolate to plot PCoA
# comapre community and clinic
# problem: some antibiotics were not tested in the clincial isolates, which will also contribute to the difference.

# load packages
pkgs<-c("viridis", "RColorBrewer", "ggplot2", "ggsci", "wesanderson","rainbow", "scales", "wesanderson",
        "readxl","dplyr","tidyr","janitor","stringr","ggpubr","ggstance","cowplot",
        "vegan", "phyloseq","openxlsx","ggthemes","pracma", "egg", "grid","dplyr","pheatmap","geomtextpath")
sapply(pkgs, require, character.only = TRUE)

# install phyloseq
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("phyloseq")

# set directory and load functions
setwd("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/")
source("~/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_intro/AMR_analysis/00_Manuscript_1/AMR_functions_20230308.R")

#============================================================================================
# Simper analysis
#============================================================================================
# input data
phy_list = readRDS ("vivli_mic_sus_fold_changes_phyloseqs.rds")
bac_level = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus")

#-------------------------------------------------------------
# add sector 1 into the metadata
#-------------------------------------------------------------
for (k in 1:length(phy_list)){
    sample_info = data.frame(sample_data(phy_list[[k]]))
    sample_info = sample_info %>% mutate(Sector1 = ifelse (Country == "Hong Kong" & Sector=="Community",  "Community (Hong Kong)",Country), .after=Sector)
    sample_info$Sector1 = factor (sample_info$Sector1, levels =c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)"))
    sample_info$Pathogen = sub("_", ' ', sample_info$Pathogen)

    # update the metadata to Phyloseq object
    sample_data(phy_list[[k]]) = sample_info
}

saveRDS (phy_list, "vivli_mic_sus_fold_changes_phyloseqs_Sector1.rds")




#---------------------
# design
#---------------------
# 1. overall simper
# 2. by clinically importamt pathogen

#---------------------
# comparing groups
#---------------------
# 1. compare between clinical data and community data in Hong Kong
# 2. compare clinical data between different contries from vivli
#---------------------
# only present the abundance data
# ---------------------
source("simper_pretty.R")
source("R_krusk.R")
source("phyloseq_sep_variable.R")


# use the original data to calculate the MIC fold changes
# input data
phy_list = readRDS ("vivli_mic_sus_fold_changes_phyloseqs_Sector1.rds")
bac_level = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus")
sector_level = c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)")
#k=1
#-------------------------------------------------------------
# split by pathogen
# Split phyloseq-class object by sample-level variable
#-------------------------------------------------------------
#install.packages("remotes")
#remotes::install_github("vmikk/metagMisc")

#-------------------------------------------------------------
# Split phyloseq-class object by pathogen
#-------------------------------------------------------------
# combine pathogen list for all species
phy_list = readRDS ("vivli_mic_sus_fold_changes_phyloseqs_Sector1.rds")

bac_level = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus")
bac_level_patho = c(rep("Escherichia coli",3),
                    rep("Klebsiella pneumoniae", 3),
                    rep("Staphylococcus aureus", 2))

sector_level = c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)")

phy_list_patho = list()


for (m in 1:length(phy_list)){
    phy_tmp = phyloseq_sep_variable(phy_list[[m]], variable = "Pathogen")
    bac_name = bac_level[m]
    phy_list_patho = c(phy_list_patho, phy_tmp)}

names(phy_list_patho) = paste0(bac_level_patho,"_",names(phy_list_patho))

saveRDS(phy_list_patho, "vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")

#-------------------------------------------------------------
phy_list_patho = readRDS("vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")

    # make class
    for ( k in 1:length(phy_list_patho)){
        asv_df = otu_table(phy_list_patho[[k]])
        sample_info = sample_data(phy_list_patho[[k]])

        patho_name = names(phy_list_patho)[k]
        #data_type = "abun" # abun occur
        output_name = paste0("SIMPER_vivli_rif_AST_mic_sus_folds_abun_", patho_name)

        ast_col=openxlsx::read.xlsx("vivli_arrow_abbr.xlsx")
        taxa_df = data.frame(Taxonomy = NA, OTU = rownames(asv_df))

        # prepare class table
        for (t in 1:nrow(taxa_df)){
            taxa_df$Taxonomy[t]=ast_col[ast_col$Antibiotic==taxa_df$OTU[t],]$Class
        }


        simper.pretty(asv_df, sample_info, c('Sector1'), perc_cutoff=0.9, low_cutoff = 'y', low_val=0.01, output_name)
        simper.results = read.csv(paste0(output_name,"_clean_simper.csv"))
        # simper.results= simper.results[,-1]
        cols= colnames(simper.results)
        simper.results[cols] <- lapply(simper.results[cols], factor)

        kruskal.pretty(asv_df, sample_info, simper.results, c('Sector1'), output_name, taxa_df) # cannot use underscore in loc-time name

        #Import
        KW.results = data.frame(read.csv(paste0(output_name,"_krusk_simper.csv")))

        #Remove non-significant
        KW.results.signif = KW.results[KW.results$fdr_krusk_p.val < 0.05,]
        KW.results.signif=KW.results.signif[!is.na(KW.results.signif$X),]

        #Order by OTU#
        KW.results.signif = KW.results.signif[with(KW.results.signif, order(OTU)),]
        KW.results.signif_1 = KW.results.signif[,-1]

        # add back class of abx
        for (m in 1:nrow(KW.results.signif_1)){
            KW.results.signif_1$Taxonomy[m] = taxa_df[taxa_df$OTU==KW.results.signif_1$OTU[m],"Taxonomy"]
        }

        sum(KW.results.signif_1$SIMPER)

        write.csv(KW.results.signif_1, paste0(output_name,"_krusk_simper_significant.csv"), row.names = T)

        #----------------------------
        # re-organise simper table (containing abx and class 2)
        #----------------------------
        KW.results.signif_1 = KW.results.signif_1 %>% mutate(group_1 = sub("(^.*?)_(.*$)","\\1", KW.results.signif_1$Comparison),
                                                             group_2=sub("(^.*?)_(.*$)","\\2", KW.results.signif_1$Comparison)) %>% relocate(group_1, group_2) %>%
            select(-Comparison)


        KW.results.signif_1$group_1 = factor(KW.results.signif_1$group_1, levels=sector_level)
        KW.results.signif_1$group_2 = factor(KW.results.signif_1$group_2, levels=sector_level)
        KW.results.signif_1$Taxonomy = factor(KW.results.signif_1$Taxonomy, levels=unique(ast_col$Class))
        KW.results.signif_1$OTU = factor(KW.results.signif_1$OTU, levels=unique(ast_col$Antibiotic))

        KW.results.signif_1 = KW.results.signif_1 %>% arrange(group_1,group_2, Taxonomy, OTU)
        kw.tidy = KW.results.signif_1 %>% select(c("group_1","group_2",  "Taxonomy","OTU","SIMPER", "Left.mean.abund","Right.mean.abund" ))

        write.csv(kw.tidy, paste0(output_name,"_krusk_significant_order.csv"), row.names = F)

    }


# --------------------------------------------------------
# Overall permutations
# --------------------------------------------------------
phy_list_patho = readRDS("vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")

#require(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",  force = TRUE)
library("pairwiseAdonis")

permanova_all = data.frame()
permanova_pair = c()

betadisper_all = data.frame()
betadisper_pair = c()

for (k in 1:length(phy_list_patho)){
    phy_abun = phy_list_patho[[k]]

    ps_dist_matrix <- phyloseq::distance(phy_abun, method ="bray")

    # PERMANOVA
    # all
    df1 = vegan::adonis2(ps_dist_matrix ~ phyloseq::sample_data(phy_abun)$Sector1, p.adjust.m ='bonferroni', perm=999)
    permanova_all_tmp = data.frame(r2.all = df1$R2[1], p.all = df1$`Pr(>F)`[1])
    permanova_all = rbind(permanova_all, permanova_all_tmp)

    # pairwise
    df2 = pairwise.adonis(ps_dist_matrix, phyloseq::sample_data(phy_abun)$Sector1, p.adjust.m ='bonferroni', perm=999)
    permanova_pair = c(permanova_pair, list(df2) )

    # betadispersion
    dispGP <-vegan::betadisper(ps_dist_matrix, phyloseq::sample_data(phy_abun)$Sector1, type=c("median"))
    df3 = anova(dispGP) # all
    betadisper_tmp = data.frame(p = df3$`Pr(>F)`[1])
    betadisper_all = rbind(betadisper_all, betadisper_tmp)

    dispGP_df<-TukeyHSD(dispGP)$group # pairwise
    betadisper_pair = c(betadisper_pair, list(dispGP_df) )
}

rownames(permanova_all) = names(phy_list_patho)
names(permanova_pair) = names(phy_list_patho)

rownames(betadisper_all) = names(phy_list_patho)
names(betadisper_pair) = names(phy_list_patho)

saveRDS(permanova_all, "PERMANOVA_All_mic_sus_fold_phyloseqs_vivli_rif.rds")
saveRDS(permanova_pair, "PERMANOVA_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")

saveRDS(betadisper_all, "Betadispersion_All_mic_sus_fold_phyloseqs_vivli_rif.rds")
saveRDS(betadisper_pair, "Betadispersion_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")


#----------------------------------------------------------------
# Calculating the MIC fold change means for each antibiotic in sector
#----------------------------------------------------------------
phy_list_patho = readRDS("vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")

# function to calculate the mean, median, top 5% of the MIC fold changes
asv_df = data.frame(otu_table(phy_list_patho[[k]]))
sample_info = data.frame(sample_data(phy_list_patho[[k]]))

sector_level = c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)")
level = sector_level
stat = "mean"


mid_fold_stat_vivli = function(asv_df, sample_info, level, stat){

    group_df <- data.frame(Antibiotic=rownames(asv_df))

    for ( n in 1:length(level)){

        # if use sector as grouping
        sample_info = sample_info %>% mutate(Isolate.ID = rownames(sample_info))
        id_list=sample_info %>% filter(Sector1==level[n]) %>% select(Isolate.ID)
        id_list = id_list$Isolate.ID

        if (length(id_list)==1) {
            col_names = sub("X(.*$)","\\1",colnames(asv_df))
            tmp_df = asv_df [,col_names%in% id_list]
            tmp_add=data.frame(antibiotic = rownames(asv_df), Variable = round(asv_df [,col_names%in% id_list],2))
        } else if (length(id_list)>1) {
            col_names = sub("X(.*$)","\\1",colnames(asv_df))
            tmp_df = asv_df [,col_names%in% id_list]
            if (stat == "mean"){
                tmp_add=data.frame(antibiotic = rownames(tmp_df), Variable = round(rowMeans (tmp_df),2)) # mean
            } else if (stat == "median"){
                tmp_add=data.frame(antibiotic = rownames(tmp_df), Variable = round(apply(tmp_df, 1, median, na.rm=TRUE), 2)) # median
            } else if (stat == "top5"){
                top5= function(x){quantile(as.numeric(x), probs = 0.95)}; tmp_add=data.frame(antibiotic = rownames(tmp_df), Variable = round(apply(tmp_df,1, top5),2)) # top 5 %
            }
        } else if (length(id_list)==0){
            tmp_add=data.frame(Antibiotic = rownames(asv_df), Variable = NA)
        }

        colnames(tmp_add) = c("Antibiotic",level[n])
        group_df = merge(group_df, tmp_add, by="Antibiotic", all = TRUE)
    }

    # add Class
    ast_col=openxlsx::read.xlsx("vivli_arrow_abbr.xlsx")
    ast_ref = ast_col %>% select(Antibiotic,  Class)
    group_df_tidy = merge(ast_ref, group_df)
    # order
    abx_level = ast_col$Antibiotic
    cls_level = unique(ast_col$Class)
    group_df_tidy$Class = factor (group_df_tidy$Class, levels=cls_level)
    group_df_tidy$Antibiotic = factor (group_df_tidy$Antibiotic, levels=abx_level)

    group_df_tidy=group_df_tidy%>% relocate(Class) %>% arrange(Class, Antibiotic)

    return(group_df_tidy)
}
#----------------------------------------------------------------
phy_list_patho = readRDS("vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")
patho_level = names(phy_list_patho)

    xlsx::write.xlsx (sector_level, "MIC_fold_change_table_sector_mean_vivli_rif.xlsx") # create an empty df

    mic_sum_all = c()

    for ( k in 1:length(phy_list_patho)) {
        sample_info = data.frame(sample_data(phy_list_patho[[k]]))
        fold_df = data.frame(otu_table(phy_list_patho[[k]]))
        mic_sum_tab = mid_fold_stat_vivli (fold_df, sample_info, sector_level,"mean")

        mic_sum_all = c(mic_sum_all, list(mic_sum_tab))
        xlsx::write.xlsx (mic_sum_tab, "MIC_fold_change_table_sector_mean_vivli_rif.xlsx", sheetName = patho_level[k], append=T, row.names = F)
    }

    names(mic_sum_all) = patho_level
    saveRDS(mic_sum_all, "MIC_fold_change_table_sector_mean_vivli_rif.rds")




#----------------------------------------------------------------
# re-organise simper aggregate table
#----------------------------------------------------------------

#----------------------------------------------------------------
# function to reformat
#----------------------------------------------------------------
# read in mic fold change table

raw.tab = read.csv("SIMPER_vivli_rif_AST_mic_sus_folds_abun_Klebsiella pneumoniae_Others_krusk_significant_order.csv")
perma_tab = readRDS("PERMANOVA_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")$`Klebsiella pneumoniae_Others`
betadisper = readRDS("Betadispersion_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")$`Klebsiella pneumoniae_Others`
mic_df = readRDS("MIC_fold_change_table_sector_mean_vivli_rif.rds")$`Klebsiella pneumoniae_Others`


permanova_pair = readRDS("PERMANOVA_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")
betadisper_pair = readRDS("Betadispersion_pairwise_mic_sus_fold_phyloseqs_vivli_rif.rds")
phy_list_patho = readRDS("vivli_rif_mic_sus_fold_changes_phyloseqs_3Bac_patho.rds")
mic_df_all = readRDS("MIC_fold_change_table_sector_mean_vivli_rif.rds")



for (k in 1:length(phy_list_patho)){
    patho_name = names(phy_list_patho)[k]
    output_name = paste0("SIMPER_vivli_rif_AST_mic_sus_folds_abun_", patho_name)
    raw.df=read.csv(paste0(output_name,"_krusk_significant_order.csv")) # simper result
    perma = permanova_pair[[k]]   # permanova
    betadisper = betadisper_pair[[k]]   # betadispersion
    mic_df=mic_df_all[[k]] # mean of MIC fold change per group
    simper_out=simper.transform_abx_vivli (raw.df, perma, betadisper, mic_df)
    xlsx::write.xlsx(simper_out, paste0("SIMPER_AST_mic_sus_folds_abun_reformat_vivli_rif_20240914.xlsx"),showNA = F,
                     sheetName = patho_name, append=T, row.names = F)
}





