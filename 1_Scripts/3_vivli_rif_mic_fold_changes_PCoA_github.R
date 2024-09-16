
# Use MIC fold change of all antibiotic profile of each isolate to plot PCoA
# comapre community and clinic
# Note!! Only incldue antibiotics with full data in both clinic and community

# load packages
pkgs<-c("ggplot2","dplyr","tidyr","stringr","cowplot","vegan", "phyloseq","openxlsx","dplyr", "grid", "gridExtra", "geomtextpath")
sapply(pkgs, require, character.only = TRUE)

setwd("/Users/pyhuang/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_conferences/2024_Vivli AMR data challenge/")
source("~/SPH Dropbox/Pei-Yu Huang/2022_SPH/AMR_intro/AMR_analysis/00_Manuscript_1/AMR_functions_20230308.R")

# input data
input_list = readRDS("20240913_atlas_rif_HK_mic_sus_fold_change_meta_3Bac.rds")
bac_level = c("Escherichia coli", "Klebsiella pneumoniae", "Staphylococcus aureus")

#-------------------------------------------------------------
# produce phyloseq object
#-------------------------------------------------------------
physeq_ecoli = phyloseq_create_vivli(input_list[[1]], input_list[[4]])
physeq_kpneu = phyloseq_create_vivli(input_list[[2]], input_list[[5]])
physeq_saurs = phyloseq_create_vivli(input_list[[3]], input_list[[6]])

phylosq_list = c(physeq_ecoli, physeq_kpneu, physeq_saurs)

names(phylosq_list) = bac_level
saveRDS (phylosq_list, "vivli_mic_sus_fold_changes_phyloseqs.rds")


#-------------------------------------------------------------
phy_4th_list = c()

for (k in 1:length(phylosq_list)){
    phy_4th=Phyloseq_transform_vivli (phylosq_list[[k]], "4th")
    phy_4th_list = c(phy_4th_list, list(phy_4th))
}

names(phy_4th_list) = bac_level
saveRDS (phy_4th_list, "vivli_mic_sus_fold_changes_phyloseqs_4th.rds")

#-------------------------------------------------------------
# ordination
#-------------------------------------------------------------
phy_4th_list = readRDS ("vivli_mic_sus_fold_changes_phyloseqs_4th.rds")

physeq_list = phy_4th_list
transform = "abun" #"abun" occur
dist = "bray" #  jaccard bray
method = "PCoA" # PCoA, RDA PCA


for (k in 1:length(bac_level)){
    physeq_pcoa = physeq_list[[k]]
    if(is.null(physeq_pcoa)==F){
        bac_name = bac_level[k]
        sample_info = data.frame(sample_data(physeq_pcoa))
        ordination = ordinate(physeq_pcoa, method, dist)
        saveRDS(ordination, paste0( "AST_ordination_vivli_",method,"_",transform,"_", dist, "_",bac_name,".rds"))
    }
}

#-------------------------------------------------------------
# prepare dataframe for PCoA plot
#-------------------------------------------------------------
pcoa_plot_list = c()
strivar_list = c()
ordination_list = c()

for (k in 1:length(bac_level)){
    physeq_pcoa = physeq_list[[k]]

    if (is.null(physeq_pcoa)==F){
        bac_name = bac_level[k]
        sample_info = data.frame(sample_data(physeq_pcoa))

        ordination =  readRDS(paste0( "AST_ordination_vivli_",method,"_",transform,"_", dist, "_",bac_name,".rds"))
        ordination_list = c(ordination_list, list(ordination))

        data.scores <- data.frame(ordination$vectors)
        data.scores$Isolate.ID = rownames(data.scores) # occur
        sample_info$Isolate.ID = rownames(sample_info) # occur
        pcoa_dist_merge<-merge(sample_info, data.scores, by="Isolate.ID") # occur

        pcoa_dist_merge$Sector<-factor(pcoa_dist_merge$Sector, levels=c("Clinical", "Community"))
        sector_level = c("Clinical", "Community")
        #loc <- pcoa_dist_merge[pcoa_dist_merge$a_sample_host == host_level[1], ][chull(pcoa_dist_merge[pcoa_dist_merge$a_sample_host == host_level[1], c("Axis.1", "Axis.2")]), ]
        loc <- data.frame()

        for (k in 1:length(sector_level))
        {
            loc = rbind(loc, pcoa_dist_merge[pcoa_dist_merge$Sector == sector_level[k], ][chull(pcoa_dist_merge[pcoa_dist_merge$Sector == sector_level[k], c("Axis.1", "Axis.2")]), ] )
        }

        hull.data.loctime <- loc

        # extract percentage of value explained
        axes=1:2
        extract_eigenvalue.pcoa = function(ordination) ordination$values$Relative_eig
        eigvec = extract_eigenvalue.pcoa(ordination)
        fracvar = eigvec[axes]/sum(eigvec)
        percvar = round(100 * fracvar, 1)
        strivar=c(paste0("PCo1 [", percvar[1], "%]"), paste0("PCo2 [", percvar[2], "%]"))
    }
    else{
        pcoa_dist_merge=NULL
        strivar=NULL
    }
    # Plot

    # write out
    pcoa_plot_list = c(pcoa_plot_list, list(pcoa_dist_merge))
    strivar_list = c(strivar_list, list(strivar))


}


names(pcoa_plot_list) = bac_level
names(strivar_list) = bac_level
names(ordination_list) = bac_level



# --------------------------------------------------------
# PCoA plot with ellipses
# --------------------------------------------------------
p1=pcoa_plot(physeq_list[[1]], ordination_list[[1]],pcoa_plot_list[[1]], strivar_list[[1]], bac_level[1],35)
p2=pcoa_plot(physeq_list[[2]], ordination_list[[2]],pcoa_plot_list[[2]], strivar_list[[2]], bac_level[2],35)
p3=pcoa_plot(physeq_list[[3]], ordination_list[[3]],pcoa_plot_list[[3]], strivar_list[[3]], bac_level[3],35)


combined_plot<-plot_grid(p1+theme(legend.position="none"),
                         p2+ theme(legend.position="none"),
                         p3+ theme(legend.position="none"),
                                   nrow=3, align = 'h')

# when position = "bottom," legend is the third grob in the list
shared_legend <- cowplot::get_plot_component(p1, "guide-box", return_all = TRUE)[[3]]


# Display the combined plot and the legend
pdf( "PCoA_vivli_rif_20240913.pdf", width=7, height=10) # KPNEU
plot_grid(combined_plot, shared_legend, ncol = 1, rel_heights = c(2, 0.2))
dev.off()





