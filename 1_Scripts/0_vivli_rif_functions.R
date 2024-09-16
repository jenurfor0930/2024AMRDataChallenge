
#-------------------------------------------------------------
# only keep antibiotics with data from both clinic and community
#-------------------------------------------------------------
keep_both_antibiotic = function (input){
    start_pos = which(colnames(input)=="Phenotype")+1
    end_pos = ncol(input)

    abx = colnames(input)[start_pos:end_pos]

    abx_rm = c()
    for(k in 1:length(abx)){
        tmp_tab = table(input[,"Sector"], input[,abx[k]])
        tmp_tab = as.matrix(tmp_tab)
        sum=rowSums(tmp_tab)
        if (sum[1]==0 | sum[2]==0){
            abx_rm = c(abx_rm, abx[k])
        }
    }
    input1 = input %>% select(-all_of(abx_rm))
    return(input1)
}

# ------------------------------------------------------------
# add a random number to the MIC value
# ------------------------------------------------------------
# x = MIC +/- random number from normal distribution with mean is 20%% of MIC and sd is 2% of MIC (adj=0.2)
mic_numeric <-  function (tab.new, adj){
    tab.new <- tab.new[,colSums(is.na(tab.new))<nrow(tab.new)] # remove NA columns
# transform to long
    library(tidyr)
    library(dplyr)
    # end column
    ncols = ncol(tab.new)
    tab.df = tab.new %>% pivot_longer(cols=!c(Isolate.Id:Phenotype),
                                   names_to = "antibiotic", values_to = "mic")

    # remove all space from mic
    tab.df$mic = gsub(" ", "", tab.df$mic)
    # make all NA ==0.0001
    tab.df$mic = ifelse(is.na(tab.df$mic), "0.0001", tab.df$mic)

    # separate the characters and number in mic
    library(stringr)
    tab.df = tab.df %>% mutate(mic1 = str_extract(mic, "[<=>]+"))
    tab.df = tab.df %>% mutate(mic2 = as.numeric(str_extract(mic, "\\d+\\.?\\d*")))

    #  add or deduct a randomly number from normal distribution to categorical MIC (e.g. <=2, > 128)
    tmp1=tab.df%>%filter(mic1=="<=")
    tmp1= tmp1%>% mutate(mic3 = mic2 - abs(rnorm(n(), mean=mic2*adj, sd=mic2*adj*0.1)))

    tmp2=tab.df%>%filter(mic1==">")
    tmp2= tmp2%>% mutate(mic3 = mic2 + abs(rnorm(n(), mean=mic2*adj, sd=mic2*adj*0.1)))


        tmp3 = tab.df%>%filter(!mic1%in%c(">","<=") ) %>%mutate(mic3= mic2)
        tab.df1 = rbind(tmp1,tmp2, tmp3)

        return(tab.df1)
    }

# ------------------------------------------------------------
# transform mic value to fold changes compared to mic of susceptilbe
# ------------------------------------------------------------
# change from long format to wide format and transform to MIC fold changes relative to CLSI susceptible cutoffs

convert_to_fold_changes = function (input_df){
    input_df1 = input_df %>% select(Isolate.Id, antibiotic, mic3)
    input_wide = input_df1 %>% pivot_wider (names_from = antibiotic, values_from = mic3) %>%
        tibble::column_to_rownames(var = "Isolate.Id")

    cutoff.df = openxlsx::read.xlsx ("20240911_ast_mic_cutoff_susceptible_all.xlsx") # imipenem cutoff should be 2 as it's our lower limit on plate.
    cutoff.df = cutoff.df %>% select(-abx)
    rownames(cutoff.df) = cutoff.df$Antibiotic
    cutoff.df = cutoff.df %>% select(-Antibiotic)
    # tidy up colnames
    colnames(cutoff.df) = sub("\\.", " ", colnames(cutoff.df))

    abx_lst=unique(input_df1$antibiotic) # target abx
    Species = as.character(unique(input_df$Species)) # target species

    for (i in 1:length(abx_lst)){
        row_pos = which(rownames(cutoff.df)==abx_lst[i])
        col_pos = which(colnames(cutoff.df)==Species)

        cutoff=cutoff.df[row_pos, col_pos]
        if (is.na(cutoff)==F){
            input_wide[,abx_lst[i]]=input_wide[,abx_lst[i]]/cutoff
        }
    }
    return (input_wide)
}

#-------------------------------------------------------------
# prepare otu table
#-------------------------------------------------------------
abun_tab_vivli = function (input_tmp_df)
{
    if (is.null(dim(input_tmp_df))==F){
        bac_mic_tidy=data.frame(apply(input_tmp_df, 2, function(x){ifelse (is.na(x)==T, 0, x)}),check.names = FALSE) # NA value gives 0
        bac_mic_tidy1 = bac_mic_tidy[,!colSums(bac_mic_tidy)==0]
        if(nrow(bac_mic_tidy1[rowSums(bac_mic_tidy1)==0,])>0){
            bac_mic_tidy1 = bac_mic_tidy1[!rowSums(bac_mic_tidy1)==0,]
        }
    }
    require(data.table)
    tmp.df.mx=as.matrix(bac_mic_tidy1)
    tmp.df.mx.t=t(tmp.df.mx)
    return (tmp.df.mx.t)
}

#-------------------------------------------------------------
# prepare the phyloseq object for PCoA plots
#-------------------------------------------------------------
phyloseq_create_vivli = function (input_fold, input_meta){
    abun.tab = abun_tab_vivli(input_fold)
    abun.meta = input_meta
    bac = unique(abun.meta$Species)

    if (is.null(dim(abun.tab))==F){
        abun.meta$Sector = factor (abun.meta$Sector, levels=c("Clinical", "Community"))
        # meta_data
        meta_tab_phy = phyloseq::sample_data(abun.meta)
        # asv table
        asv_tab_phy<-otu_table(abun.tab, taxa_are_rows=T)
        # combine metadata, taxa table and asv table
        physeq_obj<-phyloseq(meta_tab_phy, asv_tab_phy)
        saveRDS(physeq_obj, paste0(bac, "_mic_fold_abundance_meta_phyloseq.rds"))
        return (physeq_obj)
    }

}
#-------------------------------------------------------------
# Transform mic fold change data using 4th roots
#-------------------------------------------------------------

Phyloseq_transform_vivli = function (physeq_obj, transform){
    if (is.null(physeq_obj)==F){
        abun_tab_raw=as.data.frame(otu_table(physeq_obj))
        # make every column numeric
        abun_tab_num <- sapply(abun_tab_raw, as.numeric )
        rownames(abun_tab_num) = rownames(abun_tab_raw)
        if (transform == "sqrt"){
            # square root transformation
            abun_tab<-sqrt(abun_tab_num)
        } else if (transform == "4th"){
            # Fourth root transformation
            abun_tab<-nthroot(abun_tab_num, 4)
        } else if (transform == "occur"){
            # occurrence for non-susceptibles (convert anything <=1 as 0 as fold changes <=1 are susceptibles)
            abun_tab<-ifelse(abun_tab_num>1, 1, 0)
        }
        abun_tab_trans <-otu_table(abun_tab, taxa_are_rows=T)
        phy_abun_trans <- phyloseq(abun_tab_trans,  sample_data (physeq_obj))
        return (phy_abun_trans)
    }
}

#-----------------------------------------------------
# plot PCoA
#-----------------------------------------------------
pcoa_plot = function (physeq_obj, ordination, pcoa_df, strivar, bac_full, adjust){
    library("ggrepel")  # avoid overlap text
    host_col = c("China" = "#E74C3C",
                 "Taiwan" = "#2ECC71",
                 "Hong Kong" = "#F1C40F",
                 "Community (Hong Kong)" = "#3498DB")

    pcoa_df = pcoa_df %>% mutate(Sector1 = ifelse (Country == "Hong Kong" & Sector=="Community",  "Community (Hong Kong)",Country), .after=Sector)
    pcoa_df$Sector1 = factor (pcoa_df$Sector1, levels =c("China", "Taiwan", "Hong Kong", "Community (Hong Kong)"))
    pcoa_df$Pathogen = sub("_", ' ', pcoa_df$Pathogen)


    pcoa.p<-ggplot(pcoa_df, aes(x=Axis.1, y=Axis.2)) +
        facet_grid(~Pathogen)+
        geom_point(aes(col=Sector1), size=2, alpha=0.3)+
        stat_ellipse(aes(col=Sector1), lwd = 1.2, level=0.95, alpha=0.7, linetype=1)+ #geom = "polygon",
        scale_color_manual(values=host_col)+
        scale_fill_manual(values=host_col)+
        theme_bw()+
        theme(
            axis.text.y=element_text(size=11),
            axis.text.x=element_text(size=11),
            strip.text.x = element_text(size = 12),
            legend.title=element_blank(),
            legend.position="bottom",
            legend.text = element_text(face = "plain", size=14),
            plot.title = element_text(face="bold.italic", size=15),
            axis.title = element_text(size=14)
        )+
        labs(x=strivar[1], y=strivar[2], title=bac_full)+
        coord_equal()

    return(pcoa.p)
}

#-----------------------------------------------------
# Make the simper table
#-----------------------------------------------------
simper.transform_abx_vivli=function (raw.tab,perma_tab, betadisper_tab, mic_df ){

    perma_df = perma_tab %>% select(pairs, 'p.adjusted', sig) %>% rename(Comparison = pairs)

    # reformat betadispersion df
    betadisper_df = betadisper %>% data.frame() %>% mutate(Comparison = sub("(^.*?)-(.*$)","\\2 vs \\1", rownames(betadisper))) %>%
        relocate(Comparison) %>% select(Comparison, p.adj) %>% rename(beta_p = p.adj) %>% mutate(beta_sig = case_when(beta_p<=0.01~"*",
                                                                                                                      beta_p>0.01 & beta_p<=0.05 ~ ".",
                                                                                                                      beta_p>0.05 ~ "")) %>%
        mutate(beta_p = as.numeric(beta_p)) %>%
        mutate(beta_p = ifelse(beta_p>=0.001, round(beta_p, 4), format(beta_p, digit=3) ))
    # tidy up simper table
    raw.tab = raw.tab%>% rename(Antibiotic = OTU, Class= Taxonomy) %>% mutate(Comparison = paste0(group_1," vs ", group_2), .before =group_1)%>% mutate(SIMPER_per=round(SIMPER*100, 1)) %>% select (-SIMPER)

    mic_df_g1 = mic_df %>%pivot_longer(-c(Class,Antibiotic), names_to = "group_1", values_to = "L_Values")
    mic_df_g2 = mic_df %>%pivot_longer(-c(Class,Antibiotic), names_to = "group_2", values_to = "R_Values")

    raw.tab = merge(mic_df_g1, raw.tab) %>% relocate(L_Values, .after=SIMPER_per)
    raw.tab = merge(mic_df_g2, raw.tab)%>% relocate(R_Values, .after=SIMPER_per)

    # add annotation to larger cell
    raw.tab1 = raw.tab %>% mutate(SIMPER_per_annot = ifelse(L_Values>R_Values, paste0("*", SIMPER_per), SIMPER_per)) %>%
        select(-c(L_Values, R_Values,SIMPER_per,  Left.mean.abund, Right.mean.abund, group_2, group_1))

    # sum of simper
    simper_sum = raw.tab %>% select(Comparison, Antibiotic, SIMPER_per) %>% group_by (Comparison) %>% summarise(Sum=sum(SIMPER_per))

    ast_col=openxlsx::read.xlsx("vivli_arrow_abbr.xlsx")
    abx_level = ast_col$Antibiotic
    cls_level = unique(ast_col$Class)

    raw.tab1$Antibiotic=factor(raw.tab1$Antibiotic, levels=abx_level)
    raw.tab1$Class=factor(raw.tab1$Class, levels=cls_level)

    raw.tab1 = raw.tab1 %>% arrange(Class, Antibiotic)
    raw.tab2 = raw.tab1 %>% group_by (Comparison) %>% pivot_wider(names_from = c(Class, Antibiotic), values_from = SIMPER_per_annot)

    raw.tab3 = merge(perma_df,raw.tab2 ) %>% rename(permanova_p = p.adjusted, perma_sig = sig) # add permanova
    raw.tab3 = merge(betadisper_df, raw.tab3) # add permanova
    raw.tab4 = merge(raw.tab3,simper_sum ) # add SIMPER  sum

    cls=sub("(^.*?)_(.*$)","\\1", colnames(raw.tab4))
    abx=sub("(^.*?)_(.*$)","\\2", colnames(raw.tab4))
    raw.tab4 = rbind(cls, abx, data.frame(raw.tab4))
    return(raw.tab4)
}
