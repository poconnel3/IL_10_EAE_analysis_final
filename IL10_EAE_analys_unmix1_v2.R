#By: Patrick O'Connell
# July 10, 2021

#This analysis starts from FCS files and goes from there. W/ downsampling. 
#only ERAP1 and WT right now.
#V2 includes an extra IL-10-GFP WT EAE mouse

#Packages
library(CATALYST)
library(tidyverse)
library(flowCore)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(diffcyt)


##-------------------------Making metadata file----------------------------------##
#Load the content of the metadata.xlsx file into R (file name must match name of the fcs file you generate for each sample)
setwd("~/IL10_EAE_R_analys_unmix1_v2")
metadata_filename <- "eae_metadata_altX.xlsx"
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$condition <- factor(md$condition, levels = c("EAE", "naive"))
md$IV_pos <- factor(md$IV_pos, levels = c("yes", "no"))
md$IL_10_pos <- factor(md$IL_10_pos, levels = c("yes", "no"))
md$genotype <- factor(md$genotype, levels = c("ERAP1_KO", "WT", "SLAMF7_KO"))
head(data.frame(md))
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/IL10_EAE_R_analys_unmix1_SF7/clean_il10_EAE_fcs")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/IL10_EAE_R_analys_unmix1_v2")
panel <- "neuroimmune_panel_IL10.xlsx"
panel <- read_excel(panel) 

#check
all(panel$fcs_colname %in% colnames(fcs_raw))
##------------------making panel file-----------------------##

#see how many cells per file
fsApply(fcs_raw, nrow)

#downsample to 100,000 cells per sample
# Define a downsampling ceiling
sampling.ceiling <- 100000
# Being reproducible is a plus
set.seed(666)

# BUILD A DOWNSAMPLED FLOWSET
fcs_raw.dsamp <- fsApply(fcs_raw, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling.ceiling, nrow(ff)))
  ff[idx,]  
})

#check it
fcs_raw.dsamp
fsApply(fcs_raw.dsamp, nrow)

#make Catalyst data object and arcsinh t-for all channels by 6000

#First 3 must be present and sample_id must match condition (i think)
factors <- list(factors = c("file_name", "condition", "sample_id", "genotype", "IV_pos", "IL_10_pos"))

daf_X <- prepData(
  fcs_raw.dsamp,
  features = NULL,
  md = md,
  md_cols = factors,
  panel = panel,
  transform = TRUE,
  cofactor = 6000,
  by_time = FALSE,
  FACS = TRUE
)


# view number of events per sample
table(daf_X$condition)
table(daf_X$sample_id)
table(daf_X$file_name)
table(daf_X$genotype)
table(daf_X$IV_pos)
table(daf_X$IL_10_pos)

# view non-mass channels
names(int_colData(daf_alt_clean))

#view metadata parameters
names(colData(daf_X))

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")    

#plot global marker expression by condition
p <- plotExprs(daf_X, color_by = "condition")
p$facet$params$ncol <- 6                   
p 

#plot number of cells per condition
plotCounts(daf_X, color_by = "IV_pos")
table(daf$condition) 

#remove SF7 samples
daf_X <- filterSCE(daf_X, genotype != "SLAMF7_KO")

#Make MDS plot  
pbMDS(daf_X, by = "sample_id", color_by = "condition")


#Make MDS plot  clean 6/14/2021
pbMDS(daf_X, by = "sample_id", color_by = "genotype", shape_by = "condition", label_by = NULL, size_by = T, pal = c("#ff1a1a", "#0d0d0d"))

#Heatmap of all markers per sample_id and condition w/ heiracheal clustering.
plotExprHeatmap(daf_X, 
                bin_anno = TRUE, 
                scale = "first",
                row_anno = "condition"
                )

#Perform FlowSOM clustering
daf_X <- cluster(daf_X, 
               features = "type", 
               xdim = 10, 
               ydim = 10, 
               maxK = 22, 
               seed = 1234, 
               verbose = TRUE
               ) 

#plot cluster heatmap
plotClusterHeatmap(daf_X,                               
                   hm2 = NULL, #swap for "state" or "abundances" to get more info. 
                   k = "meta22", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
                   ) 

#plot marker expression by cluster
plotClusterExprs(daf2, k = "meta20")  

# run UMAP                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
             dr = "UMAP", 
             cells = 10000,
             features = "type" 
             )

# run TSNE                          
set.seed(777)                               
daf_X <- runDR(daf_X, 
               dr = "TSNE", 
               cells = 10000,
               features = "type" 
)

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#UMAP visualizations
plotDR(daf_X, "UMAP", color_by = "meta22")

#tSNE visualizations
plotDR(daf_X, "TSNE", color_by = "meta22")

daf_X_naive <- filterSCE(daf_X, condition == "naive")
pz <- plotDR(daf_X_naive, "UMAP", color_by = "IL_10_pos", facet_by = "genotype")               
pz + theme(axis.line = element_line(colour = NA),   #use this theme for all plots
    axis.ticks = element_line(colour = NA), 
    panel.grid.major = element_line(linetype = "blank"), 
    panel.grid.minor = element_line(linetype = "blank"), 
    axis.title = element_text(colour = NA), 
    axis.text = element_text(colour = NA), 
    legend.text = element_text(size = 10), 
    legend.title = element_text(size = 14), 
    panel.background = element_rect(fill = NA)) 

#make UMAP colored by important markers 
plotDR(daf_X, "UMAP", color_by = c("CD38", "NK1.1", "CD8", "CD4", "B220", "CCR2", "CD11c", "Ly6G", "MHCII"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#make UMAP colored by important markers 
plotDR(daf_X, "UMAP", color_by = c("IL-10"), facet_by = "genotype", a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#---------------plot IL10+ by cluster for only naive and only eae (WT only) 6/14/2021
daf_naive_wt <- filterSCE(daf_X_naive, genotype == "WT", IV_pos == "no")
#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
Il10_clust <- table(
  Il10_pos = daf_naive_wt$IL_10_pos, 
  cluster = cluster_ids(daf_naive_wt, "merging1")) 

plot_clust_table <- prop.table(Il10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
rr <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=Il10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rr + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)

daf_eae_wt <- filterSCE(daf_X, genotype == "WT", condition == "EAE", IV_pos == "no")
#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
Il10_clust <- table(
  Il10_pos = daf_eae_wt$IL_10_pos, 
  cluster = cluster_ids(daf_eae_wt, "merging1")) 

plot_clust_table <- prop.table(Il10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
rr <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=Il10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rr + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##---------------------------
  

## Facet per sample                                          
plotDR(daf, "UMAP", color_by = "meta20", facet = "sample_id")

## Facet per condition                                       
plotDR(daf_alt2, "UMAP", color_by = "meta20", facet = "IV_pos")

#plot FlowSOM codes to visualize similarity of clusters
plotCodes(daf_X, k = "meta22")

#Merge and annotate clusters
merging_tableX <- "merged_clusters_X.xlsx"
merging_table1 <- read_excel(merging_tableX)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("Microglia", "BAMs", "MdCs", "pDCs", "DCs", "CCR2+ monocytes", "Neutrophils",             
                                                "CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Debris"))        

daf_X <- mergeClusters(daf_X, k = "meta22",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster
#make a new SCE object when you do this!!
daf_X <- filterSCE(daf_X, cluster_id != "Debris", k = "merging1")

#change cluster colors
clust_color <- c("#996600", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC")

#UMAP of merged clusters
plotDR(daf_X, "UMAP", color_by = "merging1", k_pal = clust_color) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
                                                                axis.ticks = element_line(colour = NA), 
                                                                panel.grid.major = element_line(linetype = "blank"), 
                                                                panel.grid.minor = element_line(linetype = "blank"), 
                                                                axis.title = element_text(colour = NA), 
                                                                axis.text = element_text(colour = NA), 
                                                                legend.text = element_text(size = 10), 
                                                                legend.title = element_text(size = 14), 
                                                                panel.background = element_rect(fill = NA)) 

#heatmap of merged cluster markers
to_remove <- c("CD3")
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help

plotExprHeatmap(daf_X, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                hm_pal = brewer.pal(9, "Greys"),
                features = c("SLAMF7", "CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b","CD45", 
                             "B220", "CD19", "SiglecH", "CCR2", "NK1.1", "CD11c", "Ly6G", "MHCII", "Lyve1", "CD206", "IgD")
)


#ridgeplots of markers on merged clusters
plotClusterExprs(daf_alt_clean_rerun, k = "merging1")

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#-------------------------------compare groups now (w/o removing IV or IL-10 cell groups------------------------------##
library(diffcyt)

#remove all naive samples
daf_alt <- filterSCE(daf_X, condition != "naive")

#plot abundances as barplot
plotAbundances(daf_alt, k = "merging1", by = "sample_id", k_pal = clust_color)


#plot abundances as boxplot
plotAbundances(daf_alt, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
               )


##To analyze now we can generate a GLMM
ei <- metadata(daf_alt)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

#how to access merged cluster levels
levels(cluster_ids(daf_alt, k = "merging1"))

#cluster levels of FlowSOM clusters
levels(cluster_ids(daf_X))

da_res1 <- diffcyt(daf_alt,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
names(da_res1)

cluster_stats <- rowData(da_res1$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)

#Show all results
topTable(da_res1, show_props = TRUE, format_vals = TRUE, digits = 2)

#plot expression of markers by cluster (boxplot)
plotPbExprs(daf_alt,
            k = "merging1",
            features = NULL,
            facet_by = "cluster_id",
            color_by = "genotype",
            group_by = "genotype",
            geom = "both", 
            jitter = TRUE
                   )



##Statistically compare markers on clusters w/ GLMM. EAE only
#upd 6/14/2021 I am removing all IV+ cells and setting more stringent FDR
daf_no_IV_eae <- filterSCE(daf_X, IV_pos == "no", condition == "EAE")
rowData(daf_no_IV_eae)$marker_class <- "state"

ei <- metadata(daf_no_IV_eae)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_no_IV_eae), cols_design = "genotype")

da_res99 <- diffcyt(daf_no_IV_eae,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    #markers_to_test = all_mark,  #works, except I can't specify which markers using this command
                    clustering_to_use = "merging1", verbose = FALSE)  


#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)
FDR_cutoff <- 0.01
table(rowData(da_res99$res)$p_adj < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(daf_no_IV_eae, tbl_DS, fdr = 0.01, all = T,
                sort_by = "padj", col_anno = "genotype")


##Statistically compare markers on clusters w/ GLMM. naive only
#upd 6/14/2021 I am removing all IV+ cells and setting more stringent FDR
daf_no_IV_naive <- filterSCE(daf_X, IV_pos == "no", condition == "naive")
rowData(daf_no_IV_naive)$marker_class <- "state"

ei <- metadata(daf_no_IV_naive)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_no_IV_naive), cols_design = "genotype")

da_res99 <- diffcyt(daf_no_IV_naive,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    #markers_to_test = all_mark,  #works, except I can't specify which markers using this command
                    clustering_to_use = "merging1", verbose = FALSE)  


#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)
FDR_cutoff <- 0.01
table(rowData(da_res99$res)$p_adj < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(daf_no_IV_naive, tbl_DS, fdr = 0.01, all = T,
                sort_by = "padj", col_anno = "genotype")

##------------THIS VERSION OF THE CODE WORKS FOR THIS ANALYSIS (SAVE FOR FUTURE USE)--------------------------##

##Statistically compare markers on clusters w/ GLMM. 
ei <- metadata(daf_alt)$experiment_info

ds_formula2 <- createFormula(ei, cols_fixed = "genotype")  #NOTE: cannot have a rondom column variable here or else it fails for some reason. 

contrast <- createContrast(c(0, 1))

FDR_cutoff <- 0.05

ds_res2 <- diffcyt(daf_alt, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

cluster_stats4 <- rowData(ds_res2$res) %>%
  as.data.frame(.)
#export stat data
library(xlsx)
write.xlsx(cluster_stats4, "~/IL10_EAE_R_analys_unmix1_v2/EAE_clust_stats.xlsx")

##tbl_DS <- rowData(da_res99$res)
plotDiffHeatmap(daf_alt, rowData(ds_res2$res), top_n = 50, fdr = FDR_cutoff, lfc = 2, col_anno = c("genotype"), row_anno = TRUE) 


#now naive mice
##Statistically compare markers on clusters w/ GLMM. 
daf_naive <- filterSCE(daf_X, condition == "naive")
rowData(daf_naive)$marker_class <- "state"
ei <- metadata(daf_naive)$experiment_info

ds_formula2 <- createFormula(ei, cols_fixed = "genotype")  #NOTE: cannot have a rondom column variable here or else it fails for some reason. 

contrast <- createContrast(c(0, 1))

FDR_cutoff <- 0.05

ds_res2 <- diffcyt(daf_naive, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

table(rowData(ds_res2$res)$p_adj < FDR_cutoff)

cluster_stats5 <- rowData(ds_res2$res) %>%
  as.data.frame(.)
#export stat data
write.xlsx(cluster_stats5, "~/IL10_EAE_R_analys_unmix1_v2/naive_clust_stats.xlsx")

tbl_DS <- rowData(ds_res2$res)
plotDiffHeatmap(daf_naive, rowData(ds_res2$res), top_n = 50, fdr = FDR_cutoff, lfc = 2, col_anno = c("genotype"), row_anno = TRUE) 

##----------------------see what fraction of each cluster/sample are IV+-----------------------##
#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_X, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_X, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
daf_naive_WT_for_iv <- filterSCE(daf_X, genotype == "WT", condition == "naive")
IV_clust <- table(
  IV_pos = daf_naive_WT_for_iv$IV_pos, 
  cluster = cluster_ids(daf_naive_WT_for_iv, "merging1")) 

plot_clust_table <- prop.table(IV_clust, "cluster") %>%
  as.data.frame(.)

#make plot
rr <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
rr + theme(axis.ticks = element_line(colour = "black", 
    size = 1), axis.title = element_text(size = 15), 
    axis.text = element_text(colour = "black", 
        hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
    axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
    )

#now EAE
daf_eae_WT_for_iv <- filterSCE(daf_X, genotype == "WT", condition == "EAE")
#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IV_clust <- table(
  IV_pos = daf_eae_WT_for_iv$IV_pos, 
  cluster = cluster_ids(daf_eae_WT_for_iv, "merging1")) 

plot_clust_table <- prop.table(IV_clust, "cluster") %>%
  as.data.frame(.)

#make plot
rr <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
rr + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)

##---------------------------------------------------------------------------------------------##

#I'm going to leave the IV+ cells in for now. we can see what are the major immune cell and IL-10+ cell
#contributions to the genotypes in our dataset and can remove later should we choose. 

##---------------------Now lets look at the IL-10+ compartment during EAE----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_X, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_X$IL_10_pos, 
  cluster = cluster_ids(daf_X, "merging1")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cc <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
cc + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
}
##--------------------------------------------------------------------------------------------##


##---------------------Now lets reconstruct the IL-10+ compartment at steady state----------------------##
#1/25/2021 I re-did this w/ removal of IV_pos cells. 
#remove EAE samples
daf_IL10 <- filterSCE(daf_X, condition == "naive")
#remove IL_10_neg samples
daf_IL10 <- filterSCE(daf_IL10, IL_10_pos == "yes")
#remove IV_pos samples
daf_IL10 <- filterSCE(daf_IL10, IV_pos == "yes")

levels(daf_IL10$condition)
levels(daf_IL10$IL_10_pos)
levels(daf_IL10$sample_id)
levels(daf_IL10$IV_pos)

#update cluster colors since not all clusters are here
clust_color2 <- c("#996600", "#6665FD", "#00994D", "#FF97A7", "#006633", "#FF65CA", "#DD97FC", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC")

#UMAP of IL-10+ by genotype
plotDR(daf_IL10, "UMAP", color_by = "merging1", k_pal = clust_color2, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
        geom_point(size=2)

#comapre cell numbers
n_cells(daf_IL10)
n_cells(daf_X)

IL10_by_IV <- table(
  IL10_pos = daf_IL10$IL_10_pos, 
  IV_pos = daf_IL10$IV_pos, 
  genotype = daf_IL10$genotype) 


#pie chart of cluster contributions to IL-10 compartment in WT and ERAP1-KO (WT first)
to_remove <- c("BAMs", "Neutrophils", "CCR2+ monocytes", "MdCs", "DCs", "Act. Microglia")
`%notin%` <- Negate(`%in%`) #make this opposite of %in% operator to help
IL10_pie <- table(
  genotype = daf_IL10$genotype, 
  cluster = cluster_ids(daf_IL10, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT") %>%
  dplyr::filter(., cluster %notin% to_remove)

#make plot
hh <- ggplot(data=IL10_pie, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#FF97A7", "#006633", "#FF65CA", "#DD97FC", "#6665FD", "#00994D", "#FF97A7", "#006633", "#FF65CA"))
hh

#now ERAP1-KO
IL10_pie_erap <- table(
  genotype = daf_IL10$genotype, 
  cluster = cluster_ids(daf_IL10, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO") %>%
  dplyr::filter(., cluster %notin% to_remove)

#make plot
hz <- ggplot(data=IL10_pie_erap, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#FF97A7", "#006633", "#FF65CA", "#DD97FC", "#6665FD", "#00994D", "#FF97A7", "#006633", "#FF65CA"))
hz

##Now compare cluster IL-10+ cluster frequency b/w WT and ERAP1-KO
library(diffcyt)
eii <- metadata(daf_IL10)$experiment_info 
(da_formula11 <- createFormula(eii,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res11 <- diffcyt(daf_IL10,                                            
                   formula = da_formula11, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res11)

cluster_stats <- rowData(da_res11$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res11$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------------##

##---------------------Now lets reconstruct the IL-10+ compartment during EAE----------------------##
#remove naive samples
daf_IL10_eae <- filterSCE(daf_X, condition == "EAE")
#remove IL_10_neg samples and IV+  redidi this on 6/14/21 to remove IV+
daf_IL10_eae <- filterSCE(daf_IL10_eae, IL_10_pos == "yes", IV_pos == "no")

levels(daf_IL10_eae$condition)
levels(daf_IL10_eae$IL_10_pos)
levels(daf_IL10_eae$sample_id)

#UMAP of IL-10+ by genotype
plotDR(daf_IL10_eae, "UMAP", color_by = "merging1", k_pal = clust_color, facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=1)

#comapre cell numbers
n_cells(daf_IL10_eae)
n_cells(daf_X)
table(cluster_ids(daf_IL10_eae, k = "merging1"))

IL10_eae_by_IV <- table(
  IL10_pos = daf_IL10_eae$IL_10_pos, 
  IV_pos = daf_IL10_eae$IV_pos, 
  genotype = daf_IL10_eae$genotype) 

#pie chart of cluster contributions to IL-10 compartment in WT and ERAP1-KO (WT first)
IL10_pie_eae <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT") 

#make plot
hf <- ggplot(data=IL10_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
hf

#now ERAP1-KO
{IL10_pie_eae_erap <- table(
  genotype = daf_IL10_eae$genotype, 
  cluster = cluster_ids(daf_IL10_eae, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO")}

#make plot
fg <- ggplot(data=IL10_pie_eae_erap, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
fg

#-------------------now WT and erap during EAE no IV 6/14/2021
#remove naive samples
daf_IL10_eae_2 <- filterSCE(daf_X, condition == "EAE", IV_pos == "no", IL_10_pos == "yes")

#pie chart of cluster contributions to IL-10 compartment in WT and ERAP1-KO (WT first)
IL10_eae_by_IV <- table(
  genotype = daf_IL10_eae_2$genotype, 
  cluster = cluster_ids(daf_IL10_eae_2, "merging1")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO") 

#make plot
hf <- ggplot(data=IL10_eae_by_IV, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() + 
  scale_fill_manual(values = c("#996600", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#FF9864", "#FF97A7", "#006633", "#FF65CA", "#DD97FC"))
hf

#compare clust freq b/w steady state and eae (WT only)
daf_IL10_eae_wt <- filterSCE(daf_X, genotype == "WT", IV_pos == "no", IL_10_pos == "yes")

eii <- metadata(daf_IL10_eae_wt)$experiment_info 
(da_formula11 <- createFormula(eii,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res11 <- diffcyt(daf_IL10_eae_wt,                                            
                    formula = da_formula11, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res11)

cluster_stats434 <- rowData(da_res11$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res11$res)$p_adj < FDR_cutoff)

#plot abundances as boxplot
obj66 <- plotAbundances(daf_IL10_eae_wt, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "condition",
               shape = NULL
) 
obj66 + scale_color_manual(values = c("#a6a6a6", "#000000")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)
#-------------------

##Now compare cluster IL-10+ cluster frequency b/w WT and ERAP1-KO
eie <- metadata(daf_IL10_eae)$experiment_info 
(da_formula12 <- createFormula(eie,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res12 <- diffcyt(daf_IL10_eae,                                            
                    formula = da_formula12, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "merging1", verbose = FALSE)    

#examine output
names(da_res12)

cluster_stats2 <- rowData(da_res12$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res12$res)$p_adj < FDR_cutoff)


##-----------------------------------------------------------------------------------##

##---------------------Now subset NK cells from complete exp and analyze----------------------##
#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#subset NK cells
daf_NK <- filterSCE(daf_X, cluster_id == "NK cells", k = "merging1")

levels(daf_NK$condition)
levels(daf_NK$IL_10_pos)
levels(daf_NK$sample_id)

#Make MDS plot  
pbMDS(daf_NK, by = "sample_id", color_by = "condition", shape_by = "genotype", size_by = TRUE)

#re-cluster just the NK cells
daf_NK <- cluster(daf_NK, 
                 features = NULL, 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 8, 
                 seed = 7689, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_NK,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta8", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(777)                               
daf_NK <- runDR(daf_NK, 
               dr = "UMAP", 
               cells = 10000,
               features = NULL 
)

#UMAP NK cells
plotDR(daf_NK, "UMAP", color_by = "meta8", facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#Merge and annotate NK clusters
merging_table_NK <- "merged_clusters_NK.xlsx"
merging_table2 <- read_excel(merging_table_NK)
head(data.frame(merging_table2))  

# convert to factor with merged clusters in desired order                
merging_table2new_cluster <- factor(merging_table2$new_cluster,         
                                     levels = c("NK_1", "NK_2", "NK_3", "NK_4", "NK_5", "NK_6", "NK_7"))        

daf_NK <- mergeClusters(daf_NK, k = "meta8",                                  
                       table = merging_table2, 
                       id = "NK_clusters",
                       overwrite = TRUE
)   


#save objects
saveRDS(daf_NK, file = "./daf_NK.rds")

#change cluster colors
clust_color4 <- c("#996600", "#00CCAA", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276")

#UMAP NK cells
plotDR(daf_NK, "UMAP", color_by = "NK_clusters", facet_by = "IL_10_pos") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged NK cluster markers
plotExprHeatmap(daf_NK, 
                bin_anno = FALSE, 
                scale = "first",
                k = "NK_clusters",
                by = "cluster_id", 
                features = c("SLAMF7", "CD38", "CD49b", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD11c", "AF", "NK1.1", "CD3"),
                hm_pal = brewer.pal(9, "Greys")
)


#make UMAP colored by important markers on just NK cells
plotDR(daf_NK, "UMAP", color_by = c("CD38", "CD49b", "Ly6C", "CD90"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#------------------compare NK cluster frequencies by genotype (naive first)-----------------------------#

#remove EAE cells
daf_NK_naive <- filterSCE(daf_NK, condition == "naive")

#plot abundances as barplot
plotAbundances(daf_NK_naive, k = "NK_clusters", by = "sample_id")


#plot abundances as boxplot
plotAbundances(daf_NK_naive, 
               k = "NK_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##To analyze now we can generate a GLMM
eil <- metadata(daf_NK_naive)$experiment_info 
(da_formula13 <- createFormula(eil,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res13 <- diffcyt(daf_NK_naive,                                            
                   formula = da_formula13, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "NK_clusters", verbose = FALSE)  

#examine output
cluster_stats3 <- rowData(da_res13$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res13$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------##


#------------------compare NK cluster frequencies by genotype (EAE)-----------------------------#
#remove naive cells
daf_NK_eae <- filterSCE(daf_NK, condition == "EAE")

#plot abundances as barplot
plotAbundances(daf_NK_eae, k = "NK_clusters", by = "sample_id")


#plot abundances as boxplot
plotAbundances(daf_NK_eae, 
               k = "NK_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##To analyze now we can generate a GLMM
eim <- metadata(daf_NK_eae)$experiment_info 
(da_formula14 <- createFormula(eim,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res14 <- diffcyt(daf_NK_eae,                                            
                    formula = da_formula14, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "NK_clusters", verbose = FALSE)  

#examine output
cluster_stats4 <- rowData(da_res14$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res14$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##

##----------------------see what fraction of NK cells from each cluster/sample are IV+-----------------------##
#load main analysis object back in
daf_NK <- readRDS(file = "./daf_NK.rds")        

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_NK, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IV_clust_nk <- table(
  IV_pos = daf_NK$IV_pos, 
  cluster = cluster_ids(daf_NK, "NK_clusters")) 

plot_clust_table <- prop.table(IV_clust_nk, "cluster") %>%
  as.data.frame(.)

#make plot
rg <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rg + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##---------------------------------------------------------------------------------------------##


##---------------------Now lets see what fraction of each cluster are IL_10+----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_NK, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
IL10_clust <- table(
  IL10_pos = daf_NK$IL_10_pos, 
  cluster = cluster_ids(daf_NK, "NK_clusters")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cn <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#FF9864", "#000000"))

#make purrty
cn + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)
##--------------------------------------------------------------------------------------------##


##----------------------Compare IL-10+ NK cells b/w genotypes by IV+ and cluster-----------------------##
#load main analysis object back in
daf_NK <- readRDS(file = "./daf_NK.rds")  

#start w/ naive
#remove IL-10 neg cells
daf_NK_IL10 <- filterSCE(daf_NK, IL_10_pos == "yes")

#remove EAE cells
daf_NK_IL10 <- filterSCE(daf_NK_IL10, condition == "naive")

##NOTE: there are hardly any IL-10+ NK cells in naive. skip this. 

#now during EAE
#remove IL-10 neg cells
daf_NK_IL10 <- filterSCE(daf_NK, IL_10_pos == "yes")

#remove naive cells
daf_NK_IL10_eae <- filterSCE(daf_NK_IL10, condition == "EAE")

#pie chart of cluster contributions to IL-10+ NK cells in WT and ERAP1-KO (naive) (wt first)
IL10_NK_pie_eae <- table(
  genotype = daf_NK_IL10_eae$genotype, 
  cluster = cluster_ids(daf_NK_IL10_eae, "NK_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hk <- ggplot(data=IL10_NK_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CC0000", "#FF9999", "#0044CC", "#66CBFD", "#660099", "#BB97FC", "#FF7632"))
hk

#pie chart of cluster contributions to IL-10+ NK cells in WT and ERAP1-KO (naive) (erap1 now)
IL10_NK_pie_eae2 <- table(
  genotype = daf_NK_IL10_eae$genotype, 
  cluster = cluster_ids(daf_NK_IL10_eae, "NK_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO")

#make plot
hl <- ggplot(data=IL10_NK_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CC0000", "#FF9999", "#0044CC", "#66CBFD", "#660099", "#BB97FC", "#FF7632"))
hl

#plot abundances as boxplot
plotAbundances(daf_NK_IL10_eae, 
               k = "NK_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##Now compare cluster NK IL-10+ cluster frequency b/w WT and ERAP1-KO
library(diffcyt)
eiz <- metadata(daf_NK_IL10_eae)$experiment_info 
(da_formula17 <- createFormula(eiz,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res17 <- diffcyt(daf_NK_IL10_eae,                                            
                    formula = da_formula17, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "NK_clusters", verbose = FALSE)    

#examine output
cluster_stats7 <- rowData(da_res17$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res17$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------------##



##---------------------Now subset T cells from complete exp and analyze----------------------##
#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#subset T cells
daf_T <- filterSCE(daf_X, cluster_id %in% c("CD8+ T cells", "CD4+ T cells"), k = "merging1")

levels(daf_T$condition)
levels(daf_T$IL_10_pos)
levels(daf_T$sample_id)

#Make MDS plot  
pbMDS(daf_T, by = "sample_id", color_by = "condition", shape_by = "genotype", size_by = TRUE)

#re-cluster just the T cells     
daf_T <- cluster(daf_T, 
                  features = NULL, 
                  xdim = 10, 
                  ydim = 10, 
                  maxK = 15, 
                  seed = 7689, 
                  verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_T,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta15", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(797)                               
daf_T <- runDR(daf_T, 
                dr = "UMAP", 
                cells = 10000,
                features = NULL 
)

#UMAP T cells
plotDR(daf_T, "UMAP", color_by = "meta15", facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#Merge and annotate T clusters
merging_table_T <- "merged_clusters_T.xlsx"
merging_table3 <- read_excel(merging_table_T)
head(data.frame(merging_table3))  

# convert to factor with merged clusters in desired order                
merging_table3$new_cluster <- factor(merging_table3$new_cluster,         
                                    levels = c("Debris", "CD4 1", "CD4 2", "CD4 IL-10 1", "CD4 IL-10 2", "CD8 1", 
                                               "CD8 2", "CD8 3", "CD8 4", "NK cells", "B cells", "DN T cells"))        

daf_T <- mergeClusters(daf_T, k = "meta15",                                  
                        table = merging_table3, 
                        id = "T_clusters",
                        overwrite = TRUE
)   

#remove debris cluster and B cells and NK cells
#make a new SCE object when you do this!!
to_remove2 <- c("Debris", "B cells", "NK cells")
daf_T_clean <- filterSCE(daf_T, cluster_id %notin% to_remove2, k = "T_clusters")

#save objects
saveRDS(daf_T, file = "./daf_T.rds")
saveRDS(daf_T_clean, file = "./daf_T_clean.rds")
daf_T_clean <-  readRDS("./daf_T_clean.rds")

#change cluster colors
clust_color5 <- c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                  "#CC0022", "#FF9999", "#CC00CC", "#636363")

#UMAP T cells w/ color by markers
plotDR(daf_T_clean, "UMAP", color_by = c("CD90", "Ly6C", "CD4", "CD8", "CD38", "CD11b", "SLAMF7", "IL-10"), k_pal = clust_color5, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#UMAP T cells 
plotDR(daf_T_clean, "UMAP", color_by = "T_clusters", facet_by = "IL_10_pos") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged T cluster markers
plotExprHeatmap(daf_T_clean, 
                bin_anno = FALSE, 
                scale = "first",
                k = "T_clusters",
                by = "cluster_id", 
                k_pal = clust_color5,
                features = c("SLAMF7", "CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b"),
                hm_pal = brewer.pal(9, "Greys")
)

#------------------compare T cluster frequencies by genotype (naive first)-----------------------------#
#remove EAE cells    
daf_T_clean_naive <- filterSCE(daf_T_clean, condition == "naive")

#plot abundances as barplot
plotAbundances(daf_T_clean_naive, k = "T_clusters", by = "sample_id", k_pal = clust_color5)


#plot abundances as boxplot
plotAbundances(daf_T_clean_naive, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color5,
               shape = NULL
)


##To analyze now we can generate a GLMM
library(diffcyt)
eik <- metadata(daf_T_clean_naive)$experiment_info 
(da_formula19 <- createFormula(eik,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res19 <- diffcyt(daf_T_clean_naive,                                            
                    formula = da_formula19, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE)  

#examine output
cluster_stats9 <- rowData(da_res19$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res19$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------##


#------------------compare T cluster frequencies by genotype (EAE)-----------------------------#
daf_T_clean <- readRDS(file = "./daf_T_clean.rds")
#remove naive cells
daf_T_clean_eae <- filterSCE(daf_T_clean, condition == "EAE")
#and IV+
daf_T_clean_eae_IV <- filterSCE(daf_T_clean, condition == "EAE", IV_pos == "no")


#plot abundances as barplot
plotAbundances(daf_T_clean_eae, k = "T_clusters", by = "sample_id", k_pal = clust_color5)


#plot abundances as boxplot
plotAbundances(daf_T_clean_eae_IV, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color5,
               shape = NULL
) + scale_color_manual(values = c("#0099cc", "#ffbf00")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA) 


##To analyze now we can generate a GLMM
eib <- metadata(daf_T_clean_eae_IV)$experiment_info 
(da_formula20 <- createFormula(eib,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res20 <- diffcyt(daf_T_clean_eae_IV,                                            
                    formula = da_formula20, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE)  

#examine output
cluster_stats40 <- rowData(da_res20$res) %>%
  as.data.frame(.)
FDR_cutoff <- 0.05
table(rowData(da_res20$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##

##----------------------see what fraction of T cells from each cluster/sample are IV+-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_T_clean, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IV_clust_T <- table(
  IV_pos = daf_T_clean$IV_pos, 
  cluster = cluster_ids(daf_T_clean, "T_clusters")) 

plot_clust_table <- prop.table(IV_clust_T, "cluster") %>%
  as.data.frame(.)

#make plot
rn <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
rn + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)}
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_T_clean, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_T_clean$IL_10_pos, 
  cluster = cluster_ids(daf_T_clean, "T_clusters")) 

plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
  as.data.frame(.)

#make plot
cv <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("#00CC22", "#000000"))

#make purrty
cv + theme(axis.ticks = element_line(colour = "black", 
                                     size = 1), axis.title = element_text(size = 15), 
           axis.text = element_text(colour = "black", 
                                    hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
           axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
)}
##--------------------------------------------------------------------------------------------##


##----------------------Compare IL-10+ T cells b/w genotypes by IV+ and cluster-----------------------##
#start w/ naive
#remove IL-10 neg cells
daf_T_clean_IL10 <- filterSCE(daf_T_clean, IL_10_pos == "yes")

#remove EAE cells
daf_T_clean_IL10 <- filterSCE(daf_T_clean_IL10, condition == "naive")

#pie chart of cluster contributions to IL-10+ T cells in WT and ERAP1-KO (naive) (wt first)
IL10_T_pie_naive <- table(
  genotype = daf_T_clean_IL10$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hj <- ggplot(data=IL10_T_pie_naive, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#FF9999", "#CC00CC", "#636363", "#BB32FE", 
                               "#CC0022", "#FF9999"))
hj

#now ERAP1-KO
IL10_T_pie_naive <- table(
  genotype = daf_T_clean_IL10$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO")

#make plot
hp <- ggplot(data=IL10_T_pie_naive, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#FF9999", "#CC00CC", "#636363", "#BB32FE", 
                               "#CC0022", "#FF9999"))
hp

#plot abundances as boxplot
plotAbundances(daf_T_clean_IL10, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)

##Now compare cluster T IL-10+ cluster frequency b/w WT and ERAP1-KO (naive)
eiq <- metadata(daf_T_clean_IL10)$experiment_info 
(da_formula21 <- createFormula(eiq,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res21 <- diffcyt(daf_T_clean_IL10,                                            
                    formula = da_formula21, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE)    

#examine output
cluster_stats21 <- rowData(da_res21$res) %>%
  as.data.frame(.)

table(rowData(da_res21$res)$p_adj < FDR_cutoff)

#now during EAE
#remove IL-10 neg cells
daf_T_clean_IL10_eae <- filterSCE(daf_T_clean, IL_10_pos == "yes")

#remove naive cells
daf_T_clean_IL10_eae <- filterSCE(daf_T_clean_IL10_eae, condition == "EAE")

#pie chart of cluster contributions to IL-10+ T cells in WT and ERAP1-KO (eae) (wt first)
IL10_T_pie_eae <- table(
  genotype = daf_T_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10_eae, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hf <- ggplot(data=IL10_T_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                               "#CC0022", "#FF9999", "#CC00CC", "#636363"))
hf

#pie chart of cluster contributions to IL-10+ T cells in WT and ERAP1-KO (eae) (erap1 now)
IL10_T_pie_eae2 <- table(
  genotype = daf_T_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_T_clean_IL10_eae, "T_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO")

#make plot
hz <- ggplot(data=IL10_T_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#CCC000", "#00997F", "#00D4FF", "#3376FE", "#BB32FE", 
                               "#CC0022", "#FF9999", "#CC00CC", "#636363"))
hz

#plot abundances as boxplot
plotAbundances(daf_T_clean_IL10_eae, 
               k = "T_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##Now compare cluster T IL-10+ cluster frequency b/w WT and ERAP1-KO    
eif <- metadata(daf_T_clean_IL10_eae)$experiment_info 
(da_formula39 <- createFormula(eif,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res39 <- diffcyt(daf_T_clean_IL10_eae,                                            
                    formula = da_formula39, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "T_clusters", verbose = FALSE, )  

#examine output
cluster_stats39 <- rowData(da_res39$res) %>%
  as.data.frame(.)

table(rowData(da_res39$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------------##






##------------------Now subset microglia-----------------------------------------------------##

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#subset microglia cells
daf_micro <- filterSCE(daf_X, cluster_id %in% c("Microglia"), k = "merging1")

levels(daf_micro$condition)
levels(daf_micro$IL_10_pos)
levels(daf_micro$sample_id)

#Make MDS plot  
pbMDS(daf_micro, by = "sample_id", color_by = "condition", shape_by = "genotype", size_by = TRUE)

#re-cluster just the microglia     
daf_micro <- cluster(daf_micro, 
                 features = NULL, 
                 xdim = 10, 
                 ydim = 10, 
                 maxK = 6, 
                 seed = 7089, 
                 verbose = TRUE
) 

#plot cluster heatmap
plotClusterHeatmap(daf_micro,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta6", 
                   m = NULL,               
                   cluster_anno = TRUE, 
                   draw_freqs = TRUE
) 

# run UMAP                          
set.seed(790)                               
daf_micro <- runDR(daf_micro, 
               dr = "UMAP", 
               cells = 10000,
               features = NULL 
)

#UMAP microglia
plotDR(daf_micro, "UMAP", color_by = "meta6", facet_by = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#Merge and annotate microglia clusters
merging_table_micro <- "merged_clusters_micro.xlsx"
merging_table4 <- read_excel(merging_table_micro)
head(data.frame(merging_table4))  

# convert to factor with merged clusters in desired order                
merging_table4$new_cluster <- factor(merging_table4$new_cluster,         
                                     levels = c("Debris", "BAMs", "Microglia 1", "Microglia 2"))        

daf_micro <- mergeClusters(daf_micro, k = "meta6",                                  
                       table = merging_table4, 
                       id = "micro_clusters",
                       overwrite = TRUE
)   

#remove debris  and BAMs 
#make a new SCE object when you do this!!
daf_micro_clean <- filterSCE(daf_micro, cluster_id %in% c("Microglia 1", "Microglia 2"), k = "micro_clusters")

#save objects
saveRDS(daf_micro, file = "./daf_micro.rds")
saveRDS(daf_micro_clean, file = "./daf_micro_clean.rds")

#change cluster colors
clust_color6 <- c("#66E4FD", "#FF97A7")

#UMAP microglia
plotDR(daf_micro_clean, "UMAP", color_by = "micro_clusters", facet_by = "condition", k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)

#heatmap of merged micro cluster markers
plotExprHeatmap(daf_micro_clean, 
                bin_anno = FALSE, 
                scale = "first",
                k = "micro_clusters",
                by = "cluster_id", 
                k_pal = clust_color6,
                features = c("CD11b", "CD4", "Ly6C", "IL-10", "SiglecH", "CD11c", "B220", "CD38", "AF", "CD206", "Tim3", "Lyve1", "NKG2D", "CD90"),
                hm_pal = brewer.pal(9, "Greys")
)

#UMAP microglia (by important markers)
plotDR(daf_micro_clean, "UMAP", color_by = c("B220", "CD90", "NKG2D", "CD4", "CD38", "Ly6C", "SiglecH", "AF", "CD11c", "Tim3"), k_pal = clust_color6, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.5)


#------------------compare microglia cluster frequencies by genotype (naive first)-----------------------------#
#remove EAE cells    
daf_micro_clean_naive <- filterSCE(daf_micro_clean, condition == "naive")

#plot abundances as barplot
plotAbundances(daf_micro_clean_naive, k = "micro_clusters", by = "sample_id", k_pal = clust_color6)


#plot abundances as boxplot
plotAbundances(daf_micro_clean_naive, 
               k = "micro_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color6,
               shape = NULL
)


##To analyze now we can generate a GLMM
library(diffcyt)
eiw <- metadata(daf_micro_clean_naive)$experiment_info 
(da_formula30 <- createFormula(eiw,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

da_res30 <- diffcyt(daf_micro_clean_naive,                                            
                    formula = da_formula30, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "micro_clusters", verbose = FALSE)  

#examine output
cluster_stats30 <- rowData(da_res30$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res30$res)$p_adj < FDR_cutoff)
##--------------------------------------------------------------------------------------##


#------------------compare microglia cluster frequencies by genotype (EAE)-----------------------------#
#remove naive cells
daf_micro_clean_eae <- filterSCE(daf_micro_clean, condition == "EAE")

#plot abundances as barplot
plotAbundances(daf_micro_clean_eae, k = "micro_clusters", by = "sample_id", k_pal = clust_color6)


#plot abundances as boxplot
plotAbundances(daf_micro_clean_eae, 
               k = "micro_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               k_pal = clust_color6,
               shape = NULL
)


##To analyze now we can generate a GLMM
eiv <- metadata(daf_micro_clean_eae)$experiment_info 
(da_formula31 <- createFormula(eiv,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res31 <- diffcyt(daf_micro_clean_eae,                                            
                    formula = da_formula31, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "micro_clusters", verbose = FALSE)  

#examine output
cluster_stats31 <- rowData(da_res31$res) %>%
  as.data.frame(.)

table(rowData(da_res31$res)$p_adj < FDR_cutoff)
##-----------------------------------------------------------------------------------##

##----------------------see what fraction of microglia from each cluster/sample are IV+ (naive)-----------------------##
    
#barplot of contribution of IV+ cells to each sample
plotCounts(daf_micro_clean_naive, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IV_clust_micro <- table(
  IV_pos = daf_micro_clean_naive$IV_pos, 
  cluster = cluster_ids(daf_micro_clean_naive, "micro_clusters")) 
  
  plot_clust_table <- prop.table(IV_clust_micro, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  rp <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#00CC22", "#000000"))
  
  #make purrty
  rp + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+ (naive)----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_micro_clean_naive, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_micro_clean_naive$IL_10_pos, 
  cluster = cluster_ids(daf_micro_clean_naive, "micro_clusters")) 
  
  plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  ct <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#FF9864", "#000000"))
  
  #make purrty
  ct + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}


##----------------------see what fraction of microglia from each cluster/sample are IV+ (eae)-----------------------##

#barplot of contribution of IV+ cells to each sample
plotCounts(daf_micro_clean_eae, group_by = "sample_id", color_by = "IV_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#00CC22", "#000000"))

#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IV_clust_micro <- table(
  IV_pos = daf_micro_clean_eae$IV_pos, 
  cluster = cluster_ids(daf_micro_clean_eae, "micro_clusters")) 
  
  plot_clust_table <- prop.table(IV_clust_micro, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  rp <- ggplot(data=plot_clust_table, aes(x=cluster, y=Freq, fill=IV_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#00CC22", "#000000"))
  
  #make purrty
  rp + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}
##---------------------------------------------------------------------------------------------##

##---------------------Now lets see what fraction of each cluster are IL_10+ (eae)----------------------##
#barplot of contribution of IL-10+ cells to each sample
plotCounts(daf_micro_clean_eae, group_by = "sample_id", color_by = "IL_10_pos", prop = TRUE) + 
  scale_fill_manual(values = c("#FF9864", "#000000"))


#now by cluster   (must do this manually)
#manually extract df containing merging1 cluster ids matched to IV_pos
{IL10_clust <- table(
  IL10_pos = daf_micro_clean_eae$IL_10_pos, 
  cluster = cluster_ids(daf_micro_clean_eae, "micro_clusters")) 
  
  plot_clust_table2 <- prop.table(IL10_clust, "cluster") %>%
    as.data.frame(.)
  
  #make plot
  ct <- ggplot(data=plot_clust_table2, aes(x=cluster, y=Freq, fill=IL10_pos)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = c("#FF9864", "#000000"))
  
  #make purrty
  ct + theme(axis.ticks = element_line(colour = "black", 
                                       size = 1), axis.title = element_text(size = 15), 
             axis.text = element_text(colour = "black", 
                                      hjust = 1, angle = 0, size = 12), panel.background = element_rect(fill = NA),
             axis.text.x = element_text(angle = 45), axis.ticks.length = unit(7, "pt")
  )}

##--------------------------------------------------------------------------------------------##


##----------------------Compare IL-10+ T cells b/w genotypes by IV+ and cluster (eae)-----------------------##
# during EAE
#remove IL-10 neg cells
daf_micro_clean_IL10_eae <- filterSCE(daf_micro_clean_eae, IL_10_pos == "yes")

#pie chart of cluster contributions to IL-10+ microglia in WT and ERAP1-KO (eae) (wt first)
{IL10_micro_pie_eae <- table(
  genotype = daf_micro_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_micro_clean_IL10_eae, "micro_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "WT")

#make plot
hd <- ggplot(data=IL10_micro_pie_eae, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#66E4FD", "#FF97A7", "#9965FD"))
hd}

#pie chart of cluster contributions to IL-10+ T cells in WT and ERAP1-KO (eae) (erap1 now)
{IL10_micro_pie_eae2 <- table(
  genotype = daf_micro_clean_IL10_eae$genotype, 
  cluster = cluster_ids(daf_micro_clean_IL10_eae, "micro_clusters")) %>%
  as.data.frame(.) %>%
  dplyr::filter(., genotype == "ERAP1_KO")

#make plot
hz <- ggplot(data=IL10_micro_pie_eae2, aes(x="", y=Freq, fill=cluster)) +
  geom_bar(stat="identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("#66E4FD", "#FF97A7", "#9965FD"))
hz}

#plot abundances as boxplot
plotAbundances(daf_micro_clean_IL10_eae, 
               k = "micro_clusters", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL
)


##Now compare cluster microglia IL-10+ cluster frequency b/w WT and ERAP1-KO   
eig <- metadata(daf_micro_clean_IL10_eae)$experiment_info 
(da_formula34 <- createFormula(eig,                     
                               cols_fixed = "genotype",      
                               cols_random = "sample_id")) 

da_res34 <- diffcyt(daf_micro_clean_IL10_eae,                                            
                    formula = da_formula34, contrast = contrast,                    
                    analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                    clustering_to_use = "micro_clusters", verbose = FALSE, transform = FALSE)  

#examine output
cluster_stats34 <- rowData(da_res34$res) %>%
  as.data.frame(.)

table(rowData(da_res34$res)$p_adj < FDR_cutoff)

##--------------------------------------------------------------------------------------------##




##-------------------------Compare marker expression on IL-10+ cells b/w naive and EAE-------##
#only running this on WT, since it is more for general reference 

#load main dataset back
setwd("~/IL10_EAE_R_analys_unmix1")
daf_X <- readRDS("./daf_X.rds")

#remove IL-10 neg cells
daf_IL10_only <- filterSCE(daf_X, IL_10_pos == "yes")

#remove ERAP1-KO samples
daf_IL10_only <- filterSCE(daf_IL10_only, genotype == "WT")

##Now compare differential states b/w EAE and naive   ##this only tests state markers; need to force it to use all

#now trying to manually set all markers to be "state" markers on the raw SCE object as a work around
rowData(daf_IL10_only)$marker_class <- "state"

ei <- metadata(daf_IL10_only)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(daf_IL10_only), cols_design = "condition")

da_res99 <- diffcyt(daf_IL10_only,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    #markers_to_test = all_mark,  #works, except I can't specify which markers using this command
                    clustering_to_use = "merging1", verbose = FALSE)  


#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)

table(rowData(da_res99$res)$p_adj < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(daf_IL10_only, tbl_DS, fdr = 0.05, 
                sort_by = "lfc", col_anno = "condition")

#make specific for NK, CD8, CD4 T, B cells
k <- metadata(da_res99$res)$clustering_name
sub <- filterSCE(daf_IL10_only, cluster_id %in% c("NK cells", "CD8+ T cells", "CD4+ T cells", "B cells"), k = k)
plotDiffHeatmap(sub, tbl_DS, all = TRUE, normalize = FALSE)


##--------same as above, but removing IV+ cells  6/14/2021
#Filter out clusters that make no IL-10
for_analy <- filterSCE(daf_IL10_eae_wt, cluster_id %in% c("NK cells", "CD8+ T cells", "CD4+ T cells", "B cells", "Microglia", "DCs", "MdCs"), k = k)

#now trying to manually set all markers to be "state" markers on the raw SCE object as a work around
rowData(for_analy)$marker_class <- "state"

ei <- metadata(for_analy)$experiment_info   ##NOTE: must use "ei" here as package only recognizes this##
(da_formula99 <- createFormula(ei,                     
                               cols_fixed = "condition",      
                               cols_random = "sample_id")) 

#DS testing needs a design matrix     
design99 <- createDesignMatrix(ei(for_analy), cols_design = "condition")

da_res99 <- diffcyt(for_analy,                                            
                    formula = da_formula99, contrast = contrast,                    
                    analysis_type = "DS", method_DS = "diffcyt-DS-limma",  
                    design = design99,
                    #markers_to_test = all_mark,  #works, except I can't specify which markers using this command
                    clustering_to_use = "merging1", verbose = FALSE)  


#examine output
cluster_stats99 <- rowData(da_res99$res) %>%
  as.data.frame(.)

table(rowData(da_res99$res)$p_val < FDR_cutoff)

tbl_DS <- rowData(da_res99$res)

#plot heatmap of results
plotDiffHeatmap(for_analy, tbl_DS, fdr = 0.05, 
                sort_by = "lfc", col_anno = "condition")

plotDiffHeatmap(for_analy, tbl_DS, all = T, 
                sort_by = "padj", col_anno = "condition")

#make specific for NK, CD8, CD4 T, B cells
k <- metadata(da_res99$res)$clustering_name
sub <- filterSCE(daf_IL10_eae_wt, cluster_id %in% c("NK cells", "CD8+ T cells", "CD4+ T cells", "B cells", "Microglia", "BAMs"), k = k)
plotDiffHeatmap(sub, tbl_DS, all = TRUE, normalize = T)


#--------------------------------------------------------------------------------------------##


##--------------------Generating force-directed scaffold map-------------------------##

#this is same as from Mrjden et al. 2017 and applies the 100 FlowSOM clusters to manually gated landmark nodes. 
#source: https://github.com/ParkerICI/flow-analysis-tutorial

setwd("~/IL10_EAE_R_analys_unmix1_v2")

library(grappolo) #for clustering; dont use
library(vite)

##--now generate equivalent matrix from my data--##
#note: only outputs the markers that were used for original FlowSOM clustering. 
SOM_expr <- daf_X@metadata$SOM_codes %>%
  as.data.frame(.) %>%
  rownames_to_column(., var = "cellType") 

#write SOM_expr table out to tab separated txt file
write.table(SOM_expr, row.names = FALSE, col.names = TRUE, 'SOM_expr.txt', sep = '\t', quote = FALSE)

#read it back to confirm is has correct format.
readback <- read.delim("SOM_expr.txt", header = TRUE, stringsAsFactors = FALSE)
#move .txt file to a new folder

#now run Vite

# List the clustered.txt files contained in the "single_samples" directory
#note: I have a single example file here from all pooled samples so I gave it a single file in pattern.
input.files <- list.files(path = "~/IL10_EAE_R_analys_unmix1_v2/Vite_files", pattern = "SOM_expr.txt$", full.names = TRUE)

# Optional: load a table of sample-level metadata. All the nodes derived from the corresponding cluster file will
# have vertex properties corresponding to this metadata
#note: of no use if you only have a single file w/ all pooled; I am not using now. 
#metadata.tab <- read.table("metadata_tab.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define which columns contain variables that are going to be used to calculate similarities between the nodes
#Note: you need to only use the same markers that original FlowSOM clustering was done on. 
col.names <- c("CD45", "Ly6G", "IgD", "CD11c", "CD19", "CD8", "Ly6C", "CD4", "CD11b", "SiglecH", "B220", 
               "CD49b", "CD38", "MHCII", "NK1.1", "CCR2", "CD206", "Lyve1", "AF", "CD3")

landmarks.data <- load_landmarks_from_dir("~/IL10_EAE_R_analys_unmix1_v2/landmarks", asinh.cofactor = 6000, transform.data = T)

# Load the data for the landmarks (make sure FCS files are named properly)
landmarks.data2 <- landmarks.data$landmarks.data %>%
  select(., -"SampleID") %>%
  rename(., "AF" = "AF-A")  #remove sampleID col and fix AF syntax

landmarks.data3 <- landmarks.data$tab.landmarks %>%
  select(., -"SampleID") %>%
  rename(., "AF" = "AF-A")  #remove sampleID col and fix AF syntax

landmarks5 <- (list(landmarks.data = landmarks.data2, tab.landmarks = landmarks.data3))

# Run the analysis. By default results will be save in a directory called "scaffold_result"
scaff <- run_scaffold_analysis(input.files, ref.file = input.files[1], 
                               landmarks.data = landmarks5, col.names = col.names, process.clusters.data = FALSE)


sessionInfo()

