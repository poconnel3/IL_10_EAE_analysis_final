#By: Patrick O'Connell
#07012021

#This analysis is of clean, living CD19+ B cells from splenocytes and CNS of EAE mice
#mice are all on IL-10-GFP background


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
setwd("~/EAE_b_cell_phenotype1_analy")
metadata_filename <- "B_cell_metadata.xlsx"
md <- read_excel(metadata_filename)

# Define condition variables as named in metadata
md$location <- factor(md$location, levels = c("spleen", "CNS"))
md$genotype <- factor(md$genotype, levels = c("ERAP1_KO", "WT"))
md$IL10_pos <- factor(md$IL10_pos, levels = c("yes", "no"))
head(data.frame(md))
##-----------------------------------------------------------##

#Read FCS files in as flowset 
setwd("~/EAE_b_cell_phenotype1_analy/clean_Bcell_fcs")
fcs_raw <- read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

##------------------making panel file-----------------------##
setwd("~/EAE_b_cell_phenotype1_analy")
panel <- "B_cell_panel.xlsx"
panel <- read_excel(panel) 

#check
all(panel$fcs_colname %in% colnames(fcs_raw))
##------------------making panel file-----------------------##

#see how many cells per file
fsApply(fcs_raw, nrow)

#downsample to 90,000 cells per sample
# Define a downsampling ceiling
sampling.ceiling <- 90000
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
factors <- list(factors = c("file_name", "location", "sample_id", "genotype", "IL10_pos"))

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
table(daf_X_clean$location)
table(daf_X$sample_id)
table(daf_X$file_name)
table(daf_X$genotype)


#view metadata parameters
names(colData(daf_X_clean))

#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")    

#plot global marker expression by condition
p <- plotExprs(daf_X, color_by = "genotype")
p$facet$params$ncol <- 6                   
p 


#Make MDS plot  
pbMDS(daf_X, by = "sample_id", color_by = "genotype", size_by = FALSE)

#plot number of cells per condition
plotCounts(daf_X, group_by = "location", color_by = "genotype", prop = FALSE)
table(daf$condition) 

#Perform FlowSOM clustering
daf_X <- cluster(daf_X, 
               features = "type", 
               xdim = 10, 
               ydim = 10, 
               maxK = 18, 
               seed = 1234, 
               verbose = TRUE
               ) 

#plot cluster heatmap
plotClusterHeatmap(daf_X,                               
                   hm2 = "state", #swap for "state" or "abundances" to get more info. 
                   k = "meta18", 
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



#save main analysis object
saveRDS(daf_X, file = "./daf_X.rds")

#load main analysis object back in
daf_X <- readRDS(file = "./daf_X.rds")

#load main analysis object back in 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#UMAP visualizations
plotDR(daf_X, "UMAP", color_by = "meta18") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) 


#tSNE visualizations
plotDR(daf_X, "TSNE", color_by = "meta14") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) 


#make UMAP colored by important markers 
plotDR(daf_X_clean, "UMAP", color_by = c("CD38", "B220", "MHC-II", "IgD", "CD19", "CD21", "CD23", "CD138", "IL-10", "GL7", "IgM", "CD1d", "CD5", "CD43", "CD40", "CD80", "CD93"), facet_by = NULL, a_pal = rev(hcl.colors(10, "Oranges"))) +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) +
  geom_point(size=0.4)


## Facet per sample                                          
plotDR(daf_X, "UMAP", color_by = "meta14", facet = "genotype") +
  theme(axis.line = element_line(colour = NA),   #use this theme for all plots
        axis.ticks = element_line(colour = NA), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        axis.title = element_text(colour = NA), 
        axis.text = element_text(colour = NA), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA)) 



#Merge and annotate clusters
merging_tableX <- "merged_clusters_X.xlsx"
merging_table1 <- read_excel(merging_tableX)
head(data.frame(merging_table1))  

# convert to factor with merged clusters in desired order                
merging_table1$new_cluster <- factor(merging_table1$new_cluster,         
                                     levels = c("MZ B cells", "FO B cells", "T1 cells", "T2 cells", "T3 cells", "B1b cells", "B1a cells",             
                                                "Plasma cells", "Debris", "B220+ B1a cells", "B220+ B1b cells", "Bin cells", "Bregs"))        

daf_X <- mergeClusters(daf_X, k = "meta18",                                  
                     table = merging_table1, 
                     id = "merging1",
                     overwrite = TRUE
                     )   

#remove debris cluster
#make a new SCE object when you do this!!
daf_X_clean <- filterSCE(daf_X, cluster_id != "Debris", k = "merging1")

#change cluster colors
clust_color <- c("#996600", "#6665FD", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#ff9900", 
                 "#b3b3b3", "#0d0d0d", "#FF65CA", "#00cc66", "#cc9900")

clust_color.alt <- c("#996600", "#006699", "#FFCB64", "#8800CC", "#00994D", "#FF3276", "#ff9900", 
                 "#b3b3b3", "#0d0d0d", "#FF65CA", "#00cc66", "#cc9900")

#UMAP of merged clusters
plotDR(daf_X_clean, "UMAP", color_by = "merging1", k_pal = clust_color.alt, facet_by = "location") +
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

plotExprHeatmap(daf_X_clean, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color.alt, 
                hm_pal = brewer.pal(9, "Greys")
                #features = c("SLAMF7", "CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b","CD45", 
                           #  "B220", "CD19", "SiglecH", "CCR2", "NK1.1", "CD11c", "Ly6G", "MHCII", "Lyve1", "CD206", "IgD")
)

#heatmap only CNS
daf_X_clean_CNS <- filterSCE(daf_X_clean, location == "CNS")
plotExprHeatmap(daf_X_clean_CNS, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color.alt, 
                hm_pal = brewer.pal(9, "Greys")
                #features = c("SLAMF7", "CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b","CD45", 
                #  "B220", "CD19", "SiglecH", "CCR2", "NK1.1", "CD11c", "Ly6G", "MHCII", "Lyve1", "CD206", "IgD")
)

#now only spleen
#heatmap only CNS
daf_X_clean_spleen <- filterSCE(daf_X_clean, location == "spleen")
plotExprHeatmap(daf_X_clean_spleen, 
                bin_anno = FALSE, 
                scale = "last",
                k = "merging1",
                by = "cluster_id", 
                k_pal = clust_color, 
                hm_pal = brewer.pal(9, "Greys")
                #features = c("SLAMF7", "CD38", "CD8", "NKG2D", "IL-10", "CD90", "CD11b", "Ly6C", "CD4", "AF", "LAG3", "Tim3", "CD49b","CD45", 
                #  "B220", "CD19", "SiglecH", "CCR2", "NK1.1", "CD11c", "Ly6G", "MHCII", "Lyve1", "CD206", "IgD")
)
#ridgeplots of markers on merged clusters
plotClusterExprs(daf_alt_clean_rerun, k = "merging1")

#save main analysis object
saveRDS(daf_X_clean, file = "./daf_X_clean.rds")

#-------------------------------compare groups now------------------------------##

#plot abundances as barplot
plotAbundances(daf_X_clean, k = "merging1", group_by = "genotype", k_pal = clust_color)


#plot abundances as boxplot
plotAbundances(daf_X_clean, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
               ) +
  ylim(0, NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_X_clean)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))

#how to access merged cluster levels
levels(cluster_ids(daf_X_clean, k = "merging1"))

#cluster levels of FlowSOM clusters
levels(cluster_ids(daf_X))

da_res1 <- diffcyt(daf_X_clean,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
cluster_stats1 <- rowData(da_res1$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res1$res)$p_adj < FDR_cutoff)



#now compare just spleen 6/15/2021
daf_spleen <- filterSCE(daf_X_clean, location == "spleen")

#plot abundances as boxplot
plotAbundances(daf_spleen, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
) + scale_color_manual(values = c("#0099cc", "#ffbf00")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)

#now compare just CNS 6/15/2021
daf_cns <- filterSCE(daf_X_clean, location == "CNS")

#plot abundances as boxplot
plotAbundances(daf_cns, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
) + scale_color_manual(values = c("#0099cc", "#ffbf00")) +
  scale_fill_manual(values = c("#ffffff", "#ffffff")) +
  ylim(0,NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_spleen)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))



da_res2 <- diffcyt(daf_spleen,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
cluster_stats2 <- rowData(da_res2$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res2$res)$p_adj < FDR_cutoff)



#now compare just CNS
daf_cns <- filterSCE(daf_X_clean, location == "CNS")

#plot abundances as boxplot
plotAbundances(daf_cns, 
               k = "merging1", 
               by = "cluster_id",
               group_by = "genotype",
               shape = NULL,
) +
  ylim(0, NA)


##To analyze now we can generate a GLMM
ei <- metadata(daf_cns)$experiment_info 
(da_formula1 <- createFormula(ei,                     
                              cols_fixed = "genotype",      
                              cols_random = "sample_id")) 

contrast <- createContrast(c(0, 1))



da_res3 <- diffcyt(daf_cns,                                            
                   formula = da_formula1, contrast = contrast,                    
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",           
                   clustering_to_use = "merging1", verbose = FALSE)  

#examine output
cluster_stats3 <- rowData(da_res3$res) %>%
  as.data.frame(.)

FDR_cutoff <- 0.05

table(rowData(da_res3$res)$p_adj < FDR_cutoff)

##----------------------IGNORE ALL BELOW THIS LINE; NOT RUN--------------------------##
#now compare WT only across locations

##Note: I have to change the sample_id names for all of these so that they have unique names. 
#ie. the paired CNS and spleen samples need separate names even though they are from the same mouse. 

#load main analysis object back in 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#edit the metadata
#newest attempt 2/24/21

# unfactor sample IDs, else overwriting won't work
daf_X_clean$sample_id <- as.character(daf_X_clean$sample_id)

#for WT
i <- grep("WT_[0-9]+$", daf_X_clean$sample_id)
cd <- data.frame(colData(daf_X_clean))
new <- with(cd, paste(sample_id, location, sep = "_"))
new2 <- as.data.frame(new)
daf_X_clean$sample_id[i] <- new2$new[i]

#for ERAP1-KO
i <- grep("ERAP1_KO_[0-9]+$", daf_X_clean$sample_id)
cd <- data.frame(colData(daf_X_clean))
new <- with(cd, paste(sample_id, location, sep = "_"))
new3 <- as.data.frame(new)
#unfactor
daf_X_clean@metadata$experiment_info <- as.character(daf_X_clean@metadata$experiment_info)
metadata(daf_X_clean)$experiment_info  $sample_id <- new3$new


#check
table(daf_X_clean$sample_id)

#save
saveRDS(daf_X_clean, file = "./daf_X_clean2.rds")  


##------------compare WT markers between spleen and CNS
#6/15/2021

daf_X_WT2 <- filterSCE(daf_X_clean2, genotype == "WT")
rowData(daf_X_WT2)$marker_class <- "state"
##Statistically compare markers on clusters w/ GLMM. (CNS)
ei <- metadata(daf_X_WT2)$experiment_info

ds_formula2 <- createFormula(ei, cols_fixed = "location")  #NOTE: cannot have a rondom column variable here or else it fails for some reason. 

contrast <- createContrast(c(0, 1))

FDR_cutoff <- 0.05

ds_res7 <- diffcyt(daf_X_WT2, 
                   formula = ds_formula2, contrast = contrast,
                   analysis_type = "DS", method_DS = "diffcyt-DS-LMM",
                   clustering_to_use = "merging1", verbose = FALSE)

table(rowData(ds_res7$res)$p_adj < FDR_cutoff)

cluster_stats7 <- rowData(ds_res7$res) %>%
  as.data.frame(.)
#export stat data
write.xlsx(cluster_stats7, "~/EAE_b_cell_phenotype1_analy/cluster_stats7.xlsx")

##tbl_DS <- rowData(da_res99$res)
plotDiffHeatmap(daf_cns, rowData(ds_res7$res), top_n = 50, fdr = FDR_cutoff, lfc = 2, col_anno = c("genotype"), row_anno = TRUE) 

##----------------------------------------------------------------##

#generate correlation matrix of frequencies of each B cell subset comparing computational approach to manual gating. 
#do this separably for CNS and spleen. 
library(xlsx)
library(corrplot)
#load main analysis object back in 
daf_X_clean <- readRDS(file = "./daf_X_clean.rds")

#get number of cells in each cluster (spleen)
daf_spleen <- filterSCE(daf_X_clean, location == "spleen") 
spleen_clust_freq <- table(cluster_ids(daf_spleen, k = "merging1"), daf_spleen$sample_id)
write.xlsx(spleen_clust_freq, "~/EAE_b_cell_phenotype1_analy/spleen_clust_freq.xlsx")

#add the above values (for conserved clusters across both approaches) to excel sheet and get 
#frequency of total for HD analysis and add in flowJo values

raw_corr <- read_excel("spleen_clust_freq.xlsx")
head(raw_corr)

# First calculate corrleation coefficients to be visualized
cor_matrix_1 = cor(raw_corr[ ,c(2:19)], method='spearman',use='pairwise.complete.obs')
# Now produce the plot
all_cor_matrix1 <- corrplot(cor_matrix_1, method='color', type='upper', order="FPC",  
                            tl.col = "black", tl.srt = 90)
all_cor_matrix1

