#By: Patrick O'Connell
#6/16/2021

## Analysis of scRNA-seq from PBMCs of healthy indivduals across many data sets w/ ERAP1 SNPs called


#Goal here is to see what B cell look like comparing indivduals w/ various ERAP1 SNPs (mainly K528R)
#Will also compare this against the murine B cell RNA-seq data
#Using SCTransform normalization approach


#packages
library(tidyverse)
library(ggplot2)
library(Seurat)
library(readxl)
library(grr)
library(harmony)

#set working directory
setwd("~/Desktop/PAT/scRNA_SNP_project")
#on my macbook
setwd("~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/scRNA_SNP_project")

#First load all data sets, combined, and add metadata. Then split by indivdual and perform integration. 

#read in metadata
md <- read_xlsx("ERAP1_snp_meta.xlsx", col_names = TRUE) 
tenx_md <- md %>% dplyr::filter(., Method == "TenX")
seq_well_md <- md %>% dplyr::filter(., Method == "seq_well")

#load data
#read all in as a loop
filez <- as.character(tenx_md$sample_ID)
for (file in c(filez)){
  seurat_data <- Read10X(data.dir = paste0("count_tables/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 100, 
                                   project = file)
  assign(file, seurat_obj)
}

#not all are loading in from above chekc to make sure they match
list.files("./count_tables") %in% tenx_md$sample_ID  #they are
#corrected file names again and now they all load. 


##-----merge seurat objects--------
new_all_but <- c(SRR10211568, SRR10211569, SRR10211570, SRR10211572, SRR10338074, SRR10338075, SRR10338076, SRR10338077, SRR10338078,
                 SRR10338079, SRR10080324, SRR10080325, SRR10080326, SRR10080327, SRR10080328, SRR10080329, SRR10080330, SRR10080331,
                 SRR10080332, SRR10080333, SRR10080334, SRR10080335, SRR3561757,  SRR3561756,  SRR3561758,  HC1, HC2, HC3, HC4,  HC5, 
                 HRR059128, HRR059129, HRR059130, HRR059131,  HRR059132)
new_all_seurat <- c(SRR10211571,SRR10211568, SRR10211569, SRR10211570, SRR10211572, SRR10338074, SRR10338075, SRR10338076, SRR10338077, SRR10338078,
                    SRR10338079, SRR10080324, SRR10080325, SRR10080326, SRR10080327, SRR10080328, SRR10080329, SRR10080330, SRR10080331,
                    SRR10080332, SRR10080333, SRR10080334, SRR10080335, SRR3561757,  SRR3561756,  SRR3561758,  HC1, HC2, HC3, HC4,  HC5, 
                    HRR059128, HRR059129, HRR059130, HRR059131,  HRR059132)

temp_merge <- merge(SRR10211571, y = new_all_but, add.cell.ids = filez, project = "ERAP1_SNP")

##------------add metadata-----------
tenx_md1 <- tenx_md %>% mutate(orig.ident = tenx_md$sample_ID)
temp_merge_metadata <- temp_merge@meta.data
temp_merge_metadata <- temp_merge_metadata[,!(colnames(temp_merge_metadata)) %in% 
                                     setdiff(colnames(tenx_md1), "orig.ident")]
tenx_md_combined <- merge(temp_merge_metadata, tenx_md1, by = "orig.ident")
rownames(tenx_md_combined) <- rownames(temp_merge@meta.data)
temp_merge@meta.data <- tenx_md_combined

#temp_merge <- readRDS("./temp_merge.rds")                                                                
saveRDS(temp_merge, "./temp_merge.rds")



##-----read in seq-well data------
path = "~/Documents/Research Projects/Tr1 project:ERAP1_MS_project/scRNA_SNP_project/count_tables_seq_well/"
cm.list = paste0(path, list.files(pattern = "*.matrices.rds", path = path))
cm.files <- lapply(cm.list, readRDS)
names(cm.files) <- sub(path,"",
                       sub("\\_cell.counts.matrices.rds", "", cm.list))

#library(devtools)
library(EpicTools)
library(grr)
cm.pp <- mapply(EpicPreHS, cm.files, orig.ident = names(cm.files), SIMPLIFY = F)
covid_combined.emat <- mergeCM(cm.pp, type = "emat")
seq_well_combined <- CreateSeuratObject(counts = covid_combined.emat, min.cells = 10, names.field = 1, names.delim = "\\.")

#add metadata to seq_well samples
covid_metadata <- read_csv("https://raw.githubusercontent.com/ajwilk/2020_Wilk_COVID/master/code/COVID-19_metadata_repo.csv")
#add ERAP1 SNPs to this
seq_well_md1 <- seq_well_md %>% mutate(orig.ident = c("SRR11804725_HIP002", "SRR11804726_HIP015", "SRR11804727_HIP023", "SRR11804728_HIP043", "SRR11804729_HIP044", "SRR11804730_HIP045"))

seurat_metadata <- seq_well_combined@meta.data
seurat_metadata <- seurat_metadata[,!(colnames(seurat_metadata)) %in% 
                                     setdiff(colnames(seq_well_md1), "orig.ident")]
metadata_combined <- merge(seurat_metadata, seq_well_md1, by = "orig.ident")
rownames(metadata_combined) <- rownames(seq_well_combined@meta.data)
seq_well_combined@meta.data <- metadata_combined

saveRDS(seq_well_combined, "./seq_well_combined.rds")
seq_well_combined <- readRDS("./seq_well_combined.rds")


#---merge 10X and Seq_well data-----
compl.seurat <- merge(temp_merge, y = seq_well_combined, add.cell.ids = c("tenX", "seq_well"), project = "ERAP1_SNP")
saveRDS(compl.seurat, "./compl.seurat.rds")
#its too large
rm(list=ls(pattern= c("HRR")))

##--check cell quality and remove low quality----
compl.seurat <- SetIdent(compl.seurat, value = "study")
compl.seurat[["percent.mt"]] <- PercentageFeatureSet(compl.seurat, pattern = "^MT-")
VlnPlot(compl.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0) 
plot1 <- FeatureScatter(compl.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(compl.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#subset out good cells
clean.seurat <- subset(compl.seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4100 & percent.mt < 20)

#check how many cells were filterd
clean.out <- table(Idents(clean.seurat))/table(Idents(compl.seurat))
ggplot(data=clean.out, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity") + theme_minimal()

##---calc cell cycle score----
clean.seurat <- CellCycleScoring(
  object = clean.seurat,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)

VlnPlot(clean.seurat, features = c("S.Score","G2M.Score"), pt.size = 0)
saveRDS(clean.seurat, "./clean.seurat.rds")
clean.seurat <- readRDS("./clean.seurat.rds")

##-------split seurat objects back out again by study and perform intrgration----
indv.list <- SplitObject(clean.seurat, split.by = "sample_ID")
library(glmGamPoi)
indv.list <- lapply(X = indv.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = indv.list, nfeatures = 3000)
indv.list <- PrepSCTIntegration(object.list = indv.list, anchor.features = features)
indv.list <- lapply(X = indv.list, FUN = RunPCA, features = features)
saveRDS(indv.list, "./indv.list.rds")

##--------perform integration-----------
immune.anchors <- FindIntegrationAnchors(object.list = indv.list, normalization.method = "SCT",
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
saveRDS(immune.anchors, "./immune.anchors.rds")
immune.anchors <- readRDS("./immune.anchors.rds")
int.seurat <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30)

##----begin standard analysis--------
int.seurat <- RunPCA(int.seurat, verbose = FALSE)
int.seurat <- RunUMAP(int.seurat, reduction = "pca", dims = 1:30)



##-------test out Harmony for dataset integration------

##--------------script to run on HPCC--------------
library(Seurat)
library(harmony)
clean.seurat <- readRDS("./clean.seurat.rds")
clean.seurat <- SCTransform(clean.seurat, vars.to.regress = "percent.mt") %>%
  RunPCA(npcs = 30, verbose = FALSE)
int.seurat <- clean.seurat %>% 
  RunHarmony(c("sample_ID", "study"), plot_convergence = FALSE, verbose = TRUE, assay.use = "SCT") #may have to try chainging this to "RNA" if results are not nice. 
saveRDS(int.seurat, "./int.seurat.rds")
##----------------------------------------



#investigating feature differences b/w datasets
clean.seurat <- SetIdent(clean.seurat, value = "Method")
VlnPlot(clean.seurat, features = "nFeature_RNA", pt.size = 0)

feat_calc <- FetchData(object = clean.seurat, vars = c("nFeature_RNA", "Method")) %>%
  group_by(Method) %>% 
  mutate(avg = median(nFeature_RNA))   #They have a similar # of features per cell across tech

#Used below to compare the gene names and they are the same. 
#both technologies. 
well_feat <- subset(x = clean.seurat, subset = Method == "seq_well") %>% rownames(.)
tex_feat <- subset(x = clean.seurat, subset = Method == "TenX") %>% rownames(.)

#can use this to see how many are different b/w technologies
length(setdiff(well_feat, tex_feat)) #=0, they have the same genes


#------picking up from B cell subsetting done on lab Mac-------
DimPlot(int.seurat_norm, reduction = "umap", label = FALSE, raster = FALSE, cols = c("#bfbfbf", "#7094db"))


#-----now subset out B cells---------
B_cell_sub <- subset(int.seurat_norm, idents = "B cells") 
B_cell_sub
saveRDS(B_cell_sub, "./B_cell_sub.rds")
#------------------------------------

#--------------re-SCT, UMAP, and re-cluster---------------------
#Will first see if cells are still grouping by study, and if so, will add a harmony integration step in here. 
#NOTE: they are!
B_cell_sub <- readRDS("./B_cell_sub.rds")

#remove Wilk et al
B_cell_sub <- subset(x = B_cell_sub, subset = Method == "TenX")
saveRDS(B_cell_sub, "./B_cell_sub.rds")

#RUNNING SCT AND HARMONY ON HPCC
int.Bcell.1 <- readRDS("./int.Bcell.rds")

int.Bcell.1 <- RunUMAP(int.Bcell.1, reduction = "harmony", dims = 1:30, seed.use = 1991) %>%
  FindNeighbors(dims = 1:30, verbose = FALSE)


#RUNNING SCT SEURAT INTEGRATION ON HPCC
#----------script to run (split by study)-----------
Bcell.list <- SplitObject(B_cell_sub, split.by = "study")
Bcell.list <- lapply(X = Bcell.list, FUN = SCTransform)
features <- SelectIntegrationFeatures(object.list = Bcell.list, nfeatures = 3000)
Bcell.list <- PrepSCTIntegration(object.list = Bcell.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = Bcell.list, normalization.method = "SCT",
                                         anchor.features = features)
Bcell_comb.sct_1 <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
#----------------------------------


#---check how the seurat integration worked
#This one integrated by study
Bcell_comb.sct_1 <- readRDS("./Bcell_comb.sct_1.rds")
Bcell_comb.sct_1 <- RunPCA(Bcell_comb.sct_1, verbose = FALSE)
Bcell_comb.sct_1 <- RunUMAP(Bcell_comb.sct_1, reduction = "pca", dims = 1:30)
DefaultAssay(Bcell_comb.sct_1) <- "integrated"
Bcell_comb.sct_1 <- FindNeighbors(Bcell_comb.sct_1, dims = 1:30, verbose = FALSE)

Bcell_comb.sct_1 <- FindClusters(Bcell_comb.sct_1, verbose = FALSE, resolution = 0.8, algorithm = 2, random.seed = 1991) #play w/ resolution until I get nice balance of clusters


DimPlot(Bcell_comb.sct_1, label = TRUE, raster=FALSE) 
DimPlot(Bcell_comb.sct_1, label = FALSE, raster=FALSE, group.by = "study") 
#Looks good! Will save and continue w/ this one. Will check the seurat by indv integration first just to make sure that is not better
saveRDS(Bcell_comb.sct_1, "./Bcell_comb.sct_1.rds")


#Find DE markers
DefaultAssay(Bcell_comb.sct_1) <- "RNA"
B.sub.markers <- FindAllMarkers(Bcell_comb.sct_1, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, test.use = "MAST", verbose = TRUE)
saveRDS(B.sub.markers, "./B.sub.markers.rds")
B.sub.markers <- readRDS("./B.sub.markers.rds")


##-----annotate new clusters---------

#save orig cell cluster IDs
Bcell_comb.sct_1[["old.ident"]] <- Idents(object = Bcell_comb.sct_1)

new_B_ids <- read_xlsx("./new_B_ids.xlsx") 

new.cluster.ids <- new_B_ids$new_Ids
names(new.cluster.ids) <- levels(Bcell_comb.sct_1)
Bcell_comb.sct_1 <- RenameIdents(Bcell_comb.sct_1, new.cluster.ids)


#-----remove non-B cells and plasma cells, re-cluster/UMAP and use SingleR to predict B cell subsets-----
upd_B <- subset(Bcell_comb.sct_1, idents = c("Plasma cells", "Put.B cells"))
DimPlot(upd_B, label = TRUE, raster=FALSE) 

DefaultAssay(upd_B) <- "integrated"
upd_B <- RunPCA(upd_B, verbose = FALSE)
upd_B <- RunUMAP(upd_B, reduction = "pca", dims = 1:30)
upd_B <- FindNeighbors(upd_B, dims = 1:30, verbose = FALSE)

upd_B <- FindClusters(upd_B, verbose = TRUE, resolution = 0.5, algorithm = 2, random.seed = 1991) #play w/ resolution until I get nice balance of clusters

DimPlot(upd_B, label = TRUE, raster=FALSE) 
DimPlot(upd_B, label = FALSE, raster=FALSE, group.by = "study") 

DefaultAssay(upd_B) <- "RNA"
upd_B.markers <- FindAllMarkers(upd_B, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, test.use = "MAST", verbose = TRUE)
saveRDS(upd_B.markers, "./upd_B.markers.rds")
saveRDS(upd_B, "./upd_B.rds")
upd_B <- readRDS("./upd_B.rds")

##-----------------get cell annotations w/ SingleR-------------
library(SingleR)
library(celldex)
#load ref datasets
hpca.se <- HumanPrimaryCellAtlasData()
blueprint.se <- BlueprintEncodeData()
monaco.se <- MonacoImmuneData()

#get approp. counts from my data
DefaultAssay(upd_B) <- "SCT"
counts <- GetAssayData(upd_B)
clust_lab <- as.character(Idents(upd_B))

#run SingleR
clust_based_class1 <- SingleR(
  test = counts,
  ref = list(hpca=hpca.se, blue=blueprint.se, monac=monaco.se),
  labels = list(hpca.se$label.fine, blueprint.se$label.fine, monaco.se$label.fine),
  clusters = clust_lab,
  genes = "de",
  sd.thresh = 1,
  de.method = "classic",
  recompute = TRUE,
  prune = TRUE,
  assay.type.test = "logcounts",
  check.missing = TRUE
)
saveRDS(clust_based_class1, "./clust_based_class1.rds")
#check results
table(clust_based_class1$labels)
plotScoreHeatmap(clust_based_class1)
clust_based_class1 <- readRDS("./clust_based_class1.rds")

#----------------------Add SingleR results back to seurat object
upd_B$cell_type_by_cluster <- factor(
  upd_B$seurat_clusters,
  levels = rownames(clust_based_class1),
  labels = clust_based_class1$labels)
#rename some overlaping identies
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
upd_B <- RenameIdents(object = upd_B, `naive B-cells` = "Naive B cells")
upd_B <- RenameIdents(object = upd_B, `Memory B-cells` = "Non-switched memory B cells")
upd_B <- RenameIdents(object = upd_B, `Plasmablasts` = "Plasmacells/blasts")
upd_B[["cell_type_by_cluster"]] <- Idents(object = upd_B)
levels(upd_B@meta.data$cell_type_by_cluster)
saveRDS(upd_B, "./upd_B.rds")
DimPlot(upd_B, label = F, raster=FALSE, cols = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53")) 
DimPlot(upd_B, label = F, raster=FALSE, group.by = "study", cols = c("#34568B", "#CD212A", "#FFA500", "#56C6A9", "#798EA4", "#FA7A35")) 

#------get cluster specific markers-------
DefaultAssay(upd_B) <- "RNA"
upd_B.markers2 <- FindAllMarkers(upd_B, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, test.use = "MAST", verbose = TRUE)
saveRDS(upd_B.markers2, "./upd_B.markers2.rds")
write_csv(upd_B.markers2, "./upd_B.markers2.csv")

#-------plot specific markers per cluster------

#heatmap
library(RColorBrewer)
mapal <- colorRampPalette(RColorBrewer::brewer.pal(9,"Greys"))(90)
DefaultAssay(upd_B) <- "SCT"
DoHeatmap(upd_B, slot = "data", features = top15$gene, group.colors = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53"),
          size = 1, group.bar.height = 0.05) + scale_fill_gradientn(colours = mapal)
#alt_color
DoHeatmap(upd_B, slot = "scale.data", features = top15$gene, group.colors = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53"),
          size = 1, group.bar.height = 0.05) #+ scale_fill_gradientn(colours = mapal)

#violin plot
VlnPlot(upd_B, features = c("XBP1", "MZB1", "IGHG1", "JCHAIN", "SSR4", "GPR183", "CLECL1", "RGS2", "CRIP1", "TCL1A", "IL4R", "IGHD", "FCER2", "YBX3", "LINC01781", "SCIMP", "ITGB1", "IGHA2", "S100A10", "MT2A", "GABARAP", "PTPRCAP", "MT-ATP8", "MT1X", "MS4A1", "CD79A"), 
        pt.size = 0, cols = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53"), assay = "RNA", slot = "data", 
        log = F, stack = T, flip = T, split.by = "cell_type_by_cluster")

#----feature plot of ERAP1 expression---
DefaultAssay(upd_B) <- "RNA"
FeaturePlot(upd_B, features = "ERAP1", slot = "data", pt.size = 1, cols = c("gray","red"))

#------------plot clust freq by study   
clust_study <- as.data.frame(prop.table(table(Idents(upd_B), upd_B$study), margin = 1))
gg <- ggplot(data=clust_study, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("#34568B", "#CD212A", "#FFA500", "#56C6A9", "#798EA4", "#FA7A35")) +theme_bw()
gg

#plot clust freq by indivdual   
clust_study <- as.data.frame(prop.table(table(Idents(upd_B), upd_B$sample_ID), margin = 1))
gr <- ggplot(data=clust_study, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_boxplot(stat="identity") + scale_fill_manual(values = c("#34568B", "#CD212A", "#FFA500", "#56C6A9", "#798EA4", "#FA7A35")) +theme_bw()
gr

#UMAP of indivduals on B cells
DimPlot(upd_B, label = F, raster=FALSE, group.by = "sample_ID") 


#plot indivdual contributions to each cluster
library(dittoSeq)
dittoBarPlot(upd_B, var = "cell_type_by_cluster", group.by = "sample_ID")
dittoBarPlot(upd_B, var = "study", group.by = "cell_type_by_cluster")

#plot clust freq by ERAP1 SNPs (boxplot)
library(plyr)
#pos
table(upd_B@meta.data$cell_type_by_cluster, upd_B@meta.data$sample_ID, upd_B@meta.data$K528R_SNP) %>% 
  as.data.frame() %>% dplyr::rename(cluster = Var1, sample_ID = Var2, K528R_pos = Var3) %>% filter(K528R_pos == "yes") %>%
  ddply(.(sample_ID),transform,prop=Freq/sum(Freq)) -> pos_props
pos_props <- pos_props %>% filter(prop != "NaN")

#neg
table(upd_B@meta.data$cell_type_by_cluster, upd_B@meta.data$sample_ID, upd_B@meta.data$K528R_SNP) %>% 
  as.data.frame() %>% dplyr::rename(cluster = Var1, sample_ID = Var2, K528R_pos = Var3) %>% filter(K528R_pos == "no") %>%
  ddply(.(sample_ID),transform,prop=Freq/sum(Freq)) -> neg_props
neg_props <- neg_props %>% filter(prop != "NaN")

comb_props <- bind_rows(neg_props, pos_props)
#plot it
ggplot(comb_props, aes(x=cluster, y=prop, fill=K528R_pos)) + stat_boxplot(geom ='errorbar') +
  geom_boxplot() + 
  scale_fill_manual(values = c("#EDD59E", "#6B5876")) +
  theme(axis.line = element_line(colour = "black"),   
        axis.ticks = element_line(colour = "black", size = 0.9), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA), 
        axis.title = element_text(size = 18,  vjust = 0), axis.text = element_text(size = 12, colour = "black")) +
  labs(y = "Proportion") +
  ylim(0,1)



#plot clust freq by ERAP1 SNPs (boxplot) [more granular clusters]
DimPlot(upd_B, label = T, raster=FALSE, group.by = "old.ident") 

#re-auto-cluster since I think I forgot to save
DefaultAssay(upd_B) <- "integrated"
upd_B <- FindNeighbors(upd_B, dims = 1:30, verbose = FALSE)
upd_B <- FindClusters(upd_B, verbose = TRUE, resolution = 0.5, algorithm = 2, random.seed = 1991) #play w/ resolution until I get nice balance of clusters
DimPlot(upd_B, label = T, raster=FALSE) 
#save
upd_B[["old.ident"]] <- Idents(object = upd_B)

DefaultAssay(upd_B) <- "RNA"
upd_B.markers.orig.clust <- FindAllMarkers(upd_B, only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, test.use = "MAST", verbose = TRUE)


#pos
table(upd_B@meta.data$old.ident, upd_B@meta.data$sample_ID, upd_B@meta.data$K528R_SNP) %>% 
  as.data.frame() %>% dplyr::rename(cluster = Var1, sample_ID = Var2, K528R_pos = Var3) %>% filter(K528R_pos == "yes") %>%
  ddply(.(sample_ID),transform,prop=Freq/sum(Freq)) -> pos_props2
pos_props2 <- pos_props2 %>% filter(prop != "NaN")

#neg
table(upd_B@meta.data$old.ident, upd_B@meta.data$sample_ID, upd_B@meta.data$K528R_SNP) %>% 
  as.data.frame() %>% dplyr::rename(cluster = Var1, sample_ID = Var2, K528R_pos = Var3) %>% filter(K528R_pos == "no") %>%
  ddply(.(sample_ID),transform,prop=Freq/sum(Freq)) -> neg_props2
neg_props2 <- neg_props2 %>% filter(prop != "NaN")

comb_props2 <- bind_rows(neg_props2, pos_props2)
write_csv(comb_props2, "./comb_props2.csv")
#plot it
ggplot(comb_props2, aes(x=cluster, y=prop, fill=K528R_pos)) + stat_boxplot(geom ='errorbar') +
  geom_boxplot() + 
  scale_fill_manual(values = c("#EDD59E", "#6B5876")) +
  theme(axis.line = element_line(colour = "black"),   
        axis.ticks = element_line(colour = "black", size = 0.9), 
        panel.grid.major = element_line(linetype = "blank"), 
        panel.grid.minor = element_line(linetype = "blank"), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 14), 
        panel.background = element_rect(fill = NA), 
        axis.title = element_text(size = 18,  vjust = 0), axis.text = element_text(size = 12, colour = "black")) +
  labs(y = "Proportion") +
  ylim(0,1)

#run 2-way ANOVA 
#export as csv and do it in Prism
write_csv(comb_props2, "./comb_props.csv") 

#plot ERAP1 expression by K528R SNP
DotPlot(upd_B, features = "ERAP1", assay = "SCT", split.by = "K528R_SNP") 
VlnPlot(upd_B, features = "ERAP1", split.by = "K528R_SNP", assay = "RNA", slot = "data")

cluster.averages <- AverageExpression(upd_B, return.seurat = TRUE, group.by = c("cell_type_by_cluster", "K528R_SNP"))
DoHeatmap(cluster.averages, features = "ERAP1", size = 3, draw.lines = T, group.by = "orig.ident", assay = "RNA", slot = "data", group.colors = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53")) #+ scale_fill_gradientn(colours = mapal)

#get expression of ERAP1 by cluster and SNP
DefaultAssay(cluster.averages) <- "RNA"
FetchData(cluster.averages, vars = "ERAP1")


#get marker differences b/w clust 0 and 2 and plot them
DefaultAssay(upd_B) <- "RNA"
upd_B <- SetIdent(upd_B, value = "old.ident")
zero_2_dif_markers <- FindMarkers(upd_B, ident.1 = "0", ident.2 = "2", only.pos = F, test.use = "MAST", verbose = TRUE, slot = "data")
saveRDS(zero_2_dif_markers, "./zero_2_dif_markers.rds")

#now Wilcoxin markers (same as MAST)
zero_2_dif_markers.W <- FindMarkers(upd_B, ident.1 = "0", ident.2 = "2", only.pos = F, verbose = TRUE, slot = "data")
saveRDS(zero_2_dif_markers.W, "./zero_2_dif_markers.W.rds")

#now DESeq2
zero_2_dif_markers.DS <- FindMarkers(upd_B, ident.1 = "0", ident.2 = "2", only.pos = F, verbose = TRUE, test.use = "DESeq2")
saveRDS(zero_2_dif_markers.DS, "./zero_2_dif_markers.DS.rds")

#subset 0 and 2
Zero2_sub <- subset(upd_B, idents = c("0", "2"))
zero_2_dif_markers %>%
  top_n(n = 30, wt = avg_log2FC) -> top30
zero_2_dif_markers %>%
  top_n(n = 30, wt = -avg_log2FC) -> bottom30
up_down <- rbind(top30, bottom30)

DoHeatmap(Zero2_sub, features = rownames(up_down), assay = "SCT", slot = "scale.data")

#this one is ok, just change colors and scale dot size and gene names
fg <- DotPlot(Zero2_sub, features = rownames(up_down), assay = "SCT", cols = c("lightgrey", "#e68a00")) + RotatedAxis()
fg + theme(axis.text.x = element_text(size =7))

DoHeatmap(upd_B, features = "ERAP1", size = 3, draw.lines = T, group.by = "orig.ident", assay = "RNA", slot = "data", group.colors = c("#926AA6", "#00A170", "#E9897E", "#0072B5", "#FDAC53")) #+ scale_fill_gradientn(colours = mapal)

#volcano plot of DE genes
library(EnhancedVolcano)
EnhancedVolcano(zero_2_dif_markers.DS,lab = rownames(zero_2_dif_markers.DS),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "cluster 0 vs. 2", pCutoff = 1e-90,FCcutoff = 0.25,labSize = 3.0,pointSize = 1,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F)


#------------plot clust freq by study   [granular clusters]
upd_B <- SetIdent(upd_B, value = "old.ident")
clust_study4 <- as.data.frame(prop.table(table(Idents(upd_B), upd_B$study), margin = 1))
ggplot(data=clust_study4, aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat="identity") + scale_fill_manual(values = c("#34568B", "#CD212A", "#FFA500", "#56C6A9", "#798EA4", "#FA7A35")) +theme_bw()

#-------------DEG across all B cells----------
#first remove Exhausted cluster since no K528R+ in there
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
upd_B_min_exh <- subset(upd_B, subset = cell_type_by_cluster == "Exhausted B cells", invert = T)
upd_B_min_exh <- SetIdent(upd_B_min_exh, value = "K528R_SNP")
DefaultAssay(upd_B_min_exh) <- "RNA"
K528R_DE.MAST <- FindMarkers(upd_B_min_exh, ident.1 = "yes", ident.2 = "no", only.pos = F, verbose = TRUE, test.use = "MAST", assay = "RNA", slot = "data")
saveRDS(K528R_DE.MAST, "./K528R_DE.MAST.rds")

#plot it (this looks nice)
EnhancedVolcano(K528R_DE.MAST,lab = rownames(K528R_DE.MAST),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "cluster 0 vs. 2", pCutoff = 1e-90,FCcutoff = 0.4,labSize = 3.0,pointSize = 2,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F, 
                drawConnectors = TRUE, arrowheads = F, maxoverlapsConnectors = 7,
                widthConnectors = 0.5)
#--------------------------------------------


#-------------DEG across Plasma cells----------
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
plasmacells <- subset(upd_B, subset = cell_type_by_cluster == "Plasmacells/blasts", invert = F)
plasmacells <- SetIdent(plasmacells, value = "K528R_SNP")
DefaultAssay(plasmacells) <- "RNA"
K528R_plasm_DE.MAST <- FindMarkers(plasmacells, ident.1 = "yes", ident.2 = "no", only.pos = F, verbose = TRUE, test.use = "MAST", assay = "RNA", slot = "data")
K528R_plasm_DE.MAST <- K528R_plasm_DE.MAST %>% rownames_to_column(var = "gene")
saveRDS(K528R_plasm_DE.MAST, "./K528R_plasm_DE.MAST.rds")
write_csv(K528R_plasm_DE.MAST, "./K528R_plasm_DE.MAST.csv")
K528R_plasm_DE.MAST <- readRDS("./K528R_plasm_DE.MAST.rds")

#plot it (this looks nice)
EnhancedVolcano(K528R_plasm_DE.MAST,lab = rownames(K528R_plasm_DE.MAST),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "Plasma cells/blasts K528R DEG", pCutoff = 5e-2,FCcutoff = 0.4,labSize = 3.0,pointSize = 2,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F, 
                drawConnectors = TRUE, arrowheads = F, maxoverlapsConnectors = 7,
                widthConnectors = 0.5,  xlim = c(-3, 3), ylim = c(-0, 7), cutoffLineWidth = 0.3, ylab = bquote(~-Log[10] ~Adj. ~italic(P)))
#--------------------------------------------



#-------------DEG across Naive B cells----------
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
naive <- subset(upd_B, subset = cell_type_by_cluster == "Naive B cells", invert = F)
naive <- SetIdent(naive, value = "K528R_SNP")
DefaultAssay(naive) <- "RNA"
K528R_naive_DE.MAST <- FindMarkers(naive, ident.1 = "yes", ident.2 = "no", only.pos = F, verbose = TRUE, test.use = "MAST", assay = "RNA", slot = "data")
saveRDS(K528R_naive_DE.MAST, "./K528R_naive_DE.MAST.rds")
write_csv(K528R_naive_DE.MAST, "./K528R_naive_DE.MAST.csv")
K528R_naive_DE.MAST <- readRDS("./K528R_naive_DE.MAST.rds")
K528R_naive_DE.MAST <- K528R_naive_DE.MAST %>% rownames_to_column(var = "gene")


#plot it (this looks nice)
EnhancedVolcano(K528R_naive_DE.MAST,lab = rownames(K528R_naive_DE.MAST),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "Naive K528R DEG", pCutoff = 1e-70,FCcutoff = 0.6,labSize = 3.0,pointSize = 2,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F, 
                drawConnectors = TRUE, arrowheads = F, maxoverlapsConnectors = 7,
                widthConnectors = 0.5, cutoffLineWidth = 0.3, xlim = c(-3, 3), ylab = bquote(~-Log[10] ~Adj. ~italic(P)))
#--------------------------------------------



#-------------DEG across class-switched B cells----------
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
class_swi <- subset(upd_B, subset = cell_type_by_cluster == "Class-switched memory B-cells", invert = F)
class_swi <- SetIdent(class_swi, value = "K528R_SNP")
DefaultAssay(class_swi) <- "RNA"
K528R_class_swi_DE.MAST <- FindMarkers(class_swi, ident.1 = "yes", ident.2 = "no", only.pos = F, verbose = TRUE, test.use = "MAST", assay = "RNA", slot = "data")
saveRDS(K528R_class_swi_DE.MAST, "./K528R_class_swi_DE.MAST.rds")
write_csv(K528R_class_swi_DE.MAST, "./K528R_class_swi_DE.MAST.csv")
K528R_class_swi_DE.MAST <- readRDS("./K528R_class_swi_DE.MAST.rds")
K528R_class_swi_DE.MAST <- K528R_class_swi_DE.MAST %>% rownames_to_column(var = "gene")


#plot it (this looks nice)
EnhancedVolcano(K528R_class_swi_DE.MAST,lab = rownames(K528R_class_swi_DE.MAST),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "Naive K528R DEG", pCutoff = 1e-10,FCcutoff = 0.6,labSize = 3.0,pointSize = 2,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F, 
                drawConnectors = TRUE, arrowheads = F, maxoverlapsConnectors = 7,
                widthConnectors = 0.5, cutoffLineWidth = 0.3, ylab = bquote(~-Log[10] ~Adj. ~italic(P)))
#--------------------------------------------



#-------------DEG across non-class-switched B cells----------
upd_B <- SetIdent(upd_B, value = "cell_type_by_cluster")
non_class_swi <- subset(upd_B, subset = cell_type_by_cluster == "Non-switched memory B cells", invert = F)
non_class_swi <- SetIdent(non_class_swi, value = "K528R_SNP")
DefaultAssay(non_class_swi) <- "RNA"
K528R_non_class_swi_DE.MAST <- FindMarkers(non_class_swi, ident.1 = "yes", ident.2 = "no", only.pos = F, verbose = TRUE, test.use = "MAST", assay = "RNA", slot = "data")
saveRDS(K528R_non_class_swi_DE.MAST, "./K528R_non_class_swi_DE.MAST.rds")
write_csv(K528R_non_class_swi_DE.MAST, "./K528R_non_class_swi_DE.MAST.csv")
K528R_non_class_swi_DE.MAST <- readRDS("./K528R_non_class_swi_DE.MAST.rds")
K528R_non_class_swi_DE.MAST <- K528R_non_class_swi_DE.MAST %>% rownames_to_column(var = "gene")


#plot it (this looks nice)
EnhancedVolcano(K528R_non_class_swi_DE.MAST,lab = rownames(K528R_non_class_swi_DE.MAST),x = 'avg_log2FC',y = 'p_val_adj', 
                title = "Naive K528R DEG", pCutoff = 1e-10,FCcutoff = 0.6,labSize = 3.0,pointSize = 2,
                col=c('black', 'black', 'black', 'red3'),colAlpha = 1, gridlines.major = F, gridlines.minor = F, 
                drawConnectors = TRUE, arrowheads = F, maxoverlapsConnectors = 7,
                widthConnectors = 0.5, cutoffLineWidth = 0.3, ylab = bquote(~-Log[10] ~Adj. ~italic(P)))
#--------------------------------------------


#------pathway anaylsis (running in IPA)----------------

#make clean pathway figs (plasma cells)
plasma_K528R_path <-  read_xls("./plasma_K528R_path.xls") %>% as.data.frame() 
plasma_K528R_path <- dplyr::rename(plasma_K528R_path, neg.Log.p_val = `-log(p-value)`, Z_score = `z-score`, Pathways = `Ingenuity Canonical Pathways`)
plasma_K528R_path <-  plasma_K528R_path %>% filter(neg.Log.p_val > 2)

#since all are NA, but one, color just the one bar
area.color <- c("#248f8f", NA, NA, NA, NA, NA)
ggplot(plasma_K528R_path, aes(x=reorder(Pathways, neg.Log.p_val), y = neg.Log.p_val, fill=as.factor(Z_score))) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_manual(values = c("#6B5876", "#8c8c8c", "#8c8c8c", "#8c8c8c", "#8c8c8c", "#8c8c8c")) +
  theme(axis.line = element_line(linetype = "solid"), 
          axis.ticks = element_line(size = 0.9), 
          axis.title = element_text(size = 15), 
          axis.text = element_text(size = 12, colour = "black"), 
          axis.text.x = element_text(size = 12), 
          panel.background = element_rect(fill = NA)) + labs( x="Pathways", y="-Log(p-value)") 


#make clean pathway figs (naive B cells)
naive_K528R_path <-  read_xls("./naive_K528R_path.xls") %>% as.data.frame() 
naive_K528R_path <- dplyr::rename(naive_K528R_path, neg.Log.p_val = `-log(p-value)`, Z_score = `z-score`, Pathways = `Ingenuity Canonical Pathways`)
naive_K528R_path <-  naive_K528R_path %>% filter(neg.Log.p_val > 6)

#area.color <- c("#248f8f", NA, NA, NA, NA, NA)
ggplot(naive_K528R_path, aes(x=reorder(Pathways, neg.Log.p_val), y = neg.Log.p_val, fill=Z_score)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_gradient2(low='#EDD59E', mid='#ffffff', high='#6B5876', space='Lab') +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(size = 0.9), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        panel.background = element_rect(fill = NA)) + labs( x="Pathways", y="-Log(p-value)") 



#make clean pathway figs (class-swiched B cells)
class_swi_K528R_path <-  read_xls("./class_swi_K528R_path.xls") %>% as.data.frame() 
class_swi_K528R_path <- dplyr::rename(class_swi_K528R_path, neg.Log.p_val = `-log(p-value)`, Z_score = `z-score`, Pathways = `Ingenuity Canonical Pathways`)
class_swi_K528R_path <-  class_swi_K528R_path %>% filter(neg.Log.p_val >= 3)

ggplot(class_swi_K528R_path, aes(x=reorder(Pathways, neg.Log.p_val), y = neg.Log.p_val, fill=Z_score)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_gradient2(low='#EDD59E', mid='#ffffff', high='#6B5876', space='Lab') +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(size = 0.9), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        panel.background = element_rect(fill = NA)) + labs( x="Pathways", y="-Log(p-value)") 


#make clean pathway figs (non-swiched B cells)
non_swi_K528R_path <-  read_xls("./non_swi_K528R_path.xls") %>% as.data.frame() 
non_swi_K528R_path <- dplyr::rename(non_swi_K528R_path, neg.Log.p_val = `-log(p-value)`, Z_score = `z-score`, Pathways = `Ingenuity Canonical Pathways`)
non_swi_K528R_path <-  non_swi_K528R_path %>% filter(neg.Log.p_val >= 4)

ggplot(non_swi_K528R_path, aes(x=reorder(Pathways, neg.Log.p_val), y = neg.Log.p_val, fill=Z_score)) + 
  geom_bar(stat = 'identity') + 
  coord_flip() +
  scale_fill_gradient2(low='#EDD59E', mid='#ffffff', high='#6B5876', space='Lab') +
  theme(axis.line = element_line(linetype = "solid"), 
        axis.ticks = element_line(size = 0.9), 
        axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(size = 12), 
        panel.background = element_rect(fill = NA)) + labs( x="Pathways", y="-Log(p-value)") 
#-------------------------------------------------------

