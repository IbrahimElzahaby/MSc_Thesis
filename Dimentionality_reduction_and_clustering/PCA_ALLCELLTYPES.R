# import libraries
library(readxl)
library(tidyverse)
library(FactoMineR)
library(dplyr)
library(factoextra)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(cluster)
###########################################################
# PCA for Hela merged sets
h4 <- read.csv("h4intersect.csv")
h5 <- read.csv("h5intersect.csv")
h5_h4_merged <- merge(h5,h4, by='Protein')
h5_h4_merged <- as.data.frame(t(h5_h4_merged[,-1]))
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)
celltype_h5 <- c(replicate('Helacells5',n=dim(h5)[2]-1))
celltype_h4 <- c(replicate('Helacells4',n=dim(h4)[2]-1))
celltype = c(celltype_h5,celltype_h4)
pca.merged <- PCA(h5_h4_merged[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
HM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
merged_helaPCA <- ggpar(HM_merged,
                        title = "PCA_HelaCellTypes",
                        xlab = "PC1", ylab = "PC2",
                        legend.title = "Cell type", legend.position = "top",
                        ggtheme = theme_minimal())
merged_helaPCA
##########################################################
# PCA for Melanoma merged sets
melanoma1 <- read.csv("m1intersect.csv")
melanoma3 <- read.csv("m3intersect.csv")
m1_m3_merged <- merge(melanoma1,melanoma3, by='proteins')
m1_m3_merged <- as.data.frame(t(m1_m3_merged[,-1]))
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)
celltype_m1 <- c(replicate('Melanoma1',n=dim(melanoma1)[2]-1))
celltype_m3 <- c(replicate('Melanoma3',n=dim(melanoma3)[2]-1))
celltype = c(celltype_m1,celltype_m3)
pca.merged <- PCA(m1_m3_merged [,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
MM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
merged_melanomaPCA <- ggpar(MM_merged,
                            title = "PCA_MelanomacellTypes",
                            xlab = "PC1", ylab = "PC2",
                            legend.title = "Cell type", legend.position = "top",
                            ggtheme = theme_minimal())
merged_melanomaPCA
##########################################################
# PCA for Monocytes merged sets
U1_cells <-read_excel("U937_1.xlsx")
U3_cells <-read_excel("U937_3.xlsx")
U4_cells <-read_excel("U937_4.xlsx")
U5_cells <-read_excel("U937_5.xlsx")
U6_cells <-read_excel("U937_6.xlsx")
monocytes1 <- read.csv("mono1intersect.csv")
monocytes3 <- read.csv("mono3intersect.csv")
monocytes4 <- read.csv("mono4intersect.csv")
monocytes5 <- read.csv("mono5intersect.csv")
monocytes6 <- read.csv("mono6intersect.csv")
U1_U3_merged <- merge(monocytes1,monocytes3, by='proteins')
U4_U5_merged <- merge(monocytes4,monocytes5, by='proteins')
U4_merged <- merge(U1_U3_merged, U4_U5_merged, by='proteins')
U_merged <- merge(U4_merged, monocytes6, by='proteins')
U_merged <- as.data.frame(t(U_merged[,-1]))
nb.cols <- 5
mycolors <- colorRampPalette(brewer.pal(5, "Set1"))(nb.cols)
celltype_U1 <- c(replicate('Monocytes1',n=dim(U1_cells)[2]-1))
celltype_U3 <- c(replicate('Monocytes3',n=dim(U3_cells)[2]-1))
celltype_U4 <- c(replicate('Monocytes4',n=dim(U4_cells)[2]-1))
celltype_U5 <- c(replicate('Monocytes5',n=dim(U5_cells)[2]-1))
celltype_U6 <- c(replicate('Monocytes6',n=dim(U6_cells)[2]-1))
celltype = c(celltype_U1,celltype_U3,celltype_U4,celltype_U5,celltype_U6)
pca.merged <- PCA(U_merged [,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
MO_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
merged_monocytesPCA <- ggpar(MO_merged,
                             title = "PCA_MonocytescellTypes",
                             xlab = "PC1", ylab = "PC2",
                             legend.title = "Cell type", legend.position = "top",
                             ggtheme = theme_minimal())

merged_monocytesPCA
##################################################################################
# PCA for PDAC merged sets
PDAC2_cells <-read_excel("pdac2.xlsx")
PDAC3_cells <-read_excel("pdac3.xlsx")              
pdac2 <- read.csv("p2intersect.csv")
col_names<- pdac2$proteins
#col_names
pdac2 <- as.data.frame(t(pdac2[,-1]))
colnames(pdac2) <- col_names
proteins <-  row.names(pdac2)
pdac2 <- cbind(proteins,pdac2)
#proteins
row.names(pdac2) <- NULL 
pdac3 <- read.csv("p3intersect.csv")
# PCA PDAC merged sets
p2_p3_merged <- merge(pdac2, pdac3, by='proteins')
p2_p3_merged <- as.data.frame(t(p2_p3_merged[,-1]))
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)
celltype_Pdac2 <- c(replicate('PDAC2',n=dim(PDAC2_cells)[2]-1))
celltype_Pdac3 <- c(replicate('PDAC3',n=dim(PDAC3_cells)[2]-1))
celltype = c(celltype_Pdac2,celltype_Pdac3)
pca.merged <- PCA(p2_p3_merged[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
PD_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
MO_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
merged_pdacPCA <- ggpar(PD_merged,
                             title = "PCA_PDACcellTypes",
                             xlab = "PC1", ylab = "PC2",
                             legend.title = "Cell type", legend.position = "top",
                             ggtheme = theme_minimal())
merged_pdacPCA
###################################################################################
# Draw the PCA diagrams
res <- 300
w <- 14
h <- 10
png("ALL_pca.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(merged_monocytesPCA, merged_helaPCA, merged_melanomaPCA, merged_pdacPCA, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()



