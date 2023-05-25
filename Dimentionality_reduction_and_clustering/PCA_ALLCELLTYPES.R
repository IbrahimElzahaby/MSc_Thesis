
# load required libraries
library(readxl)
library(tidyverse)
library(FactoMineR)
library(dplyr)
library(factoextra)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(cluster)

# Read Hela data from csv files
h4 <- read.csv("h4intersect.csv")
h5 <- read.csv("h5intersect.csv")

# Merge Hela cell datasets
h5_h4_merged <- merge(h5,h4, by='Protein')
h5_h4_merged <- as.data.frame(t(h5_h4_merged[,-1]))

# Set colors for cell types
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)

# Create labels for Hela cell types
celltype_h5 <- c(replicate('Helacells5',n=dim(h5)[2]-1))
celltype_h4 <- c(replicate('Helacells4',n=dim(h4)[2]-1))
celltype = c(celltype_h5,celltype_h4)

# Perform PCA analysis on Hela merged dataset
pca.merged <- PCA(h5_h4_merged[,-1], scale.unit = TRUE, graph = F)

# Visualize eigenvalues
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))

# Visualize PCA results with colored individual points and ellipses for Hela cell types
HM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)

# Create a combined plot with PCA results for Hela cell types
merged_helaPCA <- ggpar(HM_merged,
                        title = "PCA_HelaCellTypes",
                        xlab = "PC1", ylab = "PC2",
                        legend.title = "Cell type", legend.position = "top",
                        ggtheme = theme_minimal())

###############################################################################

# Read Melanoma data from csv files
melanoma1 <- read.csv("m1intersect.csv")
melanoma3 <- read.csv("m3intersect.csv")

# Merge Melanoma cell datasets
m1_m3_merged <- merge(melanoma1,melanoma3, by='proteins')
m1_m3_merged <- as.data.frame(t(m1_m3_merged[,-1]))

# Set colors for Melanoma cell types
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)

# Create labels for Melanoma cell types
celltype_m1 <- c(replicate('Melanoma1',n=dim(melanoma1)[2]-1))
celltype_m3 <- c(replicate('Melanoma3',n=dim(melanoma3)[2]-1))
celltype = c(celltype_m1,celltype_m3)

# Perform PCA analysis on Melanoma merged dataset
pca.merged <- PCA(m1_m3_merged [,-1], scale.unit = TRUE, graph = F)

# Visualize eigenvalues
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))

# Visualize PCA results with colored individual points and ellipses for Melanoma cell types
MM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)

# Create a combined plot with PCA results for Melanoma cell types
merged_melanomaPCA <- ggpar(MM_merged,
                            title = "PCA_MelanomacellTypes",
                            xlab = "PC1", ylab = "PC2",
                            legend.title = "Cell type", legend.position = "top",
                            ggtheme = theme_minimal())

###############################################################################

# Read Monocytes data from csv and excel files
monocytes1 <- read.csv("mono1intersect.csv")
U1_cells <-read_excel("U937_1.xlsx")

monocytes3 <- read.csv("mono3intersect.csv")
U3_cells <-read_excel("U937_3.xlsx")

monocytes4 <- read.csv("mono4intersect.csv")
U4_cells <-read_excel("U937_4.xlsx")

monocytes5 <- read.csv("mono5intersect.csv")
U5_cells <-read_excel("U937_5.xlsx")

monocytes6 <- read.csv("mono6intersect.csv")
U6_cells <-read_excel("U937_6.xlsx")

# Merge Monocyte cell datasets
U1_U3_merged <- merge(monocytes1,monocytes3, by='proteins')
U4_U5_merged <- merge(monocytes4,monocytes5, by='proteins')
U4_merged <- merge(U1_U3_merged, U4_U5_merged, by='proteins')
U_merged <- merge(U4_merged, monocytes6, by='proteins')
U_merged <- as.data.frame(t(U_merged[,-1]))

# Set colors for Monocyte cell types
nb.cols <- 5
mycolors <- colorRampPalette(brewer.pal(5, "Set1"))(nb.cols)

# Create labels for Monocyte cell types
celltype_U1 <- c(replicate('Monocytes1',n=dim(U1_cells)[2]-1))
celltype_U3 <- c(replicate('Monocytes3',n=dim(U3_cells)[2]-1))
celltype_U4 <- c(replicate('Monocytes4',n=dim(U4_cells)[2]-1))
celltype_U5 <- c(replicate('Monocytes5',n=dim(U5_cells)[2]-1))
celltype_U6 <- c(replicate('Monocytes6',n=dim(U6_cells)[2]-1))
celltype = c(celltype_U1,celltype_U3,celltype_U4,celltype_U5,celltype_U6)

# Perform PCA analysis on Monocyte merged dataset
pca.merged <- PCA(U_merged [,-1], scale.unit = TRUE, graph = F)

# Visualize eigenvalues
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))

# Visualize PCA results with colored individual points and ellipses for Monocyte cell types
MO_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)

# Create a combined plot with PCA results for Monocyte cell types
merged_monocytesPCA <- ggpar(MO_merged,
                             title = "PCA_MonocytescellTypes",
                             xlab = "PC1", ylab = "PC2",
                             legend.title = "Cell type", legend.position = "top",
                             ggtheme = theme_minimal())

###############################################################################

# Read PDAC data from csv files
pdac2 <- read.csv("p2intersect.csv")
# Extract column names from the 'proteins' column
col_names<- pdac2$proteins
# Transpose the data frame and set column names
pdac2 <- as.data.frame(t(pdac2[,-1]))
# Set the column names using the extracted names
colnames(pdac2) <- col_names
# Extract the protein names as row names
proteins <-  row.names(pdac2)
# Add the protein names as a new column
pdac2 <- cbind(proteins,pdac2)
# Reset the row names to default values
row.names(pdac2) <- NULL
PDAC2_cells <-read_excel("pdac2.xlsx")
pdac3 <- read.csv("p3intersect.csv")
PDAC3_cells <-read_excel("pdac3.xlsx")

# Merge PDAC cell datasets
p2_p3_merged <- merge(pdac2, pdac3, by='proteins')
p2_p3_merged <- as.data.frame(t(p2_p3_merged[,-1]))

# Set colors for PDAC cell types
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)

# Create labels for PDAC cell types
celltype_Pdac2 <- c(replicate('PDAC2',n=dim(PDAC2_cells)[2]-1))
celltype_Pdac3 <- c(replicate('PDAC3',n=dim(PDAC3_cells)[2]-1))
celltype = c(celltype_Pdac2,celltype_Pdac3)

# Perform PCA analysis on the PDAC merged dataset
pca.merged <- PCA(p2_p3_merged[,-1], scale.unit = TRUE, graph = F)

# Visualize eigenvalues
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))

# Visualize PCA results with colored individual points and ellipses for PDAC cell types
PD_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)

# Create a combined plot with PCA results for PDAC cell types
merged_pdacPCA <- ggpar(PD_merged,
                             title = "PCA_PDACcellTypes",
                             xlab = "PC1", ylab = "PC2",
                             legend.title = "Cell type", legend.position = "top",
                             ggtheme = theme_minimal())

###############################################################################

# Save PCA plots to png file
res <- 300
w <- 14
h <- 10
png("ALL_pca.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(merged_monocytesPCA, merged_helaPCA, 
                   merged_melanomaPCA, merged_pdacPCA, ncol = 2,
                   labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()



