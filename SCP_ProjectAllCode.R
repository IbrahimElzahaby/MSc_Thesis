# import libraries
library(readxl)
library(VennDiagram)
library(clusterProfiler)
library(enrichplot)
library(VennDiagram)
library(tidyverse)
library(FactoMineR)
library(dplyr)
library(factoextra)
library(ggpubr)
library(ggplot2)
library(DOSE)
library(BiocManager)
library(base)
library(corrplot)
library(RColorBrewer)
library(visNetwork)
library(igraph)
library(corrr)
# git add SCP_ProjectAllCode.R
# git commit -m '...'
# git push origin master
# Hela intersected proteins
helacells4 <-read_excel("hc4.xlsx")
helalist4 <- helacells4$protein
helacells5 <-read_excel("hc5.xlsx")              
helalist5 <- helacells5$Protein
venn.diagram(list("Helacells4"=helalist4, "Helacells5"=helalist5),
             fill = c("cyan","red"),
             cex= 0.5, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.3, filename = "Hela.png",
             print.mode =c("raw", "percent"))
# Hela overlapped protein names
helaoverlap <- intersect(helalist4,helalist5)
hela4_unique <- setdiff(helalist4, helalist5)
hela5_unique <- setdiff(helalist5, helalist4)
write.table(helaoverlap, file = "Hela_intersected.txt", sep = "")
write.table(hela4_unique, file = "hela4_unique.txt", sep = "")
write.table(hela5_unique, file = "hela5_unique.txt", sep = "")
# Hela_G0terms_KEGG pathways
hela <- read.table(file = "Hela_intersected.txt")
hela_prots <- write.table(hela$x, file = "helacells_proteins.txt")
helaGO_unique4 <- read.table(file = "hela4_unique.txt")
helaGO_unique5 <- read.table(file = "hela5_unique.txt")
hela_prots4 <- write.table(helaGO_unique4$x, file = "hela_unique4.txt",
                          quote=F, row.names=F, col.names = F)
hela_prots5 <- write.table(helaGO_unique5$x, file = "hela_unique5.txt",
                          quote=F, row.names=F, col.names = F)
# Hela_GOterms
hela_genes <- read_excel("hela_idmap.xlsx")
hela_go <- enrichGO(gene = hela_genes$entrezgene, ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    readable=TRUE,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.3)
corex <- read_excel("trial.xlsx")
?read_excel
m <- t(corex)
colnames(m) <- m[1,]
rownames(m)
m <- as.data.frame(as.numeric(m[-1,]))
class(m[1,2])
n <- cor(corex[,2:50])
corrplot(n, method = 'circle', order = "AOE", type = 'upper')
k <- cor(h5_h4_merged[,2:50])
corrplot(k, method = 'circle', order = "AOE", type = 'upper')



h4 <- read.delim2("h4intersect.txt", sep = "\t", stringsAsFactors = F)
h5 <- read.delim2("h5intersect.txt", sep = "\t", stringsAsFactors = F)
h5_h4_merged <- merge(h5,h4, by='Proteins')
hela_cor <- as.data.frame(h4)
hela_cor <- t(hela_cor)
colnames(hela_cor) <- hela_cor[1,]
hela_cor <- hela_cor[-1,]
#hela_cor <- as.data.frame(sapply(hela_cor, as.numeric))
#library(utility)
hj <- type.convert(hela_cor)
class(hj[2,2])
hela_cor <- as.numeric(hela_cor)
library(GGally)
library(qgraph)
library(ggraph)
library(tidyr)
library(dplyr)
library(devtools)
install_github("VanderJag/anvis")
library(anvis)
net1 <- adjToNetwork(adj_mats = shit_happen,
                     node_attrs = "all",
                     edge_attrs = "all",
                     output_as = "igraph")
anvis(net1)
shit_happen <- cor(hj, method = "spearman")
shit_happen <- as.data.frame(shit_happen)
# Reshape the data frame into a long format for plotting
cor_df_long <- shit_happen %>%
  tibble::rownames_to_column("from") %>%
  gather(to, value, -from)
# Filter for correlations above a threshold
cor_df_long_filtered <- cor_df_long %>%
  filter(abs(value) > 0.7)
# Plot the filtered correlations as a network graph
library(ggraph)
library(gridExtra)
ggraph(cor_df_long_filtered, layout = "fr") +
  geom_edge_link(aes(x = 0, y = 0, xend = value, yend = 0, alpha = abs(value)), 
                 start_cap = circle(0.15, "mm"),
                 end_cap = circle(0.15, "mm")) +
  geom_node_point() +
  scale_alpha_continuous(range = c(0.1, 1)) + 
  
  theme_void()
  
# Plot 
ggraph(shit_happen, layout="fr") + 
  geom_edge_link(aes(c, width = log(weight)), alpha = 0.5, 
                 start_cap = circle(2, 'mm'), end_cap = circle(2, 'mm')) +
  scale_edge_width(range = c(0.5, 2.5)) + 
  geom_node_point(color = (shit_happen)$color, size = 5, alpha = 0.5) +
  geom_node_text(aes(label = name), repel = TRUE) +
  theme_void() + 
  theme(legend.position = "none")


ig_loc <- shit_happen %>% 
  select()



hj %>% correlate() %>% 
  network_plot(hj)


shit_happen[] <- ifelse(abs(shit_happen) >0.5, 1, 0)
View(shit_happen)


cor(df$Count, df$qvalue, method = "spearman")
cor(df$Count, df$pvalue, method = "spearman")
ggcorr(df)
ggcorr(df, nbreaks = 6,
       low = "steelblue",
       mid = "white",
       high = "darkred",
       geom = "circle")
ggcorr(n,
       nbreaks = 6,
       label = TRUE,
       label_size = 3,
       color = "grey50")
h4 <- read.csv("cor.xlsx")
k <- cor(h4)
corrplot(helacells4, method = 'number')
corrplot()
hg <- pairwise_termsim(hela_go)
emapplot(hg)
hgcenter <- setReadable(hela_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(hgcenter)
cnetplot(hgcenter,  foldChange=hela_go, circular = TRUE, colorEdge = TRUE)
hela4_genes <- read_excel("hela4uni_idmap.xlsx")
hela4_go <- enrichGO(gene = hela4_genes$entrezgene, ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    readable=TRUE,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(hela4_go)
hela5_genes <- read_excel("hela5uni_idmap.xlsx")
hela5_go <- enrichGO(gene = hela5_genes$entrezgene, ont = "BP",
                     OrgDb ="org.Hs.eg.db",
                     readable=TRUE,
                     pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(hela5_go)
# Hela_KEGG pathways
#data(hela_genes$entrezgene, package='DOSE')
hela_kegg <- enrichKEGG(gene = hela_genes$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3, keyType = "kegg")
dotplot(hela_kegg)
he <- pairwise_termsim(hela_kegg)
emapplot(he)
hecenter <- setReadable(hela_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(hecenter)
cnetplot(hecenter,  foldChange=hela_kegg, circular = TRUE, colorEdge = TRUE)
hela4_kegg <- enrichKEGG(gene = hela4_genes$entrezgene,
                        organism = 'hsa',
                        keyType = "kegg",
                        use_internal_data = F, pAdjustMethod = "none",
                        qvalueCutoff = 0.3, pvalueCutoff = 0.05)
# , qvalueCutoff = 0.3
dotplot(hela4_kegg)
hela5_kegg <- enrichKEGG(gene = hela5_genes$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)
dotplot(hela5_kegg)
###########################################3
# PCA for hela
# Hela intersected proteins
helapca4 <-read_excel("hc4.xlsx")
helalistpca4 <- helacells4$protein
helapca5 <-read_excel("hc5.xlsx")              
helalistpca5 <- helacells5$Protein
# Hela overlapped protein names
helapca_overlap <- intersect(helalistpca4,helalistpca5)
hela4_df4 <- filter(helapca4, protein %in% helapca_overlap)
hela5_df5 <- filter(helapca5, Protein %in% helapca_overlap)
write.csv(hela4_df4, "h4intersect.csv", row.names = F)
write.csv(hela5_df5, "h5intersect.csv", row.names = F)
h4 <- read.csv2("h4intersect.csv", sep = ",")
h5 <- read.csv2("h5intersect.csv", sep = ",")
h5_h4_merged <- merge(h5,h4, by='Protein')
#h5_h4_merged <- as.data.frame(t(h5_h4_merged[,-1]))

nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)
celltype_h5 <- c(replicate('Helacells5',n=dim(h5)[2]-1))
celltype_h4 <- c(replicate('Helacells4',n=dim(h4)[2]-1))
celltype = c(celltype_h5,celltype_h4)
pca.merged <- PCA(h5_h4_merged [,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
#merged_pca <- fviz_pca_var(pca.merged, col.var = "Cell Type",
#            gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
#           repel = TRUE)
HM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                   palette = mycolors)
merged_helaPCA <- ggpar(HM_merged,
      title = "PCA_HelaCellTypes",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cell type", legend.position = "top",
      ggtheme = theme_minimal())

merged_helaPCA
class(h5_h4_merged$A6NFK2)
pca_plot

ggpar(pca_plot,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cell type", legend.position = "top",
      ggtheme = theme_minimal())



###
pca.h5 <- PCA(h5[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.h5, addlabels = TRUE, ylim = c(0, 40))
h4_pca <- fviz_pca_var(pca.h5, col.var = "cos2",
                       gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
                       repel = TRUE)
HM_h5 <- fviz_pca_ind(pca.h5, col.ind = "cos2", 
                      gradient.cols = c("#64C8E9", "#96FA64", "#C86496", "#FAC896"), 
                      repel = TRUE)
pca_plot +  ggpar(HM_h5,
                  title = "PCA_pdacCells_set3",
                  xlab = "PC1", ylab = "PC2",
                  legend.title = "cos2", legend.position = "top",
                  ggtheme = theme_minimal())
pca_plot
ggpar(HM_h5,
      title = "PCA_pdacCells_set3",
      xlab = "PC1", ylab = "PC2",
      legend.title = "cos2", legend.position = "top",
      ggtheme = theme_minimal())

hela4 <- read.csv("h4intersect.csv")
pca.hela4 <- PCA(hela4[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.hela4, addlabels = TRUE, ylim = c(0, 40))
fviz_pca_var(pca.hela4, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)
pca.hela4 <- PCA(t(hela4[,-1]), scale.unit = TRUE, graph = F)
h4 <- fviz_pca_ind(pca.hela4, col.ind = "cos2", 
                   gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                   repel = TRUE)
ggpar(h4,
      title = "PCA_helaCells_set4",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
hela5 <- read.csv("h5intersect.csv")
pca.hela5 <- PCA(t(hela5[,-1]), scale.unit = TRUE, graph = T)
fviz_eig(pca.hela5, addlabels = TRUE, ylim = c(0, 40))
fviz_pca_var(pca.hela5, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)
pca.hela5 <- PCA(t(hela5[,-1]), scale.unit = TRUE, graph = F)
h5 <- fviz_pca_ind(pca.hela5, col.ind = "cos2", 
                   gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                   repel = TRUE)
ggpar(h5,
      title = "PCA_helaCells_set5",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
########################################################
# PDAC intersected proteins
PDAC2_cells <-read_excel("pdac2.xlsx")
PDAC_list2 <- PDAC2_cells$proteins
PDAC3_cells <-read_excel("pdac3.xlsx")              
PDAC_list3 <- PDAC3_cells$prot
venn.diagram(list("PDAC2"=PDAC_list2, "PDAC3"=PDAC_list3),
             fill = c("cyan","red"),
             cex= 0.5, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.3, filename = "PDAC.png",
             print.mode =c("raw", "percent"))
# PDAC overlapped protein names
PDACoverlap <- intersect(PDAC_list2,PDAC_list3)
pdac2_unique <- setdiff(PDAC_list2, PDAC_list3)
pdac3_unique <- setdiff(PDAC_list3, PDAC_list2)
write.table(PDACoverlap, file = "PDAC_intersected.txt", sep = "")
write.table(pdac2_unique, file = "pdac2_unique.txt", sep = "")
write.table(pdac3_unique, file = "pdac3_unique.txt", sep = "")
# PDAC GO terms
pdacGO_inetrsect <- read.table(file = "PDAC_intersected.txt")
pdacGO_unique2 <- read.table(file = "pdac2_unique.txt")
pdacGO_unique3 <- read.table(file = "pdac3_unique.txt")
pdac_intersect <- write.table(pdacGO_inetrsect$x, file = "pdac_intersect.txt",
                           quote=F, row.names=F, col.names = F)
pdac_prots2 <- write.table(pdacGO_unique2$x, file = "pdac_unique2.txt",
                           quote=F, row.names=F, col.names = F)
pdac_prots3 <- write.table(pdacGO_unique3$x, file = "pdac_unique3.txt",
                           quote=F, row.names=F, col.names = F)
# PDAC set2 GO terms
pdac <- read_excel("idmap.xlsx")
pdacGO <- enrichGO(gene = pdac$entrezgene, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(pdacGO)
ridgeplot(pdacGO)
pd <- pairwise_termsim(pdacGO)
emapplot(pd)
pdcenter <- setReadable(pdacGO, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(pdcenter)
cnetplot(pdcenter,  foldChange=pdacGO, circular = TRUE, colorEdge = TRUE)
pdac2_genes <- read_excel("pdac2uni_idmap.xlsx")
pdac2_GO <- enrichGO(gene = pdac2_genes$entrezgene, ont = "ALL",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(pdac2_GO)
# PDAC set3 GO terms
pdac3_genes <- read_excel("pdac3uni_idmap.xlsx")
pdac3_GO <- enrichGO(gene = pdac3_genes$entrezgene, ont = "ALL",
                     OrgDb ="org.Hs.eg.db",
                     readable=TRUE,
                     pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(pdac3_GO)
# PDAC KEGG pathways set 2
pdac_kegg <- enrichKEGG(gene = pdac$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(pdac_kegg)
pdkegg <- pairwise_termsim(pdac_kegg)
emapplot(pdkegg)
pkcenter <- setReadable(pdac_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(pkcenter)
cnetplot(pkcenter,  foldChange=pdac_kegg, circular = TRUE, colorEdge = TRUE)
pdac2_KEGG <- enrichKEGG(gene = pdac2_genes$entrezgene,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
dotplot(pdac2_KEGG) #+ ggtitle("KEGG")
# PDAC KEGG pathways set3
pdac3_KEGG <- enrichKEGG(gene = pdac3_genes$entrezgene,
                       organism = 'hsa',
                       pvalueCutoff = 0.05)
dotplot(pdac3_KEGG) #+ ggtitle("KEGG")
##############################
# PCA for PDAC
# PDAC intersected proteins
PDAC2_cells <-read_excel("pdac2.xlsx")
PDAC_list2 <- PDAC2_cells$proteins
PDAC3_cells <-read_excel("pdac3.xlsx")              
PDAC_list3 <- PDAC3_cells$prot
# PDAC overlapped protein names
PDACoverlap <- intersect(PDAC_list2,PDAC_list3)
pdac2_df2 <- filter(PDAC2_cells, proteins %in% PDACoverlap)
pdac3_df3 <- filter(PDAC3_cells, prot %in% PDACoverlap)
write.csv(pdac2_df2, "p2intersect.csv", row.names = F)
write.csv(pdac3_df3, "p3intersect.csv", row.names = F)
pdac2 <- read.csv("p2intersect.csv")
pdac3 <- read.csv("p3intersect.csv")
# PCA PDAC merged sets
p2_p3_merged <- merge(pdac2,pdac3, by='proteins')
p2_p3_merged <- as.data.frame(t(p2_p3_merged[,-1]))
nb.cols <- 2
mycolors <- colorRampPalette(brewer.pal(2, "Set1"))(nb.cols)
###########################UMAP############
library(tidyverse)
library(palmerpenguins)
library(umap)
theme_set(theme_bw(18))
pdac2 <- read.csv("p2intersect.csv")
hela5 <- read.csv("h5intersect.csv")
View(hela5)
View(pdac2)
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
library(STRINGdb)
string_db <- STRINGdb$new( version="11.5", species=9606,
                           score_threshold=400, input_directory="")
c1 <- as.data.frame(read.csv("avgC2_pdac.csv"))
mappedc1 <- string_db$map(c1, 'proteins')
string_db$plot_network(mappedc1$STRING_id)
enrichmentGOc1 <- string_db$get_enrichment(mappedc1$STRING_id)
enc1 <- data.frame(fd= enrichmentGOc1$fdr, p_val= enrichmentGOc1$p_value, 
                 description=enrichmentGOc1$description)
#en3 <- en[order(p_value),]
#dotplot(en)
enc1 <- subset(enc1, fd <0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001)
dotchart(enc1$p_val, labels = enc1$description, bg = '#56B4E9', xlab = "p_value")
########
sce <- SingleCellExperiment(assays = (counts = as.matrix(hela5)))
hela5_umap <- umap(hela5[,-1])
head(hela5_umap$layout)
umap_df_hela5 <- data.frame(umap1=hela5_umap$layout[,1],
                      
                      umap2=hela5_umap$layout[,2],
                      cell=hela5$Protein
                      #labels= labels
                      )

#library(lattice)
dotplot(data=en,  description ~ ratio.se | number_of_genes)
head(en)
dotchart(x=en$ratio.se, y=en$description)

View(umap_df_hela5)
pdac2_umap$knn
#ggplot(data= umap_df, aes(x=umap1,y=umap2 #color=label, 
#                          ))+
#  geom_point()
dist = dist(hela5_umap$layout)
helaclust <- hclust(dist, method = "complete")
plot(helaclust)
helaclusters <- cutree(helaclust, k=2) 


umap_df_hela5$hcluster_label <- as.factor(helaclusters)
head(umap_df_hela5)
ggplot(data= umap_df_hela5,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="Hela subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
#ggsave("PDAC_Clusters")

cluster_1 <- subset(umap_df,cluster_label== '1')
cluster_2 <- subset(umap_df, cluster_label=='2')
cluster_3 <- subset(umap_df, cluster_label=='3')
View(cluster_1)
plot.PCA(pdac2, x = 1, y = 2, n = 206, point_size = 3)
clust1 <- write.table(cluster_1, file = "pdac_clust1.txt",
                           quote=F, row.names=F, col.names = F)
clust2 <- write.table(cluster_2, file = "pdac_clust2.txt",
                      quote=F, row.names=F, col.names = F)
clust3 <- write.table(cluster_3, file = "pdac_clust3.txt",
                      quote=F, row.names=F, col.names = F)




d <- read.table("pdac_clust1.txt")

BiocManager::available("DEP")
library(DEP)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("STRINGdb")
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
################
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

celltype_p2 <- c(replicate('PDAC2',n=dim(pdac2)[2]-1))

celltype_p3 <- c(replicate('PDAC3',n=dim(pdac3)[2]-1))
celltype = c(celltype_p2,celltype_p3)
pca.merged <- PCA(celltype_p2[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))

pca.merged <- PCA(p2_p3_merged [,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.merged, addlabels = TRUE, ylim = c(0, 20))
#merged_pca <- fviz_pca_var(pca.merged, col.var = "Cell Type",
#            gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
#           repel = TRUE)
PM_merged <- fviz_pca_ind(pca.merged , col.ind = celltype, addEllipses = T,
                          palette = mycolors)
merged_pdacPCA <- ggpar(PM_merged,
                        title = "PCA_PDACcellTypes",
                        xlab = "PC1", ylab = "PC2",
                        legend.title = "Cell type", legend.position = "top",
                        ggtheme = theme_minimal())

merged_pdacPCA
pca.p2_p3_merged <- PCA(p2_p3_merged[,-1], scale.unit = TRUE, graph = F)
fviz_eig(pca.p2_p3_merged, addlabels = TRUE, ylim = c(0, 40))
fviz_pca_var(pca.p2_p3_merged, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

pdac2 <- read.csv("p2intersect.csv")
pca.pdac2 <- PCA(pdac2[,-1], scale.unit = TRUE, graph = T)
fviz_eig(pca.pdac2, addlabels = TRUE, ylim = c(0, 40))
fviz_pca_var(pca.pdac2, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)
pca.pdac2 <- PCA(t(pdac2[,-1]), scale.unit = TRUE, graph = T)
p2 <- fviz_pca_ind(pca.pdac2, col.ind = "cos2", 
                   gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                   repel = TRUE)
ggpar(p2,
      title = "PCA_pdacCells_set2",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
###
pdac3 <- read.csv("p3intersect.csv")
pca.pdac3 <- PCA(pdac3[,-1], scale.unit = TRUE, graph = T)
fviz_eig(pca.pdac3, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(pca.pdac3, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)
pca.pdac3 <- PCA(t(pdac3[,-1]), scale.unit = TRUE, graph = F)
p3 <- fviz_pca_ind(pca.pdac3, col.ind = "cos2", 
                   gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                   repel = TRUE)
ggpar(p3,
      title = "PCA_pdacCells_set3",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
##########################################################
# melanoma intersected proteins
melanoma3_cells <- read_excel("melanoma3.xlsx")
melanoma1_cells <- read_excel("melanoma1.xlsx")
melanoma_list3 <- melanoma3_cells$proteins
melanoma_list1 <- melanoma1_cells$proteins
venn.diagram(list("melanoma3"=melanoma_list3, "melanoma1"=melanoma_list1),
             fill = c("cyan","red"),
             cex= 0.5, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.3, filename = "melanoma.png",
             print.mode =c("raw", "percent"))
# melanoma overlapped protein names
melanomaoverlap <- intersect(melanoma_list3,melanoma_list1)
write.table(melanomaoverlap, file = "melanoma_intersected.txt", sep = "")
mela1_unique <- setdiff(melanoma_list1, melanoma_list3)
mela3_unique <- setdiff(melanoma_list3, melanoma_list1)
write.table(mela1_unique, file = "mela1_unique.txt", sep = "")
write.table(mela3_unique, file = "mela3_unique.txt", sep = "")
melanoma <- read.table(file = "melanoma_intersected.txt")
melanoma_prots <- write.table(melanoma$x, file = "melanoma_proteins.txt",
                              quote=F, row.names=F, col.names = F)
mela1 <- read.table(file = "mela1_unique.txt")
mela3 <- read.table(file = "mela3_unique.txt")
melanoma_uni1 <- write.table(mela1$x, file = "mela1_unique.txt",
                              quote=F, row.names=F, col.names = F)
melanoma_uni3 <- write.table(mela3$x, file = "mela3_unique.txt",
                              quote=F, row.names=F, col.names = F)
# melanoma_go intersected cutoff = 0.3
melanoma_genes <- read_excel("melanoma_idmap.xlsx")
melanoma_go <- enrichGO(gene = melanoma_genes$entrezgene, ont = "BP",
                        OrgDb ="org.Hs.eg.db",
                        readable=TRUE,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)

dotplot(melanoma_go)
mg <- pairwise_termsim(melanoma_go)
emapplot(mg)
mgcenter <- setReadable(melanoma_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(mgcenter)
cnetplot(mgcenter,  foldChange=melanoma_go, circular = TRUE, colorEdge = TRUE)
# melanoma set 1 cutoff = 0.2
melanoma1_genes <- read_excel("mela1uni_idmap.xlsx")
melanoma1_go <- enrichGO(gene = melanoma1_genes$entrezgene, ont = "ALL",
                        OrgDb ="org.Hs.eg.db",
                        readable=TRUE,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(melanoma1_go)
# melanoma set 3 cutoff = 0.2
melanoma3_genes <- read_excel("mela3uni_idmap.xlsx")
melanoma3_go <- enrichGO(gene = melanoma3_genes$entrezgene, ont = "BP",
                         OrgDb ="org.Hs.eg.db",
                         readable=TRUE,
                         pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(melanoma3_go)
# melanoma KEGG cutoff = 0.2
melanoma_kegg <- enrichKEGG(gene = melanoma_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(melanoma_kegg)
mk <- pairwise_termsim(melanoma_kegg)
emapplot(mk)
mkcenter <- setReadable(melanoma_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(mkcenter)
cnetplot(mkcenter,  foldChange=melanoma_kegg, circular = TRUE, colorEdge = TRUE)
# melanoma KEGG set1 cutoff = 0.2
melanoma1_kegg <- enrichKEGG(gene = melanoma1_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(melanoma1_kegg)
# melanoma KEGG set3 cutoff = 0.2
melanoma3_kegg <- enrichKEGG(gene = melanoma3_genes$entrezgene,
                             organism = 'hsa',
                             pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(melanoma3_kegg)
##########################
# PCA for melanoma
# melanoma intersected proteins
melanoma3_cells <- read_excel("melanoma3.xlsx")
melanoma1_cells <- read_excel("melanoma1.xlsx")
melanoma_list3 <- melanoma3_cells$proteins
melanoma_list1 <- melanoma1_cells$proteins
melanomaoverlap <- intersect(melanoma_list3,melanoma_list1)
melanoma_df1 <- filter(melanoma1_cells, proteins %in% melanomaoverlap)
melanoma3_df3 <- filter(melanoma3_cells, proteins %in% melanomaoverlap)
write.csv(melanoma_df1, "m1intersect.csv", row.names = F)
write.csv(melanoma3_df3, "m3intersect.csv", row.names = F)

# PCA Melanoma merged sets
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


pca.m1_m3_merged <- PCA(t(m1_m3_merged[,-1]), scale.unit = TRUE, graph = F)
fviz_eig(pca.m1_m3_merged, addlabels = TRUE, ylim = c(0, 40))
fviz_pca_ind(pca.m1_m3_merged, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)


melanoma1 <- read.csv("m1intersect.csv")
melanoma3 <- read.csv("m3intersect.csv")
pca.melanoma3 <- PCA(melanoma3[,-1], scale.unit = TRUE, graph = T)
fviz_eig(pca.melanoma3, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(pca.melanoma3, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)
pca.melanoma3 <- PCA(t(melanoma3[,-1]), scale.unit = TRUE, graph = F)
m3 <- fviz_pca_ind(pca.melanoma3, col.ind = "cos2", 
                   gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                   repel = TRUE)
ggpar(m3,
      title = "PCA_MelanomaCells_set3",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
##############
#####################
# Monocytes (U937) intersected proteins
U1_cells <-read_excel("U937_1.xlsx")
U3_cells <-read_excel("U937_3.xlsx")
U4_cells <-read_excel("U937_4.xlsx")
U5_cells <-read_excel("U937_5.xlsx")
U6_cells <-read_excel("U937_6.xlsx")
U937_list1 <- U1_cells$proteins
U937_list3 <- U3_cells$proteins
U937_list4 <- U4_cells$proteins
U937_list5 <- U5_cells$proteins
U937_list6 <- U6_cells$proteins
venn.diagram(list("U1"=U937_list1, "U3"=U937_list3, "U4"=U937_list4, "U5"=U937_list5, "U6"=U937_list6),
             fill = c("cyan","red","green","blue","yellow"),
             cex= 0.5, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.3, filename = "U937.png",
             print.mode =c("raw", "percent"))
# Monocytes (U937) overlapped protein names
monocytesoverlap <- reduce(U937_list1, U937_list3, U937_list4, U937_list5, U937_list6)
View(monocytesoverlap)
monocytesoverlap <- intersect(U937_list1, intersect(U937_list3,U937_list4))
monooverlap <- intersect(monocytesoverlap, intersect(U937_list5,U937_list6))
write.table(monooverlap, file = "U937_intersected.txt", sep = "")
monocytes <- read.table(file = "U937_intersected.txt")
monocytes_prots <- write.table(monocytes$x, file = "monocytes_proteins.txt",
                              quote=F, row.names=F, col.names = F)
# monocytes_go intersected cutoff = 0.3
monocytes_genes <- read_excel("monocytes_idmap.xlsx")
monocytes_go <- enrichGO(gene = monocytes_genes$entrezgene, ont = "BP",
                        OrgDb ="org.Hs.eg.db",
                        readable=TRUE,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(monocytes_go)
# All the plots in one GO_terms
h <- dotplot(hela_go) + ggtitle("Hela")
p <- dotplot(pdacGO) + ggtitle("PDAC")
m <- dotplot(melanoma_go) + ggtitle("Melanoma")
u <- dotplot(monocytes_go) + ggtitle("Monocytes")
plot_list(h, p, m, u, ncol=4)
ug <- pairwise_termsim(monocytes_go)
emapplot(ug)
ugcenter <- setReadable(monocytes_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(ugcenter)
cnetplot(ugcenter,  foldChange=monocytes_go, circular = TRUE, colorEdge = TRUE)
# monocytes KEGG cutoff = 0.2
monocytes_kegg <- enrichKEGG(gene = monocytes_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, qvalueCutoff = 0.3)
dotplot(monocytes_kegg)
# All the plots in one KEGG
hk <- dotplot(hela_kegg) + ggtitle("Hela")
pk <- dotplot(pdac_kegg) + ggtitle("PDAC")
mk <- dotplot(melanoma_kegg) + ggtitle("Melanoma")
uk <- dotplot(monocytes_kegg) + ggtitle("Monocytes")
plot_list(hk, pk, mk, uk, ncol=4)
uk <- pairwise_termsim(monocytes_kegg)
emapplot(uk)
ukcenter <- setReadable(monocytes_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(ukcenter)
cnetplot(ukcenter,  foldChange=monocytes_kegg, circular = TRUE, colorEdge = TRUE)
#####################################
# PCA Melanoma merged sets
mono_df1 <- filter(U1_cells, proteins %in% monooverlap)
write.csv(mono_df1, "mono1intersect.csv", row.names = F)
monocytes1 <- read.csv("mono1intersect.csv")
mono_df3 <- filter(U3_cells, proteins %in% monooverlap)
write.csv(mono_df3, "mono3intersect.csv", row.names = F)
monocytes3 <- read.csv("mono3intersect.csv")
mono_df4 <- filter(U4_cells, proteins %in% monooverlap)
write.csv(mono_df4, "mono4intersect.csv", row.names = F)
monocytes4 <- read.csv("mono4intersect.csv")
mono_df5 <- filter(U5_cells, proteins %in% monooverlap)
write.csv(mono_df5, "mono5intersect.csv", row.names = F)
monocytes5 <- read.csv("mono5intersect.csv")
mono_df6 <- filter(U6_cells, proteins %in% monooverlap)
write.csv(mono_df6, "mono6intersect.csv", row.names = F)
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
#####################################
derek <- as.data.frame(`GSE47552_Multiple myeloma`, row.names = NULL, check.rows = FALSE, check.names = TRUE)
View(derek)

derek <- as.data.frame(`GSE47552_Multiple myeloma`[[1
                             ]])
View(derek)
class(`GSE47552_Multiple myeloma`)
derek2 <- as.data.frame(`GSE47552_Multiple myeloma`[[2
]])
View(derek2)
df <- as.data.frame(hela_go)
View(df)
df_kegg <- as.data.frame(monocytes_kegg)
View(df_kegg)
df_mono <- as.data.frame(monocytes_go)
View(df_mono)
merged_go <- merge(df,df_mono, by="ID")
View(merged_go)
go <- df$Description 
go_mono <- df_mono$Description
go_mela <- as.data.frame(melanoma_go)$Description
go_pdac <- pdacGO$Description
View(as.data.frame(go))
View(as.data.frame(go_mono))
View(as.data.frame(go_mela))
View(as.data.frame(go_pdac))
go_overlap <- intersect(go, go_mono)
mela_pdac <- intersect(go_mela, go_pdac)
all_overlap <- intersect(go_overlap, mela_pdac)
go_uni <- setdiff(go, all_overlap)
mono_uni <- setdiff(go_mono, all_overlap)
melauni <- setdiff(go_mela, all_overlap)
pdacuni <- setdiff(go_pdac, all_overlap)
View(as.data.frame(all_overlap))
test_df <- data.frame(Celltype = c("Helacells", "Monocytes(U937)", "Melanoma", "PDAC"))
View(test_df)
test_df[c(all_overlap)] <- 1
View(test_df)
test_df[c(go_uni)] <- c(1, 0, 0, 0)
test_df[c(mono_uni)] <- c(0, 1, 0, 0)
test_df[c(melauni)] <- c(0, 0, 1, 0)
test_df[c(pdacuni)] <- c(0, 0, 0, 1)
write.table(test_df, file = "test.csv")
read.csv("test.csv")
celltype = test_df$Celltype
test_df <- as.data.frame(t(test_df))
colnames(test_df) <- celltype
View(test_df)
pca.test_df <- PCA(test_df, scale.unit = TRUE, graph = F)
fviz_eig(pca.test_df, addlabels = TRUE, ylim = c(0, 80))
allcells <- fviz_pca_var(pca.test_df, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

nb.cols <- 4
mycolors <- colorRampPalette(brewer.pal(4, "Set1"))(nb.cols)
all_GO <- fviz_pca_ind(pca.test_df , col.ind = test_df$Celltype, addEllipses = T,
                          palette = mycolors)
ggpar(allcells,
      title = "PCA_Celltypes_GO",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())
########$#############$###############$##################$############$##################
?ggpar
View(as.data.frame(hela_go))
hexcols = cancer_colors[match(strsplit2(bi$name, ' ')[,1], word(names(cancer_colors),1))]
hexcols[is.na(hexcols)] = '#000000'
####################################
hkegg <- hela_kegg$Description 
kegg_mono <- monocytes_kegg$Description
kegg_mela <- melanoma_kegg$Description
kegg_pdac <- pdac_kegg$Description
View(as.data.frame(hkegg))
View(as.data.frame(kegg_mono))
View(as.data.frame(kegg_mela))
View(as.data.frame(kegg_pdac))
kegg_overlap <- intersect(hkegg, kegg_mono)
mela_pdac_kegg <- intersect(kegg_mela, kegg_pdac)
allkegg_overlap <- intersect(kegg_overlap, mela_pdac_kegg)
hkegg_uni <- setdiff(go, allkegg_overlap)
keggmono_uni <- setdiff(kegg_mono, all_overlap)
keggmela_uni <- setdiff(kegg_mela, all_overlap)
keggpdac_uni <- setdiff(kegg_pdac, all_overlap)
testkegg_df <- data.frame(Celltype = c("Hela cells", "Monocytes (U937)", "Melanoma", "PDAC"))
View(testkegg_df)
testkegg_df[c(allkegg_overlap)] <- 1
View(testkegg_df)
testkegg_df[c(hkegg_uni)] <- c(1, 0, 0, 0)
testkegg_df[c(keggmono_uni)] <- c(0, 1, 0, 0)
testkegg_df[c(keggmela_uni)] <- c(0, 0, 1, 0)
testkegg_df[c(keggpdac_uni)] <- c(0, 0, 0, 1)
write.table(testkegg_df, file = "testkegg.csv")
read.csv("testkegg.csv")
celltype = testkegg_df$Celltype
testkegg_df <- as.data.frame(t(testkegg_df[,-1]))
colnames(testkegg_df) <- celltype
pca.testkegg_df <- PCA(testkegg_df, scale.unit = TRUE, graph = F)
fviz_eig(pca.testkegg_df, addlabels = TRUE, ylim = c(0, 80))
allkeggcells <- fviz_pca_var(pca.testkegg_df, col.var = "cos2",
                         gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
                         repel = TRUE)
ggpar(allkeggcells,
      title = "PCA_Celltype_KEGG",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())


## ---- Correlation plot --------------------



corrplot(All_Plot_Data$Cancer_clustering_similarity,
         p.mat = All_Plot_Data$Cancer_clustering_similarity,
         sig.level = -1,
         insig = "p-value",
         order = "hclust",
         method = "circle", type = "upper",
         pch.cex = .9,
         pch.col = "white", tl.col = "black", tl.srt = 45)

dev.print(pdf, 'Cancer_Type_Clustering.pdf', width = 10, height = 10)
