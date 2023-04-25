library(tidyverse)
library(palmerpenguins)
library(umap)
library(lattice)
theme_set(theme_bw(18))
###########################################################
# UMAP for hela cells
hela5 <- read.csv("h5intersect.csv")
Protein <- hela5$Protein
hela5 <- as.data.frame(t(hela5[,-1]))
colnames(hela5) <- Protein
Protein <- row.names(hela5)
hela5 <- cbind(Protein, hela5)
row.names(hela5) <- NULL
hela5_umap <- umap(hela5[,-1])
umap_df_hela5 <- data.frame(umap1=hela5_umap$layout[,1],
                            umap2=hela5_umap$layout[,2],
                            cell=hela5$Protein)
dist = dist(hela5_umap$layout)
helaclust <- hclust(dist, method = "complete")
helaclusters <- cutree(helaclust, k=2) 
umap_df_hela5$hcluster_label <- as.factor(helaclusters)
h5_umap <- ggplot(data= umap_df_hela5,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="Hela5 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Hela5_cluster1
hcluster_1 <- subset(umap_df_hela5,hcluster_label== '1')
hcluster_1$cell
row.names(hela5) <- hela5$Protein
clust_abun1 <- hela5[c(hcluster_1$cell), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(hela5), clust_abun1)
hele_cluster1 <- write_csv(clust_abun1[-1,], file = "hela_cluster1.csv", col_names = F)
hela1_cluster1 <- read.csv("hela_cluster1.csv")
# Hela5_cluster2
hcluster_2 <- subset(umap_df_hela5, hcluster_label=='2')
hcluster_2$cell
row.names(hela5) <- hela5$Protein
clust_abun2 <- hela5[c(hcluster_2$cell), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(hela5), clust_abun2)
hele_cluster2 <- write_csv(clust_abun2[-1,], file = "hela_cluster2.csv", col_names = F)
######
hela4 <- read.csv("h4intersect.csv")
Protein <- hela4$Protein
hela4 <- as.data.frame(t(hela4[,-1]))
colnames(hela4) <- Protein
Protein <- row.names(hela4)
hela4 <- cbind(Protein, hela4)
row.names(hela4) <- NULL
hela4_umap <- umap(hela4[,-1])
umap_df_hela4 <- data.frame(umap1=hela4_umap$layout[,1],
                            umap2=hela4_umap$layout[,2],
                            cell=hela4$Protein)
hela4_dist = dist(hela4_umap$layout)
helaclust <- hclust(hela4_dist, method = "complete")
helaclusters <- cutree(helaclust, k=4)
umap_df_hela4$hcluster_label <- as.factor(helaclusters)
h4_umap <- ggplot(data= umap_df_hela4,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="Hela4 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
#Hela4_cluster1
h4cluster_1 <- subset(umap_df_hela4, hcluster_label=='1')
h4cluster_1$cell
row.names(hela4) <- hela4$Protein
clust_abun1 <- hela4[c(h4cluster_1$cell), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(hela4), clust_abun1)
hela4_cluster1 <- write_csv(clust_abun1[-1,], file = "hela4_cluster1.csv", col_names = F)
# Hela4_cluster2
h4cluster_2 <- subset(umap_df_hela4, hcluster_label=='2')
h4cluster_2$cell
row.names(hela4) <- hela4$Protein
clust_abun2 <- hela4[c(h4cluster_2$cell), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(hela4), clust_abun2)
hela4_cluster2 <- write_csv(clust_abun2[-1,], file = "hela4_cluster2.csv", col_names = F)
# Hela4_cluster3
h4cluster_3 <- subset(umap_df_hela4, hcluster_label=='3')
h4cluster_3$cell
row.names(hela4) <- hela4$Protein
clust_abun3 <- hela4[c(h4cluster_3$cell), ]
clust_abun3 <- t(clust_abun3)
clust_abun3 <- as.data.frame(clust_abun3)
clust_abun3 <- cbind(colnames(hela4), clust_abun3)
hela4_cluster3 <- write_csv(clust_abun3[-1,], file = "hela4_cluster3.csv", col_names = F)
# Hela4_cluster4
h4cluster_4 <- subset(umap_df_hela4, hcluster_label=='4')
h4cluster_4$cell
row.names(hela4) <- hela4$Protein
clust_abun4 <- hela4[c(h4cluster_4$cell), ]
clust_abun4 <- t(clust_abun4)
clust_abun4 <- as.data.frame(clust_abun4)
clust_abun4 <- cbind(colnames(hela4), clust_abun4)
hela4_cluster4 <- write_csv(clust_abun4[-1,], file = "hela4_cluster4.csv", col_names = F)
######
# Draw Hela UMAP diagrams
res <- 300
w <- 9
h <- 6
png("Hela_umap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(h4_umap, h5_umap, ncol = 2, labels = c('A', 'B'), label_size = 16)
dev.off()
###############################################################
# UMAP for PDAC cells
pdac2 <- read.csv("p2intersect.csv")
pdac2_umap <- umap(pdac2[,-1])
umap_df_pdac2 <- data.frame(umap1=pdac2_umap$layout[,1],
                            umap2=pdac2_umap$layout[,2],
                            cell=pdac2$proteins)
pdac2_dist = dist(pdac2_umap$layout)
pdac2clust <- hclust(pdac2_dist, method = "complete")
pdac2clusters <- cutree(pdac2clust, k=3) 
umap_df_pdac2$hcluster_label <- as.factor(pdac2clusters)
pd2_umap <- ggplot(data= umap_df_pdac2,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="PDAC2 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# PDAC2_cluster1
pd2cluster_1 <- subset(umap_df_pdac2, hcluster_label=='1')
pd2cluster_1$cell
row.names(pdac2) <- pdac2$proteins
clust_abun1 <- pdac2[c(pd2cluster_1$cell), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(pdac2), clust_abun1)
pdac2_cluster1 <- write_csv(clust_abun1[-1,], file = "pdac2_cluster1.csv", col_names = F)
# PDAC2_cluster2
pd2cluster_2 <- subset(umap_df_pdac2, hcluster_label=='2')
pd2cluster_2$cell
row.names(pdac2) <- pdac2$proteins
clust_abun2 <- pdac2[c(pd2cluster_2$cell), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(pdac2), clust_abun2)
pdac2_cluster2 <- write_csv(clust_abun2[-1,], file = "pdac2_cluster2.csv", col_names = F)
# PDAC2_cluster3
pd2cluster_3 <- subset(umap_df_pdac2, hcluster_label=='3')
pd2cluster_3$cell
row.names(pdac2) <- pdac2$proteins
clust_abun3 <- pdac2[c(pd2cluster_3$cell), ]
clust_abun3 <- t(clust_abun3)
clust_abun3 <- as.data.frame(clust_abun3)
clust_abun3 <- cbind(colnames(pdac2), clust_abun3)
pdac2_cluster3 <- write_csv(clust_abun3[-1,], file = "pdac2_cluster3.csv", col_names = F)
######
pdac3 <- read.csv("p3intersect.csv")
proteins <- pdac3$proteins
pdac3 <- as.data.frame(t(pdac3[,-1]))
colnames(pdac3) <- proteins
proteins <- row.names(pdac3)
pdac3 <- cbind(proteins, pdac3)
row.names(pdac3) <- NULL
pdac3_umap <- umap(pdac3[,-1])
umap_df_pdac3 <- data.frame(umap1=pdac3_umap$layout[,1],
                            umap2=pdac3_umap$layout[,2],
                            cell=pdac3$proteins)
pdac3_dist = dist(pdac3_umap$layout)
pdac3clust <- hclust(pdac3_dist, method = "complete")
pdac3clusters <- cutree(pdac3clust, k=2) 
umap_df_pdac3$hcluster_label <- as.factor(pdac3clusters)
pd3_umap <- ggplot(data= umap_df_pdac3,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="PDAC3 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
pd3_umap
# PDAC3_cluster1
pd3cluster_1 <- subset(umap_df_pdac3, hcluster_label=='1')
pd3cluster_1$cell
row.names(pdac3) <- pdac3$proteins
clust_abun1 <- pdac3[c(pd3cluster_1$cell), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(pdac3), clust_abun1)
pdac3_cluster1 <- write_csv(clust_abun1[-1,], file = "pdac3_cluster1.csv", col_names = F)
# PDAC3_cluster2
pd3cluster_2 <- subset(umap_df_pdac3, hcluster_label=='2')
pd3cluster_2$cell
row.names(pdac3) <- pdac3$proteins
clust_abun2 <- pdac3[c(pd3cluster_2$cell), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(pdac3), clust_abun2)
pdac3_cluster2 <- write_csv(clust_abun2[-1,], file = "pdac3_cluster2.csv", col_names = F)
######
res <- 300
w <- 9
h <- 6
png("PDAC_umap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(pd2_umap, pd3_umap, ncol = 2, labels = c('A', 'B'), label_size = 16)
dev.off()
##############################################################
# UMAP for Melanoma cells
melanoma1 <- read.csv("m1intersect.csv")
proteins <- melanoma1$proteins
melanoma1 <- as.data.frame(t(melanoma1[,-1]))
colnames(melanoma1) <- proteins
proteins <- row.names(melanoma1)
melanoma1 <- cbind(proteins, melanoma1)
row.names(melanoma1) <- NULL
melanoma1_umap <- umap(melanoma1[,-1])
umap_df_melanoma1 <- data.frame(umap1=melanoma1_umap$layout[,1],
                            umap2=melanoma1_umap$layout[,2],
                            cell=melanoma1)
melanoma1_dist = dist(melanoma1_umap$layout)
melanoma1clust <- hclust(melanoma1_dist, method = "complete")
melanoma1clusters <- cutree(melanoma1clust, k=3) 
umap_df_melanoma1$hcluster_label <- as.factor(melanoma1clusters)
m1_umap <- ggplot(data= umap_df_melanoma1,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="Melanoma1 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Melanoma1_cluster1
me1cluster_1 <- subset(umap_df_melanoma1, hcluster_label=='1')
me1cluster_1$cell.proteins
row.names(melanoma1) <- melanoma1$proteins
clust_abun1 <- melanoma1[c(me1cluster_1$cell.proteins), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(melanoma1), clust_abun1)
melanoma1_cluster1 <- write_csv(clust_abun1[-1,], file = "melanoma1_cluster1.csv", col_names = F)
# Melanoma1_cluster2
me1cluster_2 <- subset(umap_df_melanoma1, hcluster_label=='2')
me1cluster_2$cell.proteins
row.names(melanoma1) <- melanoma1$proteins
clust_abun2 <- melanoma1[c(me1cluster_2$cell.proteins), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(melanoma1), clust_abun2)
melanoma1_cluster2 <- write_csv(clust_abun2[-1,], file = "melanoma1_cluster2.csv", col_names = F)
# Melanoma1_cluster3
me1cluster_3 <- subset(umap_df_melanoma1, hcluster_label=='3')
me1cluster_3$cell.proteins
row.names(melanoma1) <- melanoma1$proteins
clust_abun3 <- melanoma1[c(me1cluster_3$cell.proteins), ]
clust_abun3 <- t(clust_abun3)
clust_abun3 <- as.data.frame(clust_abun3)
clust_abun3 <- cbind(colnames(melanoma1), clust_abun3)
melanoma1_cluster3 <- write_csv(clust_abun3[-1,], file = "melanoma1_cluster3.csv", col_names = F)
######
melanoma3 <- read.csv("m3intersect.csv")
proteins <- melanoma3$proteins
melanoma3 <- as.data.frame(t(melanoma3[,-1]))
colnames(melanoma3) <- proteins
proteins <- row.names(melanoma3)
melanoma3 <- cbind(proteins, melanoma3)
row.names(melanoma3) <- NULL
melanoma3_umap <- umap(melanoma3[,-1])
umap_df_melanoma3 <- data.frame(umap1=melanoma3_umap$layout[,1],
                                umap2=melanoma3_umap$layout[,2],
                                cell=melanoma3)
melanoma3_dist = dist(melanoma3_umap$layout)
melanoma3clust <- hclust(melanoma3_dist, method = "complete")
melanoma3clusters <- cutree(melanoma3clust, k=3) 
umap_df_melanoma3$hcluster_label <- as.factor(melanoma3clusters)
m3_umap <- ggplot(data= umap_df_melanoma3,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="Melanoma3 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Melanoma3_cluster1
me3cluster_1 <- subset(umap_df_melanoma3, hcluster_label=='1')
me3cluster_1$cell.proteins
row.names(melanoma3) <- melanoma3$proteins
clust_abun1 <- melanoma3[c(me3cluster_1$cell.proteins), ]
clust_abun1 <- t(clust_abun1)
clust_abun1 <- as.data.frame(clust_abun1)
clust_abun1 <- cbind(colnames(melanoma3), clust_abun1)
melanoma3_cluster1 <- write_csv(clust_abun1[-1,], file = "melanoma3_cluster1.csv", col_names = F)
# Melanoma3_cluster2
me3cluster_2 <- subset(umap_df_melanoma3, hcluster_label=='2')
me3cluster_2$cell.proteins
row.names(melanoma3) <- melanoma3$proteins
clust_abun2 <- melanoma3[c(me3cluster_2$cell.proteins), ]
clust_abun2 <- t(clust_abun2)
clust_abun2 <- as.data.frame(clust_abun2)
clust_abun2 <- cbind(colnames(melanoma3), clust_abun2)
melanoma3_cluster2 <- write_csv(clust_abun2[-1,], file = "melanoma3_cluster2.csv", col_names = F)
# Melanoma3_cluster3
me3cluster_3 <- subset(umap_df_melanoma3, hcluster_label=='3')
me3cluster_3$cell.proteins
row.names(melanoma3) <- melanoma3$proteins
clust_abun3 <- melanoma3[c(me3cluster_3$cell.proteins), ]
clust_abun3 <- t(clust_abun3)
clust_abun3 <- as.data.frame(clust_abun3)
clust_abun3 <- cbind(colnames(melanoma3), clust_abun3)
melanoma3_cluster3 <- write_csv(clust_abun3[-1,], file = "melanoma3_cluster3.csv", col_names = F)
######
res <- 300
w <- 9
h <- 6
png("Melanoma_umap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(m1_umap, m3_umap, ncol = 2, labels = c('A', 'B'), label_size = 16)
dev.off()
##################################################################
# UMAP for Monocytes cells
monocytes1 <- read.csv("mono1intersect.csv")
proteins <- monocytes1$proteins
monocytes1 <- as.data.frame(t(monocytes1[,-1]))
colnames(monocytes1) <- proteins
proteins <- row.names(monocytes1)
monocytes1 <- cbind(proteins, monocytes1)
row.names(monocytes1) <- NULL
monocytes1_umap <- umap(monocytes1[,-1])
umap_df_monocytes1 <- data.frame(umap1=monocytes1_umap$layout[,1],
                                umap2=monocytes1_umap$layout[,2],
                                cell=monocytes1)
monocytes1_dist = dist(monocytes1_umap$layout)
monocytes1clust <- hclust(monocytes1_dist, method = "complete")
monocytes1clusters <- cutree(monocytes1clust, k=2) 
umap_df_monocytes1$hcluster_label <- as.factor(monocytes1clusters)
mo1_umap <- ggplot(data= umap_df_monocytes1,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="U1 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Monocytes1_cluster1
mo1cluster_1 <- subset(umap_df_monocytes1, hcluster_label=='1')
mo1cluster_1$cell.proteins
row.names(monocytes1) <- monocytes1$proteins
clustmo_abun1 <- monocytes1[c(mo1cluster_1$cell.proteins), ]
clustmo_abun1 <- t(clustmo_abun1)
clustmo_abun1 <- as.data.frame(clustmo_abun1)
clustmo_abun1 <- cbind(colnames(monocytes1), clustmo_abun1)
monocytes1_cluster1 <- write_csv(clustmo_abun1[-1,], file = "monocytes1_cluster1.csv", col_names = F)
# Monocytes1_cluster2
mo1cluster_2 <- subset(umap_df_monocytes1, hcluster_label=='2')
mo1cluster_2$cell.proteins
row.names(monocytes1) <- monocytes1$proteins
clustmo_abun2 <- monocytes1[c(mo1cluster_2$cell.proteins), ]
clustmo_abun2 <- t(clustmo_abun2)
clustmo_abun2 <- as.data.frame(clustmo_abun2)
clustmo_abun2 <- cbind(colnames(monocytes1), clustmo_abun2)
monocytes1_cluster2 <- write_csv(clustmo_abun2[-1,], file = "monocytes1_cluster2.csv", col_names = F)
######
monocytes3 <- read.csv("mono3intersect.csv")
proteins <- monocytes3$proteins
monocytes3 <- as.data.frame(t(monocytes3[,-1]))
colnames(monocytes3) <- proteins
proteins <- row.names(monocytes3)
monocytes3 <- cbind(proteins, monocytes3)
row.names(monocytes3) <- NULL
monocytes3_umap <- umap(monocytes3[,-1])
umap_df_monocytes3 <- data.frame(umap1=monocytes3_umap$layout[,1],
                                 umap2=monocytes3_umap$layout[,2],
                                 cell=monocytes3)
monocytes3_dist = dist(monocytes3_umap$layout)
monocytes3clust <- hclust(monocytes3_dist, method = "complete")
monocytes3clusters <- cutree(monocytes3clust, k=3) 
umap_df_monocytes3$hcluster_label <- as.factor(monocytes3clusters)
mo3_umap <- ggplot(data= umap_df_monocytes3,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="U3 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Monocytes3_cluster1
mo3cluster_1 <- subset(umap_df_monocytes3, hcluster_label=='1')
mo3cluster_1$cell.proteins
row.names(monocytes3) <- monocytes3$proteins
clustmo3_abun1 <- monocytes3[c(mo3cluster_1$cell.proteins), ]
clustmo3_abun1 <- t(clustmo3_abun1)
clustmo3_abun1 <- as.data.frame(clustmo3_abun1)
clustmo3_abun1 <- cbind(colnames(monocytes3), clustmo3_abun1)
monocytes3_cluster1 <- write_csv(clustmo3_abun1[-1,], file = "monocytes3_cluster1.csv", col_names = F)
# Monocytes3_cluster2
mo3cluster_2 <- subset(umap_df_monocytes3, hcluster_label=='2')
mo3cluster_2$cell.proteins
row.names(monocytes3) <- monocytes3$proteins
clustmo3_abun2 <- monocytes3[c(mo3cluster_2$cell.proteins), ]
clustmo3_abun2 <- t(clustmo3_abun2)
clustmo3_abun2 <- as.data.frame(clustmo3_abun2)
clustmo3_abun2 <- cbind(colnames(monocytes3), clustmo3_abun2)
monocytes3_cluster2 <- write_csv(clustmo3_abun2[-1,], file = "monocytes3_cluster2.csv", col_names = F)
# Monocytes3_cluster3
mo3cluster_3 <- subset(umap_df_monocytes3, hcluster_label=='3')
mo3cluster_3$cell.proteins
row.names(monocytes3) <- monocytes3$proteins
clustmo3_abun3 <- monocytes3[c(mo3cluster_3$cell.proteins), ]
clustmo3_abun3 <- t(clustmo3_abun3)
clustmo3_abun3 <- as.data.frame(clustmo3_abun3)
clustmo3_abun3 <- cbind(colnames(monocytes3), clustmo3_abun3)
monocytes3_cluster3 <- write_csv(clustmo3_abun3[-1,], file = "monocytes3_cluster3.csv", col_names = F)
######
monocytes4 <- read.csv("mono4intersect.csv")
proteins <- monocytes4$proteins
monocytes4 <- as.data.frame(t(monocytes4[,-1]))
colnames(monocytes4) <- proteins
proteins <- row.names(monocytes4)
monocytes4 <- cbind(proteins, monocytes4)
row.names(monocytes4) <- NULL
monocytes4_umap <- umap(monocytes4[,-1])
umap_df_monocytes4 <- data.frame(umap1=monocytes4_umap$layout[,1],
                                 umap2=monocytes4_umap$layout[,2],
                                 cell=monocytes4)
monocytes4_dist = dist(monocytes4_umap$layout)
monocytes4clust <- hclust(monocytes4_dist, method = "complete")
monocytes4clusters <- cutree(monocytes4clust, k=2) 
umap_df_monocytes4$hcluster_label <- as.factor(monocytes4clusters)
mo4_umap <- ggplot(data= umap_df_monocytes4,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="U4 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Monocytes4_cluster1
mo4cluster_1 <- subset(umap_df_monocytes4, hcluster_label=='1')
mo4cluster_1$cell.proteins
row.names(monocytes4) <- monocytes4$proteins
clustmo4_abun1 <- monocytes4[c(mo4cluster_1$cell.proteins), ]
clustmo4_abun1 <- t(clustmo4_abun1)
clustmo4_abun1 <- as.data.frame(clustmo4_abun1)
clustmo4_abun1 <- cbind(colnames(monocytes4), clustmo4_abun1)
monocytes4_cluster1 <- write_csv(clustmo4_abun1[-1,], file = "monocytes4_cluster1.csv", col_names = F)
# Monocytes4_cluster2
mo4cluster_2 <- subset(umap_df_monocytes4, hcluster_label=='2')
mo4cluster_2$cell.proteins
row.names(monocytes4) <- monocytes4$proteins
clustmo4_abun2 <- monocytes4[c(mo4cluster_2$cell.proteins), ]
clustmo4_abun2 <- t(clustmo4_abun2)
clustmo4_abun2 <- as.data.frame(clustmo4_abun2)
clustmo4_abun2 <- cbind(colnames(monocytes4), clustmo4_abun2)
monocytes4_cluster2 <- write_csv(clustmo4_abun2[-1,], file = "monocytes4_cluster2.csv", col_names = F)
######
monocytes5 <- read.csv("mono5intersect.csv")
proteins <- monocytes5$proteins
monocytes5 <- as.data.frame(t(monocytes5[,-1]))
colnames(monocytes5) <- proteins
proteins <- row.names(monocytes5)
monocytes5 <- cbind(proteins, monocytes5)
row.names(monocytes5) <- NULL
monocytes5_umap <- umap(monocytes5[,-1])
umap_df_monocytes5 <- data.frame(umap1=monocytes5_umap$layout[,1],
                                 umap2=monocytes5_umap$layout[,2],
                                 cell=monocytes5)
monocytes5_dist = dist(monocytes5_umap$layout)
monocytes5clust <- hclust(monocytes5_dist, method = "complete")
monocytes5clusters <- cutree(monocytes5clust, k=4) 
umap_df_monocytes5$hcluster_label <- as.factor(monocytes5clusters)
mo5_umap <- ggplot(data= umap_df_monocytes5,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="U5 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Monocytes5_cluster1
mo5cluster_1 <- subset(umap_df_monocytes5, hcluster_label=='1')
mo5cluster_1$cell.proteins
row.names(monocytes5) <- monocytes5$proteins
clustmo5_abun1 <- monocytes5[c(mo5cluster_1$cell.proteins), ]
clustmo5_abun1 <- t(clustmo5_abun1)
clustmo5_abun1 <- as.data.frame(clustmo5_abun1)
clustmo5_abun1 <- cbind(colnames(monocytes5), clustmo5_abun1)
monocytes5_cluster1 <- write_csv(clustmo5_abun1[-1,], file = "monocytes5_cluster1.csv", col_names = F)
# Monocytes5_cluster2
mo5cluster_2 <- subset(umap_df_monocytes5, hcluster_label=='2')
mo5cluster_2$cell.proteins
row.names(monocytes5) <- monocytes5$proteins
clustmo5_abun2 <- monocytes5[c(mo5cluster_2$cell.proteins), ]
clustmo5_abun2 <- t(clustmo5_abun2)
clustmo5_abun2 <- as.data.frame(clustmo5_abun2)
clustmo5_abun2 <- cbind(colnames(monocytes5), clustmo5_abun2)
monocytes5_cluster2 <- write_csv(clustmo5_abun2[-1,], file = "monocytes5_cluster2.csv", col_names = F)
# Monocytes5_cluster3
mo5cluster_3 <- subset(umap_df_monocytes5, hcluster_label=='3')
mo5cluster_3$cell.proteins
row.names(monocytes5) <- monocytes5$proteins
clustmo5_abun3 <- monocytes5[c(mo5cluster_3$cell.proteins), ]
clustmo5_abun3 <- t(clustmo5_abun3)
clustmo5_abun3 <- as.data.frame(clustmo5_abun3)
clustmo5_abun3 <- cbind(colnames(monocytes5), clustmo5_abun3)
monocytes5_cluster3 <- write_csv(clustmo5_abun3[-1,], file = "monocytes5_cluster3.csv", col_names = F)
# Monocytes5_cluster4
mo5cluster_4 <- subset(umap_df_monocytes5, hcluster_label=='4')
mo5cluster_4$cell.proteins
row.names(monocytes5) <- monocytes5$proteins
clustmo5_abun4 <- monocytes5[c(mo5cluster_4$cell.proteins), ]
clustmo5_abun4 <- t(clustmo5_abun4)
clustmo5_abun4 <- as.data.frame(clustmo5_abun4)
clustmo5_abun4 <- cbind(colnames(monocytes5), clustmo5_abun4)
monocytes5_cluster4 <- write_csv(clustmo5_abun4[-1,], file = "monocytes5_cluster4.csv", col_names = F)
######
monocytes6 <- read.csv("mono6intersect.csv")
proteins <- monocytes6$proteins
monocytes6 <- as.data.frame(t(monocytes6[,-1]))
colnames(monocytes6) <- proteins
proteins <- row.names(monocytes6)
monocytes6 <- cbind(proteins, monocytes6)
row.names(monocytes6) <- NULL
monocytes6_umap <- umap(monocytes6[,-1])
umap_df_monocytes6 <- data.frame(umap1=monocytes6_umap$layout[,1],
                                 umap2=monocytes6_umap$layout[,2],
                                 cell=monocytes6)
monocytes6_dist = dist(monocytes6_umap$layout)
monocytes6clust <- hclust(monocytes6_dist, method = "complete")
monocytes6clusters <- cutree(monocytes6clust, k=3) 
umap_df_monocytes6$hcluster_label <- as.factor(monocytes6clusters)
mo6_umap <- ggplot(data= umap_df_monocytes6,
       aes(x=umap1,y=umap2, color=hcluster_label))+
  geom_point(size = 3, alpha=0.5)+
  labs(x = "UMAP1", y = "UMAP2", colour="U6 subgroups")+
  scale_y_continuous()+scale_x_continuous()+
  theme(legend.position = "bottom")
# Monocytes6_cluster1
mo6cluster_1 <- subset(umap_df_monocytes6, hcluster_label=='1')
mo6cluster_1$cell.proteins
row.names(monocytes6) <- monocytes6$proteins
clustmo6_abun1 <- monocytes6[c(mo6cluster_1$cell.proteins), ]
clustmo6_abun1 <- t(clustmo6_abun1)
clustmo6_abun1 <- as.data.frame(clustmo6_abun1)
clustmo6_abun1 <- cbind(colnames(monocytes6), clustmo6_abun1)
monocytes6_cluster1 <- write_csv(clustmo6_abun1[-1,], file = "monocytes6_cluster1.csv", col_names = F)
# Monocytes6_cluster2
mo6cluster_2 <- subset(umap_df_monocytes6, hcluster_label=='2')
mo6cluster_2$cell.proteins
row.names(monocytes6) <- monocytes6$proteins
clustmo6_abun2 <- monocytes6[c(mo6cluster_2$cell.proteins), ]
clustmo6_abun2 <- t(clustmo6_abun2)
clustmo6_abun2 <- as.data.frame(clustmo6_abun2)
clustmo6_abun2 <- cbind(colnames(monocytes6), clustmo6_abun2)
monocytes6_cluster2 <- write_csv(clustmo6_abun2[-1,], file = "monocytes6_cluster2.csv", col_names = F)
# Monocytes6_cluster1
mo6cluster_3 <- subset(umap_df_monocytes6, hcluster_label=='3')
mo6cluster_3$cell.proteins
row.names(monocytes6) <- monocytes6$proteins
clustmo6_abun3 <- monocytes6[c(mo6cluster_3$cell.proteins), ]
clustmo6_abun3 <- t(clustmo6_abun3)
clustmo6_abun3 <- as.data.frame(clustmo6_abun3)
clustmo6_abun3 <- cbind(colnames(monocytes6), clustmo6_abun3)
monocytes6_cluster3 <- write_csv(clustmo6_abun3[-1,], file = "monocytes6_cluster3.csv", col_names = F)
######
res <- 300
w <- 18
h <- 6
png("Monocytes_umap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mo1_umap, mo3_umap, mo4_umap, mo5_umap, mo6_umap, ncol = 5, labels = c('A', 'B', 'C', 'D', 'E'), label_size = 16)
dev.off()














