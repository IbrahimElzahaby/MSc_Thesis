# import libraries
library(readxl)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(corrplot)
library(visNetwork)
library(igraph)
library(corrr)
library(GGally)
library(qgraph)
library(ggraph)
library(gridExtra)
library(anvis)
library(pals)
cyto_message <- testthat::capture_condition(RCy3::cytoscapePing())$message
cyto_active <- cyto_message == "You are connected to Cytoscape!\n"
cyto_active
layout_func <- function(graph) igraph::layout_on_grid(graph, width = 4, height = 9)
##########
h4 <- read.delim2("h4intersect.txt", sep = "\t", stringsAsFactors = F)
hela_cor <- as.data.frame(h4)
hela_cor <- t(hela_cor)
colnames(hela_cor) <- hela_cor[1,]
hela_cor <- hela_cor[-1,]
hj <- type.convert(hela_cor)
hela_cor <- as.numeric(hela_cor)
cor_h4 <- cor(hj, method = "spearman")
cor_h4 <- as.data.frame(cor_h4)
cor_h4[abs(cor_h4) < 0.75] <- 0
diag(cor_h4) <- 0
cor_net <- cor_h4[,colSums(abs(cor_h4), na.rm = T) > 0]
cor_net <- cor_net[rowSums(abs(cor_net), na.rm = T) > 0,]
custom_cols <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2")
rev_pal <- function (n) rev(pals::ocean.oxy(n))
pal_gwo <- colorRampPalette(colors = c("darkgreen", "white", "darkorange"))
net0_h4 <- adjToNetwork(cor_net, directed = F, self_loops = F,
            node_attrs = "none", group_colors = custom_cols,
            edge_color_func = rev_pal,
            edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
            group_vec = NULL, size_type = "cytoscape",
            width_type = "partcor",output_as = 'list')
A <- anvis(net0_h4, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T
           )
A_h4 <- anvis(net0_h4,
      igr_plot_opts = list(layout = layout_func, 
                           vertex.label.cex = 0.6, 
                           vertex.label.color = "black", 
                           vertex.label.family = "sans"),
      vis_radial_labs = FALSE, cyto_close_session = F)
######
h5 <- read.delim2("h5intersect.txt", sep = "\t", stringsAsFactors = F)
hela5_cor <- as.data.frame(h5)
hela5_cor <- t(hela5_cor)
colnames(hela5_cor) <- hela5_cor[1,]
hela5_cor <- hela5_cor[-1,]
hj5 <- type.convert(hela5_cor)
hela5_cor <- as.numeric(hela5_cor)
cor_h5 <- cor(hj5, method = "spearman")
cor_h5 <- as.data.frame(cor_h5)
cor_h5[abs(cor_h5) < 0.87] <- 0
diag(cor_h5) <- 0
cor_net5 <- cor_h5[,colSums(abs(cor_h5), na.rm = T) > 0]
cor_net5 <- cor_net5[rowSums(abs(cor_net5), na.rm = T) > 0,]
net0_h5 <- adjToNetwork(cor_net5, directed = F, self_loops = F,
                        node_attrs = "none", group_colors = custom_cols,
                        edge_color_func = rev_pal,
                        edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                        group_vec = NULL, size_type = "cytoscape",
                        width_type = "partcor",output_as = 'list')
A <- anvis(net0_h5, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
A_h4 <- anvis(net0_h5,
              igr_plot_opts = list(layout = layout_func, 
                                   vertex.label.cex = 0.6, 
                                   vertex.label.color = "black", 
                                   vertex.label.family = "sans"),
              vis_radial_labs = FALSE, cyto_close_session = F)
###
net0_h5 <- adjToNetwork(adj_mats = cor_net5, directed = FALSE, self_loops = F, 
                     node_attrs = "all",
                     edge_attrs = "all",
                     group_vec = NULL, size_type = "igraph", width_type = "percentile")
B <- anvis(net0_h5, directed = FALSE, output_type = "igraph")
##########
grid_hela <- c("Hela4", "Hela5")
G_H <- anvis(net0_h4,net0_h5, igr_grid = F,
             igr_grid_names = grid_hela,
             igr_par_opts = list(mar=c(4,4,4,4)))
res <- 300
w <- 14
h <- 10
png("Hela_network.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(A, B, ncol = 2, labels = c('A', 'B'), label_size = 16)
dev.off()
###########################################################################################
mela1 <- read.csv("m1intersect.csv")
mela1_cor <- as.data.frame(mela1)
mela1_cor <- t(mela1_cor)
colnames(mela1_cor) <- mela1_cor[1,]
mela1_cor <- mela1_cor[-1,]
m1 <- type.convert(mela1_cor)
mela1_cor <- as.numeric(mela1_cor)
cor_m1 <- cor(m1, method = "spearman")
cor_m1 <- as.data.frame(cor_m1)
cor_m1[abs(cor_m1) < 0.45] <- 0
diag(cor_m1) <- 0
cor_net_m1 <- cor_m1[,colSums(abs(cor_m1), na.rm = T) > 0]
cor_net_m1 <- cor_net_m1[rowSums(abs(cor_net_m1), na.rm = T) > 0,]
net0_m1 <- adjToNetwork(cor_net_m1, directed = F, self_loops = F,
                        node_attrs = "none", group_colors = custom_cols,
                        edge_color_func = rev_pal,
                        edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                        group_vec = NULL, size_type = "cytoscape",
                        width_type = "partcor",output_as = 'list')
C <- anvis(net0_m1, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T
)
C_m1 <- anvis(net0_m1,
              igr_plot_opts = list(layout = layout_func, 
                                   vertex.label.cex = 0.6, 
                                   vertex.label.color = "black", 
                                   vertex.label.family = "sans"),
              vis_radial_labs = FALSE, cyto_close_session = F)


net0_m1 <- adjToNetwork(adj_mats = cor_net_m1, directed = FALSE, self_loops = F, 
                        node_attrs = "all",
                        edge_attrs = "all",
                        group_vec = NULL, size_type = "igraph", width_type = "percentile")
C <- anvis(net0_m1, directed = FALSE, output_type = "igraph")
######
mela3 <- read.csv("m3intersect.csv")
mela3_cor <- as.data.frame(mela3)
mela3_cor <- t(mela3_cor)
colnames(mela3_cor) <- mela3_cor[1,]
mela3_cor <- mela3_cor[-1,]
m3 <- type.convert(mela3_cor)
mela3_cor <- as.numeric(mela3_cor)
cor_m3 <- cor(m3, method = "spearman")
cor_m3 <- as.data.frame(cor_m3)
cor_m3[abs(cor_m3) < 0.6] <- 0
diag(cor_m3) <- 0
cor_net_m3 <- cor_m3[,colSums(abs(cor_m3), na.rm = F) > 0]
cor_net_m3 <- cor_net_m3[rowSums(abs(cor_net_m3), na.rm = F) > 0,]

#net0_m3 <- adjToNetwork(adj_mats = cor_net_m3, directed = FALSE, self_loops = F, 
 #                       node_attrs = "all",
  #                      edge_attrs = "all",
   #                     group_vec = NULL, size_type = "igraph", width_type = "percentile")

net0_m3 <- adjToNetwork(cor_net_m3, directed = F, self_loops = F,
                        node_attrs = "none", group_colors = custom_cols,
                        edge_color_func = rev_pal,
                        edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                        group_vec = NULL, size_type = "igraph",
                        width_type = "partcor",
                        output_as = 'list'
                        )
D <- anvis(net0_m3, 
           igr_plot_opts = list(main=title), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T,
)

A_h4 <- anvis(net0_h4,
              igr_plot_opts = list(layout = layout_func, 
                                   vertex.label.cex = 0.6, 
                                   vertex.label.color = "black", 
                                   vertex.label.family = "sans"),
              vis_radial_labs = FALSE, cyto_close_session = F)
net0_m3 <- adjToNetwork(adj_mats = cor_net_m3, directed = FALSE, self_loops = F, 
                        node_attrs = "all",
                        edge_attrs = "all",
                        group_vec = NULL, size_type = "igraph", width_type = "percentile")
C <- anvis(net0_m3, directed = FALSE, output_type = "igraph")
############################################################################################
pdac2 <- read.csv("p2intersect.csv")
pdac2_cor <- as.data.frame(pdac2)
colnames(pdac2_cor) <- pdac2_cor[1,]
pdac2_cor <- pdac2_cor[-1,]
pd2 <- type.convert(pdac2_cor)
pdac2_cor <- as.numeric(pdac2_cor)
cor_pd2 <- cor(pdac2_cor, method = "spearman")
cor_pd2 <- as.data.frame(cor_pd2)
cor_pd2[abs(cor_pd2) < 0.5] <- 0
diag(cor_pd2) <- 0
cor_net_m1 <- cor_m1[,colSums(abs(cor_m1), na.rm = T) > 0]
cor_net_m1 <- cor_net_m1[rowSums(abs(cor_net_m1), na.rm = T) > 0,]
net0_m1 <- adjToNetwork(cor_net_m1, directed = F, self_loops = F,
                        node_attrs = "none", group_colors = custom_cols,
                        edge_color_func = rev_pal,
                        edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                        group_vec = NULL, size_type = "cytoscape",
                        width_type = "partcor",output_as = 'list')
C <- anvis(net0_m1, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
######
pdac3 <- read.csv("p3intersect.csv") %>% as.data.frame() %>% t()
colnames(pdac3) <- pdac3[1,]
pdac3 <- pdac3[-1,] %>% as.matrix()
pdac3 <- type.convert(pdac3)
cor_pdac3 <- cor(pdac3, method = "spearman") %>% as.data.frame()
cor_pdac3[abs(cor_pdac3) < 0.9] <- 0
diag(cor_pdac3) <- 0
cor_net_pdac3 <- cor_pdac3[,colSums(abs(cor_pdac3), na.rm = T) > 0]
cor_net_pdac3 <- cor_net_pdac3[rowSums(abs(cor_net_pdac3), na.rm = T) > 0,]
cor_net_pdac3 <- as.data.frame(cor_net_pdac3)
net0_pdac3$
net0_pdac3 <- adjToNetwork(cor_net_pdac3, directed = F, self_loops = F,
                        node_attrs = "none", group_colors = custom_cols,
                        edge_color_func = rev_pal,
                        edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                        group_vec = NULL, size_type = "cytoscape",
                        width_type = "partcor",output_as = 'list')
C <- anvis(net0_pdac3, 
           igr_plot_opts = list(vertex.frame.color = c("black", "red", "blue")), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
######
monocytes1 <- read.csv("mono1intersect.csv")
mo1_cor <- as.data.frame(monocytes1)
mo1_cor <- t(mo1_cor)
colnames(mo1_cor) <- mo1_cor[1,]
mo1_cor <- mo1_cor[-1,]
mo1 <- type.convert(mo1_cor)
mo1_cor <- as.numeric(mo1_cor)
cor_mo1 <- cor(mo1, method = "spearman")
cor_mo1 <- as.data.frame(cor_mo1)
cor_mo1[abs(cor_mo1) < 0.8] <- 0
diag(cor_mo1) <- 0
cor_net_mo1 <- cor_mo1[,colSums(abs(cor_mo1), na.rm = T) > 0]
cor_net_mo1 <- cor_net_mo1[rowSums(abs(cor_net_mo1), na.rm = T) > 0,]
net0_mo1 <- adjToNetwork(cor_net_mo1, directed = F, self_loops = F,
                         node_attrs = "none", group_colors = custom_cols,
                         edge_color_func = rev_pal,
                         edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                         group_vec = NULL, size_type = "cytoscape",
                         width_type = "partcor",output_as = 'list')
C <- anvis(net0_mo1, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
######
monocytes3 <- read.csv("mono3intersect.csv")
mo3_cor <- as.data.frame(monocytes3)
mo3_cor <- t(mo3_cor)
colnames(mo3_cor) <- mo3_cor[1,]
mo3_cor <- mo3_cor[-1,]
mo3 <- type.convert(mo3_cor)
mo3_cor <- as.numeric(mo3_cor)
cor_mo3 <- cor(mo3, method = "spearman")
cor_mo3 <- as.data.frame(cor_mo3)
cor_mo3[abs(cor_mo3) < 0.9] <- 0
diag(cor_mo3) <- 0
cor_net_mo3 <- cor_mo3[,colSums(abs(cor_mo3), na.rm = T) > 0]
cor_net_mo3 <- cor_net_mo3[rowSums(abs(cor_net_mo3), na.rm = T) > 0,]
net0_mo3 <- adjToNetwork(cor_net_mo3, directed = F, self_loops = F,
                         node_attrs = "none", group_colors = custom_cols,
                         edge_color_func = rev_pal,
                         edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                         group_vec = NULL, size_type = "cytoscape",
                         width_type = "partcor",output_as = 'list')
C <- anvis(net0_mo3, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
C <- anvis(net0_mo3, directed = FALSE, output_type = "igraph")
######
monocytes4 <- read.csv("mono4intersect.csv")
mo4_cor <- as.data.frame(monocytes4)
mo4_cor <- t(mo4_cor)
colnames(mo4_cor) <- mo4_cor[1,]
mo4_cor <- mo4_cor[-1,]
mo4 <- type.convert(mo4_cor)
mo4_cor <- as.numeric(mo4_cor)
cor_mo4 <- cor(mo4, method = "spearman")
cor_mo4 <- as.data.frame(cor_mo4)
cor_mo4[abs(cor_mo4) < 0.8] <- 0
diag(cor_mo4) <- 0
cor_net_mo4 <- cor_mo4[,colSums(abs(cor_mo4), na.rm = T) > 0]
cor_net_mo4 <- cor_net_mo4[rowSums(abs(cor_net_mo4), na.rm = T) > 0,]
net0_mo4 <- adjToNetwork(cor_net_mo4, directed = F, self_loops = F,
                         node_attrs = "none", group_colors = custom_cols,
                         edge_color_func = rev_pal,
                         edge_attrs = c("width", "color"), colorblind = F,arrange_co = F,
                         group_vec = NULL, size_type = "cytoscape",
                         width_type = "partcor",output_as = 'list')
C <- anvis(net0_mo4, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)






