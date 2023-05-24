# import libraries
library(corrplot)
library(visNetwork)
library(igraph)
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
B <- anvis(net0_h5, 
           igr_plot_opts = list(vertex.frame.color = "black"), 
           directed = FALSE, output_type = "igraph",
           vis_edge_factor = 3, cyto_node_space = 2,
           igr_grid = T,  vis_radial_labs = T)
A_h5 <- anvis(net0_h5,
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


