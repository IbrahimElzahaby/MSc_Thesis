
# Load required libraries
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(pathfindR)

# Read Hela gene IDs data from Excel file
hela_genes <- read_excel("hela_idmap.xlsx")

# Read PDAC gene IDs data from Excel file
pdac <- read_excel("idmap.xlsx")

# Read Melanoma gene IDs data from Excel file
melanoma_genes <- read_excel("melanoma_idmap.xlsx")

# Read Monocytes gene IDs data from Excel file
monocytes_genes <- read_excel("monocytes_idmap.xlsx")

################################################################################

# Perform KEGG pathway enrichment analysis for Hela cells
hela_kegg <- enrichKEGG(gene = hela_genes$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.3, 
                        keyType = "kegg")

# Generate a dotplot to visualize the enriched pathways for Hela cells
hela_kegg_dot <- dotplot(hela_kegg)

# Calculate pairwise similarity between enriched terms
he <- pairwise_termsim(hela_kegg)

# Generate an enriched map plot to visualize the connectivity between enriched terms
he_map <- emapplot(he, showCategory=10)

# Convert the gene IDs in the enriched pathways to readable format using the org.Hs.eg.db package
hecenter <- setReadable(hela_kegg, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate a cnetplot to visualize the enriched pathways for Hela cells
hela_kegg_cnet <- cnetplot(hecenter,  foldChange=hela_kegg, 
                           circular = F, colorEdge = TRUE, 
                           node_lable="all", showCategory = 5)

################################################################################

# Perform KEGG pathway enrichment analysis for PDAC cells
pdac_kegg <- enrichKEGG(gene = pdac$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.3)

# Generate a dotplot to visualize the enriched pathways for PDAC cells
pdac_kegg_dot <- dotplot(pdac_kegg)

# Calculate pairwise similarity between enriched terms
pdkegg <- pairwise_termsim(pdac_kegg)

# Generate an enriched map plot to visualize the connectivity between enriched terms
pd_map <- emapplot(pdkegg, showCategory=10)

# Convert the gene IDs in the enriched pathways to readable format using the org.Hs.eg.db package
pkcenter <- setReadable(pdac_kegg, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate a cnetplot to visualize the enriched pathways for PDAC cells
pdac_kegg_cnet <- cnetplot(pkcenter,  foldChange=pdac_kegg, 
                           circular = F, colorEdge = TRUE, 
                           node_lable="all", showCategory = 5)

################################################################################

# Perform KEGG pathway enrichment analysis for Melanoma cells
melanoma_kegg <- enrichKEGG(gene = melanoma_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 0.3)

# Generate a dotplot to visualize the enriched pathways for Melanoma cells
melanoma_kegg_dot <- dotplot(melanoma_kegg)

# Calculate pairwise similarity between enriched terms
mk <- pairwise_termsim(melanoma_kegg)

# Generate an enriched map plot to visualize the connectivity between enriched terms
mela_map <- emapplot(mk, showCategory=10)

# Convert the gene IDs in the enriched pathways to readable format using the org.Hs.eg.db package
mkcenter <- setReadable(melanoma_kegg, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate a cnetplot to visualize the enriched pathways for Melanoma cells
mela_kegg_cnet <- cnetplot(mkcenter,  foldChange=melanoma_kegg, 
                           circular = F, colorEdge = TRUE, 
                           node_lable="all", showCategory = 5)

################################################################################

# Perform KEGG pathway enrichment analysis for Monocytes cells
monocytes_kegg <- enrichKEGG(gene = monocytes_genes$entrezgene,
                             organism = 'hsa',
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 0.3)

# Generate a dotplot to visualize the enriched pathways for Monocytes cells
monocytes_kegg_dot <- dotplot(monocytes_kegg)

# Calculate pairwise similarity between enriched terms
mok <- pairwise_termsim(monocytes_kegg)

# Generate an enriched map plot to visualize the connectivity between enriched terms
mok_emap <- emapplot(mok, showCategory=10)

# Convert the gene IDs in the enriched pathways to readable format using the org.Hs.eg.db package
mokcenter <- setReadable(monocytes_kegg, OrgDb = "org.Hs.eg.db", 
                         keyType = "ENTREZID")

# Generate a cnetplot to visualize the enriched pathways for Monocytes cells
mono_kegg_cnet <- cnetplot(mokcenter,  foldChange=monocytes_kegg, 
                           circular = F, colorEdge = TRUE, 
                           node_lable="all", showCategory = 5)

################################################################################

# Save dotplots KEGG pathway analysis to png file
res <- 300
w <- 16
h <- 12
png("kegg_dot.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes_kegg_dot, hela_kegg_dot, 
                   melanoma_kegg_dot, pdac_kegg_dot, 
                   ncol = 2, labels = c('A', 'B', 'C', 'D'), 
                   label_size = 24)
dev.off()

#####

# Save E-Map plots KEGG pathway analysis to png file
res <- 300
w <- 16
h <- 12
png("kegg_E-Map.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mok_emap, he_map, mela_map, pd_map, ncol = 2, 
                   labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()

#####

# Save C-net plots of KEGG pathway analysis to png file
res <- 250
w <- 20
h <- 15
png("kegg_Cnet2.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mono_kegg_cnet, hela_kegg_cnet, mela_kegg_cnet, 
                   pdac_kegg_cnet, ncol = 2, labels = c('A', 'B', 'C', 'D'), 
                   label_size = 24)
dev.off()



