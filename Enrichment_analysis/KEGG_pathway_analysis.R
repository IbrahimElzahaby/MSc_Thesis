library(readxl)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(pathfindR)
BiocManager::install("clusterProfiler", version = "3.16")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
hela_genes <- read_excel("hela_idmap.xlsx")
pdac <- read_excel("idmap.xlsx")
melanoma_genes <- read_excel("melanoma_idmap.xlsx")
monocytes_genes <- read_excel("monocytes_idmap.xlsx")
# Hela_KEGG pathways
hela_kegg <- enrichKEGG(gene = hela_genes$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3, keyType = "kegg")
hela_kegg_dot <- dotplot(hela_kegg)
he <- pairwise_termsim(hela_kegg)
he_map <- emapplot(he, showCategory=10)
hecenter <- setReadable(hela_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
hela_kegg_cnet <- cnetplot(hecenter,  foldChange=hela_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
# PDAC KEGG pathways
pdac_kegg <- enrichKEGG(gene = pdac$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)
pdac_kegg_dot <- dotplot(pdac_kegg)
pdkegg <- pairwise_termsim(pdac_kegg)
pd_map <- emapplot(pdkegg, showCategory=10)
pkcenter <- setReadable(pdac_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
pdac_kegg_cnet <- cnetplot(pkcenter,  foldChange=pdac_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
# melanoma KEGG
melanoma_kegg <- enrichKEGG(gene = melanoma_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, qvalueCutoff = 0.3)
melanoma_kegg_dot <- dotplot(melanoma_kegg)
mk <- pairwise_termsim(melanoma_kegg)
mela_map <- emapplot(mk, showCategory=10)
mkcenter <- setReadable(melanoma_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
mela_kegg_cnet <- cnetplot(mkcenter,  foldChange=melanoma_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
# monocytes KEGG
monocytes_kegg <- enrichKEGG(gene = monocytes_genes$entrezgene,
                             organism = 'hsa',
                             pvalueCutoff = 0.05, qvalueCutoff = 0.3)
monocytes_kegg_dot <- dotplot(monocytes_kegg)
mok <- pairwise_termsim(monocytes_kegg)
mok_emap <- emapplot(mok, showCategory=10)
mokcenter <- setReadable(monocytes_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
mono_kegg_cnet <- cnetplot(mokcenter,  foldChange=monocytes_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
#####
# Draw the dot KEGG diagrams
res <- 300
w <- 16
h <- 12
png("kegg_dot.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes_kegg_dot, hela_kegg_dot, melanoma_kegg_dot, pdac_kegg_dot, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()
#####
# Draw the E-Map KEGG diagrams
res <- 300
w <- 16
h <- 12
png("kegg_E-Map.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mok_emap, he_map, mela_map, pd_map, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()
#####
# Draw the Cnet KEGG diagrams
res <- 250
w <- 20
h <- 15
png("kegg_Cnet2.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mono_kegg_cnet, hela_kegg_cnet, mela_kegg_cnet, pdac_kegg_cnet, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()



