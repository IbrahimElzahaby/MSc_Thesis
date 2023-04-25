library(readxl)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(pathfindR)
BiocManager::install("clusterProfiler", version = "3.16")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
# Hela_GOterms
hela_genes <- read_excel("hela_idmap.xlsx")
hela_go <- enrichGO(gene = hela_genes$entrezgene, ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    readable=TRUE,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.3)
hela_go_dot <- dotplot(hela_go)
hg <- pairwise_termsim(hela_go)
hela_emap <- emapplot(hg, showCategory=10)
hgcenter <- setReadable(hela_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
helanet <- cnetplot(hgcenter,  foldChange=hela_go, circular = F, colorEdge = TRUE, showCategory=5)
helacnet <- cnetplot(hela_go,colorEdge = TRUE,node_lable = "all",circular = TRUE,showCategory = 3)
cnetplot(hgcenter, color_gene='steelblue', showCategory=3,
         foldChange=hela_genes$symbol, node_label="all", circular = TRUE, colorEdge = TRUE)
# Hela_KEGG pathways
hela_kegg <- enrichKEGG(gene = hela_genes$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3, keyType = "kegg")
hela_kegg_dot <- dotplot(hela_kegg)
he <- pairwise_termsim(hela_kegg)
he_map <- emapplot(he, showCategory=10)
hecenter <- setReadable(hela_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
hela_kegg_cnet <- cnetplot(hecenter,  foldChange=hela_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
#####
# PDAC set2 GO terms
pdac <- read_excel("idmap.xlsx")
pdacGO <- enrichGO(gene = pdac$entrezgene, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.3)
pdac_go_dot <- dotplot(pdacGO)
ridgeplot(pdacGO)
pd <- pairwise_termsim(pdacGO)
pdac_emap <- emapplot(pd, showCategory=10)
pdcenter <- setReadable(pd, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
pdacnet <- cnetplot(pdcenter,  foldChange=pdacGO, circular = F, colorEdge = TRUE, showCategory=5)
pdaccnet <- cnetplot(pdcenter,  foldChange=pdacGO, circular = TRUE, colorEdge = TRUE, showCategory=3)
# PDAC KEGG pathways set 2
pdac_kegg <- enrichKEGG(gene = pdac$entrezgene,
                        organism = 'hsa',
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)
pdac_kegg_dot <- dotplot(pdac_kegg)
pdkegg <- pairwise_termsim(pdac_kegg)
pd_map <- emapplot(pdkegg, showCategory=10)
pkcenter <- setReadable(pdac_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
pdac_kegg_cnet <- cnetplot(pkcenter,  foldChange=pdac_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
#####
# melanoma_go intersected cutoff = 0.3
melanoma_genes <- read_excel("melanoma_idmap.xlsx")
melanoma_go <- enrichGO(gene = melanoma_genes$entrezgene, ont = "BP",
                        OrgDb ="org.Hs.eg.db",
                        readable=TRUE,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)

melanoma_go_dot <- dotplot(melanoma_go)
mg <- pairwise_termsim(melanoma_go)
mela_emap <- emapplot(mg, showCategory=10)
mgcenter <- setReadable(melanoma_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
melanet <- cnetplot(mgcenter,  foldChange=melanoma_go, circular = F, colorEdge = TRUE, showCategory=5)

# melanoma KEGG cutoff = 0.2
melanoma_kegg <- enrichKEGG(gene = melanoma_genes$entrezgene,
                            organism = 'hsa',
                            pvalueCutoff = 0.05, qvalueCutoff = 0.3)
melanoma_kegg_dot <- dotplot(melanoma_kegg)
mk <- pairwise_termsim(melanoma_kegg)
mela_map <- emapplot(mk, showCategory=10)
mkcenter <- setReadable(melanoma_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
mela_kegg_cnet <- cnetplot(mkcenter,  foldChange=melanoma_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
#####
# monocytes_go intersected cutoff = 0.3
monocytes_genes <- read_excel("monocytes_idmap.xlsx")
monocytes_go <- enrichGO(gene = monocytes_genes$entrezgene, ont = "BP",
                         OrgDb ="org.Hs.eg.db",
                         readable=TRUE,
                         pvalueCutoff = 0.05, qvalueCutoff = 0.3)
monocytes_go_dot <- dotplot(monocytes_go)
mg <- pairwise_termsim(monocytes_go)
mono_emap <- emapplot(mg, showCategory=10)
mkcenter <- setReadable(monocytes_go, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
mononet <- cnetplot(mkcenter,  foldChange=monocytes_go, circular = F, colorEdge = TRUE, showCategory=5)
monocnet <- cnetplot(mkcenter,  foldChange=monocytes_go, circular = TRUE, colorEdge = TRUE, showCategory=3)
# monocytes KEGG cutoff = 0.2
monocytes_kegg <- enrichKEGG(gene = monocytes_genes$entrezgene,
                             organism = 'hsa',
                             pvalueCutoff = 0.05, qvalueCutoff = 0.3)
monocytes_kegg_dot <- dotplot(monocytes_kegg)
mok <- pairwise_termsim(monocytes_kegg)
mok_emap <- emapplot(mok, showCategory=10)
mokcenter <- setReadable(monocytes_kegg, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
mono_kegg_cnet <- cnetplot(mokcenter,  foldChange=monocytes_kegg, circular = F, colorEdge = TRUE, node_lable="all", showCategory = 5)
#################################
# Draw the GO diagrams
res <- 300
w <- 16
h <- 12
png("GO.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes_go_dot, hela_go_dot, melanoma_go_dot, pdac_go_dot, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()
#####
# Draw the E-map GO diagrams
res <- 300
w <- 16
h <- 12
png("emap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mono_emap, hela_emap, mela_emap, pdac_emap, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()
#####
# Draw the Cnet GO diagrams
res <- 300
w <- 16
h <- 14
png("Cnet.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mononet, helanet, melanet, pdacnet, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()
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











