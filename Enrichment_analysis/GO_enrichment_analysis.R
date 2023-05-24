
# Load required libraries
library(readxl)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Hs.eg.db)
library(pathfindR)

# Read Hela gene IDs from Excel file
hela_genes <- read_excel("hela_idmap.xlsx")

# Perform Gene Ontology (GO) enrichment analysis for Biological Processes (BP)
hela_go <- enrichGO(gene = hela_genes$entrezgene, ont = "BP",
                    OrgDb ="org.Hs.eg.db",
                    readable=TRUE,
                    pvalueCutoff = 0.05, qvalueCutoff = 0.3)

# Generate dotplot for enriched GO terms
hela_go_dot <- dotplot(hela_go)

# Calculate pairwise term similarities
hg <- pairwise_termsim(hela_go)

# Generate enrichment map plot
hela_emap <- emapplot(hg, showCategory=10)

# Set gene ID names using org.Hs.eg.db
hgcenter <- setReadable(hela_go, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate circular network plot
helanet <- cnetplot(hgcenter,  foldChange=hela_go, circular = F, 
                    colorEdge = TRUE, showCategory=5)

###############################################################

# Read PDAC gene IDs from Excel file
pdac <- read_excel("idmap.xlsx")

# Perform Gene Ontology (GO) enrichment analysis for Biological Processes (BP)
pdacGO <- enrichGO(gene = pdac$entrezgene, ont = "BP",
                   OrgDb ="org.Hs.eg.db",
                   readable=TRUE,
                   pvalueCutoff = 0.05, qvalueCutoff = 0.3)

# Generate dotplot for enriched GO terms
pdac_go_dot <- dotplot(pdacGO)

# Calculate pairwise term similarities
pd <- pairwise_termsim(pdacGO)

# Generate enrichment map plot
pdac_emap <- emapplot(pd, showCategory=10)

# Set gene ID names using org.Hs.eg.db
pdcenter <- setReadable(pd, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate circular network plot
pdacnet <- cnetplot(pdcenter,  foldChange=pdacGO, circular = F, 
                    colorEdge = TRUE, showCategory=5)

###############################################################

# Read Melanoma gene IDs from Excel file
melanoma_genes <- read_excel("melanoma_idmap.xlsx")

# Perform Gene Ontology (GO) enrichment analysis for Biological Processes (BP)
melanoma_go <- enrichGO(gene = melanoma_genes$entrezgene, ont = "BP",
                        OrgDb ="org.Hs.eg.db",
                        readable=TRUE,
                        pvalueCutoff = 0.05, qvalueCutoff = 0.3)

# Generate dotplot for enriched GO terms
melanoma_go_dot <- dotplot(melanoma_go)

# Calculate pairwise term similarities
mg <- pairwise_termsim(melanoma_go)

# Generate enrichment map plot
mela_emap <- emapplot(mg, showCategory=10)

# Set gene ID names using org.Hs.eg.db
mgcenter <- setReadable(melanoma_go, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate circular network plot
melanet <- cnetplot(mgcenter,  foldChange=melanoma_go, circular = F, 
                    colorEdge = TRUE, showCategory=5)

###############################################################

# Read Monocyte gene IDs from Excel file
monocytes_genes <- read_excel("monocytes_idmap.xlsx")

# Perform Gene Ontology (GO) enrichment analysis for Biological Processes (BP)
monocytes_go <- enrichGO(gene = monocytes_genes$entrezgene, ont = "BP",
                         OrgDb ="org.Hs.eg.db",
                         readable=TRUE,
                         pvalueCutoff = 0.05, qvalueCutoff = 0.3)

# Generate dotplot for enriched GO terms
monocytes_go_dot <- dotplot(monocytes_go)

# Calculate pairwise term similarities
mg <- pairwise_termsim(monocytes_go)

# Generate enrichment map plot
mono_emap <- emapplot(mg, showCategory=10)

# Set gene ID names using org.Hs.eg.db
mkcenter <- setReadable(monocytes_go, OrgDb = "org.Hs.eg.db", 
                        keyType = "ENTREZID")

# Generate circular network plot
mononet <- cnetplot(mkcenter,  foldChange=monocytes_go, circular = F, 
                    colorEdge = TRUE, showCategory=5)

###############################################################

# Save Gene Ontology (GO) dot plots to png file
res <- 300
w <- 16
h <- 12
png("GO.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes_go_dot, hela_go_dot, melanoma_go_dot, pdac_go_dot, 
                   ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()

#####

# Save the E-map plots to png file
res <- 300
w <- 16
h <- 12
png("emap.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mono_emap, hela_emap, mela_emap, pdac_emap, ncol = 2, 
                   labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()

#####

# Save the C-net plots to png file
res <- 300
w <- 16
h <- 14
png("Cnet.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(mononet, helanet, melanet, pdacnet, ncol = 2, 
                   labels = c('A', 'B', 'C', 'D'), label_size = 24)
dev.off()


