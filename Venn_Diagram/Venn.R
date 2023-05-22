
# Load required libraries
library(readxl)
library(VennDiagram)

# Read Hela proteins from Excel files
Hela_Set1 <-read_excel("hc4.xlsx")
helalist1 <- Hela_Set1$protein

Hela_Set2 <-read_excel("hc5.xlsx")
helalist2 <- Hela_Set2$Protein

# Create Venn diagram for Hela proteins
hela <- venn.diagram(list("Hela1"=helalist1, "Hela2"=helalist2),
             fill = c("cyan","red"), height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))

# Obtain Hela overlapped protein names
helaoverlap <- intersect(helalist4,helalist5)

# Write Hela overlapped protein names to a text file
write.table(helaoverlap, file = "Hela_intersected.txt", sep = "",
            quote = F, col.names = F, row.names = F)

# Read PDAC proteins from Excel files
PDAC_Set1 <-read_excel("PDAC2.xlsx")
PDAC_list1 <- PDAC_Set1$proteins

PDAC_Set2 <-read_excel("PDAC3.xlsx")              
PDAC_list2 <- PDAC_Set2$prot

# Create Venn diagram for PDAC proteins
pdac <- venn.diagram(list("PDAC1"=PDAC_list1, "PDAC2"=PDAC_list2),
             fill = c("cyan","red"), height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))

# Obtain PDAC overlapped protein names
PDACoverlap <- intersect(PDAC_list2,PDAC_list3)

# Write PDAC overlapped protein names to a text file
write.table(PDACoverlap, file = "PDAC_intersected.txt", sep = "",
            quote = F, col.names = F, row.names = F)

# Read Monocytes (U937) proteins from Excel files
U1_cells <-read_excel("U937_1.xlsx")
U937_list1 <- U1_cells$proteins

U2_cells <-read_excel("U937_3.xlsx")
U937_list2 <- U2_cells$proteins

U3_cells <-read_excel("U937_4.xlsx")
U937_list3 <- U3_cells$proteins

U4_cells <-read_excel("U937_5.xlsx")
U937_list4 <- U4_cells$proteins

U5_cells <-read_excel("U937_6.xlsx")
U937_list5 <- U5_cells$proteins

# Create Venn diagram for Monocytes (U937) proteins
monocytes <- venn.diagram(list("U1"=U937_list1, "U2"=U937_list2, "U3"=U937_list3, "U4"=U937_list4, "U5"=U937_list5),
             fill = c("cyan","red","green","blue","yellow"),  height = 500, width = 500, imagetype = "png",
             cex= 1, col="Red", lty="blank", lwd=0.8, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))

# Obtain Monocytes (U937) overlapped protein names
monocytesoverlap <- intersect(U937_list1, intersect(U937_list3,U937_list4))
monooverlap <- intersect(monocytesoverlap, intersect(U937_list5,U937_list6))

# Write Monocytes (U937) overlapped protein names to a text file
write.table(monooverlap, file = "U937_intersected.txt", sep = "",
            quote = F, col.names = F, row.names = F)

# Read Melanoma proteins from Excel files
Melanoma_Set2 <- read_excel("melanoma3.xlsx")
melanoma_list2 <- Melanoma_Set2$proteins

Melanoma_Set1 <- read_excel("melanoma1.xlsx")
melanoma_list1 <- Melanoma_Set1$proteins

# Create Venn diagram for Melanoma proteins
melanoma <- venn.diagram(list("Mela2"=melanoma_list2, "Mela1"=melanoma_list1),
             fill = c("cyan","red"),main.cex = 2, height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))

# Obtain Melanoma overlapped protein names
melanomaoverlap <- intersect(melanoma_list3,melanoma_list1)

# Write Melanoma overlapped protein names to a text file
write.table(melanomaoverlap, file = "melanoma_intersected.txt", sep = "",
            quote = F, col.names = F, row.names = F)

# Save Venn diagrams to png file
res <- 300
w <- 8
h <- 8
png("Venn.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes, hela, melanoma, pdac, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 16)
dev.off()


