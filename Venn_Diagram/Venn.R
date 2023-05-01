# import libraries
library(readxl)
library(VennDiagram)
library(ggplot2)
# Hela intersected proteins
Hela_Set1 <-read_excel("hc4.xlsx")
helalist1 <- Hela_Set1$protein
Hela_Set2 <-read_excel("hc5.xlsx")              
helalist2 <- Hela_Set2$Protein
hela <- venn.diagram(list("Hela1"=helalist1, "Hela2"=helalist2),
             fill = c("cyan","red"), height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))
# Hela overlapped protein names
helaoverlap <- intersect(helalist4,helalist5)
setdiff()
write.table(helaoverlap, file = "Hela_intersected.txt", sep = "")
# PDAC intersected proteins
PDAC_Set1 <-read_excel("PDAC2.xlsx")
PDAC_list1 <- PDAC_Set1$proteins
PDAC_Set2 <-read_excel("PDAC3.xlsx")              
PDAC_list2 <- PDAC_Set2$prot
pdac <- venn.diagram(list("PDAC1"=PDAC_list1, "PDAC2"=PDAC_list2),
             fill = c("cyan","red"), height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))
# PDAC overlapped protein names
PDACoverlap <- intersect(PDAC_list2,PDAC_list3)
write.table(PDACoverlap, file = "PDAC_intersected.txt", sep = "")
# Monocytes (U937) intersected proteins
U1_cells <-read_excel("U937_1.xlsx")
U2_cells <-read_excel("U937_3.xlsx")
U3_cells <-read_excel("U937_4.xlsx")
U4_cells <-read_excel("U937_5.xlsx")
U5_cells <-read_excel("U937_6.xlsx")
U937_list1 <- U1_cells$proteins
U937_list2 <- U2_cells$proteins
U937_list3 <- U3_cells$proteins
U937_list4 <- U4_cells$proteins
U937_list5 <- U5_cells$proteins
monocytes <- venn.diagram(list("U1"=U937_list1, "U2"=U937_list2, "U3"=U937_list3, "U4"=U937_list4, "U5"=U937_list5),
             fill = c("cyan","red","green","blue","yellow"),  height = 500, width = 500, imagetype = "png",
             cex= 1, col="Red", lty="blank", lwd=0.8, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))
# Monocytes (U937) overlapped protein names
monocytesoverlap <- intersect(U937_list1,U937_list3)
monocyteoverlap <- intersect(U937_list4,U937_list5)

write.table(monocytesoverlap, file = "U937_intersected.txt", sep = "")
# melanoma intersected proteins
Melanoma_Set2 <- read_excel("melanoma3.xlsx")
Melanoma_Set1 <- read_excel("melanoma1.xlsx")
melanoma_list2 <- Melanoma_Set2$proteins
melanoma_list1 <- Melanoma_Set1$proteins
melanoma <- venn.diagram(list("Mela2"=melanoma_list2, "Mela1"=melanoma_list1),
             fill = c("cyan","red"),main.cex = 2, height = 2500, width = 2500, imagetype = "png",
             cex= 1.2, col="Red", lty="blank", lwd=1, fontface="bold",
             fontfamily="sans", alpha=0.5, filename = NULL,
             print.mode =c("raw"))
# melanoma overlapped protein names
melanomaoverlap <- intersect(melanoma_list3,melanoma_list1)
write.table(melanomaoverlap, file = "melanoma_intersected.txt", sep = "")
##################3
# Draw the diagrams
res <- 300
w <- 8
h <- 8
png("Venn.png", width = w*res, height = h*res, res = res)
cowplot::plot_grid(monocytes, hela, melanoma, pdac, ncol = 2, labels = c('A', 'B', 'C', 'D'), label_size = 16)
dev.off()




