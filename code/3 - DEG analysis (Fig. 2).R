#This script will analyse the differential expressed genes (DEGs), creating the Venn diagram and 
#KEGG enrichment plots shown in Figure 2.
#DEG lists are created in the `1 - RNAseq analysis.R` scripts:
#DEG_24h = MDA-MB-468 cells, 24 hours
#DEG_6h = MDA-MB-468 cells, 6 hours
#DEG_H70 = HCC70 cells, 24 hours

#-------------------------------------------------------
# Variables definitions
#-------------------------------------------------------

#colors definition
colJB2 <- "#afc7e3"
colJB3 <- "#c9ddb4"
colJB6 <- "#e9cac7"
colJ <- "#cfcecb"
colB2 <- "#c7a9d9"
colB3 <- "#d3cdb1"
colB6 <- "#f9e1b7"
colJB2_H70 <- "#81d8f8"
colJB2_6h <- "#c8e7a7"

#------------------------------------------------------------------------------
# Venn diagram of JB2, JB3 and JB6 vs. DMSO in 24 hours (Fig. 2A, S2A)
#------------------------------------------------------------------------------

#upregulated DEGs in MDA468
VennA_MDA468_list <- list("JB2" = DEG_24h$genes[DEG_24h$`JB2_24h_fdr` == 1],
                          "JB3" = DEG_24h$genes[DEG_24h$`JB3_24h_fdr` == 1],
                          "JB6" = DEG_24h$genes[DEG_24h$`JB6_24h_fdr` == 1])
VennA_MDA468_partitions <- get.venn.partitions(VennA_MDA468_list)

tiff(filename = "../results/Figure 2A (Venn).tiff", width = 800, height = 800)
plot(euler(VennA_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2, colJB3, colJB6))
dev.off()

#upregulated DEGs in HCC70
VennA_HCC70_list <- list("JB2" = DEG_H70$genes[DEG_H70$`JB2_H70_fdr` == 1],
                         "JB3" = DEG_H70$genes[DEG_H70$`JB3_H70_fdr` == 1],
                         "JB6" = DEG_H70$genes[DEG_H70$`JB6_H70_fdr` == 1])
VennA_HCC70_partitions <- get.venn.partitions(VennA_HCC70_list)

tiff(filename = "../results/Figure S2A (HCC70 upregulated).tiff", width = 800, height = 800)
plot(euler(VennA_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2_H70, colJB3, colJB6))
dev.off()

#Downregulated DEGs in MDA468
VennB_MDA468_list <- list("JB2" = DEG_24h$genes[DEG_24h$`JB2_24h_fdr` == -1],
                          "JB3" = DEG_24h$genes[DEG_24h$`JB3_24h_fdr` == -1],
                          "JB6" = DEG_24h$genes[DEG_24h$`JB6_24h_fdr` == -1])
VennB_MDA468_partitions <- get.venn.partitions(VennB_MDA468_list)

tiff(filename = "../results/Figure S2A (MDA468 downregulated).tiff", width = 800, height = 800)
plot(euler(VennB_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2, colJB3, colJB6))
dev.off()

#Downregulated DEGs in HCC70
VennB_HCC70_list <- list("JB2" = DEG_H70$genes[DEG_H70$`JB2_H70_fdr` == -1],
                         "JB3" = DEG_H70$genes[DEG_H70$`JB3_H70_fdr` == -1],
                         "JB6" = DEG_H70$genes[DEG_H70$`JB6_H70_fdr` == -1])
VennB_HCC70_partitions <- get.venn.partitions(VennB_HCC70_list)

tiff(filename = "../results/Figure S2A (HCC70 downregulated).tiff", width = 800, height = 800)
plot(euler(VennB_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2_H70, colJB3, colJB6))
dev.off()

#Intersection of both cell lines - DEGs in all combinations (JB2, JB3 and JB6 vs. DMSO) [combined_genes]
combined_genes_468 <- VennA_MDA468_partitions %>% 
  filter(JB2 == TRUE & JB3 == TRUE & JB6 == TRUE) %>% 
  pull(..values..) %>% unlist()
combined_genes_H70 <- VennA_HCC70_partitions %>% 
  filter(JB2 == TRUE & JB3 == TRUE & JB6 == TRUE) %>% 
  pull(..values..) %>% unlist()
combined_genes <- intersect(combined_genes_468, combined_genes_H70)

#Intersection of both cell lines - DEGs only increased by JB2 (and not JB3 or JB6) vs. DMSO [JB2_unique]
JB2_unique_468 <- VennA_MDA468_partitions %>% 
  filter(JB2 == TRUE & JB3 == FALSE & JB6 == FALSE) %>% 
  pull(..values..) %>% unlist()
JB2_unique_H70 <- VennA_HCC70_partitions %>% 
  filter(JB2 == TRUE & JB3 == FALSE & JB6 == FALSE) %>% 
  pull(..values..) %>% unlist()
JB2_unique <- intersect(JB2_unique_468, JB2_unique_H70)

#perform KEGG enrichment of JB2_unique and combined_genes in enrichr
KEGG_combined <- EnrichrFun(combined_genes, library = "KEGG_2019_Human", adjust_pval = TRUE)
KEGG_JB2 <- EnrichrFun(JB2_unique, library = "KEGG_2019_Human", adjust_pval = TRUE)

thisFig <- KEGG_combined[["results"]] %>% filter(logP > 1.3) %>% slice_max(logP, n = 10) %>% 
  mutate(Term = str_replace(Term, "endoplasmic reticulum", "ER"))
p1 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = "gray") + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "-log(p-value)", title = "JB2, JB3, JB6 common genes") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

thisFig <- KEGG_JB2[["results"]] %>% filter(logP > 1) %>% slice_max(logP, n = 10)
p2 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB2) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "-log(p-value)", title = "JB2 only genes") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

ggarrange(p1, p2, ncol = 1, align = "hv", heights = c(1, 1))
ggsave("../results/Figure 2A (KEGG).tiff", dpi = 300, height = 8, width = 8)

#------------------------------------------------------------------------------
# KEGG enrichment of JB2 vs. JB6 in 24 hours (Fig. 2B, S2B)
#------------------------------------------------------------------------------

#NOTE - only in here, I don't do the enrichment of the intersection of the two cell lines (there are too few)
#but for each cell line individualy.
JB6_vs_JB2_468_up <- DEG_24h$genes[DEG_24h$JB6_to_JB2_24h_fdr == 1]
JB6_vs_JB2_468_down <- DEG_24h$genes[DEG_24h$JB6_to_JB2_24h_fdr == -1]
JB6_vs_JB2_H70_up <- DEG_H70$genes[DEG_H70$JB6_to_JB2_H70_fdr == 1]
JB6_vs_JB2_H70_down <- DEG_H70$genes[DEG_H70$JB6_to_JB2_H70_fdr == -1]

KEGG_JB6_vs_JB2_up <- EnrichrFun(JB6_vs_JB2_468_up, library = "KEGG_2019_Human", adjust_pval = TRUE)
KEGG_JB6_vs_JB2_down <- EnrichrFun(JB6_vs_JB2_468_down, library = "KEGG_2019_Human", adjust_pval = TRUE)
KEGG_JB6_vs_JB2_up_H70 <- EnrichrFun(JB6_vs_JB2_H70_up, library = "KEGG_2019_Human", adjust_pval = TRUE)
KEGG_JB6_vs_JB2_down_H70 <- EnrichrFun(JB6_vs_JB2_H70_down, library = "KEGG_2019_Human", adjust_pval = TRUE)

#KEGG enrichment plot for MDA468 (Fig. 2B)
thisFig <- KEGG_JB6_vs_JB2_down[["results"]] %>% filter(logP > 1.3) %>% slice_max(logP, n = 10)
p1 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB2) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "-log(p-value)", title = "Up in JB2 vs. JB6") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

thisFig <- KEGG_JB6_vs_JB2_up[["results"]] %>% filter(logP > 1.3) %>% slice_max(logP, n = 10) %>% 
  mutate(Term = str_replace(Term, "endoplasmic reticulum", "ER"))
p2 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB6) + 
  coord_flip() +
  labs(x = "", y = "-log(p-value)", title = "Up in JB6 vs. JB2") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

ggarrange(p1, p2, ncol = 1, align = "hv")
ggsave("../results/Figure 2B.tiff", dpi = 300, height = 8, width = 8)

#KEGG enrichment plot for HCC70 (Fig. S2B)
thisFig <- KEGG_JB6_vs_JB2_down_H70[["results"]] %>% filter(logP > 1.3) %>% slice_max(logP, n = 10)
p1 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB2) + 
  coord_flip() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  labs(x = "", y = "-log(p-value)", title = "Up in JB2 vs. JB6 (HCC70)") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

thisFig <- KEGG_JB6_vs_JB2_up_H70[["results"]] %>% filter(logP > 1.3) %>% slice_max(logP, n = 10) %>% 
  mutate(Term = str_replace(Term, "endoplasmic reticulum", "ER"))
p2 <- ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB6) + 
  coord_flip() +
  labs(x = "", y = "-log(p-value)", title = "Up in JB6 vs. JB2 (HCC70)") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black"))

ggarrange(p1, p2, ncol = 1, align = "hv")
ggsave("../results/Figure S2B.tiff", dpi = 300, height = 8, width = 8)

#-------------------------------------------------------------------------------
# Venn diagrams for unregulated synergistic genes (Fig. 2C)
# Venns for MDA468 are shown in the manuscript
# Venns for HCC70 are not shown, but are prepared and used in KEGG enrichment
#-------------------------------------------------------------------------------

#synergistic genes in MDA468 Upregulated JB2, J and B2 vs. DMSO
Venn1_MDA468_list <- list("JB2" = DEG_24h$genes[DEG_24h$`JB2_24h_fdr` == 1],
                          "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == 1],
                          "B2" = DEG_24h$genes[DEG_24h$`B2_24h_fdr` == 1])

Venn1_MDA468_partitions <- get.venn.partitions(Venn1_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure 2C (JB2).tiff", width = 800, height = 800)
plot(euler(Venn1_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2, colJ, colB2))
dev.off()

#synergistic genes in MDA468 Upregulated JB2, J and B2 vs. DMSO
Venn1_HCC70_list <- list("JB2" = DEG_H70$genes[DEG_H70$`JB2_H70_fdr` == 1],
                         "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == 1],
                         "B2" = DEG_H70$genes[DEG_H70$`B2_H70_fdr` == 1])

Venn1_HCC70_partitions <- get.venn.partitions(Venn1_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn1_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE, 
     fill = c(colJB2_H70, colJ, colB2))
dev.off()

Venn1_combined <- Venn1_MDA468_partitions %>% left_join(Venn1_HCC70_partitions)

#synergistic genes in MDA468 Upregulated JB3, J and B3 vs. DMSO
Venn2_MDA468_list <- list("JB3" = DEG_24h$genes[DEG_24h$`JB3_24h_fdr` == 1],
                          "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == 1],
                          "B3" = DEG_24h$genes[DEG_24h$`B3_24h_fdr` == 1])

Venn2_MDA468_partitions <- get.venn.partitions(Venn2_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure 2C (JB3).tiff", width = 800, height = 800)
plot(euler(Venn2_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB3, colJ, colB3))
dev.off()

#synergistic genes in HCC70 Upregulated JB3, J and B3 vs. DMSO
Venn2_HCC70_list <- list("JB3" = DEG_H70$genes[DEG_H70$`JB3_H70_fdr` == 1],
                         "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == 1],
                         "B3" = DEG_H70$genes[DEG_H70$`B3_H70_fdr` == 1])

Venn2_HCC70_partitions <- get.venn.partitions(Venn2_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn2_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB3, colJ, colB3))
dev.off()

Venn2_combined <- Venn2_MDA468_partitions %>% left_join(Venn2_HCC70_partitions)

#synergistic genes in MDA468 Upregulated JB6, J and B6 vs. DMSO
Venn3_MDA468_list <- list("JB6" = DEG_24h$genes[DEG_24h$`JB6_24h_fdr` == 1],
                          "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == 1],
                          "B6" = DEG_24h$genes[DEG_24h$`B6_24h_fdr` == 1])

Venn3_MDA468_partitions <- get.venn.partitions(Venn3_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure 2C (JB6).tiff", width = 800, height = 800)
plot(euler(Venn3_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB6, colJ, colB6))
dev.off()

#synergistic genes in HCC70 Upregulated JB6, J and B6 vs. DMSO
Venn3_HCC70_list <- list("JB6" = DEG_H70$genes[DEG_H70$`JB6_H70_fdr` == 1],
                         "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == 1],
                         "B6" = DEG_H70$genes[DEG_H70$`B6_H70_fdr` == 1])

Venn3_HCC70_partitions <- get.venn.partitions(Venn3_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn3_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB6, colJ, colB6))
dev.off()

Venn3_combined <- Venn3_MDA468_partitions %>% left_join(Venn3_HCC70_partitions)

#-------------------------------------------------------------------------------
# Venn diagrams for downregulated synergistic genes (Fig. S2C)
# Venns for MDA468 are shown in the manuscript
# Venns for HCC70 are not shown, but are prepared and used in KEGG enrichment
#-------------------------------------------------------------------------------

#synergistic genes in MDA468 downregulated JB2, J and B2 vs. DMSO
Venn1D_MDA468_list <- list("JB2" = DEG_24h$genes[DEG_24h$`JB2_24h_fdr` == -1],
                           "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == -1],
                           "B2" = DEG_24h$genes[DEG_24h$`B2_24h_fdr` == -1])

Venn1D_MDA468_partitions <- get.venn.partitions(Venn1D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure S2C (JB2).tiff", width = 800, height = 800)
plot(euler(Venn1D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2, colJ, colB2))
dev.off()

#synergistic genes in HCC70 downregulated JB2, J and B2 vs. DMSO
Venn1D_HCC70_list <- list("JB2" = DEG_H70$genes[DEG_H70$`JB2_H70_fdr` == -1],
                          "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == -1],
                          "B2" = DEG_H70$genes[DEG_H70$`B2_H70_fdr` == -1])

Venn1D_HCC70_partitions <- get.venn.partitions(Venn1D_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn1D_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE, 
     fill = c(colJB2_H70, colJ, colB2))
dev.off()

Venn1D_combined <- Venn1D_MDA468_partitions %>% left_join(Venn1D_HCC70_partitions)

#synergistic genes in MDA468 downregulated JB3, J and B3 vs. DMSO
Venn2D_MDA468_list <- list("JB3" = DEG_24h$genes[DEG_24h$`JB3_24h_fdr` == -1],
                           "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == -1],
                           "B3" = DEG_24h$genes[DEG_24h$`B3_24h_fdr` == -1])

Venn2D_MDA468_partitions <- get.venn.partitions(Venn2D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure S2C (JB3).tiff", width = 800, height = 800)
plot(euler(Venn2D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB3, colJ, colB3))
dev.off()

#synergistic genes in HCC70 downregulated JB3, J and B3 vs. DMSO
Venn2D_HCC70_list <- list("JB3" = DEG_H70$genes[DEG_H70$`JB3_H70_fdr` == -1],
                          "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == -1],
                          "B3" = DEG_H70$genes[DEG_H70$`B3_H70_fdr` == -1])

Venn2D_HCC70_partitions <- get.venn.partitions(Venn2D_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn2D_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB3, colJ, colB3))
dev.off()

Venn2D_combined <- Venn2D_MDA468_partitions %>% left_join(Venn2D_HCC70_partitions)

#synergistic genes in MDA468 downregulated JB6, J and B6 vs. DMSO
Venn3D_MDA468_list <- list("JB6" = DEG_24h$genes[DEG_24h$`JB6_24h_fdr` == -1],
                           "J" = DEG_24h$genes[DEG_24h$`J_24h_fdr` == -1],
                           "B6" = DEG_24h$genes[DEG_24h$`B6_24h_fdr` == -1])

Venn3D_MDA468_partitions <- get.venn.partitions(Venn3D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure S2C (JB6).tiff", width = 800, height = 800)
plot(euler(Venn3D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB6, colJ, colB6))
dev.off()

#synergistic genes in HCC70 downregulated JB6, J and B6 vs. DMSO
Venn3D_HCC70_list <- list("JB6" = DEG_H70$genes[DEG_H70$`JB6_H70_fdr` == -1],
                          "J" = DEG_H70$genes[DEG_H70$`J_H70_fdr` == -1],
                          "B6" = DEG_H70$genes[DEG_H70$`B6_H70_fdr` == -1])

Venn3D_HCC70_partitions <- get.venn.partitions(Venn3D_HCC70_list) %>% 
  dplyr::rename("HCC70_values" = "..values..", "HCC70_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn3D_HCC70_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB6, colJ, colB6))
dev.off()

Venn3D_combined <- Venn3D_MDA468_partitions %>% left_join(Venn3D_HCC70_partitions)

#-------------------------------------------------------------------------------
# Venn diagram of the synergistic genes (24h) in both cell lines (Fig. S2D)
#-------------------------------------------------------------------------------

#get the synergistic genes from the venns for JB2, J and B2 vs. DMSO from above
SynGenes_468 <- Venn1_MDA468_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(MDA468_values) %>% unlist()
SynGenes_H70 <- Venn1_HCC70_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(HCC70_values) %>% unlist()
SynGenes_both <- intersect(SynGenes_468, SynGenes_H70)

SynGenes_468_down <- Venn1D_MDA468_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(MDA468_values) %>% unlist()
SynGenes_H70_down <- Venn1D_HCC70_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(HCC70_values) %>% unlist()
SynGenes_both_down <- intersect(SynGenes_468_down, SynGenes_H70_down)

#plot the Venn diagram
Venn_SynGenes_list <- list("MDA468" = SynGenes_468, "HCC70" = SynGenes_H70)

tiff(filename = "../results/Figure S2D.tiff", width = 800, height = 800)
plot(euler(Venn_SynGenes_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = FALSE,
     fill = c(colJB2, colJB2_H70))
dev.off()

#-------------------------------------------------------------------------------
# KEGG enrichment plot for JB2 synergistic genes in 24h (Fig. 2E)
#-------------------------------------------------------------------------------

KEGG_JB2_syn_up <- EnrichrFun(SynGenes_both, library = "KEGG_2019_Human", adjust_pval = TRUE)

thisFig <- KEGG_JB2_syn_up[["results"]] %>% 
  filter(logP > 0.1) %>% 
  slice_max(logP, n = 20) %>% 
  mutate(Term = str_replace(Term, "endoplasmic reticulum", "ER"))

ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB2) + 
  coord_flip() +
  labs(x = "", y = "-log(p-value)") +
  theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 18, color = "black"),
        axis.title = element_text(size = 18, color = "black"))

ggsave("../results/Figure 2E.tiff", dpi = 300, height = 5.5, width = 7)

#-------------------------------------------------------------------------------
# Venn diagram of synergistic genes in 6h (not shown in the manuscript, 
# but used in the KEGG enrichment
#-------------------------------------------------------------------------------

#synergistic genes in MDA468 6h up-regulated JB2, J and B2 vs. DMSO
Venn6_MDA468_list <- list("JB2" = DEG_6h$genes[DEG_6h$`JB2_6h_fdr` == 1],
                          "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == 1],
                          "B2" = DEG_6h$genes[DEG_6h$`B2_6h_fdr` == 1])

Venn6_MDA468_partitions <- get.venn.partitions(Venn6_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn6_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB2_6h, colJ, colB2))
dev.off()

#synergistic genes in MDA468 6h up-regulated JB3, J and B3 vs. DMSO
Venn7_MDA468_list <- list("JB3" = DEG_6h$genes[DEG_6h$`JB3_6h_fdr` == 1],
                          "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == 1],
                          "B3" = DEG_6h$genes[DEG_6h$`B3_6h_fdr` == 1])

Venn7_MDA468_partitions <- get.venn.partitions(Venn7_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn7_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB3, colJ, colB3))
dev.off()

#synergistic genes in MDA468 6h up-regulated JB6, J and B6 vs. DMSO
Venn8_MDA468_list <- list("JB6" = DEG_6h$genes[DEG_6h$`JB6_6h_fdr` == 1],
                          "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == 1],
                          "B6" = DEG_6h$genes[DEG_6h$`B6_6h_fdr` == 1])

Venn8_MDA468_partitions <- get.venn.partitions(Venn8_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn8_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB6, colJ, colB6))
dev.off()

#synergistic genes in MDA468 6h down-regulated JB2, J and B2 vs. DMSO
Venn6D_MDA468_list <- list("JB2" = DEG_6h$genes[DEG_6h$`JB2_6h_fdr` == -1],
                           "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == -1],
                           "B2" = DEG_6h$genes[DEG_6h$`B2_6h_fdr` == -1])

Venn6D_MDA468_partitions <- get.venn.partitions(Venn6D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn6D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB2_6h, colJ, colB2))
dev.off()

#synergistic genes in MDA468 6h down-regulated JB3, J and B3 vs. DMSO
Venn7D_MDA468_list <- list("JB3" = DEG_6h$genes[DEG_6h$`JB3_6h_fdr` == -1],
                           "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == -1],
                           "B3" = DEG_6h$genes[DEG_6h$`B3_6h_fdr` == -1])

Venn7D_MDA468_partitions <- get.venn.partitions(Venn7D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn7D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB3, colJ, colB3))
dev.off()

#synergistic genes in MDA468 6h down-regulated JB6, J and B6 vs. DMSO
Venn8D_MDA468_list <- list("JB6" = DEG_6h$genes[DEG_6h$`JB6_6h_fdr` == -1],
                           "J" = DEG_6h$genes[DEG_6h$`J_6h_fdr` == -1],
                           "B6" = DEG_6h$genes[DEG_6h$`B6_6h_fdr` == -1])

Venn8D_MDA468_partitions <- get.venn.partitions(Venn8D_MDA468_list) %>% 
  dplyr::rename("MDA468_values" = "..values..", "MDA468_count" = "..count..")

tiff(filename = "../results/Figure not shown.tiff", width = 800, height = 800)
plot(euler(Venn8D_MDA468_list), quantities = list(cex = 5), lty = 1, lwd = 3, labels = FALSE, 
     fill = c(colJB6, colJ, colB6))
dev.off()

#-------------------------------------------------------------------------------
# Extract the up- and down-regulated JB2 synergistic genes (6h)
#-------------------------------------------------------------------------------

SynGenes_468_6h_up <- Venn6_MDA468_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(MDA468_values) %>% unlist()

SynGenes_468_6h_down <- Venn6D_MDA468_partitions %>% 
  filter(JB2 == TRUE & J == FALSE & B2 == FALSE) %>% 
  pull(MDA468_values) %>% unlist()

#combine all JB2 genes (24 hours up and down, 6 hours up and down)
#this partitions are visualized in Figure 2D (created in excel based on this data)
Venn_all_synergistic_list <- list("24h_up" = SynGenes_both, "24h_down" = SynGenes_both_down,
                    "6h_up" = SynGenes_468_6h_up, "6h_down" = SynGenes_468_6h_down)
Venn_all_synergistic_partitions <- get.venn.partitions(Venn_all_synergistic_list)

#-------------------------------------------------------------------------------
# KEGG enrichment of these genes (Fig S2E)
#-------------------------------------------------------------------------------

KEGG_6h <- EnrichrFun(SynGenes_468_6h_up, library = "KEGG_2019_Human", adjust_pval = TRUE)

thisFig <- KEGG_6h[["results"]] %>% 
  filter(logP > 1.3) %>% 
  slice_max(logP, n = 15) 

ggplot(thisFig, aes(x = reorder(Term, logP), y = logP)) +
  geom_col(fill = colJB2) + 
  coord_flip() +
  labs(x = "", y = "-log(p-value)") +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(size = 16, color = "black"))

ggsave("../results/Figure S2E.tiff", dpi = 300, height = 5, width = 10)

#-------------------------------------------------------------------------------
# Enrichment of ATF4-only targets (Fig. S4F)
#-------------------------------------------------------------------------------

#Enrichment analysis is performed using the Camera method in the Limma package.
#Helper function for the Camera method:
CameraFunction <- function(contr.matrix, V, design, idx, sort = TRUE) {
  cam.list <- list()
  
  for(n in 1:8) {
    name <- colnames(contr.matrix)[n]
    
    cam.list[[n]] <- camera(V, idx, design, contrast = contr.matrix[,n], sort = sort) %>% 
      as_tibble(rownames = "pathway")   %>% 
      mutate(result = ifelse(Direction == "Up", PValue, -PValue)) %>% 
      dplyr::select(pathway, result) 
    
    names(cam.list[[n]])[names(cam.list[[n]]) == "result"] <- name
  }
  cam.object <- Reduce(full_join, cam.list)
  
  cam.object.log <- cam.object %>% 
    mutate(across(where(is.numeric), ~ifelse(.x>0, -log10(.x), log10(abs(.x)))))
  
  return(cam.object.log)
}

#get genesets (taken from this paper: PMID: 23624402)
genesets <- gmt2list("../data/genesets used in the paper.gmt")
ATF4_genesets <- genesets[c("ATF4 only", "ATF4 and CHOP")]

#perform Camera enrichment
idx <- ids2indices(ATF4_genesets, identifiers = V_24h$genes$genes)
cam.object_24h <- CameraFunction(contr.matrix_24h, V_24h, design_24h, idx)
idx <- ids2indices(ATF4_genesets, identifiers = V_6h$genes$genes)
cam.object_6h <- CameraFunction(contr.matrix_6h, V_6h, design_6h, idx)
cam.object <- cam.object_6h %>% full_join(cam.object_24h)

thisFig <- cam.object %>% 
  filter(pathway %in% c("ATF4 and CHOP", "ATF4 only")) %>% 
  gather(-pathway, key = "treat", value = "value") %>% 
  mutate(cell = case_when(str_detect(treat, "24h") ~ "MDA468 24h",
                          str_detect(treat, "6h") ~ "MDA468 6h")) %>% 
  mutate(treat = str_replace(treat, "_6h|_24h", "")) %>% 
  filter(treat != "JB6_to_JB2") %>% 
  mutate(treat = factor(treat, levels = c("JB2", "JB3", "JB6", "J", "B2", "B3", "B6"))) %>% 
  mutate(cell = factor(cell, levels = c("MDA468 6h", "MDA468 24h"))) %>% 
  mutate(pathway = ifelse(pathway == "ATF4 and CHOP", "ATF4 + CHOP", pathway))

ggplot(thisFig, aes(x = treat, y = value, fill = cell)) +
  geom_col() + 
  scale_fill_manual(values = c("black", "gray50", "gray70")) +
  geom_hline(yintercept = 1.33, linetype = "dashed") +
  facet_grid(pathway~cell) +
  labs(x = "", y = "Enrichment of targets (-log[p-value])") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave("../results/Figure S4F.tiff", dpi = 200, width = 7, height = 5)
