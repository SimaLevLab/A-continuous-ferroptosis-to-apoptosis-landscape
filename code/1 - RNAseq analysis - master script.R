#Scripts written by Yaron Vinik

#The main script takes the count matrix as input, and calculate the differential expressed genes
#Further, this script will create the Gradient Gene Signature (GGS), needed for figure 5.
#This script will also define the libraries and helper functions used in all other scripts.
#This script should be run before all others.

#libraries for all the scripts:
library(tidyverse)
library(limma)
library(edgeR)
library(factoextra)
library(ggrepel)
library(ggpubr)
library(ggforce)
library(VennDiagram)
library(eulerr)
library(enrichR)
library(cogena)
library(ggpp)
library(umap)
library(tidymodels)
library(pracma)
library(simputation)
library(GSVA)
library(rstatix)

#-------------------------------------------------------------------------------
# helper functions
#-------------------------------------------------------------------------------

#DGE function to create DEGList object, normalize counts, filter low expressing genes, measure log CPM, perform voom method
DEGCreate <- function (cts, name, plot = FALSE) {
  
  dge <- DGEList(counts = cts, genes = rownames(cts))
  groups <- colnames(dge)
  groups <- as.factor(str_sub(groups, 1, str_length(groups)-2))
  dge$samples$group <- groups
  
  library(org.Hs.eg.db)
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    m <- match(dge$genes$genes, egSYMBOL$symbol)
    dge$genes$gene_id <- egSYMBOL$gene_id[m]
  detach(package:org.Hs.eg.db)
  detach(package:AnnotationDbi)
  
  keep <- filterByExpr(dge, group = groups)
  dge <- dge[keep,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  CPM <- cpm(dge)
  LogCPM <- cpm(dge, log = TRUE)

  if(plot) {
    par(mfrow=c(1,2))
    col.group <- groups
    levels(col.group) <-  brewer.pal(nlevels(col.group), "Set2")
    col.group <- as.character(col.group)
    plotMDS(LogCPM, labels=groups, col=col.group)
    title(main="MDS dim 1 & 2")
  }
  
  designMatrix <- model.matrix(~0+groups)
  colnames(designMatrix) <- gsub("groups", "", colnames(designMatrix))

  V <- voom(dge, designMatrix, plot = plot)
}

#Function used to perform KEGG enrichment using the enrichR site
EnrichrFun <- function(genes, library, adjust_pval = TRUE) {
  x <- enrichr(genes, library)
  EnrichrResult <- x[[1]] %>% 
    as_tibble() %>% 
    dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Combined.Score, Genes) %>% 
    rowwise() %>% 
    mutate(logP = ifelse(adjust_pval, -log10(Adjusted.P.value), -log10(P.value))) %>% 
    ungroup()
  
  for(i in 1:length(genes)) {
    gene_name = genes[i]
    
    EnrichrResult <- EnrichrResult %>% 
      mutate(!!gene_name := ifelse(str_detect(EnrichrResult$Genes, gene_name), 1, ""))
  }
  
  EnrichrResult_Long <- EnrichrResult %>% 
    dplyr::select(-c(Overlap, P.value, Adjusted.P.value, Combined.Score, Genes, logP)) %>% 
    gather(-Term, key = "Gene", value = "value") %>% 
    spread(key = "Term", value = "value") %>%  
    mutate(across(-Gene, ~ifelse(.x == 1, cur_column(), NA))) %>% 
    unite(col = "Terms", -Gene, sep = " | ", na.rm = TRUE)
  
  EnrichrResult_List <- list(
    "results" = EnrichrResult, 
    "long results" = EnrichrResult_Long)
  
  return(EnrichrResult_List)
}

#-------------------------------------------------------------------------------
# count matrix analysis using Limma
#-------------------------------------------------------------------------------

#Analyze MDA468 in 24hrs
cts_24h <- read.csv("../data/RNAseq count matrix for combination - MDA468 24h.csv", row.names = "Gene")

V_24h <- DEGCreate(cts_24h, name = "MDA468 24hr 1st set")

design_24h <- V_24h$design
vfit_24h <- lmFit(V_24h, design_24h)
contr.matrix_24h <- makeContrasts(
  JB2_24h = MDA468_JB2_24h - MDA468_D_24h,
  JB3_24h = MDA468_JB3_24h - MDA468_D_24h,
  JB6_24h = MDA468_JB6_24h - MDA468_D_24h,
  J_24h = MDA468_J_24h - MDA468_D_24h,
  B2_24h = MDA468_B2_24h - MDA468_D_24h,
  B3_24h = MDA468_B3_24h - MDA468_D_24h,
  B6_24h = MDA468_B6_24h - MDA468_D_24h,
  JB6_to_JB2_24h = MDA468_JB6_24h - MDA468_JB2_24h,
  B6_to_B2_24h = MDA468_B6_24h - MDA468_B2_24h,
  levels = colnames(design_24h))

vfit2_24h <- contrasts.fit(vfit_24h, contr.matrix_24h)
efit_24h <- eBayes(vfit2_24h, 0.01)
tT_24h <- topTable(efit_24h, adjust.method = "fdr", sort.by="B", n = Inf)

decide_24h <- decideTests(efit_24h) 
decide_24h <- as.data.frame(decide_24h)
decide_24h <- as_tibble(decide_24h, rownames = "genes") %>% 
  rename_with(~paste(.x, "_fdr", sep = ""), .cols = -genes)

DEG_24h <- tT_24h %>% left_join(decide_24h, by = "genes")

#Analyze HCC70 in 24hrs
cts_H70 <- read.csv("../data/RNAseq count matrix for combination - HCC70 24h.csv", row.names = "Gene")

V_H70 <- DEGCreate(cts_H70, name = "HCC70 24hr 2nd set")

design_H70 <- V_H70$design
vfit_H70 <- lmFit(V_H70, design_H70)
contr.matrix_H70 <- makeContrasts(
  JB2_H70 = HCC70_JB2 - HCC70_D,
  JB3_H70 = HCC70_JB3 - HCC70_D,
  JB6_H70 = HCC70_JB6 - HCC70_D,
  J_H70 = HCC70_J - HCC70_D,
  B2_H70 = HCC70_B2 - HCC70_D,
  B3_H70 = HCC70_B3 - HCC70_D,
  B6_H70 = HCC70_B6 - HCC70_D,
  JB6_to_JB2_H70 = HCC70_JB6 - HCC70_JB2,
  B6_to_B2_H70 = HCC70_B6 - HCC70_B2,
  levels = colnames(design_H70))

vfit2_H70 <- contrasts.fit(vfit_H70, contr.matrix_H70)
efit_H70 <- eBayes(vfit2_H70, 0.01)
tT_H70 <- topTable(efit_H70, adjust.method = "fdr", sort.by="B", n = Inf)

decide_H70 <- decideTests(efit_H70) 
decide_H70 <- as.data.frame(decide_H70)
decide_H70 <- as_tibble(decide_H70, rownames = "genes") %>% 
  rename_with(~paste(.x, "_fdr", sep = ""), .cols = -genes)

DEG_H70 <- tT_H70 %>% left_join(decide_H70, by = "genes")

# Analyze MDA468 in 6 hrs
cts_6h <- read.csv("../data/RNAseq count matrix for combination - MDA468 6hr.csv", row.names = "Gene")

V_6h <- DEGCreate(cts_6h, name = "MDA468 6hr 2nd set")

design_6h <- V_6h$design
vfit_6h <- lmFit(V_6h, design_6h)
contr.matrix_6h <- makeContrasts(
  JB2_6h = MDA468_JB2_6h - MDA468_D_6h,
  JB3_6h = MDA468_JB3_6h - MDA468_D_6h,
  JB6_6h = MDA468_JB6_6h - MDA468_D_6h,
  J_6h = MDA468_J_6h - MDA468_D_6h,
  B2_6h = MDA468_B2_6h - MDA468_D_6h,
  B3_6h = MDA468_B3_6h - MDA468_D_6h,
  B6_6h = MDA468_B6_6h - MDA468_D_6h,
  JB6_to_JB2_6h = MDA468_JB6_6h - MDA468_JB2_6h,
  B6_to_B2_6h = MDA468_B6_6h - MDA468_B2_6h,
  levels = colnames(design_6h))

vfit2_6h <- contrasts.fit(vfit_6h, contr.matrix_6h)
efit_6h <- eBayes(vfit2_6h, 0.01)
tT_6h <- topTable(efit_6h, adjust.method = "fdr", sort.by="B", n = Inf)

decide_6h <- decideTests(efit_6h) 
decide_6h <- as.data.frame(decide_6h)
decide_6h <- as_tibble(decide_6h, rownames = "genes") %>% 
  rename_with(~paste(.x, "_fdr", sep = ""), .cols = -genes)

DEG_6h <- tT_6h %>% left_join(decide_6h, by = "genes")

#combine all DEG tables
DEG_combined <- DEG_24h %>% 
  full_join(DEG_H70, by = "genes") %>% full_join(DEG_6h, by = "genes") %>% as_tibble()

#----------------------------------------------------------------------
# creating the "Gradient Gene Set" (GGS) subsets
#----------------------------------------------------------------------

# See 'supplemental text' for details
# In each subset, either JB2 or JB6 must be significantly up (fdr == 1) or down regulated (fdr == -1)
# JB2 and JB6 can't be both significant vs. DMSO. 
# After this filteration, the genes are sorted by their gradient value (JB6 to JB2 fold change),
# and a cutoff is set to take the genes with the highest gradient.
# created here are 2 version of the GGS, with cutoff = 100 and cutoff = 300
# NOTE: due to this definition, there are 3 genes overlap between JB2U and JB6D groups, those were added to the JB2U group,
# and 1 gene overlap between JB2D and JB6U, this gene was added to the JB2D group.

cutoff <- 300

#create the 300 GGS in MDA468 cells
gradient_468 <- DEG_24h %>%       
  mutate(gradient = case_when(JB2_24h_fdr == 1 & JB6_24h_fdr != 1 ~ "JB2_up",
                              JB6_24h_fdr == -1 & JB2_24h_fdr != -1 ~ "JB6_down",
                              JB6_24h_fdr == 1 & JB2_24h_fdr != 1 ~ "JB6_up",
                              JB2_24h_fdr == -1 & JB6_24h_fdr != -1 ~"JB2_down")) 
group_JB2U_468_300 <- gradient_468 %>% filter(gradient == "JB2_up") %>% 
  slice_min(JB6_to_JB2_24h, n = cutoff)
group_JB6D_468_300 <- gradient_468 %>% filter(gradient == "JB6_down") %>% 
  slice_min(JB6_to_JB2_24h, n = cutoff)
group_JB6U_468_300 <- gradient_468 %>% filter(gradient == "JB6_up") %>% 
  slice_max(JB6_to_JB2_24h, n = cutoff)
group_JB2D_468_300 <- gradient_468 %>% filter(gradient == "JB2_down") %>% 
  slice_max(JB6_to_JB2_24h, n = cutoff)

#create the 300 GGS in HCC70 cells
gradient_H70 <- DEG_H70 %>%       
  mutate(gradient = case_when(JB2_H70_fdr == 1 & JB6_H70_fdr != 1 ~ "JB2_up",
                              JB6_H70_fdr == -1 & JB2_H70_fdr != -1 ~ "JB6_down",
                              JB6_H70_fdr == 1 & JB2_H70_fdr != 1 ~ "JB6_up",
                              JB2_H70_fdr == -1 & JB6_H70_fdr != -1 ~"JB2_down")) 
group_JB2U_H70_300 <- gradient_H70 %>% filter(gradient == "JB2_up") %>% 
  slice_min(JB6_to_JB2_H70, n = cutoff)
group_JB6D_H70_300 <- gradient_H70 %>% filter(gradient == "JB6_down") %>% 
  slice_min(JB6_to_JB2_H70, n = cutoff)
group_JB6U_H70_300 <- gradient_H70 %>% filter(gradient == "JB6_up") %>% 
  slice_max(JB6_to_JB2_H70, n = cutoff)
group_JB2D_H70_300 <- gradient_H70 %>% filter(gradient == "JB2_down") %>% 
  slice_max(JB6_to_JB2_H70, n = cutoff)

#combine the two cell lines:
#the condition - the gene follows the gradient rules in one cell line, and has a trend in the other.
group_JB2U_C_300 <- DEG_combined %>% 
  filter(genes %in% c(group_JB2U_468_300$genes, group_JB2U_H70_300$genes)) %>% 
  filter(JB2_24h > JB6_24h, JB2_H70 > JB6_H70, JB2_24h > 0, JB2_H70 > 0)

group_JB2D_C_300 <- DEG_combined %>% 
  filter(genes %in% c(group_JB2D_468_300$genes, group_JB2D_H70_300$genes)) %>% 
  filter(JB2_24h < JB6_24h, JB2_H70 < JB6_H70, JB2_24h < 0, JB2_H70 < 0)

group_JB6U_C_300 <- DEG_combined %>% 
  filter(genes %in% c(group_JB6U_468_300$genes, group_JB6U_H70_300$genes)) %>% 
  filter(JB2_24h < JB6_24h, JB2_H70 < JB6_H70, JB6_24h > 0, JB6_H70 > 0)

group_JB6D_C_300 <- DEG_combined %>% 
  filter(genes %in% c(group_JB6D_468_300$genes, group_JB6D_H70_300$genes)) %>% 
  filter(JB2_24h > JB6_24h, JB2_H70 > JB6_H70, JB6_24h < 0, JB6_H70 < 0)

cutoff <- 100

#create the 100 GGS in MDA468 cells
gradient_468 <- DEG_24h %>%       
  mutate(gradient = case_when(JB2_24h_fdr == 1 & JB6_24h_fdr != 1 ~ "JB2_up",
                              JB6_24h_fdr == -1 & JB2_24h_fdr != -1 ~ "JB6_down",
                              JB6_24h_fdr == 1 & JB2_24h_fdr != 1 ~ "JB6_up",
                              JB2_24h_fdr == -1 & JB6_24h_fdr != -1 ~"JB2_down")) 
group_JB2U_468_100 <- gradient_468 %>% filter(gradient == "JB2_up") %>% 
  slice_min(JB6_to_JB2_24h, n = cutoff)
group_JB6D_468_100 <- gradient_468 %>% filter(gradient == "JB6_down") %>% 
  slice_min(JB6_to_JB2_24h, n = cutoff)
group_JB6U_468_100 <- gradient_468 %>% filter(gradient == "JB6_up") %>% 
  slice_max(JB6_to_JB2_24h, n = cutoff)
group_JB2D_468_100 <- gradient_468 %>% filter(gradient == "JB2_down") %>% 
  slice_max(JB6_to_JB2_24h, n = cutoff)

#create the 100 GGS in HCC70 cells
gradient_H70 <- DEG_H70 %>%       
  mutate(gradient = case_when(JB2_H70_fdr == 1 & JB6_H70_fdr != 1 ~ "JB2_up",
                              JB6_H70_fdr == -1 & JB2_H70_fdr != -1 ~ "JB6_down",
                              JB6_H70_fdr == 1 & JB2_H70_fdr != 1 ~ "JB6_up",
                              JB2_H70_fdr == -1 & JB6_H70_fdr != -1 ~"JB2_down")) 
group_JB2U_H70_100 <- gradient_H70 %>% filter(gradient == "JB2_up") %>% 
  slice_min(JB6_to_JB2_H70, n = cutoff)
group_JB6D_H70_100 <- gradient_H70 %>% filter(gradient == "JB6_down") %>% 
  slice_min(JB6_to_JB2_H70, n = cutoff)
group_JB6U_H70_100 <- gradient_H70 %>% filter(gradient == "JB6_up") %>% 
  slice_max(JB6_to_JB2_H70, n = cutoff)
group_JB2D_H70_100 <- gradient_H70 %>% filter(gradient == "JB2_down") %>% 
  slice_max(JB6_to_JB2_H70, n = cutoff)

#combine the two cell lines:
group_JB2U_C_100 <- DEG_combined %>% 
  filter(genes %in% c(group_JB2U_468_100$genes, group_JB2U_H70_100$genes)) %>% 
  filter(JB2_24h > JB6_24h, JB2_H70 > JB6_H70, JB2_24h > 0, JB2_H70 > 0)

group_JB2D_C_100 <- DEG_combined %>% 
  filter(genes %in% c(group_JB2D_468_100$genes, group_JB2D_H70_100$genes)) %>% 
  filter(JB2_24h < JB6_24h, JB2_H70 < JB6_H70, JB2_24h < 0, JB2_H70 < 0)

group_JB6U_C_100 <- DEG_combined %>% 
  filter(genes %in% c(group_JB6U_468_100$genes, group_JB6U_H70_100$genes)) %>% 
  filter(JB2_24h < JB6_24h, JB2_H70 < JB6_H70, JB6_24h > 0, JB6_H70 > 0)

group_JB6D_C_100 <- DEG_combined %>% 
  filter(genes %in% c(group_JB6D_468_100$genes, group_JB6D_H70_100$genes)) %>% 
  filter(JB2_24h > JB6_24h, JB2_H70 > JB6_H70, JB6_24h < 0, JB6_H70 < 0)
