#This script will create the univesal set and the specific ferroptosis biomarkers set, as shown in figure 6E-F

set.seed(376)

#--------------------------------------------------------------------------------------------
# Create the universal ferroptosis set
#--------------------------------------------------------------------------------------------

#Recreate the ferroptosis-apoptosis universal set, the set of genes most differentiating between the 18 FIN datasets
#and 18 AINs datasets. This set is biased to those datasets. In order to reduce the bias, the geneset was determined by 
#randomly selecting 10 FINs and 10 AINs datasets (with replacement, so the number of datasets can be less than 10) 
#and extract the top 400 differential genes.
#This was repeated 1000 times. The genes were ranked according to their frequency of appearing among those 1000 iterations.
#See supplemental text for further details. Run-time is ~10 minutes.

#load the 36 datasets:
Ferroptosis_combined <- read_csv("../data/Ferroptosis inducers datasets.csv")
Apoptosis_combined <- read_csv("../data/Apoptosis inducers datasets.csv")

Lineages <- tribble(~"dataset", ~"lineage", "F1", "Fibrosarcoma", "F2", "Liver", "F3", "Kidney", "F4", "Neuroblastoma",
                    "F5", "Mouse emb fibro", "F6", "Colorectal", "F7", "Normal epithelial", "F8", "Lymphoblastoid",
                    "F9", "Lymphoblastoid", "F10", "Pancreas", "F11", "Choriocarcinoma", "F12", "Tregs",
                    "F13", "Myeloma", "F14", "Oligodendrocyte progenitor", "F15", "Glioblastoma", "F16", "TNBC",
                    "F17", "MDSC-like", "F18", "TNBC", "A1", "TNBC", "A2", "Pancreas",
                    "A3", "TNBC", "A4", "TNBC", "A5", "Breast", "A6", "Breast",
                    "A7", "TNBC", "A8", "Breast", "A9", "TNBC", "A10", "Colorectal",
                    "A11", "Prostate", "A12", "Prostate", "A13", "Esophageal", "A14", "Liver",
                    "A15", "Neuroblastoma", "A16", "Skin", "A17", "NK/T-cell lymphoma", "A18", "Myeloma")

#Add the combinations from this project - this will create similar table as above, for the JB2-JB3-JB6 combinations
combinations_t <- efit_24h[["t"]] %>% 
  as_tibble(rownames = "ID") %>% 
  dplyr::select(-contains("to")) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("t_", .x)) %>% 
  dplyr::select(ID, t_JB2_24h:t_JB6_24h)

combinations_FC <- DEG_24h %>% 
  select("ID" = "genes", JB2_24h:JB6_24h) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("logFC_", .x))
len <- dim(combinations_t)[1]

combinations_t2 <- combinations_t %>% 
  mutate(across(-ID, ~min_rank(.x))) %>% 
  mutate(across(-ID, ~.x/len)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "t_", "Rank_t_"))
combinations_FC2 <- combinations_FC %>% 
  mutate(across(-ID, ~min_rank(.x))) %>% 
  mutate(across(-ID, ~.x/len)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "logFC_", "Rank_FC_"))

combinations <- left_join(combinations_FC, combinations_t, by = "ID") %>% 
  left_join(combinations_t2, by = "ID") %>% 
  left_join(combinations_FC2, by = "ID")

#create the univesal set (this part takes a few minutes)
take <- "Rank_t"    #perform this process on the ranked t-statistics.
iterations <- 1000  #repeat the process 1000 times
Ferroptosis_consensus_keep <- list()
Apoptosis_consensus_keep <- list()

for(n in 1:iterations) {
  cat(paste0("Iteration: ", n, " / ", iterations, "\r"))
  ferroptosis_dataset_select <- round(runif(n = 10, min = 2, max = 19), 0)
  apoptosis_dataset_select <- round(runif(n = 10, min = 2, max = 19), 0)
  
  Ferroptosis_DF <- Ferroptosis_combined %>% 
    select(ID, contains(take)) %>%   #taking Rank_t as the metric by which the universal set will be built
    select(ID, all_of(ferroptosis_dataset_select)) %>% 
    mutate(count_na = rowSums(is.na(.))) %>% 
    filter(count_na < 4) %>% #remove genes with more than 4 missing values across the 10 datasets
    select(-count_na) %>% 
    gather(-ID, key = "dataset", value = "value") %>% 
    group_by(ID) %>% 
    mutate(mean_ferroptosis = mean(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(-dataset, -value) %>% 
    distinct(ID, .keep_all = TRUE)
  
  Apoptosis_DF <- Apoptosis_combined %>% 
    select(ID, contains(take)) %>% 
    select(ID, all_of(apoptosis_dataset_select)) %>% 
    mutate(count_na = rowSums(is.na(.))) %>% 
    filter(count_na < 4) %>% 
    select(-count_na) %>%   
    gather(-ID, key = "dataset", value = "value") %>% 
    group_by(ID) %>% 
    mutate(mean_apoptosis = mean(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    select(-dataset, -value) %>% 
    distinct(ID, .keep_all = TRUE)
  
  thisDF <- left_join(Ferroptosis_DF, Apoptosis_DF, by = "ID") %>% 
    na.omit() %>% 
    mutate(mean_difference = mean_ferroptosis - mean_apoptosis)
  
  Ferroptosis_consensus_keep[[n]] <- thisDF %>% slice_max(mean_difference, n = 400) %>% pull(ID)
  Apoptosis_consensus_keep[[n]] <- thisDF %>% slice_min(mean_difference, n = 400) %>% pull(ID)
}

#Create the universal set - this is a list of all genes, ranked by their specific ferroptosis score
#as defined above.
Ferroptosis_signature <- purrr::reduce(Ferroptosis_consensus_keep, c)
Ferroptosis_signature <- as_tibble(Ferroptosis_signature) %>% group_by(value) %>% summarise(ferroptosis_count = n()) 
Apoptosis_signature <- purrr::reduce(Apoptosis_consensus_keep, c)
Apoptosis_signature <- as_tibble(Apoptosis_signature) %>% group_by(value) %>% summarise(apoptosis_count = n()) 

Universal_set <- full_join(Ferroptosis_signature, Apoptosis_signature, by = "value") %>% 
  replace_na(replace = list(ferroptosis_count = 0, apoptosis_count = 0)) %>% 
  mutate(difference = ferroptosis_count - apoptosis_count) 

#plot the universal set (the top 80 genes in the ranked list)
thisFig2 <- Universal_set %>% 
  left_join(DEG_24h, by = c("value" = "genes")) %>% 
  left_join(DEG_6h, by = c("value" = "genes")) %>% 
  select(value, difference, 
         "JB2_24h_exp" = JB2_24h, "JB3_24h_exp" = JB3_24h, "JB6_24h_exp" = JB6_24h, 
         "JB2_24h_fdr" = JB2_24h_fdr, "JB3_24h_fdr" = JB3_24h_fdr, "JB6_24h_fdr" = JB6_24h_fdr,
         "JB2_6h_exp" = JB2_6h, "JB3_6h_exp" = JB3_6h, "JB6_6h_exp" = JB6_6h, 
         "JB2_6h_fdr" = JB2_6h_fdr, "JB3_6h_fdr" = JB3_6h_fdr, "JB6_6h_fdr" = JB6_6h_fdr) %>% 
  na.omit() %>% 
  mutate(JB2_greater = JB2_24h_exp > 0 & JB2_24h_exp > JB6_24h_exp) %>% 
  mutate(JB6_greater = JB6_24h_exp > 0 & JB6_24h_exp > JB2_24h_exp)

thisFig_ferro <- thisFig2 %>% 
  arrange(desc(difference)) %>% 
  slice_max(difference, n = 80) %>% 
  pivot_longer(-c(value, difference, JB2_greater, JB6_greater), names_to = c("treat", "time", ".value"), names_sep = "_") %>% 
  na.omit() %>% 
  mutate(value = fct_inorder(value)) %>% 
  mutate(fdr_text = ifelse(abs(fdr) == 1, "*", NA))

ggplot(thisFig_ferro, aes(x = treat, y = value, fill = exp)) +
  geom_tile(color = "gray") +
  geom_text(aes(label = fdr_text)) +
  geom_point(data = thisFig_ferro %>% filter(JB2_greater), x = -5.1, color = "red") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = "Treatments", y = "Genes (ordered by best ferroptosis indicator)", fill = "Log2 FC") +
  facet_wrap(~time) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 10, color = "black"))

ggsave("../results/Figure S6D.tiff", dpi = 300, height = 10, width = 5)

#--------------------------------------------------------------------------------------------
# Create the selective ferroptosis biomarkers set
#--------------------------------------------------------------------------------------------

#This is a mostly unbiased process to generate the best biomarkers from the GGS.
#We hypothesize, that the best biomarkers will be a part of the universal set, created above, and
#will be elevated in JB2 vs. JB6. We wanted to identify how the universal set genes which follow those rules behave 
#in our RNAseq, and use that info to build a set of biomarkers in an unbiased manner (without training on the 36 datasets)

#First, determine the universal set genes that are JB2 > JB6:
cutoff <- 100
Ferroptosis_finetune <- thisFig2 %>% slice_max(difference, n = cutoff) %>% filter(JB2_greater) %>% 
  filter(JB2_24h_fdr == 1) %>% pull(value)

#The folowing plot will show how the above 18 genes (`Ferroptosis_finetune`) behave in our RNAseq.
#The general rule for these genes is positive gradient in 6h (JB6>JB2) and negative in 24h (JB2>JB6)
#Notice, the aes_x and aes_y are minus the gradient - to look at JB2 vs. JB6 (rather than JB6 vs. JB2)
Inverse_gradient_t <- as_tibble(efit_24h[["t"]], rownames = "genes") %>% 
  left_join(as_tibble(efit_6h[["t"]], rownames = "genes"), by = "genes") %>% 
  mutate(in_selective = genes %in% Ferroptosis_finetune, .after = genes) 

ggplot(Inverse_gradient_t, aes(x = -JB6_to_JB2_24h, y = -JB6_to_JB2_6h)) +
  geom_point(data = Inverse_gradient_t %>% filter(!in_selective), color = "gray", alpha = 0.8) +
  geom_point(data = Inverse_gradient_t %>% filter(in_selective), color = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(data = Inverse_gradient_t %>% filter(in_selective), aes(label = genes)) +
  stat_cor() +
  labs(x = "JB2 to JB6 gradient, 24h", y = "JB2 to JB6 gradient, 6h") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.2, 0.8))

ggsave("../results/Supplementary text figure 1A.tiff", dpi = 300, height = 6, width = 6)

#From the "Inverse_gradient_t" made above, try to build a decision tree to see how the selective signature 
#might be derived unbiased from the JB2/JB6 data. To generate this, randomly selected genes are used, creating different trees each time.
this_data <- Inverse_gradient_t %>% 
  group_by(in_selective) %>% 
  slice_sample(n = 100, replace = TRUE) %>% 
  na.omit() %>% 
  select(-genes)
dt1 <- rpart(in_selective ~ ., data = this_data, method = "class", control = rpart.control(cp = 0))
rpart.plot(dt1, type = 5, branch.lwd = 7)

#using these decision trees, we generate the selective_gradient set of genes, which include all genes from our RNA-seq,
#that pass the following criteria on their gradient with the following cutoffs, determined based on the desicion trees.
selective_gradient <- Inverse_gradient_t %>% 
  filter(JB2_24h > 2.4) %>%           #first filter
  filter(JB6_to_JB2_24h < -0.7) %>%   #second filter
  filter(B6_to_B2_6h > 1) %>%         #third filter
  #here the selective genes must be part of the original GGS (created using cutoff = 300, to be as broad as possible):
  filter(genes %in% c(group_JB2U_C_300$genes, group_JB2D_C_300$genes, group_JB6U_C_300$genes, group_JB6D_C_300$genes))

#We are left with 32 genes, selected from the GGS (cutoff = 300 version).
#We rank the set by three different ways, and validate in the 36 FINs and AINs datasets by GSVA. 
#this analysis is compared to a list of 33 genes upregulated by erastin in HT-1080 cells, 
#taken from a paper by Dixon et al (2014) (dataset F1 in our ferroptosis_combined dataset):
Dixon_33 <- Ferroptosis_combined %>% slice_max(logFC_F1, n = 36) %>% filter(!str_detect(ID, "LOC100")) %>% pull(ID)
cutoffs <- seq(10, 32, 1)   #fine-tune the signature sizes from 10 to 32
take <- c("^Rank_t")        #the analysis will be done on the ranked t-statistics.

SFB_GSVA_results <- list()
for(n in seq_along(cutoffs)) {
  cutoff <- cutoffs[n]
  print(paste0("working on: ", n, ", cutoff: ", cutoff))
  gene_sets <- list("Dixon_33" = Dixon_33,
                    "universal_set" = Universal_set %>% slice_max(difference, n = cutoff) %>% pull(value),
                    "SFB_sorted_by_JB6vJB2_24h" = selective_gradient %>% slice_min(JB6_to_JB2_24h, n = cutoff) %>% pull(genes),
                    "SFB_sorted_by_JB6vJB2_6h" = selective_gradient %>% slice_max(JB6_to_JB2_6h, n = cutoff) %>% pull(genes),
                    "SFB_sorted_by_B6vB2_6h" = selective_gradient %>% slice_max(B6_to_B2_6h, n = cutoff) %>% pull(genes))
  
  Datasets_combined <- Ferroptosis_combined %>% 
    full_join(Apoptosis_combined, by = "ID") %>% 
    full_join(combinations, by = "ID") %>%
    select(ID, matches(take)) %>% 
    rename_with(.cols = -ID, .fn = ~str_replace(.x, paste0(take, "_"), ""))
    
    #perform GSVA analysis on the combined dataset
    this_GSVA <- Datasets_combined %>% 
      distinct(ID, .keep_all = TRUE) %>% 
      column_to_rownames(var = "ID") %>% 
      as.matrix()
    gsva_result <- gsva(this_GSVA, gene_sets, method = "gsva", verbose = TRUE)
    result_GSVA <- gsva_result %>% as_tibble(rownames = "pathway") %>% 
      gather(-pathway, key = "dataset", value = "score") %>% 
      mutate(method = paste0("GSVA_", str_remove(take, "\\^"))) %>% 
      mutate(cutoff = cutoff)
    
    SFB_GSVA_results[[n]]  <- bind_rows(result_GSVA) %>% 
      mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                              str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                              str_detect(dataset, "JB") ~ "JB combinations", 
                              TRUE ~ "NA"))
}

SFB_GSVA_results_tib <- purrr::reduce(SFB_GSVA_results, bind_rows)

#calculate the ROC-AUCs for each group
AUCs <- SFB_GSVA_results_tib %>% 
  filter(str_detect(method, "GSVA")) %>% 
  filter(!str_detect(dataset, "JB")) %>% 
  filter(score != Inf) %>% 
  select(-dataset) %>% 
  group_by(pathway, method, cutoff) %>% 
  dplyr::summarise(auc_total = as.numeric(pROC::auc(type, score)))

#draw plot of the AUCs for each geneset (unpublished). This plot was used to determine the best cutoff for the SFB size = 16 genes
ggplot(AUCs, aes(x = cutoff, y = auc_total, color = pathway)) +
  geom_line(size = 2) +
  scale_x_continuous(breaks = seq(10, 32, 1)) +
  labs(x = "Cutoff", y = "AUC") +
  facet_wrap(~method) +
  theme_gray(base_size = 16)

#quantify the p-values of these AUCs by 1000 random permutations of GSVA using 20 random genes
#create 1000 sets of 20 randomly selected genes each
set.seed(123)
iterations <- 1000
cutoff <- 20
take <- c("^Rank_t")
random_sets <- list()
for(n in 1:iterations) {
  random_set <- Universal_set %>% slice_sample(n = cutoff) %>% pull(value)        
  cat(paste0("working on: ", n, " / ", iterations, ": " , paste0(random_set, collapse = " ") ,"\r"))
  random_sets[n] <- list("random" = random_set)
}
  
Datasets_combined <- Ferroptosis_combined %>% 
  full_join(Apoptosis_combined, by = "ID") %>% 
  full_join(combinations, by = "ID") %>%
  select(ID, matches(take)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, paste0(take, "_"), ""))
  
#perform GSVA analysis of the 1000 random genesets on the combined dataset
this_GSVA_random <- Datasets_combined %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix()
  
gsva_result_random <- gsva(this_GSVA_random, random_sets, method = "gsva", verbose = TRUE)
  
result_GSVA_random <- gsva_result_random %>% 
  as_tibble(rownames = "pathway") %>% 
  gather(-pathway, key = "dataset", value = "score") %>% 
  mutate(method = paste0("GSVA_", str_remove(take, "\\^"))) %>% 
  mutate(cutoff = cutoff) %>% 
  mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                          str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                          str_detect(dataset, "JB") ~ "JB combinations", 
                          TRUE ~ "NA"))

#calculate AUCs for the random genesets
random_AUCs <- result_GSVA_random %>% 
  filter(str_detect(method, "GSVA")) %>% 
  filter(!str_detect(dataset, "JB")) %>% 
  filter(score != Inf) %>% 
  select(-dataset) %>% 
  group_by(pathway, method, cutoff) %>% 
  dplyr::summarise(auc_total = as.numeric(pROC::auc(type, score))) %>% 
  mutate(pathway = "random")
  
#p-value calculations of our AUCs
random_AUCs_c <- random_AUCs %>% arrange(auc_total) %>% pull(auc_total)
AUCs_p <- AUCs %>% 
  mutate(rank_total = map_dbl(auc_total, ~sum(random_AUCs_c > .x))) %>% 
  mutate(pval_total = rank_total / length(random_AUCs_c))

#plot the best ferroptosis selective signature (Figure 6F)
cutoff_value <- 16 #the optimal size of genes for the selective_GGS set, as determined above
thisFig5 <- SFB_GSVA_results_tib %>% 
  left_join(AUCs_p) %>% 
  left_join(Lineages, by = "dataset") %>% 
  filter(method == "GSVA_Rank_t") %>% 
  filter(cutoff == cutoff_value) %>% 
  replace_na(list(lineage = "TNBC")) %>% 
  filter(score != Inf) %>% 
  filter(pathway %in% c("universal_set", "SFB_sorted_by_JB6vJB2_6h", "Dixon_33")) %>% 
  mutate(pathway = case_when(pathway == "Dixon_33" ~ "Dixon set",
                             pathway == "universal_set" ~ "Universal set",
                             pathway == "SFB_sorted_by_JB6vJB2_6h" ~ "Ferroptosis \nBiomarkers set")) %>% 
  filter(!str_detect(dataset, "JB")) %>% 
  ungroup() %>% 
  mutate(pathway = factor(pathway, levels = c("Ferroptosis \nBiomarkers set", "Universal set", "Dixon set"), ordered = TRUE))

stat.test <- thisFig5 %>% 
  group_by(pathway, method, cutoff) %>% 
  rstatix::t_test(score ~ type, p.adjust.method = "none", var.equal = TRUE, 
                  comparisons = list(c("Ferroptosis", "Apoptosis"))) %>%
  add_xy_position() %>% 
  add_significance() 

this_table <- thisFig5 %>% 
  distinct(pathway, method, cutoff, .keep_all = TRUE) %>% 
  mutate(auc_total = paste0("ROC-AUC = ", round(auc_total, 3))) %>% 
  mutate(pval_total = paste0("AUC p-value = ", round(pval_total, 3))) %>% 
  mutate(across(contains("pval"), ~ifelse(.x == "AUC p-value = 0", "AUC p-value < 0.001", .x))) %>% 
  select(pathway, auc_total, pval_total) %>% 
  gather(-pathway, key = "metric", value = "value") %>% 
  spread(key = "pathway", value = "value") %>% 
  mutate(metric = factor(metric, levels = c("auc_total", "pval_total"))) %>% 
  arrange(metric) %>% 
  select(-metric) %>% 
  rename_with(.cols = everything(), .fn = ~str_wrap(.x, width = 20))

p1 <- ggtexttable(this_table[1:2,], rows = NULL, 
                  theme = ttheme(
                    colnames.style = colnames_style(color = "black", fill = "white", size = 15),
                    tbody.style = tbody_style(color = "black", 
                                              fill = c("gray90", "gray90", "gray80", "gray80")))) +
  theme(plot.margin=unit(c(1,0,0,0.8), "cm"))

p2 <- ggplot(thisFig5 %>% filter(cutoff == cutoff_value, pathway != "Selective Ferroptosis set"), 
             aes(x = type, y = score)) +
  geom_point(aes(color = str_detect(lineage, "TNBC"))) +
  geom_text_repel(data = thisFig5 %>% filter(cutoff == cutoff_value, pathway != "Selective Ferroptosis set", lineage == "TNBC"), 
                  aes(label = dataset), direction = "y", nudge_x = 0.5, size = 4.5) +
  stat_pvalue_manual(stat.test %>% filter(cutoff == cutoff_value, pathway != "Selective Ferroptosis set"), 
                     hide.ns = TRUE, size = 6) +
  scale_color_manual(values = c("grey", "red")) +
  facet_wrap(~pathway, labeller = label_wrap_gen(width = 20)) +
  labs(x = "", y = "Enrichment score", title = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", 
        plot.margin=unit(c(-1,1,0.5,0), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 45, hjust = 0.8, vjust = 0.75, color = "black"))

ggarrange(p1, p2, ncol = 1, heights = c(1, 2))
ggsave("../results/Figure 6F.tiff", dpi = 300, height = 6, width = 6)

#plot a heatmap of the genes in the SFB set (figure 6E)
SFB <- selective_gradient %>% slice_max(JB6_to_JB2_6h, n = 16) %>% pull(genes)

JB_genes <- DEG_24h %>% 
  left_join(DEG_6h, by = "genes") %>% 
  filter(genes %in% SFB) %>% 
  mutate(genes = factor(genes, levels = SFB, ordered = TRUE)) %>% 
  arrange(genes) %>% 
  select(genes, "JB2_24h_exp" = JB2_24h, "JB3_24h_exp" = JB3_24h, "JB6_24h_exp" = JB6_24h, "Gradient_24h_exp" = JB6_to_JB2_24h,
         "JB2_24h_fdr" = JB2_24h_fdr, "JB3_24h_fdr" = JB3_24h_fdr, "JB6_24h_fdr" = JB6_24h_fdr, 
         "Gradient_24h_fdr" = JB6_to_JB2_24h_fdr,
         "JB2_6h_exp" = JB2_6h, "JB3_6h_exp" = JB3_6h, "JB6_6h_exp" = JB6_6h, "Gradient_6h_exp" = JB6_to_JB2_6h,
         "JB2_6h_fdr" = JB2_6h_fdr, "JB3_6h_fdr" = JB3_6h_fdr, "JB6_6h_fdr" = JB6_6h_fdr, 
         "Gradient_6h_fdr" = JB6_to_JB2_6h_fdr) %>% 
  mutate(Gradient_24h_exp = -Gradient_24h_exp) %>%  #NOTE: gradient is JB2 vs JB6 (i.e. minus the previous gradient)
  mutate(Gradient_6h_exp = -Gradient_6h_exp) %>%    
  na.omit() %>% 
  pivot_longer(-genes, names_to = c("treat", "time", ".value"), names_sep = "_") %>% 
  na.omit() %>% 
  mutate(fdr_text = ifelse(abs(fdr) == 1, "*", NA)) %>% 
  mutate(treat = case_when(treat == "Gradient" ~ "JB2-JB6 Gradient",
                           TRUE ~ str_c(treat, " vs DMSO")))

gene_sorter <- JB_genes %>% 
  filter(treat == "JB2 vs DMSO", time == "24h") %>% 
  arrange(desc(exp))
JB_genes$genes <- factor(JB_genes$genes, levels = gene_sorter$genes, ordered = TRUE)
JB_genes$time <- factor(JB_genes$time, levels = c("6h", "24h"))
JB_genes$treat <- factor(JB_genes$treat, levels = c("JB2-JB6 Gradient", "JB2 vs DMSO", "JB3 vs DMSO", "JB6 vs DMSO"))

ggplot(JB_genes, aes(x = treat, y = genes, fill = exp)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = fdr_text)) +
  geom_rect(xmin = 1.5, xmax = 4.5, ymin = 0.5, ymax = Inf, fill = NA, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = "", y = "Genes", fill = "Log2 FC") +
  facet_wrap(~time) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right", 
        axis.text = element_text(color = "black", size = 14),
        strip.text = element_text(color = "black", size = 14),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("../results/Figure 6E.tiff", dpi = 300, height = 5, width = 6)

#Venn diagram for the universal set, FSB and Dixon sets
#NOTE - The Venn diagram shown in Fig. 6G was made in powerpoint, based on the Venn diagram created here.
Ferroptosis_list <- list("Universal set" = Universal_set %>% slice_max(difference, n = 20) %>% pull(value), 
                         "Dixon set" = Dixon_33,
                         "Selective Ferroptosis Biomarkers set" = SFB)
Ferro_partitions <- get.venn.partitions(Ferroptosis_list) 

tiff(filename = "../results/Figure 6G.tiff", width = 800, height = 800)
plot(euler(Ferroptosis_list), quantities = list(cex = 5), lty = 1, lwd = 2, labels = TRUE,
     fill = c("salmon", "lightblue", "yellowgreen"))
dev.off()

Genn_names_per_partition <- Ferro_partitions %>% 
  rowwise() %>% 
  mutate(..values.. = paste0(..values.., collapse = " ")) %>% 
  write_csv("../results/Figure 6G gene lists.csv")

#calculate the mean fold change of expression of the 33 Dixon genes in the 36 datasets
Dixon_mean <- Ferroptosis_combined %>% 
  full_join(Apoptosis_combined) %>% 
  select(ID, contains("logFC")) %>% 
  filter(ID %in% Dixon_33) %>% 
  summarise(across(-ID, mean, na.rm = TRUE)) %>% 
  gather(key = "dataset", value = "score") %>% 
  mutate(dataset = str_remove(dataset, "logFC_")) %>% 
  arrange(`score`) %>% 
  mutate(dataset = fct_inorder(dataset)) %>% 
  mutate(type = ifelse(str_detect(dataset, "F"), "Ferroptosis", "Apoptosis")) 

ggplot(Dixon_mean, aes(x = reorder(dataset, -score) , y = score, fill = type)) +
    geom_col() +
    geom_text(aes(label = dataset), y = 0.05, angle = 90) +
    scale_fill_manual(values = c("#79A7D0", "#70AD47")) +
    labs(x = "", y = "Median fold change (log2) of Dixon genes") +
    theme_bw(base_size = 12) +
    theme(legend.position = "non",
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          axis.text.x = element_blank())

ggsave("../results/Figure S6E.tiff", dpi = 300, height = 4, width = 6)

#prepare gene set data for supplementary table S4 - this code will compare the gene sets we discovered (GGS, universal, SFB)
#to known ferroptosis related gene sets - Dixon, FerrDB, Wikipathways,
#and will calculate the intersection betwween all gene sets.
other_genesets <- gmt2list("../data/genesets used in the paper.gmt")

Genesets_list <- list("Gradient genes" = rbind(group_JB2U_C_100, group_JB2D_C_100, group_JB6U_C_100, group_JB6D_C_100) %>% 
                        distinct(genes) %>% pull(genes),
                      "Universal top 100" = Universal_set %>% slice_max(difference, n = 100) %>% pull(value), 
                      "Dixon set" = Dixon_33,
                      "Selective set SFB" = selective_gradient %>% slice_max(JB6_to_JB2_6h, n = 16) %>% pull(genes),
                      "FerrDB Drivers" = other_genesets[["FerrDB Drivers"]],
                      "FerrDB Suppresor" = other_genesets[["FerrDB Suppresor"]],
                      "FerrDB Marker" = other_genesets[["FerrDB Marker"]],
                      "WP_Ferroptosis" = other_genesets[["WP_Ferroptosis"]])

analyse_geneset <- tribble(~"Group 1", ~"Group 2", ~"Intersect", ~"n")
for(i in 1:length(Genesets_list)) {
  for(j in 1:length(Genesets_list)) {
    group1 <- names(Genesets_list)[i]
    group2 <- names(Genesets_list)[j]
    intersect <- intersect(Genesets_list[[i]], Genesets_list[[j]])
    n <- length(intersect)
    analyse_geneset <- bind_rows(analyse_geneset, tibble("Group 1" = group1, "Group 2" = group2, 
                                                         "Intersect" = paste0(intersect, collapse = ", "),
                                                         "n" = n))
  }
}

thisFig_genesets <- analyse_geneset %>% 
  select(-Intersect) %>% 
  spread(key = "Group 2", value = "n") %>% 
  write_csv("../results/prepared genesets for supplemental table.csv")
