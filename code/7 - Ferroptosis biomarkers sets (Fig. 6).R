#This script will create the ferroptosis vs. apoptosis specific biomarkers set, as shown in figure 6 and S7

#helper function to draw green-blue backgrounds in plots
make_gradient <- function(deg = 90, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(data = rep(seq(0, 1, length.out = n) * cos(rad), n), byrow = TRUE, ncol = n) +
    matrix(data = rep(seq(0, 1, length.out = n) * sin(rad), n), byrow = FALSE, ncol = n)
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(image = mat, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = TRUE)
}
cols <- colorRampPalette(c("#79A7D0", "#A2C4E3", "#BCDAA7", "#70AD47"))(100)

#--------------------------------------------------------------------------------------------
# Create the dataset-derived biomarkers set
#--------------------------------------------------------------------------------------------

#the dataset-derived biomarkers set are the genes most differentiating between the 19 ferroptosis inducers (FIN) datasets
#and 26 apoptosis inducers (AINs) datasets. This set is biased to those datasets. 
#In order to reduce the bias, the geneset was determined by randomly selecting 12 FINs and 12 AINs datasets 
#(with replacement, so the number of datasets can be less than 10) and extract the top 150 differential genes.
#This was repeated 1000 times. The genes were ranked according to their frequency of appearing among those 1000 iterations.
#See supplemental text for further details.

#load the ferrotposis and apoptosis inducers datasets
#see supplemental table S3 for list of datasets used. The datasets are annotated F1-F19 (for ferroptosis inducers)
#and A1-A26 (for apoptosis inducers). These files contain 4 column for each dataset:
#logFC_X - log2 fold change for each gene of the specific inducer vs. its control
#t_X - the t-statistics for the above FC, calculated by Limma.
#Rank_FC - the fold change data, ranked on a scale from 0 to 1
#Rank_t - the t-statistic data, ranked on a scale from 0 to 1
#note: each dataset was analyzed separately, these tables (already analysed and ranked) are provided for reproducability.
#For more details on these datasets please see the supplemental text.

#load the ferroptosis/apoptosis datasets
Ferroptosis1 <- read_csv("../Data/Datasets/Ferroptosis dataset - Dixon (erastin in HT1080).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F1"))
Ferroptosis2 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE104462 (erastin in HepG2).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F2"))
Ferroptosis3 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE121689 (Erastin in RCC4).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F3"))
Ferroptosis4 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE112384 (withaferin A in neuroblastoma).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F4"))
Ferroptosis5 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE131444 (erastin in MEF).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F5"))
Ferroptosis6 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE126868 (ferroptocide in HT29).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F6"))
Ferroptosis7 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE142591 (NaIO3 in Arpe-19).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F7"))
Ferroptosis8 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE96760 (ML210 in H-STS).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F8"))
Ferroptosis9 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE96760 (Erastin in H-STS).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F9"))
Ferroptosis10 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE119628 (SLC7A11 KO in mice).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F10"))
Ferroptosis11 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE147625 (shGPX4 in BeWo).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F11"))
Ferroptosis12 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE160338 (GPX4 KO in Treg).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F12"))
Ferroptosis13 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE182638 (RSL3 in MM1).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F13"))
Ferroptosis14 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE197104 (hemin in OPC).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F14"))
Ferroptosis15 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE193295 (sevoflurane in U251).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F15"))
Ferroptosis16 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE162069 (ML162 in MB231).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F16"))
Ferroptosis17 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE163399 (NC06 in J774M).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F17"))
Ferroptosis18 <- read_csv("../Data/Datasets/Ferroptosis dataset - GSE154425 (Erastin in HCC38).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F18"))
Ferroptosis19 <- read_csv("../Data/Datasets/Ferroptosis dataset - NASH-liver - in-vivo ferroptosis model.csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_F19"))

Apoptosis1 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE70115 (EC-70124 in MDA231).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A1"))
Apoptosis2 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE161803 (Etoposide in Capan2).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A2"))
Apoptosis3 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE207382 (AZD5153 and Entospletinib in SU-DHL-4).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A3"))  
Apoptosis4 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE124715 (5FU in 231).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A4"))
Apoptosis5 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE129012 (AZD4573 in MCF7).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A5"))
Apoptosis6 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE140758 (CGM097 in MCF7).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A6"))
Apoptosis7 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE121093 (VitC_in_MDAMB231).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A7"))
Apoptosis8 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE124597 (DOX+RM in MCF7).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A8"))
Apoptosis9 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE95375 (BETd246 in MDA468).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A9"))
Apoptosis10 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE154754 (Paeonol in HCT116).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A10"))
Apoptosis11 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE48056 (Roscovitin in LNCaP).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A11"))
Apoptosis12 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE48056 (Bortezomib in LNCaP).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A12"))
Apoptosis13 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE171167 (Elaiophylin in KYSE450).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A13"))
Apoptosis14 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE193660 (NaHS in HepG2).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A14"))
Apoptosis15 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE197320 (LY500307 in SKOV3).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A15")) 
Apoptosis16 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE179079 (TP472 in A375).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A16"))
Apoptosis17 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE201959 (Icaritin  in SNT8).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A17"))
Apoptosis18 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE200385 (IT848 in MM1S).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A18"))
Apoptosis19 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE217070 (Erianin in A375).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A19")) 
Apoptosis20 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE199088 (SMIP34 in ZR-75).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A20")) 
Apoptosis21 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE205596 (2c in SKUT1).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A21"))   
Apoptosis22 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE137467 (Docetaxel in MDA-MB-231 tumors).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A22"))
Apoptosis23 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE163909 (epigambogic acid A in PC9).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A23"))
Apoptosis24 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE200748 (cisplatin in HeLa).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A24"))
Apoptosis25 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE203203 (DSF in THP-1).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A25"))
Apoptosis26 <- read_csv("../Data/Datasets/Apoptosis dataset - GSE143713 (EPZ-6438 in SK-N-SH).csv") %>% 
  rename_with(.cols = -ID, .fn = ~paste0(.x, "_A26"))

Ferroptosis_combined <- Ferroptosis1 %>% 
  full_join(Ferroptosis2) %>% full_join(Ferroptosis3) %>% full_join(Ferroptosis4) %>% 
  full_join(Ferroptosis5) %>% full_join(Ferroptosis6) %>% full_join(Ferroptosis7) %>% 
  full_join(Ferroptosis8) %>% full_join(Ferroptosis9) %>% full_join(Ferroptosis10) %>%   
  full_join(Ferroptosis11) %>% full_join(Ferroptosis12) %>% full_join(Ferroptosis13) %>% 
  full_join(Ferroptosis14) %>% full_join(Ferroptosis15) %>% full_join(Ferroptosis16) %>% 
  full_join(Ferroptosis17) %>% full_join(Ferroptosis18) %>% full_join(Ferroptosis19)

Apoptosis_combined <- Apoptosis1 %>% full_join(Apoptosis2) %>% 
  full_join(Apoptosis3) %>% full_join(Apoptosis4) %>%
  full_join(Apoptosis5) %>% full_join(Apoptosis6) %>% full_join(Apoptosis7) %>% 
  full_join(Apoptosis8) %>% full_join(Apoptosis9) %>% full_join(Apoptosis10) %>% 
  full_join(Apoptosis11) %>% full_join(Apoptosis12) %>% full_join(Apoptosis13) %>% 
  full_join(Apoptosis14) %>% full_join(Apoptosis15) %>% full_join(Apoptosis16) %>% 
  full_join(Apoptosis17) %>% full_join(Apoptosis18) %>% full_join(Apoptosis19) %>% 
  full_join(Apoptosis20) %>% full_join(Apoptosis21) %>% full_join(Apoptosis22) %>% 
  full_join(Apoptosis23) %>% full_join(Apoptosis24) %>% full_join(Apoptosis25) %>% 
  full_join(Apoptosis26)

#For further validation, 5 additional datasets were created in the lab. Analysis was perfomred similarly to the public datasets above,
#and is given here in the same format:
Validation_sets <- read_csv("../data/Validation datasets.csv")
Datasets_combined <- Ferroptosis_combined %>% 
  full_join(Apoptosis_combined, by = "ID") %>% 
  full_join(Validation_sets, by = "ID")

#Add the JB2, JB3, JB6 combinations from this project - 
#this will create similar table as above, for the JB2-JB3-JB6 combinations vs. DMSO
combinations_t <- efit_24h[["t"]] %>% 
  as_tibble(rownames = "ID") %>% 
  mutate(JB2_to_JB6_24h = -JB6_to_JB2_24h) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("t_", .x)) %>% 
  dplyr::select(ID, t_JB2_24h:t_JB6_24h, t_JB6_to_JB2_24h, t_JB2_to_JB6_24h)

combinations_FC <- DEG_24h %>% 
  select("ID" = "genes", JB2_24h:JB6_24h, JB6_to_JB2_24h) %>% 
  mutate(JB2_to_JB6_24h = -JB6_to_JB2_24h) %>% 
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

#Develop the dataset-derived biomarkers set by bootstrapping aggregation
set.seed(123)
take_possibilites <- "^logFC"     #for this analysis, the actual logFC values of each gene will be used 
                                  #(since logFC is important for biomarkers, which needs to be identified later by qRT-PCR)
iterations <- 1000
cutoff <- 150                     #how many genes to take from the top-listed genes in each iteration?
dataset_selected <- 12            #how many FINs and AINs datasets (with repeats) to select in each iteration?

Datasets_combined <- full_join(Ferroptosis_combined, Apoptosis_combined, by = "ID") %>% 
  select(-contains("Normed"))
FvA_mean <- AvF_mean <- list()

#run time of this for loop - ~90 minutes
##option: To save time, run the following line in order to load the objects created by this for loop:
#load("../data/Datasets derived biomarkers analysis object.RData")
for(n in 1:iterations) {
  cat(paste0("Iteration: ", n, " / ", iterations, "\r"))
  ferro_dataset_select <- unique(str_c("F", floor(runif(n = dataset_selected, min = 1, max = 20)), "$"))
  apo_dataset_select <- unique(str_c("A", floor(runif(n = dataset_selected, min = 1, max = 27)), "$"))
  
  temp_DF <- Datasets_combined %>% 
    select(ID, matches(ferro_dataset_select), matches(apo_dataset_select)) %>% 
    filter(rowSums(is.na(.)) < 13) %>% 
    pivot_longer(-ID, names_to = c("stat", "dataset"), values_to = "value", names_pattern = "(.*)_(.*)$") %>% 
    mutate(type = ifelse(str_detect(dataset, "F"), "Ferroptosis", "Apoptosis")) %>% 
    group_by(ID, type, stat) %>% 
    summarise(mean = mean(value, na.rm = TRUE), median = median(value, na.rm = TRUE), .groups = "drop") %>% 
    pivot_wider(names_from = type, values_from = c("mean", "median")) %>% 
    mutate(iteration = n, .before = "ID") %>% 
    mutate(mean_diff = mean_Ferroptosis - mean_Apoptosis,
           median_diff = median_Ferroptosis - median_Apoptosis) %>% 
    ungroup()
  
  FvA_mean[[n]] <- temp_DF %>% group_by(stat) %>% slice_max(mean_diff, n = cutoff)
  AvF_mean[[n]] <- temp_DF %>% group_by(stat) %>% slice_min(mean_diff, n = cutoff)
}

FvA_mean_tb <- purrr::reduce(FvA_mean, bind_rows) %>% group_by(stat, ID) %>% 
  summarise(FvA_mean_count = n())
AvF_mean_tb <- purrr::reduce(AvF_mean, bind_rows) %>% group_by(stat, ID) %>% 
  summarise(AvF_mean_count = n())

Datasets_derived_sig <- full_join(FvA_mean_tb, AvF_mean_tb, by = c("stat", "ID")) %>%  
  replace_na(replace = list(FvA_mean_count = 0, AvF_mean_count = 0)) %>% 
  mutate(mean_difference = FvA_mean_count - AvF_mean_count) 

#1D plot of the ferroptosis vs. apoptosis scores from the above analysis (Fig. S7B)
set.seed(123)
g <- make_gradient(deg = 0, n = 500, cols = cols)
y <- runif(n = dim(Datasets_derived_sig)[1])
thisFig <- bind_cols(Datasets_derived_sig, "y" = y) %>% filter(stat == "logFC")

ggplot(thisFig, aes(x = mean_difference, y = y)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -0.15, ymax = 1.05) + 
  geom_point(color = "gray20", alpha = 0.25) +
  geom_text_repel(data = thisFig %>% slice_max(mean_difference, n = 15), aes(label = ID), color = "black", min.segment.length = 0) +
  geom_text_repel(data = thisFig %>% slice_min(mean_difference, n = 15), aes(label = ID), color = "black", min.segment.length = 0) +
  annotate(geom = "text", x = 550, y = 1.12, label = "Ferroptosis", color = "darkgreen", size = 5, fontface = "bold") +
  annotate(geom = "text", x = -550, y = 1.12, label = "Apoptosis", color = "blue", size = 5, fontface = "bold") +
  scale_x_continuous(breaks = seq(-1000, 1000, 100), limits = c(-1000, 1000)) +
  labs(x = "Ferroptosis - Apoptosis Score", y = "") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.text.x = element_text(color = "black"))

ggsave("../results/Figure S7B - Dataset derived set visualized in 1D.tiff", dpi = 300, height = 2.5, width = 9)

#Score the dataset derived set in the FINs and AINs datasets using ssGSEA, for different set sizes
#The genes are ranked by their score (see above), and cutoffs are set to generate sets with 15-400 genes
set.seed(123)
cutoffs <- c(seq(15, 19, 1), seq(20, 90, 10), seq(100, 400, 20))
take <- "^logFC"
result_combined_list <- list()
for(n in seq_along(cutoffs)) {
  cutoff <- cutoffs[n]
  signatures <- list("Datasets_derived_Ferro_mean_logFC" = Datasets_derived_sig %>% 
                       filter(stat == "logFC") %>% slice_max(mean_difference, n = cutoff) %>% pull(ID), 
                     "Datasets_derived_Apop_mean_logFC" = Datasets_derived_sig %>% 
                       filter(stat == "logFC") %>% slice_min(mean_difference, n = cutoff) %>% pull(ID))
  print(paste0("working on: ", n, ", cutoff: ", cutoff))
  Datasets_combined <- Ferroptosis_combined %>% 
    full_join(Apoptosis_combined, by = "ID") %>%
    full_join(Validation_sets, by = "ID") %>%
    full_join(combinations, by = "ID") %>%
    filter(rowSums(is.na(.)) < 18) %>% 
    select(ID, matches(take)) %>% 
    rename_with(.cols = -ID, .fn = ~str_replace(.x, paste0(take, "_"), "")) %>% 
    distinct(ID, .keep_all = TRUE) %>% 
    column_to_rownames(var = "ID") %>% 
    as.matrix() %>% 
    na.omit()
  
  ssgsea_result <- gsva(Datasets_combined, signatures, method = "ssgsea", ssgsea.norm = FALSE)
    
  result_ssGSEA <- ssgsea_result %>% t() %>% 
      as_tibble(rownames = "dataset") %>% 
      gather(-dataset, key = "pathway", value = "score") %>% 
      mutate(method = paste0("SSGSEA_", str_remove(take, "\\^"))) %>%
      mutate(cutoff = cutoff)
    
  result_combined_list[[n]] <- result_ssGSEA %>% 
      mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                              str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                              str_detect(dataset, "JB") ~ "JB combinations", 
                              TRUE ~ "Validation"))
}

Dataset_derived_scoring_results <- purrr::reduce(result_combined_list, bind_rows)

#calculate ROC-AUCs for classification between ferroptosis and apoptosis datasets using the scores calculated above
Dataset_derived_AUCs <- Dataset_derived_scoring_results %>% 
  filter(score != Inf) %>% 
  select(-dataset) %>% 
  group_by(pathway, method, cutoff) %>% 
  dplyr::summarise(auc_total = as.numeric(pROC::auc(response = type, predictor = score, levels = c("Ferroptosis", "Apoptosis"))))

#plot the AUCs (Fig. S7C)
thisFig <- Dataset_derived_AUCs %>% 
  ungroup() %>% 
  mutate(pathway = fct_recode(pathway,
                              "Apoptosis genes ranked by logFC" = "Datasets_derived_Apop_mean_logFC",
                              "Ferroptosis genes ranked by logFC" = "Datasets_derived_Ferro_mean_logFC")) 

ggplot(thisFig, aes(x = cutoff, y = auc_total, color = pathway)) +
  geom_point(size = 2) +
  geom_line(size = 1) +
  scale_x_continuous(trans = "log10", breaks = c(15, 30, 60, 120, 240, 400)) +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  scale_color_manual(values = c("blue", "darkgreen")) +
  labs(x = "Number of genes in signature", y = "ROC-AUC of ferroptosis\nto apoptosis classification", color = "",
       title = "Datasets-derived signature") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.4, 0.15), legend.direction = "vertical",
        legend.background = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"))

ggsave("../results/Figure S7C (left) - AUCs of datasets derived biomarkers.tiff", dpi = 300, height = 5, width = 5)

#--------------------------------------------------------------------------------------------
# Create the gradient-derived biomarkers set
#--------------------------------------------------------------------------------------------

#This is an unbiased process to generate the best biomarkers from the GGS (gradient gene set).
#We have noticed that many of the dataset-derived biomarkers are upregulated in JB2 vs. JB6 in 24h,
#and downregulated in JB2 vs. JB6 in 6h. Therefore, we take all the genes from our RNAseq who show this behavior,
#and these are the gradient-derived biomarkers. In total, there are 26 such genes.

hallmark <- gmt2list("../data/MSigDB_Hallmark_2020.gmt") #will be used to score the "Hallmark_apoptosis" set
Gradient_derived <- DEG_24h %>% 
  left_join(DEG_6h, by = "genes") %>%
  as_tibble() %>% 
  filter(JB6_to_JB2_24h_fdr == -1) %>% 
  filter(JB6_to_JB2_6h_fdr == 1) %>% 
  mutate(rank_by_JB2vJB6_24h = min_rank(JB6_to_JB2_24h)) %>% 
  left_join(Datasets_derived_sig %>% ungroup() %>% filter(stat == "logFC") %>% 
            mutate(universal_rank = min_rank(-mean_difference)) %>% select("genes" = ID, universal_rank))

signatures <- list(#our sets - the gradient derived and datasets derived biomarkers
                   "Gradient_derived_15" = Gradient_derived %>% slice_min(rank_by_JB2vJB6_24h, n = 15) %>% pull(genes),
                   "Datasets_derived_100" = Datasets_derived_sig %>% filter(stat == "logFC") %>% 
                        slice_max(mean_difference, n = 100, with_ties = FALSE) %>% pull(ID), 
                   "Datasets_derived_15" = Datasets_derived_sig %>% filter(stat == "logFC") %>% 
                        slice_max(mean_difference, n = 15, with_ties = FALSE) %>% pull(ID),
                   #combining the top ranking biomarkers in both sets (these are the biomarkers validated by qRT-PCR,
                   #except for SACS which was not validated in the qRT-PCR)
                   "optimized_biomarkers" = unique(c(Gradient_derived %>% slice_min(universal_rank, n = 12) %>% pull(genes),
                                                     Datasets_derived_sig %>% filter(stat == "logFC") %>% 
                                                     slice_max(mean_difference, n = 15, with_ties = FALSE) %>% pull(ID))),
                   #compare to other sets - genes upregulated by erastin in HT1080 (DOI: 10.7554/eLife.02523)
                   "F1_set" = Ferroptosis1 %>% slice_max(logFC_F1, n = 36) %>% filter(!str_detect(ID, "LOC100")) %>% pull(ID),
                   #Hallmark apoptosis genes
                   "Hallmark_Apoptosis" = unique(hallmark[["Hallmark Apoptosis"]]))

#prepare datasets for scoring of the signatures enrichment
take <- c("^logFC")
result_combined_list <- list()
Datasets_combined <- Ferroptosis_combined %>% 
  full_join(Apoptosis_combined, by = "ID") %>%
  full_join(Validation_sets, by = "ID") %>%
  full_join(combinations, by = "ID") %>%
  filter(rowSums(is.na(.)) < 18) %>% 
  select(ID, matches(take)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, paste0(take, "_"), "")) %>% 
  distinct(ID, .keep_all = TRUE) 

Datasets_matrix <- Datasets_combined %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix() 

#Add 1000 random genesets, to calculate later the p-values of the AUCs
set.seed(1234)
random_sets <- list()
for(n in 1:1000) {
  random_set <- Datasets_combined %>% slice_sample(n = 20) %>% pull(ID)        
  cat(paste0("working on: ", n, " / 1000 : " , paste0(random_set, collapse = " ") ,"\r"))
  random_sets[n] <- list("random" = random_set)
  names(random_sets)[n] <- str_c("Random_", n)
}
signatures_with_random <- c(signatures, random_sets)

#score the enrichment of all signatures using ssGSEA
ssgsea_result <- gsva(Datasets_matrix %>% na.omit(), signatures_with_random, method = "ssgsea", ssgsea.norm = FALSE) %>% 
  t() %>% 
  as_tibble(rownames = "dataset") %>% 
  gather(-dataset, key = "pathway", value = "score") %>% 
  mutate(method = paste0("SSGSEA_", str_remove(take, "\\^"))) %>% 
  mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                          str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                          str_detect(dataset, "JB") ~ "JB combinations", 
                          TRUE ~ "Validation"))

#calculate ROC-AUCs and p-values for classification between ferroptosis and apoptosis datasets
#(note: p-values might be slightly different than those reported in the paper due to different random seed)
AUCs_all <- ssgsea_result %>% 
  filter(score != Inf) %>% 
  mutate(type = ifelse(type == "Validation", "Ferroptosis", type)) %>% 
  filter(type %in% c("Ferroptosis", "Apoptosis")) %>% 
  select(-dataset) %>% 
  group_by(pathway, method) %>% 
  dplyr::summarise(auc_total = as.numeric(pROC::auc(response = type, predictor = score, levels = c("Ferroptosis", "Apoptosis")))) %>% 
  group_by(method) %>% 
  mutate(rank_total = min_rank(-auc_total)) %>% 
  mutate(pval_total = rank_total / 1000) %>% 
  ungroup() %>% 
  filter(!str_detect(pathway, "Random")) %>% 
  mutate(datasets_taken = "all")

#plot a scatterplot of the gradient-derived and datasets-derived scores for all datasets (Fig. S7D-E)
thisFig <- ssgsea_result %>% 
  left_join(AUCs_all, by = c("pathway", "method")) %>% 
  filter(pathway %in% c("Datasets_derived_15", "Gradient_derived_15")) %>% 
  mutate(dataset = str_remove(dataset, "_24h")) %>% 
  filter(score != Inf) %>% 
  select(pathway, type, dataset, score) %>% 
  pivot_wider(names_from = "pathway", values_from = "score") 

#create contours for the scatterplots, by k-nearest neighbors (see explanation in the `6 - GGS validation (Fig. 5C).R` script)
set.seed(123)
thisDF1 <- thisFig %>% filter(type %in% c("Ferroptosis", "Apoptosis"))
kknn_recipe <- recipe(formula = type ~ Gradient_derived_15 + Datasets_derived_15, data = thisDF1)
kknn_spec <- nearest_neighbor(neighbors = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("kknn") 
kknn_workflow <- workflow() %>% 
  add_recipe(kknn_recipe) %>% 
  add_model(kknn_spec) 
kknn_grid <- tribble(~neighbors, 15, 17, 19, 21)
dataset_folds <- vfold_cv(data = thisDF1, v = 5, repeats = 3, strata = type)
kknn_tune <- tune_grid(kknn_workflow, resamples = dataset_folds, grid = kknn_grid)
knn_final <- finalize_workflow(kknn_workflow, parameters = select_best(kknn_tune)) %>% 
  fit(data = thisDF1)

mesh_size = 10
margin = 50
x_min = min(thisFig[, 4]) - margin
x_max = max(thisFig[, 4]) + margin + 700
y_min = min(thisFig[, 3]) - margin
y_max = max(thisFig[, 3]) + margin + 700
xrange <- seq(x_min, x_max, mesh_size)
yrange <- seq(y_min, y_max, mesh_size)
xy <- meshgrid(x = xrange, y = yrange)
xx <- xy$X
yy <- xy$Y
dim_val <- dim(xx)
xx1 <- matrix(xx, length(xx), 1)
yy1 <- matrix(yy, length(yy), 1)
final <- data.frame(xx1, yy1)
colnames(final) <- c("Datasets_derived_15", "Gradient_derived_15")
pred <- knn_final %>% predict(final, type = 'prob')
grid <- cbind(final, pred)

#plot the FINs, AINs, and validation datasets
ggplot(thisFig %>% filter(type != "JB combinations"), aes(x = Datasets_derived_15, y = Gradient_derived_15)) +
  geom_raster(data = grid, aes(fill = `.pred_Ferroptosis`)) +
  geom_point(aes(color = type), size = 3) +
  geom_text_repel(aes(label = dataset), color = "black", size = 4, segment.colour = "black", show.legend = FALSE,
                  max.overlaps = 30, nudge_x = 30, nudge_y = 30) +
  scale_fill_gradient2(low = "#79A7D0", mid = "white", high = "#70AD47", midpoint = 0.5, breaks = seq(0,1,0.25),
                       guide = "none") +
  scale_color_manual(values = c("blue", "darkgreen", "red")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "ssGSEA score of 15 dataset-derived genes\n(AUC = 0.90, p-value = 0.004)",
       y = "ssGSEA score of 15 gradient-derived genes\n(AUC = 0.88, p-value = 0.005)", 
       color = "", fill = "Ferroptosis\nprobability") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.16, 0.87),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 12),
        legend.margin = margin(0,5,5,5),
        axis.text = element_text(color = "black"))

ggsave("../results/Figure S7D - Contour plot of biomarkers in both sets (showing datasets).tiff", dpi = 300, height = 6, width = 7)

#plot the same contour plot, showing only the JB2/3/6 combinations
thisFig2 <- thisFig %>% 
  filter(type == "JB combinations") %>% 
  mutate(dataset = ifelse(!str_detect(dataset, "to"), str_c(dataset, "_to_DMSO"), dataset))

ggplot(thisFig2, aes(x = Datasets_derived_15, y = Gradient_derived_15)) +
  geom_raster(data = grid, aes(fill = `.pred_Ferroptosis`)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = dataset), color = "black", size = 4, segment.colour = "black", show.legend = FALSE,
                  max.overlaps = 30, nudge_x = 30, nudge_y = 30) +
  scale_fill_gradient2(low = "#79A7D0", mid = "white", high = "#70AD47", midpoint = 0.5, breaks = seq(0,1,0.25),
                       guide = "none") +
  scale_color_manual(values = c("purple")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "ssGSEA score of 15 dataset-derived genes",
       y = "ssGSEA score of 15 gradient-derived genes", 
       color = "", fill = "Ferroptosis\nprobability") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.22, 0.90),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black"),
        legend.text = element_text(size = 12),
        legend.margin = margin(0,5,5,5),
        axis.text = element_text(color = "black"))

ggsave("../results/Figure S7E - Contour plot of biomarkers in both sets (showing JB combinations).tiff", 
       dpi = 300, height = 6, width = 7)

#box-plots of the signatures (Figure 6A)
thisFig <- ssgsea_result %>% 
  left_join(AUCs_all, by = c("pathway", "method")) %>% 
  filter(type != "JB combinations") %>% 
  filter(pathway %in% c("Datasets_derived_15", "Gradient_derived_15", "F1_set", "Hallmark_Apoptosis")) %>% 
  filter(score != Inf) %>% 
  mutate(group_1 = case_when(pathway == "Datasets_derived_15" ~ "Datasets-derived genes (top 15)",
                             pathway == "Gradient_derived_15" ~ "Gradient-derived genes (top 15)",
                             pathway == "F1_set" ~ "Erastin in HT1080 (F1, 33 genes)",
                             pathway == "Hallmark_Apoptosis" ~ "Hallmark_Apoptosis (154 genes)")) %>% 
  na.omit() %>% 
  mutate(auc_label = str_c("AUC = ", format(round(auc_total, digits = 2)), "\np-value ", pval_total)) %>% 
  mutate(group_1 = factor(group_1, levels = c("Datasets-derived genes (top 15)", "Gradient-derived genes (top 15)",
                                              "Erastin in HT1080 (F1, 33 genes)", "Hallmark_Apoptosis (154 genes)")))

stat.test <- thisFig %>%
  group_by(group_1) %>%
  rstatix::t_test(score ~ type, p.adjust.method = "none", var.equal = TRUE,
                  comparisons = list(c("Ferroptosis", "Apoptosis"), c("Validation", "Apoptosis"))) %>%
  add_xy_position() %>%
  add_significance()

stat.test_onesample <- thisFig %>%
  group_by(group_1, type) %>%
  rstatix::t_test(score ~ 1, p.adjust.method = "none", var.equal = TRUE, mu = 0) %>%
  add_y_position() %>% 
  add_significance() %>% 
  mutate(p.lab = ifelse(p > 0.05, "#", NA)) %>% #the "#" symbol will mark non-significance compared to the baseline in the one-sample t-test
  filter(type == "Apoptosis")

ggplot(thisFig, aes(x = type, y = score)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(aes(fill = type)) +
  geom_point(aes(fill = type), size = 3, shape = 21, color = "black") +
  geom_text(data = thisFig %>% distinct(group_1, .keep_all = TRUE), aes(label = auc_label), y = -1600, size = 5) +
  geom_text(data = stat.test_onesample, aes(label = p.lab, y = y.position + 100), size = 5) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE, size = 5) +
  scale_fill_manual(values = c("#79A7D0", "#70AD47", "salmon")) +
  scale_y_continuous(expand = c(0, 500)) +
  facet_wrap(~ group_1, nrow = 1, labeller = label_wrap_gen(width = 20)) +
  labs(x = "", y = "ssGSEA score", fill = "") +
  theme_bw(base_size = 16) +
  theme(legend.position = "top", 
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_blank())

ggsave("../results/Figure 6A - ssGSEA scores of the biomarker sets for datasets.tiff", dpi = 300, height = 5, width = 8)

#plot enrichment of the optimized biomarkers of all datasets + the validation sets
thisFig <- ssgsea_result %>% 
  filter(type != "JB combinations") %>% 
  filter(pathway == "optimized_biomarkers") %>% 
  mutate(type = ifelse(type == "Validation", dataset, type)) %>% 
  group_by(type) %>% 
  get_summary_stats() %>% 
  mutate(name = case_when(str_detect(type, "_E") ~ str_replace(type, "_E", " + Erastin"),
                          str_detect(type, "_R") ~ str_replace(type, "_R", " + RSL3"),
                          TRUE ~ str_c("Public ", type))) %>% 
  mutate(type = ifelse(str_detect(type, "_E|_R"), "Validation", type)) %>% 
  arrange(type, mean) %>% 
  mutate(name = fct_inorder(name))

ggplot(thisFig, aes(x = name, y = mean)) +
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se), width = 0.2) +
  geom_col(aes(fill = type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#79A7D0", "#70AD47", "gray50")) +
  labs(x = "", y = "ssGSEA score", fill = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", 
        legend.direction = "vertical",
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("../results/Figure 6G - Final biomarkers (25) enrichment in all datasets (aggregated).tiff", dpi = 300, height = 6, width = 2.5)

#plot the AUCs of the gradient-derived set by using different geneset size (10-26 genes) (Figure S7C, right)
set.seed(123)
take <- c("^logFC")
cutoffs <- seq(10, 26, 1)
gradient_derived_sets <- list()

for(i in seq_along(cutoffs)) {
  gradient_derived_sets[[i]] <- Gradient_derived %>% slice_min(rank_by_JB2vJB6_24h , n = cutoffs[i]) %>% pull(genes)
  names(gradient_derived_sets)[i] <- str_c("cutoff_", cutoffs[i])
}

gradient_derived_sets <- compact(gradient_derived_sets)
genesets <- c(gradient_derived_sets, random_sets)
ssgsea_result_all <- gsva(Datasets_matrix %>% na.omit(), genesets, method = "ssgsea", ssgsea.norm = FALSE) %>% 
  t() %>% 
  as_tibble(rownames = "dataset") %>% 
  gather(-dataset, key = "pathway", value = "score") %>% 
  mutate(method = paste0("SSGSEA_", str_remove(take, "\\^"))) %>% 
  mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                          str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                          str_detect(dataset, "JB") ~ "JB combinations", 
                          TRUE ~ "Validation")) 

combined_AUCs <- ssgsea_result_all %>% 
  mutate(type = ifelse(type == "Validation", "Ferroptosis", type)) %>% 
  select(-dataset) %>% 
  group_by(pathway, method) %>% 
  dplyr::summarise(auc_total = as.numeric(pROC::auc(response = type, predictor = score, levels = c("Ferroptosis", "Apoptosis")))) %>% 
  group_by(method) %>% 
  mutate(rank_total = min_rank(-auc_total)) %>% 
  mutate(pval_total = rank_total / 1000) %>% 
  ungroup() %>% 
  filter(!str_detect(pathway, "Random"))

thisFig <- combined_AUCs %>% 
  mutate(cutoff = as.numeric(str_extract(pathway, "[:digit:][:digit:]$"))) %>% 
  mutate(method = str_remove(pathway, "_[:digit:][:digit:]$"))

ggplot(thisFig, aes(x = cutoff, y = auc_total, color = method)) +
  geom_point(size = 2, color = "black") +
  geom_line(size = 1, color = "black") +
  scale_x_continuous(breaks = seq(10, 26, 4)) +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  labs(x = "Number of genes in signature", y = "ROC-AUC of ferroptosis\nto apoptosis classification", color = "",
       title = "Gradient-derived signature") +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        legend.text = element_text(size = 16),
        plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"))

ggsave("../results/Figure S7C (right) - AUCs of gradient derived biomarkers.tiff", dpi = 300, height = 5, width = 5)

#venn diagram depicting the inverse gradient idea (Figure S7A)
InvGrad_Venn <- list("JB2 up vs. JB6 in 24h" = DEG_24h %>% 
                       filter(JB6_to_JB2_24h_fdr == -1) %>% 
                       pull(genes),
                     "JB2 down vs. JB6 in 6h" = DEG_24h %>% 
                       left_join(DEG_6h, by = "genes") %>%
                       filter(JB6_to_JB2_6h_fdr == 1) %>% 
                       pull(genes))

InvGrad_partitions <- get.venn.partitions(InvGrad_Venn)
tiff(filename = "../results/Figure S7A - Euler plot for inverse gradient genes.tiff", 
     width = 600, height = 600)
plot(euler(InvGrad_Venn), quantities = list(cex = 2), lty = 1, lwd = 2)
dev.off()

#venn diagram of all the biomarkers (Fig 6B)
this_venn <- list("Gradient_derived_26" = Gradient_derived %>% pull(genes),
                  "Datasets_derived_15" = Datasets_derived_sig %>% filter(stat == "logFC") %>% 
                    slice_max(mean_difference, n = 15, with_ties = FALSE) %>% pull(ID), 
                  "F1_set" = Ferroptosis1 %>% slice_max(logFC_F1, n = 36) %>% filter(!str_detect(ID, "LOC100")) %>% pull(ID))

this_venn_partitions <- get.venn.partitions(this_venn)
venn_partitions <- this_venn_partitions %>% 
  rowwise() %>% 
  mutate(..values.. = paste0(..values.., collapse = " ")) #actual Venn shown in Fig. 6B was drawn in powerpoint based on this object
