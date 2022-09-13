#This script will create Figure 1H, showing the JB2-JB3-JB6 transcriptomics on a ferroptosis-apoptosis plane.

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

#make similar dataset out of the combinations (our RNAseq). For that, I'm extracting the t-statistics of the treatments vs. DMSO
#from the efit object created in the main script.
#for MDA468 cells, 24 hours:
combinations_t_24h <- efit_24h[["t"]] %>% 
  as_tibble(rownames = "ID") %>% 
  dplyr::select(-contains("to")) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("t_", .x)) %>% 
  dplyr::select(ID, t_JB2_24h:t_JB6_24h)
len <- dim(combinations_t_24h)[1]
combinations_t2_24h <- combinations_t %>% 
  mutate(across(-ID, ~min_rank(.x))) %>% 
  mutate(across(-ID, ~.x/len)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "t_", "Rank_t_"))

#get similar info for HCC70 (just t-statistics)
combinations_t_H70 <- efit_H70[["t"]] %>% 
  as_tibble(rownames = "ID") %>% 
  dplyr::select(-contains("to")) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("t_", .x)) %>% 
  dplyr::select(ID, t_JB2_H70:t_JB6_H70)
len <- dim(combinations_t_H70)[1]
combinations_t2_H70 <- combinations_t_H70 %>% 
  mutate(across(-ID, ~min_rank(.x))) %>% 
  mutate(across(-ID, ~.x/len)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "t_", "Rank_t_"))

#get similar info for MDA468 6hrs (just t-statistics)
combinations_t_6h <- efit_6h[["t"]] %>% 
  as_tibble(rownames = "ID") %>% 
  dplyr::select(-contains("to")) %>% 
  rename_with(.cols = -ID, .fn = ~str_c("t_", .x)) %>% 
  dplyr::select(ID, t_JB2_6h:t_JB6_6h)
len <- dim(combinations_t_6h)[1]
combinations_t2_6h <- combinations_t_6h %>% 
  mutate(across(-ID, ~min_rank(.x))) %>% 
  mutate(across(-ID, ~.x/len)) %>% 
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "t_", "Rank_t_"))

#------------------------------------------------------------------------------------------------
# Create apoptosis-ferroptosis plane 
# This plot will show the combination RNAseq results on a plane based on a ferroptosis score
# in the x-axis, and apoptosis score on the y-axis. 
#------------------------------------------------------------------------------------------------

#datasets used to generate the x and y axis - ferroptosis and apoptosis inducers in TNBC
Ferroptotis_datasets <- c("F16$", "F18$")
Apoptosis_datasets <- c("A4$")

#combine for match() function
Selected_ferroptosis_datasets <- paste0(Ferroptotis_datasets, collapse = "|")
Selected_apoptosis_datasets <- paste0(Apoptosis_datasets, collapse = "|")
NA_cutoff <- 1

#create the ferroptosis and apoptosis consensus signatures based on these sets
Ferroptosis_DF <- Ferroptosis_combined %>% 
  select(ID, contains("Rank_t")) %>% 
  select(ID, matches(Selected_ferroptosis_datasets)) %>% 
  mutate(count_na = rowSums(is.na(.))) %>% 
  filter(count_na < NA_cutoff) %>% 
  select(-count_na) %>% 
  gather(-ID, key = "dataset", value = "value") %>% 
  group_by(ID) %>% 
  mutate(mean_F = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(-dataset, -value) %>% 
  distinct(ID, .keep_all = TRUE)

Apoptosis_DF <- Apoptosis_combined %>% 
  select(ID, contains("Rank_t")) %>% 
  select(ID, matches(Selected_apoptosis_datasets)) %>% 
  mutate(count_na = rowSums(is.na(.))) %>% 
  filter(count_na < NA_cutoff) %>% 
  select(-count_na) %>%   gather(-ID, key = "dataset", value = "value") %>% 
  group_by(ID) %>% 
  mutate(mean_A = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(-dataset, -value) %>% 
  distinct(ID, .keep_all = TRUE)

thisDF_FAP <- left_join(Ferroptosis_DF, Apoptosis_DF, by = "ID") %>% 
  na.omit() %>% 
  mutate(mean_dif = mean_F - mean_A)

sig_len <- 100

Ferro_consensus <- thisDF_FAP %>% slice_max(mean_dif, n = sig_len) %>% pull(ID)
Apop_consensus <- thisDF_FAP %>% slice_min(mean_dif, n = sig_len) %>% pull(ID)
this_geneset <- list("Ferroptosis" = Ferro_consensus, "Apoptosis" = Apop_consensus)

#Quantify the enrichemnt of these signatures in the JB2-JB3-JB6 combinations, and in the datasets used to create the signatures, by GSVA
this_GSVA <- combinations_t2_24h %>% 
  left_join(combinations_t2_H70, by = "ID") %>% 
  left_join(combinations_t2_6h, by = "ID") %>% 
  left_join(Ferroptosis_combined, by = "ID") %>% 
  left_join(Apoptosis_combined, by = "ID") %>% 
  select(ID, contains(c("24h", "Rank_t"))) %>% 
  select(ID, contains("JB"), matches(Selected_ferroptosis_datasets), matches(Selected_apoptosis_datasets)) %>% 
  na.omit() %>% 
  distinct(ID, .keep_all = TRUE) %>% 
  column_to_rownames(var = "ID") %>% 
  as.matrix()
  
gsva_result <- gsva(this_GSVA, this_geneset, method = "gsva")

FAP_GSVA <- gsva_result %>% t() %>% 
  as_tibble(rownames = "dataset") %>% 
  mutate(method = "GSVA") %>% 
  mutate(dataset = str_remove_all(dataset, "Rank_t_"))
  
#helper function that will create the backgroud color of the plot
make_gradient <- function(deg = 45, n = 100, cols = blues9) {
  cols <- colorRampPalette(cols)(n + 1)
  rad <- deg / (180 / pi)
  mat <- matrix(data = rep(seq(0, 1, length.out = n) * cos(rad), n),
                byrow = TRUE,
                ncol = n) +
    matrix(data = rep(seq(0, 1, length.out = n) * sin(rad), n),
           byrow = FALSE,
           ncol = n)
  mat <- mat - min(mat)
  mat <- mat / max(mat)
  mat <- 1 + mat * n
  mat <- matrix(data = cols[round(mat)], ncol = n)
  grid::rasterGrob(image = mat,
                   width = unit(1, "npc"),
                   height = unit(1, "npc"), 
                   interpolate = TRUE
  )
}

#create the plot
Inducers <- tribble(~"dataset", ~"inducer", "A4", "5FU", "F16", "ML162", "F18", "Erastin")
  
thisFig <- FAP_GSVA %>% 
  mutate(group = case_when(str_detect(dataset, "F") ~ "Ferroptosis",
                           str_detect(dataset, "A") ~ "Apoptosis",
                           str_detect(dataset, "H70") ~ "HCC70 24h",
                           str_detect(dataset, "24h") ~ "MDA468 24h",
                           str_detect(dataset, "6h") ~ "MDA468 6h",
                           TRUE ~ "other")) %>% 
  left_join(Inducers, by = "dataset") %>% 
  mutate(inducer = ifelse(is.na(inducer), dataset, inducer)) %>% 
  mutate(inducer = str_remove(inducer, "_24h|_6h|_H70")) %>% 
  filter(group != "MDA468 6h")
  
cols <- colorRampPalette(c("#79A7D0", "#A2C4E3", "#BCDAA7", "#70AD47"))(100)
g <- make_gradient(deg = 45, n = 500, cols = cols)
  
ggplot(thisFig, aes(x = `Ferroptosis`, y = `Apoptosis`, color = group)) +
  annotation_custom(grob = g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
  geom_segment(aes(x = -Inf, y = -Inf, xend = Inf, yend = Inf), 
                 color = "white", linetype = "dashed", size = 1.5) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = inducer), size = 6, show.legend = FALSE) +
  scale_color_manual(values = c("blue", "darkgreen", "black", "gray40", "yellowgreen")) +
  labs(x = "Ferroptosis score", y = "Apoptosis score") +
  theme_bw(base_size = 14) + 
  theme(legend.position = c(0.2, 0.2),
        legend.background = element_rect(color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, color = "black"))
    
ggsave("../results/Figure 1H.tiff", dpi = 200, width = 6, height = 6)
