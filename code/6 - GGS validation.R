#Scripts written by Yaron Vinik
#This script will validate the Gradient Gene Set, as shown in figure 3C

#load the ferrotposis and apoptosis inducers datasets
#see supplemental table S1 for list of datasets used. The datasets are annotated F1-F19 (for ferroptosis inducers)
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

#Function definition: contour_plot
#This function accepts a list of genes (such as the GGS genes, random genes, etc.),
#and output the UMAP plot of the 19 FINs and 26 AINs datasets.
#This function also performs k-nearest neighbors (kNN) in order to compute classification metrics.
#Arguments:
#geneset = the genes used to create the UMAP
#plot = should plots be produced?
#title = a title identifing the geneset
#See supplemental text for details.

set.seed(123)
Contour_plot <- function(geneset, plot = TRUE, title = "test") {
  geneset <- unique(geneset)
  imputation_cutoff = 12 #genes with more than this ammont of missing data will be droped; those with lower amount of missing data
                         #will be imputaed (see supplemental text)
  Datasets_combined_2 <- Datasets_combined %>% 
  select(ID, matches("^Rank_t_")) %>% #this analysis is done on the ranked t-statistics
  rename_with(.cols = -ID, .fn = ~str_replace(.x, "^Rank_t_", ""))
  
  this_dataset <- Datasets_combined_2 %>% 
    filter(ID %in% geneset) %>% 
    gather(-ID, key = "dataset", value = "value") %>% 
    distinct(ID, dataset, .keep_all = TRUE) %>% 
    spread(key = "ID", value = "value") %>% 
    mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                            str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis", 
                            TRUE ~ "Ferroptosis"), .after = dataset) %>% 
    select_if(~!sum(is.na(.)) > imputation_cutoff) %>%  
    group_by(type) %>%                  
    mutate(across(-dataset, ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) %>% #missing data imputation using median
    ungroup() %>% 
    column_to_rownames(var = "dataset") %>% 
    select(-type)
  
  thisDF.map <- umap(this_dataset, preserve.seed = TRUE, random_state = 173)
  result <- thisDF.map$layout %>% 
    as_tibble(rownames = "dataset") %>% 
    mutate(group = case_when(str_detect(dataset, "^A[:digit:]") ~ "Apoptosis",
                             str_detect(dataset, "^F[:digit:]") ~ "Ferroptosis",
                             TRUE ~ "Ferroptosis"))
  
  #perform kNN on the UMAP x- and y- coordinates: this model will be used to generate the contours in the countour plot
  if(plot) {
    thisDF1 <- result
    knn_recipe <- recipe(formula = group ~ V1 + V2, data = thisDF1)
    knn_spec <- nearest_neighbor(neighbors = tune()) %>% 
      set_mode("classification") %>% 
      set_engine("kknn") 
    knn_workflow <- workflow() %>% 
      add_recipe(knn_recipe) %>% 
      add_model(knn_spec) 
    knn_grid <- tribble(~neighbors, 5, 7, 9, 11, 13, 15)
    dataset_folds <- vfold_cv(data = thisDF1, v = 5, repeats = 1, strata = group)
    knn_tune <- tune_grid(knn_workflow, resamples = dataset_folds, grid = knn_grid)
    knn_final <- finalize_workflow(knn_workflow, parameters = select_best(knn_tune, metric = "roc_auc")) %>% 
      fit(data = thisDF1)
    knn_table <- knn_tune %>% 
      collect_metrics() %>% 
      mutate(parameter = str_c("k = ", neighbors, " ", .metric)) %>%
      select(parameter, "value" = mean)
  }

  #perform kNN on all the selected genes: this model will be used to determine the AUC of classification by the GGS
  #(so this AUC is not dependent on the UMAP visualization, which is sensitive to random state)
  this_dataset_for_kNN <- this_dataset %>% 
    as_tibble(rownames = "dataset") %>% 
    mutate(group = case_when(str_detect(dataset, "^A[:digit:]") ~ "Apoptosis",
                             str_detect(dataset, "^F[:digit:]") ~ "Ferroptosis",
                             TRUE ~ "Ferroptosis"), .after = dataset) %>% 
    select(-dataset)
  knnAll_recipe <- recipe(formula = group ~ ., data = this_dataset_for_kNN) %>% 
    step_zv(all_predictors()) %>% 
    step_normalize(all_predictors())
  knn_spec <- nearest_neighbor(neighbors = tune()) %>% 
    set_mode("classification") %>% 
    set_engine("kknn") 
  knnAll_workflow <- workflow() %>% 
    add_recipe(knnAll_recipe) %>% 
    add_model(knn_spec)
  knn_grid <- tribble(~neighbors, 5, 7, 9, 11, 13, 15)
  dataset_folds <- vfold_cv(data = this_dataset_for_kNN, v = 5, repeats = 1, strata = group)
  knnAll_tune <- tune_grid(knnAll_workflow, resamples = dataset_folds, grid = knn_grid)
  knnAll_table <- knnAll_tune %>% 
    collect_metrics() %>% 
    mutate(parameter = str_c("k_all = ", neighbors, " ", .metric)) %>%
    select(parameter, "value" = mean)
  
  #create the contour plot using the kNN model built using the UMAP coordiantes
  if(plot) {
    margin = 0.1
    mesh_size = 0.02
    x_min = min(result[, 2]) - margin
    x_max = max(result[, 2]) + margin
    y_min = min(result[, 3]) - margin
    y_max = max(result[, 3]) + margin
    xrange <- seq(x_min, x_max, mesh_size)
    yrange <- seq(y_min, y_max, mesh_size)
    xy <- meshgrid(x = xrange, y = yrange)
    xx <- xy$X
    yy <- xy$Y
    dim_val <- dim(xx)
    xx1 <- matrix(xx, length(xx), 1)
    yy1 <- matrix(yy, length(yy), 1)
    final <- data.frame(xx1, yy1)
    colnames(final) <- c("V1", "V2")
    pred <- knn_final %>% predict(final, type = 'prob')
    grid <- cbind(final, pred)
  
    p1 <- ggplot(grid, aes(x = V1, y = V2)) +
      geom_raster(aes(fill = `.pred_Ferroptosis`)) +
      geom_point(data = result %>% mutate(group_1 = ifelse(str_detect(dataset, "_E|_R"), "Validation", group)),
                 aes(x = V1, y = V2, color = group_1), size = 3) +
      geom_text_repel(data = result, aes(label = dataset), size = 4) +
      scale_fill_gradient2(low = "#79A7D0", mid = "white", high = "#70AD47", midpoint = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(values = c("blue", "darkgreen", "red")) +
      labs(x = "", y = "", fill = "Ferroptosis\nprobability", color = "Dataset") +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(),
            legend.justification = "bottom") 
    
    p1
    ggsave(paste0("../results/Figure 2C - Contour plot (UMAP method) for ", title, ".tiff"), dpi = 300, height = 5, width = 7)
  
    #calculate the probability of ferroptosis for the datasets
    datasets_ferroptosis_prob <- grid %>% 
      mutate(group = "grid") %>% 
      bind_rows(result) %>% 
      select(V1, V2, .pred_Ferroptosis, group, dataset) %>% 
      impute_knn(.pred_Ferroptosis ~ V1 + V2) %>% 
      filter(group != "grid") %>% 
      mutate(dataset = ifelse(group == "validation", str_c(dataset, "_V"), dataset))
  
    #return accuracy calculations (treat ferroptosis as positive and apoptosis as negative)
    Accuracies <- datasets_ferroptosis_prob %>% 
      mutate(predict = ifelse(.pred_Ferroptosis > 0.5, "Ferroptosis", "Apoptosis")) %>% 
      mutate(parameter = case_when(group == "Ferroptosis" & predict == "Ferroptosis" ~ "TP",
                                  group == "Ferroptosis" & predict == "Apoptosis" ~ "FN",
                                  group == "Apoptosis" & predict == "Ferroptosis" ~ "FP",
                                  group == "Apoptosis" & predict == "Apoptosis" ~ "TN")) %>% 
      group_by(parameter) %>% 
      summarise(value = n()) %>% 
      na.omit() 
  
    tempDF2 <- datasets_ferroptosis_prob %>% select("parameter" = dataset, "value" = .pred_Ferroptosis)
    final_result <- bind_rows(Accuracies, tempDF2, knnAll_table)
    return(final_result)
  }
  
  final_result <- bind_rows(knnAll_table)
  return(final_result)
}

#plot contour for the GGS genes (the collected 4 genesets prepared in the `1 - RNAseq analysis - master script.R`)
set.seed(123)
GGS <- rbind(group_JB2U_C_100, group_JB2D_C_100, group_JB6U_C_100, group_JB6D_C_100) %>% 
  distinct(genes) %>% pull(genes)

CP_for_GGS <- Contour_plot(geneset = GGS, plot = TRUE, title = "GGS") %>% 
  rename_with(.cols = value, .fn = ~"SM100")

