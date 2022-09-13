#This script will validate the Gradient Gene Set, as shown in figure 6C, D

#NOTES: 
#1. Run-time of this script is ~6 hours, due to the need to analyze 1000 iterations of random and datasets-derived genesets.
#   For fast reproduction of the code, we offer two saved objects replacing two "for" loops (see later in the code)
#2. The resulting UMAP might be different than the one appearing in the manuscripot (which was prepared in R version 3.5),
#   however the results of the classification [Fig. 6D] should be almost similar (since the random seed was changed in R version 4.2).

#load the ferrotposis and apoptosis inducers datasets
#see supplemental table S3 for list of datasets used. The datasets are annotated F1-F18 (for ferroptosis inducers)
#and A1-A18 (for apoptosis inducers). These files contains 4 column for each dataset:
#logFC_X - log2 fold change for each gene of the specific inducer vs. its control
#t_X - the t-statistics for the above FC, calculated by Limma.
#Rank_FC - the fold change data, ranked on a scale from 0 to 1
#Rank_t - the t-statistic data, ranked on a scale from 0 to 1 (this is the metric used mostely)
#note: each dataset was analyzed separately, these tables (already analysed and ranked) are provided for reproducability.
#For more details on these datasets please see the supplemental text.

Ferroptosis_combined <- read_csv("../data/Ferroptosis inducers datasets.csv")
Apoptosis_combined <- read_csv("../data/Apoptosis inducers datasets.csv")

#Function definition: contour_plot
#This function accepts a list of genes (such as the GGS genes, random genes, etc.),
#and output the UMAP plot of the 18 FINs and 18 AINs datasets.
#This function also performs k-nearest neighbors (kNN) in order to compute classification metrics.
#Arguments:
#geneset = the genes used to create the UMAP
#plot = should plots be produced?
#title = a title identifing the geneset
#See supplemental text for details.

set.seed(123) #a rnadom seed is set for reproducability (for the kNN method; the umap might still be different)

Contour_plot <- function(geneset, plot = TRUE, title = "") {
  geneset <- unique(geneset)
  print(paste0("Working on ", length(geneset), " genes"))
  normed_method <- "Rank_t_"
  
  #combine the ferroptosis and apoptosis datasets
  Datasets_combined <- Ferroptosis_combined %>% 
    full_join(Apoptosis_combined, by = "ID") %>% 
    select(ID, matches(normed_method)) %>% 
    rename_with(.cols = -ID, .fn = ~str_replace(.x, normed_method, ""))
  
  #filter to include only the genes supplied to the function, and deal with missing values
  this_dataset <- Datasets_combined %>% 
    filter(ID %in% geneset) %>% 
    gather(-ID, key = "dataset", value = "value") %>% 
    distinct(ID, dataset, .keep_all = TRUE) %>% 
    spread(key = "ID", value = "value") %>% 
    mutate(type = case_when(str_detect(dataset, "^F[:digit:]*") ~ "Ferroptosis",
                            str_detect(dataset, "^A[:digit:]*") ~ "Apoptosis"), .after = dataset) %>% 
    select_if(~!sum(is.na(.)) > 4) %>%  #remove coloumn with more than 4 missing values, and then..
    group_by(type) %>%                  # ..impute missing values by median for each type
    mutate(across(-dataset, ~ifelse(is.na(.x), median(.x, na.rm = TRUE), .x))) %>% 
    column_to_rownames(var = "dataset") %>% 
    select(-type)
  
  thisDF.map <- umap(this_dataset, preserve.seed = TRUE, random_state = 84517, min_dist = 0.8)
  result <- thisDF.map$layout %>% 
    as_tibble(rownames = "treat") %>% 
    mutate(group = case_when(str_detect(treat, "A[:digit:]") ~ "Apoptosis",
                             str_detect(treat, "F[:digit:]") ~ "Ferroptosis",
                             TRUE ~ "NA"))
  
  #perform kNN to create the contour plot, and collect the ROC-AUCs of classification
  thisDF1 <- result
  
  kknn_recipe <- recipe(formula = group ~ V1 + V2, data = thisDF1)
  kknn_spec <- nearest_neighbor(neighbors = tune()) %>% 
    set_mode("classification") %>% 
    set_engine("kknn") 
  kknn_workflow <- workflow() %>% 
    add_recipe(kknn_recipe) %>% 
    add_model(kknn_spec) 
  kknn_grid <- tribble(~neighbors, 9, 11, 13, 15)
  dataset_folds <- vfold_cv(v = 5, thisDF1, repeats = 3, strata = group)
  kknn_tune <- tune_grid(kknn_workflow, resamples = dataset_folds, grid = kknn_grid)
  knn_final <- finalize_workflow(kknn_workflow, parameters = select_best(kknn_tune)) %>% 
    fit(data = thisDF1)
  knn_table <- kknn_tune %>% show_best(n = 10) %>% 
    arrange(neighbors) %>% select("parameter" = neighbors, "value" = mean) %>% 
    mutate(parameter = str_c("k = ", parameter)) %>% 
    mutate(value = round(value, 3))
  
  #for the contour plot, create a grid of points across the plane and predict ferroptosis probability
  #in each point using the kNN model trained above
  mesh_size = .02
  margin = 0.1
  x_min =  min(result[, 2]) - margin
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
  
  #plot the contour plot
  if(plot) {
    ggplot(grid, aes(x = V1, y = V2)) +
      geom_raster(aes(fill = `.pred_Ferroptosis`)) +
      geom_point(data = result, aes(x = V1, y = V2, color = group), size = 3) +
      geom_text_repel(data = result, aes(label = treat), size = 4) +
      scale_fill_gradient2(low = "#79A7D0", mid = "white", high = "#70AD47", midpoint = 0.5) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_color_manual(values = c("blue", "darkgreen", "red")) +
      labs(x = "", y = "", fill = "Ferroptosis\nprobability", color = "Dataset") +
      theme_minimal(base_size = 14) +
      theme(panel.grid = element_blank(),
            legend.justification = "bottom") 
    
    ggsave(paste0("../results/Figure 6C - representing ", title, ".tiff"), dpi = 300, height = 5, width = 7)
  }
  
  #calculate the probability of ferroptosis for each of the 36 inducers using kNN imputation
  grid_with_results <- grid %>% 
    mutate(group = "grid") %>% 
    bind_rows(result) %>% 
    select(V1, V2, .pred_Ferroptosis, group, treat)
  transcript_comp <- impute_knn(grid_with_results, .pred_Ferroptosis ~ V1 + V2)
  thisFig_2 <- transcript_comp %>% 
    filter(group != "grid")
  
  if(plot) {
    ggplot(thisFig_2, aes(x = reorder(treat, -.pred_Ferroptosis) , y = .pred_Ferroptosis, fill = group)) +
      geom_col() +
      geom_text(aes(label = treat), y = 0.05, angle = 90) +
      geom_hline(yintercept = 0.5) +
      scale_fill_manual(values = c("#79A7D0", "#70AD47", "salmon")) +
      labs(x = "", y = "Ferroptosis probability") +
      theme_bw(base_size = 14) +
      theme(legend.position = "non",
            panel.grid = element_blank(),
            axis.text.x = element_blank())
    
    ggsave(paste0("../results/Figure S6B - representing ", title, ".tiff"), dpi = 300, height = 4, width = 6)
  }
  
  #return accuracy calculations (treat ferroptosis as positive and apoptosis as negative)
  thisFig_3 <- thisFig_2 %>% 
    mutate(predict = ifelse(.pred_Ferroptosis > 0.5, "Ferroptosis", "Apoptosis")) %>% 
    mutate(parameter = case_when(group == "Ferroptosis" & predict == "Ferroptosis" ~ "TP",
                                 group == "Ferroptosis" & predict == "Apoptosis" ~ "FN",
                                 group == "Apoptosis" & predict == "Ferroptosis" ~ "FP",
                                 group == "Apoptosis" & predict == "Apoptosis" ~ "TN")) %>% 
    group_by(parameter) %>% 
    summarise(value = n())
  
  thisFig_4 <- thisFig_2 %>% 
    select("parameter" = treat, "value" = .pred_Ferroptosis)
  
  result <- bind_rows(thisFig_3, knn_table, thisFig_4)
  return(result)
}

#plot contour for the GGS genes (the collected 4 genesets prepared in the `5 - Gradient Gene Set creation.R`)
#FAP = "Ferroptosis Apoptosis Plane"
GGS <- rbind(group_JB2U_C_100, group_JB2D_C_100, group_JB6U_C_100, group_JB6D_C_100) %>% 
  distinct(genes) %>% pull(genes)

FAP_for_GGS <- Contour_plot(geneset = GGS, plot = TRUE, title = "GGS") %>% 
  dplyr::rename("GGS" = "value")

#predictions by random genesets, In each iteration pick a random set of genes to plot the contour plot on.
#this is done in order to compute the p-values of the ROC=AUCs
set.seed(123)
iterations <- 1000
set_size <- 306
random_results <- list()

#NOTE - the following "for" loop will run for ~3 hours. To reproduce the results, we offer two options:
#option 1: uncomment this "for" loop to run it:

  # for(n in 1:iterations) {
  #   print(paste0("Iteration: ", n, " / ", iterations))
  #   Random_set <- Datasets_combined %>% slice_sample(n = set_size) %>% pull(ID)
  #   random_results[[n]] <- Contour_plot(geneset = Random_set, plot = TRUE, title = "Random_set") %>% 
  #     rename_with(.cols = -parameter, .fn = ~str_c("random_", n))
  # }

#option 2: load the object created by this "for" loop ("random_results"):
load("../data/Response signature results for random genesets (1000 permutations, seed 84517, R 4.2).Rdata")

random_results_tib <- purrr::reduce(random_results, left_join, by = "parameter") 

#predictions based on the 36 datasets: the universal set.
#In each iterations, 10 ferroptosis and 10 apoptosis datasets will be picked at random, and 
#a signature differentiating these 20 will be made. This signature will be sent to the contour plot function to determine how well
#the other 8 datasets of each type is predicted.
set.seed(123)
iterations <- 1000
set_size <- 306
FAP_datasets_results <- list()
Ferro_consensus_keep <- list()
Apo_consensus_keep <- list()

#NOTE - the following for loop will run for ~3 hours. To reproduce the results, we offer two options:
#option 1: uncomment the "for" loop to run it:
  # for(n in 1:iterations) {
  #   print(paste0("Iteration: ", n, " / ", iterations))
  # 
  #   #randomly select 10 FINs and 10 AINs datasets (with repetitions)
  #   ferro_dataset_select <- round(runif(n = 10, min = 2, max = 19), 0)
  #   appo_dataset_select <- round(runif(n = 10, min = 2, max = 19), 0)
  #   
  #   Ferroptosis_DF <- Ferroptosis_combined %>% 
  #     select(ID, contains("Rank_t")) %>% 
  #     select(ID, all_of(ferro_dataset_select)) %>% 
  #     mutate(count_na = rowSums(is.na(.))) %>% 
  #     filter(count_na < 4) %>% 
  #     select(-count_na) %>% 
  #     gather(-ID, key = "dataset", value = "value") %>% 
  #     group_by(ID) %>% 
  #     mutate(mean_F = mean(value, na.rm = TRUE)) %>% 
  #     ungroup() %>% 
  #     select(-dataset, -value) %>% 
  #     distinct(ID, .keep_all = TRUE)
  #   
  #   Apoptosis_DF <- Apoptosis_combined %>% 
  #     select(ID, contains("Rank_t")) %>% 
  #     select(ID, all_of(appo_dataset_select)) %>% 
  #     mutate(count_na = rowSums(is.na(.))) %>% 
  #     filter(count_na < 4) %>% 
  #     select(-count_na) %>%   gather(-ID, key = "dataset", value = "value") %>% 
  #     group_by(ID) %>% 
  #     mutate(mean_A = mean(value, na.rm = TRUE)) %>% 
  #     ungroup() %>% 
  #     select(-dataset, -value) %>% 
  #     distinct(ID, .keep_all = TRUE)
  #   
  #   thisDF <- left_join(Ferroptosis_DF, Apoptosis_DF, by = "ID") %>% 
  #     na.omit() %>% 
  #     mutate(mean_dif = mean_F - mean_A)
  #   
  #   Ferro_consensus_keep[[n]] <- thisDF %>% slice_max(mean_dif, n = set_size/2) %>% pull(ID)
  #   Apo_consensus_keep[[n]] <- thisDF %>% slice_min(mean_dif, n = set_size/2) %>% pull(ID)
  #   
  #   FAP_geneset <- c(Ferro_consensus_keep[[n]], Apo_consensus_keep[[n]])
  # 
  #   FAP_datasets_results[[n]] <- Contour_plot(geneset = FAP_geneset, plot = TRUE, title = "universal_set") %>% 
  #     rename_with(.cols = -parameter, .fn = ~str_c("FAP_", n))
  # }

#option 2: load the object created by this "for" loop ("FAP_datasets_results"):
load("../data/Response signature results for FAP datasets genesets (1000 permutations, seed 84517, R 4.2).Rdata")

FAP_results_tib <- purrr::reduce(FAP_datasets_results, left_join, by = "parameter") 

#plot the classification metrics (Fig. S63)
thisFig <- FAP_for_GGS %>% 
  full_join(random_results_tib, by = "parameter") %>% 
  full_join(FAP_results_tib, by = "parameter") %>% 
  gather(-parameter, key = "geneset", value = "value") %>% 
  spread(key = "parameter", value = "value") %>% 
  mutate(FN = ifelse(is.na(FN), 0, FN)) %>%   
  mutate(FP = ifelse(is.na(FP), 0, FP)) %>% 
  mutate(`Ferroptosis accuracy` = TP/(TP+FN)) %>% 
  mutate(`Apoptosis accuracy` = TN/(TN+FP)) %>% 
  rowwise() %>% 
  mutate(`ROC-AUC of prediction` = max(`k = 9`, `k = 11`, `k = 13`, `k = 15`)) %>% 
  ungroup() %>% 
  mutate(group = case_when(str_detect(geneset, "random") ~ "Random",
                           str_detect(geneset, "FAP") ~ "Datasets derived",
                           str_detect(geneset, "GGS") ~ "GGS"))

thisFig2 <- thisFig %>% 
  select(group, geneset, `Ferroptosis accuracy`, `Apoptosis accuracy`, `ROC-AUC of prediction`) %>% 
  gather(`Ferroptosis accuracy`:`ROC-AUC of prediction`, key = "metric", value = "metric_value") %>% 
  mutate(metric = factor(metric, levels = c("ROC-AUC of prediction", "Ferroptosis accuracy", "Apoptosis accuracy")))

mean_types <- thisFig2 %>% 
  group_by(group, metric) %>% 
  summarise(mean = mean(metric_value)) %>% 
  ungroup()

ggplot(data = thisFig2, aes(x = metric_value)) +
  geom_vline(data = mean_types %>% filter(group == "Random"), aes(xintercept = mean), 
             color = "salmon", linetype = "dashed", size = 1) +
  geom_vline(data = mean_types %>% filter(group == "Datasets derived"), aes(xintercept = mean), 
             color = "blue", linetype = "dashed", size = 1) +
  geom_vline(data = mean_types %>% filter(group == "GGS"), aes(xintercept = mean), 
             color = "black", linetype = "dashed", size = 1) +
  geom_density(data = thisFig2 %>% filter(group == "Random"), alpha = .3, fill="salmon") +
  geom_density(data = thisFig2 %>% filter(group == "Datasets derived"), alpha = .3, fill="lightblue") +
  facet_wrap(~metric, ncol = 1, scales = "free_y") +
  labs(x = "Metric value", y = "Density") +
  theme_bw(base_size = 14)

ggsave("../results/Figure S6C.tiff", dpi = 300, height = 7, width = 6)

#Calculate p-values based on random distribution
Random_distribution <- thisFig2 %>% 
  filter(group == "Random") %>% 
  filter(metric == "ROC-AUC of prediction") %>% 
  select(group, metric, "mean" = metric_value)

random_length <- dim(Random_distribution)[1] + 4

pvalues <- mean_types %>% 
  filter(metric == "ROC-AUC of prediction") %>% 
  filter(group != "Random") %>% 
  bind_rows(Random_distribution) %>% 
  mutate(rank = min_rank(-mean)) %>% 
  mutate(pval = rank / random_length) %>% 
  filter(group != "Random") %>% 
  mutate(`Gene set` = case_when(str_detect(group, "GGS") ~ "Gradient genes", 
                                str_detect(group, "Datasets derived") ~ "Universal set")) %>% 
  select(`Gene set`, "p-value of AUC" = pval)

#prepare table of metrics
this_table <- thisFig2 %>% 
  mutate(`Gene set` = case_when(str_detect(geneset, "GGS") ~ "Gradient genes", 
                                str_detect(geneset, "random") ~ "Random genes",
                                str_detect(geneset, "FAP") ~ "Universal set")) %>% 
  group_by(`Gene set`, metric) %>% 
  summarise(mean = round(mean(metric_value), 3)) %>% 
  spread(key = "metric", value = "mean") %>% 
  left_join(pvalues, by = "Gene set") %>%
  mutate(`p-value of AUC` = round(`p-value of AUC`, 3)) %>% 
  arrange(desc(`ROC-AUC of prediction`)) %>% 
  rename_with(.cols = everything(), .fn = ~str_wrap(.x, width = 10)) %>% 
  filter(`Gene set` %in% c("Gradient genes", "Universal set", "Random genes"))

ggtexttable(this_table, rows = NULL, theme = ttheme(
  colnames.style = colnames_style(color = "black", size = 15),
  tbody.style = tbody_style(color = "black", size = 15)))

ggsave("../results/Figure 6D.tiff", dpi = 300, height = 4, width = 8)