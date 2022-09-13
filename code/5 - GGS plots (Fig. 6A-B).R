#This script will create the plots of the Gradient Gene Set (Fig. 6A, B)

#define gradient colors for the plots
JB2Ucol <- brewer.pal(n=8, name="Set1")[1]
JB2Dcol <- brewer.pal(n=8, name="Set1")[2]
JB6Ucol <- brewer.pal(n=8, name="Set1")[3]
JB6Dcol <- brewer.pal(n=8, name="Set1")[4]
colors <- c(JB2Ucol, JB2Dcol, JB6Ucol, JB6Dcol)

#plot representative genes, one for each gradient group (Fig. 6A)
choose <- c("PCK2", "REST", "KDM6A", "BZW2")

thisFig <- DEG_24h %>% 
  gather(-genes, key = "treat", value = "value") %>% 
  spread(key = "genes", value = "value") %>% 
  select(treat, any_of(choose)) %>% 
  filter(str_detect(treat, "24h$")) %>% 
  filter(!str_detect(treat, "to")) %>% 
  gather(-treat, key = "gene", value = "value") %>% 
  mutate(treat = factor(treat, levels = c("JB2_24h", "JB3_24h", "JB6_24h"))) %>% 
  mutate(name = case_when(gene == "PCK2" ~ "JB2 up\n(PCK2)",
                          gene == "REST" ~ "JB6 down\n(REST)",
                          gene == "KDM6A" ~ "JB6 up\n(KDM6A)",
                          gene == "BZW2" ~ "JB2 down\n(BZW2)")) %>% 
  mutate(name = fct_relevel(name, "JB2 up\n(PCK2)", "JB2 down\n(BZW2)",
                            "JB6 up\n(KDM6A)", "JB6 down\n(REST)")) %>% 
  mutate(value = as.numeric(value)) %>% 
  na.omit()

ggplot(thisFig, aes(x = treat, y = value, fill = name)) +
  geom_col(alpha = 0.75) +
  scale_fill_manual(values = c(JB2Ucol, JB2Dcol, JB6Ucol, JB6Dcol)) +
  scale_x_discrete(labels = function(x) str_replace(x, "_24h", "")) +
  labs(x = "", y = "Fold change\nfrom DMSO (Log2)") +
  facet_wrap(~name, nrow = 1) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave("../results/Figure 6A.tiff", dpi = 300, height = 2, width = 4)

#plot the GGS scatter plots for MDA468 (Fig. 6B)
gradient <- DEG_24h %>% 
  mutate(gradient = case_when(genes %in% group_JB2U_C_100$genes ~ "JB2 up",
                              genes %in% group_JB2D_C_100$genes ~ "JB2 down",
                              genes %in% group_JB6U_C_100$genes ~ "JB6 up",
                              genes %in% group_JB6D_C_100$genes ~ "JB6 down",
                              TRUE ~ "Non-gradient"))

thisFig <- gradient %>% 
  select(gradient, genes, JB2_24h, JB6_24h) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient")))

#inset table showing the number of genes in the GGS
sum_tibble <- thisFig %>% group_by(gradient) %>% count() %>% dplyr::rename("Group" = "gradient")
forChart <- tibble(x = -3.5, y = 10, tb = list(sum_tibble))

#inset plot showing JB6/DMSO (y-axis) and JB2/DMSO(x-axis).
#(due to font issues, the axes titles were removed and added outside R using photoshop)
p1 <- ggplot(thisFig, aes(x = JB2_24h, y = JB6_24h, fill = gradient)) +
  geom_point(data = thisFig %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 2, shape = 21, stroke = 0) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
  guides(alpha = "none", colour = guide_legend(override.aes = list(size=4))) +
  scale_fill_manual(values = c("JB2 up" =  JB2Ucol, "JB2 down" = JB2Dcol, 
                               "JB6 up" = JB6Ucol, "JB6 down" = JB6Dcol, "Non-gradient" = "grey70")) +
  scale_x_continuous(limits = c(-6, 6)) +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(x = "JB2 vs DMSO", y = "JB6 vs DMSO", title = "", color = "") +
  guides(size = FALSE) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        text = element_blank(),
        panel.grid = element_blank())

#volcano plot for MDA468, with annotations for selected KEGG pathways (using enrichr)
take <- c("Arginine biosynthesis", "TCA cycle", "Lysosome", "Alanine, aspartate and glutamate metabolism",
          "Arginine and proline metabolism", "Oxidative phosphorylation", "Glutathione metabolism", "Autophagy")
take_string <- paste0(take, collapse = "|")
library <- "KEGG_2021_Human"
JB2U_KEGG <- EnrichrFun(group_JB2U_C_100$genes, library)
JB2D_KEGG <- EnrichrFun(group_JB2D_C_100$genes, library)
GeneList <- list(JB2U_KEGG$`long results` %>% mutate(group = "JB2U"), 
                 JB2D_KEGG$`long results` %>% mutate(group = "JB2D"))
annotations <- purrr::reduce(GeneList, bind_rows) %>% 
  filter(str_detect(Terms, take_string)) %>% 
  mutate(Terms = str_extract(Terms, take_string))

fdr <- as_tibble(efit_24h[["p.value"]], rownames = "genes") %>% 
  select(genes, JB6_to_JB2_24h) %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pval")) %>% 
  mutate(JB6_to_JB2_24h_fdr = p.adjust(JB6_to_JB2_24h_pval, method = "fdr"))

thisFig2 <- gradient %>% 
  dplyr::select(gradient, genes, JB6_to_JB2_24h, JB2_24h) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient"))) %>% 
  left_join(fdr, by = "genes") %>% 
  mutate(mlogp = -log10(JB6_to_JB2_24h_fdr)) %>% 
  left_join(annotations, by = c("genes" = "Gene")) 

ggplot(thisFig2, aes(x = -JB6_to_JB2_24h, y = mlogp)) +
  geom_point(data = thisFig2 %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig2 %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 2, shape = 21, stroke = 0) +
  geom_table(data = forChart, aes(x, y, label = tb),
             stat = "fmt_tb", size = 5, parse = TRUE, 
             table.theme = ttheme_gtplain(base_size = 13, core=list(fg_params=list(hjust=0, x=0.1)))) +
  geom_text_repel(data = thisFig2 %>% filter(!is.na(Terms)), aes(label = genes, color = Terms), 
                  nudge_x = 1, nudge_y = 1, size = 5) +
  geom_text_repel(data = thisFig2 %>% filter(genes == "PDAP1"), aes(label = genes), 
                  nudge_x = 3, nudge_y = 1, size = 5, color = "black") +
  annotation_custom(ggplotGrob(p1), xmin = -3.8, xmax = -0.5, ymin = 3, ymax = 6.8) +
  guides(alpha = "none", color = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  scale_fill_manual(values = c(JB2Ucol, JB2Dcol, JB6Ucol, JB6Dcol, "grey70")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Fold change in JB2 vs. JB6 (Log2)", y = "P-value of JB2 vs. JB6 (-log10)", 
       title = "", color = "", fill = "") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.09,0.84),
        legend.background = element_blank(),
        legend.key.size = unit(0.60, "cm"),
        legend.text = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave("../results/Figure 6B.tiff", dpi = 300, height = 6, width = 6)

#same scatter plot as above, for HCC70 (Fig. S6A)
gradient <- DEG_H70 %>% 
  mutate(gradient = case_when(genes %in% group_JB2U_C_100$genes ~ "JB2 up",
                              genes %in% group_JB2D_C_100$genes ~ "JB2 down",
                              genes %in% group_JB6U_C_100$genes ~ "JB6 up",
                              genes %in% group_JB6D_C_100$genes ~ "JB6 down",
                              TRUE ~ "Non-gradient"))

thisFig <- gradient %>% 
  select(gradient, genes, JB2_H70, JB6_H70) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient")))

sum_tibble <- thisFig %>% group_by(gradient) %>% count() %>% dplyr::rename("Group" = "gradient")
forChart <- tibble(x = 6, y = -6, tb = list(sum_tibble))

#inset plot - HCC70
p1 <- ggplot(thisFig, aes(x = JB2_H70, y = JB6_H70, fill = gradient)) +
  geom_point(data = thisFig %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 2, shape = 21, stroke = 0) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", size = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
  guides(alpha = "none", colour = guide_legend(override.aes = list(size=4))) +
  scale_fill_manual(values = c("JB2 up" =  JB2Ucol, "JB2 down" = JB2Dcol, 
                               "JB6 up" = JB6Ucol, "JB6 down" = JB6Dcol, "Non-gradient" = "grey70")) +  scale_x_continuous(limits = c(-6, 6)) +
  scale_y_continuous(limits = c(-6, 6)) +
  labs(x = "JB2 vs DMSO", y = "JB6 vs DMSO", title = "", color = "") +
  guides(size = FALSE) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none",
        text = element_blank(),
        panel.grid = element_blank())

#volcano plot - HCC70
take <- c("Arginine biosynthesis", "TCA cycle", "Lysosome", "Alanine, aspartate and glutamate metabolism",
          "Arginine and proline metabolism", "Oxidative phosphorylation", "Glutathione metabolism", "Autophagy")
take_string <- paste0(take, collapse = "|")
library <- "KEGG_2021_Human"
JB2U_KEGG <- EnrichrFun(group_JB2U_C_100$genes, library)
JB2D_KEGG <- EnrichrFun(group_JB2D_C_100$genes, library)
GeneList <- list(JB2U_KEGG$`long results` %>% mutate(group = "JB2U"), 
                 JB2D_KEGG$`long results` %>% mutate(group = "JB2D"))
annotations <- purrr::reduce(GeneList, bind_rows) %>% 
  filter(str_detect(Terms, take_string)) %>% 
  mutate(Terms = str_extract(Terms, take_string))

fdr <- as_tibble(efit_H70[["p.value"]], rownames = "genes") %>% 
  select(genes, JB6_to_JB2_H70) %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pval")) %>% 
  mutate(JB6_to_JB2_H70_fdr = p.adjust(JB6_to_JB2_H70_pval, method = "fdr"))

thisFig2 <- gradient %>% 
  dplyr::select(gradient, genes, JB6_to_JB2_H70, JB2_H70) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient"))) %>% 
  left_join(fdr, by = "genes") %>% 
  mutate(mlogp = -log10(JB6_to_JB2_H70_fdr)) %>% 
  left_join(annotations, by = c("genes" = "Gene")) 

forChart <- tibble(x = -4.3, y = 10, tb = list(sum_tibble))

ggplot(thisFig2, aes(x = -JB6_to_JB2_H70, y = mlogp)) +
  geom_point(data = thisFig2 %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig2 %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 2, shape = 21, stroke = 0) +
  geom_table(data = forChart, aes(x, y, label = tb),
             stat = "fmt_tb", size = 5, parse = TRUE, 
             table.theme = ttheme_gtplain(base_size = 13, core=list(fg_params=list(hjust=0, x=0.1)))) +
  geom_text_repel(data = thisFig2 %>% filter(!is.na(Terms)), aes(label = genes, color = Terms), 
                  nudge_x = 2, nudge_y = 0, size = 5) +
  geom_text_repel(data = thisFig2 %>% filter(genes == "PDAP1"), aes(label = genes), 
                  nudge_x = -1, nudge_y = 1, size = 5, color = "black") +
  annotation_custom(ggplotGrob(p1), xmin = 0, xmax = 3.5, ymin = 6.5, ymax = 10) +
  guides(alpha = "none", color = "none",
         fill = guide_legend(override.aes = list(size = 4))) +
  scale_fill_manual(values = c(JB2Ucol, JB2Dcol, JB6Ucol, JB6Dcol, "grey70")) +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Fold change in JB2 vs. JB6 (Log2)", y = "P-value of JB2 vs. JB6 (-log10)", 
       title = "", color = "", fill = "") +
  theme_bw(base_size = 16) +
  theme(legend.position = c(0.09,0.84),
        legend.background = element_blank(),
        legend.key.size = unit(0.60, "cm"),
        legend.text = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 16, color = "black"),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggsave("../results/Figure S6A.tiff", dpi = 300, height = 6, width = 6)
