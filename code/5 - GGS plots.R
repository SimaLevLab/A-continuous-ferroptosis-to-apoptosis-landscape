#Scripts written by Yaron Vinik
#This script will create the plots of the Gradient Gene Set (Fig. 2A, B)

#define gradient colors for the plots
JB2Ucol <- "#4DAF4A"
JB2Dcol <- "salmon"
JB6Ucol <- "#377EB8"
JB6Dcol <- "mediumorchid"
colors <- c(JB2Ucol, JB2Dcol, JB6Ucol, JB6Dcol)

#plot mean expression of all genes in the GGS
thisFig <- DEG_24h %>% 
  as_tibble() %>% 
  mutate(group = case_when(genes %in% group_JB2U_C_100$genes ~ "JB2 up",
                           genes %in% group_JB2D_C_100$genes ~ "JB2 down",
                           genes %in% group_JB6U_C_100$genes ~ "JB6 up",
                           genes %in% group_JB6D_C_100$genes ~ "JB6 down")) %>% 
  filter(!is.na(group)) %>% 
  select(group, JB2_24h, JB3_24h, JB6_24h) %>%
  gather(-group, key = "treat", value = "value") %>% 
  group_by(group, treat) %>% 
  summarise(mean = mean(value)) %>% 
  mutate(treat = factor(treat, levels = c("JB2_24h", "JB3_24h", "JB6_24h"))) %>% 
  mutate(group = factor(group, levels = c("JB6 up", "JB2 up", "JB6 down", "JB2 down"))) %>%
  mutate(group1 = str_extract(group, "JB2|JB6")) %>% 
  mutate(group2 = str_extract(group, "up|down")) %>% 
  mutate(group1 = fct_relevel(group1, "JB6", "JB2")) %>% 
  mutate(group2 = fct_relevel(group2, "up", "down"))

ggplot(thisFig, aes(x = treat, y = mean, fill = group)) +
  geom_col(alpha = 0.85) +
  geom_label(data = thisFig %>% filter(group2 == "up"), aes(label = group, color = group), 
             fill = "white", x = 2, y = 0.2, fontface = "bold", size = 5.5) +
  geom_label(data = thisFig %>% filter(group2 == "down"), aes(label = group, color = group), 
             fill = "white", x = 2, y = -0.2, fontface = "bold", size = 5.5) +
  scale_fill_manual(values = c("JB6 up" = JB6Ucol, "JB2 up" = JB2Ucol, "JB6 down" = JB6Dcol, "JB2 down" = JB2Dcol)) +
  scale_color_manual(values = c("JB6 up" = JB6Ucol, "JB2 up" = JB2Ucol, "JB6 down" = JB6Dcol, "JB2 down" = JB2Dcol)) +
  scale_x_discrete(labels = function(x) str_replace(x, "_24h", "")) +
  labs(x = "", y = "Fold change\nfrom DMSO (Log2)") +
  facet_grid(group2 ~ group1, scales = "free_y") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.text = element_text(size = 18, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("../results/Figure 2A - GGS expression plot.tiff", dpi = 300, height = 4, width = 4)

#plot the GGS scatter plots for MDA468 (Fig. 2B)
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

#volcano plot for MDA468, with annotations for selected KEGG pathways (using enrichr)
selected_genes <- c("DDIT4", "PCK2", "MT2A", "CKB", "ASL", "RRM2", "GGA3", "CHAC2", "IDS", "PTEN", "GPX8", "PDAP1")

fdr <- as_tibble(efit_24h[["p.value"]], rownames = "genes") %>% 
  select(genes, JB6_to_JB2_24h) %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pval")) %>% 
  mutate(JB6_to_JB2_24h_fdr = p.adjust(JB6_to_JB2_24h_pval, method = "fdr"))

thisFig2 <- gradient %>% 
  dplyr::select(gradient, genes, JB6_to_JB2_24h, JB2_24h) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient"))) %>% 
  left_join(fdr, by = "genes") %>% 
  mutate(mlogp = -log10(JB6_to_JB2_24h_fdr))

ggplot(thisFig2, aes(x = -JB6_to_JB2_24h, y = mlogp)) +
  geom_point(data = thisFig2 %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig2 %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 3, shape = 21, stroke = 0) +
  geom_table(data = forChart, aes(x, y, label = tb),
             stat = "fmt_tb", size = 5, parse = TRUE, 
             table.theme = ttheme_gtplain(base_size = 13, core=list(fg_params=list(hjust=0, x=0.1)))) +
  geom_text_repel(data = thisFig2 %>% filter(genes %in% selected_genes), aes(label = genes), 
                  nudge_x = 1, nudge_y = 1, size = 5, color = "black") +
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

ggsave("../results/Figure 2B - Volcano plot of GGS genes in MDA468.tiff", dpi = 300, height = 6, width = 6)

#same scatter plot as above, for HCC70 (Fig. S2A)
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

#volcano plot - HCC70
fdr <- as_tibble(efit_H70[["p.value"]], rownames = "genes") %>% 
  select(genes, JB6_to_JB2_H70) %>% 
  rename_with(.cols = -genes, .fn = ~str_c(.x, "_pval")) %>% 
  mutate(JB6_to_JB2_H70_fdr = p.adjust(JB6_to_JB2_H70_pval, method = "fdr"))

thisFig2 <- gradient %>% 
  dplyr::select(gradient, genes, JB6_to_JB2_H70, JB2_H70) %>% 
  mutate(gradient = factor(gradient, levels = c("JB2 up", "JB2 down", "JB6 up", "JB6 down",
                                                "Non-gradient"))) %>% 
  left_join(fdr, by = "genes") %>% 
  mutate(mlogp = -log10(JB6_to_JB2_H70_fdr))

forChart <- tibble(x = -4.3, y = 10, tb = list(sum_tibble))

ggplot(thisFig2, aes(x = -JB6_to_JB2_H70, y = mlogp)) +
  geom_point(data = thisFig2 %>% filter(gradient == "Non-gradient"), alpha = 0.5, color = "gray") +
  geom_point(data = thisFig2 %>% filter(gradient != "Non-gradient"), aes(fill = gradient), 
             size = 3, shape = 21, stroke = 0) +
  geom_table(data = forChart, aes(x, y, label = tb),
             stat = "fmt_tb", size = 5, parse = TRUE, 
             table.theme = ttheme_gtplain(base_size = 13, core=list(fg_params=list(hjust=0, x=0.1)))) +
  geom_text_repel(data = thisFig2 %>% filter(genes %in% selected_genes), aes(label = genes), 
                  nudge_x = 1, nudge_y = 0, size = 5, color = "black") +
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

ggsave("../results/Figure S2A - Volcano plot of GGS genes in HCC70.tiff", dpi = 300, height = 6, width = 6)
