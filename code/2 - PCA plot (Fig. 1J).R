#This script will create the PCA plot, using info from all 3 count matrices.

#To create PCA plot for all the RNAseq samples, I combined all count matrix into one
cts_combined <- read.csv("../data/RNAseq count matrix combined all sets.csv", row.names = "Gene")
V_combined <- DEGCreate(cts_combined, name = "combined both sets", plot = FALSE)
normCPM <- t(V_combined$E)

#PCA plot - for MDA468 cells (Fig. 1J)
thisFig <- as_tibble(normCPM, rownames = "sample") %>% 
  filter(str_detect(sample, "MDA468")) %>% 
  mutate(sample = str_replace(sample, "MDA468_", "")) %>% 
  mutate(sample = str_replace(sample, "_a|_b", "")) %>% 
  mutate(timepoint = str_extract(sample, "[:digit:]+h$")) %>% 
  select(sample, timepoint, everything())

PCA <- prcomp(thisFig[,c(-1,-2)])

p <- fviz_pca_ind(PCA, axes = c(1, 2), geom = "point", axes.linetype = "blank") +
  geom_point(aes(color = thisFig$timepoint), size = 3) +
  geom_text_repel(aes(label = str_replace(thisFig$sample, "_[:graph:]*$", ""), color = thisFig$timepoint),
                  show.legend = FALSE, size = 5) +
  scale_color_manual(values = c("dodgerblue4", "red")) +
  labs(x = paste0("PC1 (", round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100, 1), "%)"),
       y = paste0("PC2 (", round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100, 1), "%)"),
      title = "", color = "Time point") +
  theme_classic(base_size = 16) +
  theme(legend.position = c(0.15, 0.85),
        legend.background = element_rect(color = "black"),
        axis.text = element_text(size = 16, color = "black"))
p
ggsave("../results/Figure 1J - PCA plot.tiff", dpi = 200, height = 6, width = 6)

#Measure distance in PCA (Fig. S1H)
Dim1 = (PCA$sdev[1]^2 / sum(PCA$sdev^2))*100
Dim2 = (PCA$sdev[2]^2 / sum(PCA$sdev^2))*100

distances <- p[["data"]] %>% 
  cbind("ID" = thisFig$sample) %>% 
  mutate(xw = x * Dim1, yw = y * Dim2) %>% 
  group_by(ID) %>% 
  summarise(mean_xw = mean(xw), mean_yw = mean(yw)) %>% 
  column_to_rownames(var = "ID") %>% 
  dist(method = "euclidean") %>% 
  as.matrix() %>% 
  as_tibble(rownames = "ID1") %>% 
  gather(-ID1, key = "ID2", value = "distance") %>% 
  mutate(timepoint1 = str_extract(ID1, "[:digit:]+h$")) %>% 
  mutate(timepoint2 = str_extract(ID2, "[:digit:]+h$")) %>% 
  mutate(treat1 = str_extract(ID1, "^[:alnum:]*")) %>% 
  mutate(treat2 = str_extract(ID2, "^[:alnum:]*")) %>%  
  filter(treat1 %in% c("JB2", "JB3", "JB6")) %>% 
  filter(treat2 %in% c("B2", "B3", "B6")) %>% 
  filter(timepoint1 == timepoint2) %>% 
  filter(str_extract(treat1, "B\\d$") == treat2) %>% 
  mutate(distance = distance / 1000) %>% 
  mutate(name = str_c(treat1, "\u2192", treat2))

ggplot(distances, aes(x = name, y = distance, fill = timepoint1, alpha = name)) +
  geom_col() +
  scale_alpha_manual(values = c(0.9, 0.6, 0.3)) +
  scale_fill_manual(values = c("dodgerblue4", "red")) +
  labs(x = "", y = "Distance on PCA", fill = "Time point") +
  facet_row(~timepoint1, scales = "free_x") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.text = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave("../results/Figure S1H - Distances in PCA plot.tiff", dpi = 300, height = 4, width = 2.5)
