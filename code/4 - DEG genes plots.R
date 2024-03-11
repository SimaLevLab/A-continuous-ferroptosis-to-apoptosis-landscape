#Scripts written by Yaron Vinik
#create the genes bar graphs (Figure 5H, 6G, 7I)

#load genesets used to annotate genes in these figures
C2 <- gmt2list("../data/c2.all.v7.2.symbols.gmt") #load MsigDB genesets, to assign functional groups to genes
genesets <- gmt2list("../data/genesets used in the paper.gmt") #other genesets used in this study

#-------------------------------------------------------------------------------------------
# Ferroptosis related JB2 exclusive genes (Fig. 5H)
#-------------------------------------------------------------------------------------------

#1. Ferroptosis related genes, based on the Venn diagrams:
pick <- c(SynGenes_both, SynGenes_both_down, SynGenes_468_6h_up, SynGenes_468_6h_down) #create these objects in the "3 - DEG analysis.R" script

BH4_synthesis <- tibble(group = "BH4 synthesis", genes = C2[["REACTOME_TETRAHYDROBIOPTERIN_BH4_SYNTHESIS_RECYCLING_SALVAGE_AND_REGULATION"]])
GSH_metabolism <- tibble(group = "GSH metabolism", genes = unique(c(C2[["WP_GLUTATHIONE_METABOLISM"]], C2[["REACTOME_GLUTATHIONE_SYNTHESIS_AND_RECYCLING"]], "SLC7A11")))
Iron_transport <- tibble(group = "Iron matabolism", genes = C2[["REACTOME_IRON_UPTAKE_AND_TRANSPORT"]])
Lipids_beta_oxidation <- tibble(group = "Lipids beta oxidation", genes = C2[["WP_FATTY_ACID_BETA_OXIDATION"]])
Lipids <- tibble(group = "Lipid metabolism", genes = C2[["REACTOME_METABOLISM_OF_LIPIDS"]]) %>% 
  filter(genes %in% pick) 
Ferroptosis <- tibble(group = "Other ferroptosis related", genes = c(C2[["WP_FERROPTOSIS"]], "AIFM2"))
FerrDB1 <- tibble(group = "FerrDB drivers", genes = genesets[["FerrDB Drivers"]])
FerrDB2 <- tibble(group = "FerrDB suppressors", genes = genesets[["FerrDB Suppresor"]])
FerrDB3 <- tibble(group = "FerrDB markers", genes = genesets[["FerrDB Marker"]])
FerrDB <- bind_rows(FerrDB1, FerrDB2, FerrDB3)
take <- bind_rows(FerrDB, BH4_synthesis, GSH_metabolism, Iron_transport, Ferroptosis) %>% 
  distinct(genes, .keep_all = TRUE) %>% 
  filter(genes %in% pick) %>% 
  filter(!genes %in% Lipids$genes)

#add GO annotations, using enrichR
Enrichr_result1 <- EnrichrFun(genes = take$genes, library = "GO_Molecular_Function_2021") 
Enrichr_result2 <- EnrichrFun(genes = take$genes, "GO_Biological_Process_2021")

take_ferroptosis <- take %>% 
  left_join(Enrichr_result1$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(Enrichr_result2$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(DEG_combined) %>% 
  mutate(group = ifelse(str_detect(group, "FerrDB"), 
                        case_when(str_detect(Terms.y, "GO:0034599") ~ "Oxidative stress",
                                  str_detect(Terms.y, "starvation") ~ "Response to starvation",
                                  str_detect(Terms.y, "iron") ~ "Iron matabolism",
                                  str_detect(Terms.y, "GO:0034599") ~ "Oxidative stress",
                                  TRUE ~ group),
                        group)) %>% 
  mutate(group = ifelse(str_detect(group, "FerrDB"), "Other ferroptosis related", group))

pval_ferroptosis <- take_ferroptosis %>% 
  select(group, genes, JB2_6h_fdr, JB2_24h_fdr) %>% 
  gather(-group, -genes, key = "treat", value = "pval") %>% 
  mutate(pval = ifelse(pval != 0, "*", "")) %>% 
  mutate(treat = str_replace(treat, "_fdr", "")) %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h"))

thisFig_ferroptosis <- take_ferroptosis %>% 
  select(group, genes, JB2_6h, JB2_24h) %>% 
  na.omit() %>% 
  gather(-group, -genes, key = "treat", value = "value") %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h")) %>% 
  left_join(pval_ferroptosis, by = c("group", "genes", "treat", "cell_line")) %>% 
  mutate(cell_line = factor(cell_line, levels = c("MDA468 6h", "MDA468 24h", "HCC70 24h"))) %>% 
  mutate(group = fct_relevel(group, "Other ferroptosis related", after = Inf)) 

gene_sorter_ferroptosis <- thisFig_ferroptosis %>% 
  select(group, genes, cell_line, value) %>% 
  filter(cell_line == "MDA468 24h") %>% 
  arrange(group, value) %>% 
  mutate(FerrDB = case_when(genes %in% FerrDB1$genes ~ "Ind",
                            genes %in% FerrDB2$genes ~ "Supp",
                            genes %in% FerrDB3$genes ~ "Mark"))
thisFig_ferroptosis$genes <- factor(thisFig_ferroptosis$genes, levels = gene_sorter_ferroptosis$genes, ordered = TRUE)

p1 <- ggplot(thisFig_ferroptosis, aes(x = value, y = reorder(genes, value), fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = gene_sorter_ferroptosis %>% filter(!is.na(FerrDB)), 
             aes(color = FerrDB), x = -20.5, size = 3, show.legend = FALSE) +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

#2. Lipids related genes
Enrichr_result1 <- EnrichrFun(genes = Lipids$genes, "GO_Molecular_Function_2021", "temp") 
Enrichr_result2 <- EnrichrFun(genes = Lipids$genes, "GO_Biological_Process_2021", "temp")

take_lipids <- Lipids %>% 
  left_join(Enrichr_result1$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(Enrichr_result2$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(DEG_combined) %>% 
  mutate(group = case_when(str_detect(Terms.x, "acyltransferase") ~ "Acyltransferase activity",
                           str_detect(Terms.y, "GO:0006637") ~ "Acyl-CoA metabolic process", 
                           str_detect(Terms.y, "GO:0006695") ~ "Cholesterol biosynthesis", 
                           str_detect(Terms.y, "arachidonic") ~ "Arachidonic acid", 
                           str_detect(Terms.y, "GO:0006635") ~ "Beta-oxidation", 
                           str_detect(Terms.y, "GO:0001561") ~ "Alpha-oxidation", 
                           str_detect(Terms.y, "GO:0034626") ~ "PUFA elongation", 
                           str_detect(Terms.y, "fatty-acyl-CoA") ~ "Fatty-acyl-CoA biosynthesis",
                           str_detect(Terms.y, "GO:0030148") ~ "Sphingolipid biosynthesis", 
                           str_detect(Terms.y, "GO:0046474") ~ "Glycerophospholipid biosynthesis", 
                           str_detect(Terms.y, "GO:0006665") ~ "Sphingolipid metabolic process", 
                           TRUE ~ "Other lipids related")) 

pval_lipids <- take_lipids %>% 
  select(group, genes, JB2_6h_fdr, JB2_24h_fdr) %>% 
  gather(-group, -genes, key = "treat", value = "pval") %>% 
  mutate(pval = ifelse(pval != 0, "*", "")) %>% 
  mutate(treat = str_replace(treat, "_fdr", "")) %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h"))

thisFig_lipids <- take_lipids %>% 
  select(group, genes, JB2_6h, JB2_24h) %>% 
  filter(!is.na(JB2_24h)) %>% 
  gather(-group, -genes, key = "treat", value = "value") %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h")) %>% 
  left_join(pval_lipids, by = c("group", "genes", "treat", "cell_line")) %>% 
  mutate(cell_line = factor(cell_line, levels = c("MDA468 6h", "MDA468 24h", "HCC70 24h"))) %>% 
  mutate(group = fct_relevel(group, "Other lipids related", after = Inf)) 

gene_sorter_lipids <- thisFig_lipids %>% 
  select(group, genes, cell_line, value) %>% 
  filter(cell_line == "MDA468 24h") %>% 
  arrange(group, value) %>% 
  mutate(FerrDB = case_when(genes %in% FerrDB1$genes ~ "Ind",
                            genes %in% FerrDB2$genes ~ "Supp",
                            genes %in% FerrDB3$genes ~ "Mark"))

thisFig_lipids$genes <- factor(thisFig_lipids$genes, levels = gene_sorter_lipids$genes, ordered = TRUE)
faceting_groups <- c("PUFA elongation", "Sphingolipid biosynthesis", "Sphingolipid metabolic process", "Other lipids related")

p2 <- ggplot(thisFig_lipids %>% filter(!group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = gene_sorter_lipids %>% filter(!group %in% faceting_groups) %>% filter(!is.na(FerrDB)), 
             aes(color = FerrDB), x = -11, size = 3, show.legend = FALSE) +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p3 <- ggplot(thisFig_lipids %>% filter(group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = gene_sorter_lipids %>% filter(group %in% faceting_groups) %>% filter(!is.na(FerrDB)), 
             aes(color = FerrDB), x = -15, size = 3, show.legend = FALSE) +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggarrange(p1, p2, p3, nrow = 1, align = "h")
ggsave("../results/Figure 5H - List of JB2 exclusive genes.tiff", dpi = 300, height = 8, width = 20)

#-------------------------------------------------------------------------------------------
# Lysosomal genes in JB2 (Fig. 6G)
#-------------------------------------------------------------------------------------------

LysosomeBiogenesis <- c("LGMN", "CTSV", "CTSB", "TPP1", "CTSD", "GLA", "NEU1", "GLB1", "GUSB", 
                        "NAGA", "IDS", "LYPLA2", "LIPA", "ACP2", "PSAP", "CTNS", "MFSD12", "NPC2", "LAMP1")
pick <- DEG_combined %>% filter(JB2_24h_fdr == 1) %>% pull(genes)
Lysosome <- tibble(genes = genesets[["Lysosome"]]) %>% filter(genes %in% pick)
Enrichr_result1 <- EnrichrFun(genes = Lysosome$genes, library = "GO_Molecular_Function_2021") 
Enrichr_result2 <- EnrichrFun(genes = Lysosome$genes, library = "GO_Biological_Process_2021")
Enrichr_result3 <- EnrichrFun(genes = Lysosome$genes, library = "ChEA_2016")

take <- Lysosome %>% 
  left_join(Enrichr_result1$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(Enrichr_result2$`long results`, by = c("genes" = "Gene"))  %>% 
  left_join(Enrichr_result3$`long results`, by = c("genes" = "Gene"))  %>% 
  left_join(DEG_combined) %>% 
  mutate(group = case_when(str_detect(Terms.y, "GO:0006508") ~ "Proteolysis", 
                           str_detect(Terms.x, "ATPase") ~ "ATPase binding/activity", 
                           str_detect(Terms.x, "GO:1990459") ~ "Iron metabolism",
                           str_detect(Terms.x, "GO:0015485") ~ "Cholesterol binding", 
                           str_detect(Terms.y, "GO:1990928") ~ "Response to amino acid starvation", 
                           str_detect(Terms.y, "GO:0010506") ~ "Regulation of autophagy", 
                           str_detect(Terms.y, "iron") ~ "Iron metabolism", 
                           str_detect(Terms.y, "transport") ~ "Transport", 
                           genes %in% LysosomeBiogenesis ~ "Lysosome Biogenesis",
                           TRUE ~ "Other"), .after = genes) %>% 
  mutate(MITF = ifelse(str_detect(Terms, "MITF"), "MITF", NA), .after = genes) %>% 
  mutate(TFEB = ifelse(str_detect(Terms, "TFEB"), "TFEB", NA), .after = MITF)

pval <- take %>% 
  select(group, genes, MITF, TFEB, JB2_6h_fdr, JB2_24h_fdr) %>% 
  gather(-group, -genes, -MITF, -TFEB, key = "treat", value = "pval") %>% 
  mutate(pval = ifelse(pval != 0, "*", "")) %>% 
  mutate(treat = str_replace(treat, "_fdr", "")) %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h"))

thisFig <- take %>% 
  select(group, genes, MITF, TFEB, JB2_6h, JB2_24h) %>% 
  filter(!is.na(JB2_24h)) %>% 
  gather(-group, -genes, -MITF, -TFEB, key = "treat", value = "value") %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h")) %>% 
  left_join(pval, by = c("group", "genes", "MITF", "TFEB", "treat", "cell_line")) %>% 
  mutate(cell_line = factor(cell_line, levels = c("MDA468 6h", "MDA468 24h", "HCC70 24h"))) %>% 
  mutate(group = fct_relevel(group, "Other", after = Inf)) 

gene_sorter <- thisFig %>% 
  select(group, genes, cell_line, value, MITF, TFEB) %>% 
  filter(cell_line == "MDA468 24h") %>% 
  arrange(group, value) 

thisFig$genes <- factor(thisFig$genes, levels = gene_sorter$genes, ordered = TRUE)
faceting_groups <- c("Transport", "Other")

p1 <- ggplot(thisFig %>% filter(!group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_point(data = gene_sorter %>% filter(!group %in% faceting_groups) %>% filter(!is.na(MITF)), 
             x = -10, color = "red", size = 3, show.legend = FALSE) +
  geom_point(data = gene_sorter %>% filter(!group %in% faceting_groups) %>% filter(!is.na(TFEB)), 
             x = -10.5, color = "blue", size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p2 <- ggplot(thisFig %>% filter(group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_point(data = gene_sorter %>% filter(group %in% faceting_groups) %>% filter(!is.na(MITF)), 
             x = -8, color = "red", size = 3, show.legend = FALSE) +
  geom_point(data = gene_sorter %>% filter(group %in% faceting_groups) %>% filter(!is.na(TFEB)), 
             x = -8.5, color = "blue", size = 3, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggarrange(p1, p2, align = "hv")
ggsave("../results/Figure 6G - List of lysosome related genes.tiff", dpi = 300, height = 9, width = 12)

#-------------------------------------------------------------------------------------------
# ATF4-only targets (Figure 7I)
#-------------------------------------------------------------------------------------------

pick <- DEG_combined %>% filter(JB2_24h_fdr == 1) %>% pull(genes)
ATF4 <- genesets[["ATF4 targets"]] %>% as_tibble() %>% 
  dplyr::rename("genes" = "value") %>% 
  filter(genes %in% pick)
Enrichr_result1 <- EnrichrFun(genes = ATF4$genes, "GO_Molecular_Function_2021", "temp") 
Enrichr_result2 <- EnrichrFun(genes = ATF4$genes, "GO_Biological_Process_2021", "temp")
Enrichr_result3 <- EnrichrFun(genes = ATF4$genes, "GO_Cellular_Component_2021", "temp")

take <- ATF4 %>% 
  left_join(Enrichr_result1$`long results`, by = c("genes" = "Gene")) %>% 
  left_join(Enrichr_result2$`long results`, by = c("genes" = "Gene"))  %>% 
  left_join(Enrichr_result3$`long results`, by = c("genes" = "Gene"))  %>% 
  left_join(DEG_combined) %>% 
  mutate(group = case_when(str_detect(Terms.y, "glutamine|glutamate") ~ "Glutamine/gluthamate metabolism",
                           str_detect(Terms.y, "GO:0034599") ~ "Response to Oxidative stress",
                           str_detect(Terms.y, "cysteine|serine|glutathione|aspartate") ~ "Serine / cysteine / aspartate metabolism",
                           str_detect(Terms.y, "glucose") ~ "Cellular glucose homeostasis",
                           str_detect(Terms.y, "GO:0034976") ~ "Response to ER stress",
                           str_detect(Terms.y, "(GO:0006357)") ~ "Regulation of transcription",
                           str_detect(Terms.y, "(GO:0006865)") ~ "Amino acid transport",
                           str_detect(Terms.y, "(GO:0046474)") ~ "Glycerophospholids",
                           str_detect(Terms.y, "(GO:0042304)") ~ "Glycerophospholids",
                           str_detect(Terms.y, "(GO:0008652)") ~ "Amino acid biosynthesis",
                           str_detect(Terms.y, "lipid") ~ "Lipids",
                           str_detect(Terms.y, "oxidative stress") ~ "Oxidative stress",
                           str_detect(Terms.y, "oxygen species") ~ "Oxidative stress",
                           TRUE ~ "Other"), .after = genes) 

pval <- take %>% 
  select(group, genes, JB2_6h_fdr, JB2_24h_fdr) %>% 
  gather(-group, -genes, key = "treat", value = "pval") %>% 
  mutate(pval = ifelse(pval != 0, "*", "")) %>% 
  mutate(treat = str_replace(treat, "_fdr", "")) %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h"))

thisFig <- take %>% 
  select(group, genes, JB2_6h, JB2_24h) %>% 
  filter(!is.na(JB2_24h)) %>% 
  gather(-group, -genes, key = "treat", value = "value") %>% 
  mutate(cell_line = case_when(str_detect(treat, "H70") ~ "HCC70 24h",
                               str_detect(treat, "6h") ~ "MDA468 6h",
                               str_detect(treat, "24h") ~ "MDA468 24h")) %>% 
  left_join(pval, by = c("group", "genes", "treat", "cell_line")) %>% 
  mutate(cell_line = factor(cell_line, levels = c("MDA468 6h", "MDA468 24h", "HCC70 24h"))) %>% 
  mutate(group = fct_relevel(group, "Other", after = Inf)) 

gene_sorter <- thisFig %>% 
  select(group, genes, cell_line, value) %>% 
  filter(cell_line == "MDA468 24h") %>% 
  arrange(group, value) 

thisFig$genes <- factor(thisFig$genes, levels = gene_sorter$genes, ordered = TRUE)
faceting_groups <- c("Other")

p1 <- ggplot(thisFig %>% filter(!group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

p2 <- ggplot(thisFig %>% filter(group %in% faceting_groups), aes(x = value, y = genes, fill = group)) +
  geom_col() +
  geom_text(aes(label = pval), size = 5) +
  scale_fill_brewer(palette = "Set2") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_grid(group~cell_line, scales = "free_y", space = "free", labeller = label_wrap_gen(width = 20)) +
  labs(x = "Fold change vs. DMSO (Log2)", y = "", fill = "Function") +
  coord_cartesian(clip = "off") +
  theme_gray(base_size = 16) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle = 0, vjust = 1, hjust = 0),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

ggarrange(p1, p2, align = "hv")
ggsave("../results/Figure 7I - List of ATF4 target genes.tiff", dpi = 300, height = 8, width = 12)
