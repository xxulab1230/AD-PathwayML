# Load necessary libraries for data processing and visualization
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggplot2)
library(gghalves)
library(geomtextpath)

# Filter genes with a count above the threshold and subset the main dataset accordingly
gene <- result_summary[result_summary$gene_count >= 8,]$Gene
datExpr0sub <- datExpr0[rownames(datExpr0) %in% gene,]

# Transform the dataset to a longer format and filter for specific types
df <- datExpr0sub %>%
  mutate("ID" = rownames(.)) %>%
  pivot_longer(-ID) %>%
  mutate(type = map_chr(name, ~ str_split(.x, "@")[[1]][length(str_split(.x, "@")[[1]])])) %>%
  filter(type %in% c("Dementia", "CN")) %>%
  mutate(type = as.factor(type)) %>%
  mutate(value = as.numeric(value))

# Perform Wilcoxon test for each gene and adjust p-values using the Bonferroni method
df_pvalue <- df %>%
  group_by(ID) %>%
  wilcox_test(value ~ type) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  select(ID, p.adj)

# (Alternative) No adjustment for ADNI dataset
df_pvalue <- df %>%
  group_by(ID) %>%
  wilcox_test(value ~ type) %>%
  adjust_pvalue(p.col = "p", method = "none") %>%
  add_significance(p.col = "p.adj") %>%
  select(ID, p.adj)

# Join p-values back to the main data, label differentially expressed genes, and add p-value labels
dff <- df %>%
  left_join(df_pvalue, by = "ID") %>%
  mutate(group = case_when(p.adj < 0.05 ~ "Differentially-expressed", TRUE ~ "NS")) %>%
  mutate(p.adj = paste("p=", p.adj, sep = ""))

# Plot boxplot for each gene
pdf("./ADNI-degboxplot.pdf", height = 8, width = 8)
col <- c("Dementia" = "#5785C1", "CN" = "#00A08A")
dff %>%
  ggplot(aes(value, ID)) +
  geom_stripped_rows(odd = "grey90", xfrom = -0.5, xto = 15) +
  geom_boxplot(aes(color = type), position = position_dodge(0.5), width = 0.5, outliers = FALSE, staplewidth = 0.5, key_glyph = draw_key_smooth) +
  geom_text(data = dff %>% select(1, p.adj) %>% distinct(), aes(x = 13.8, y = ID, label = p.adj), size = 2.5, color = "black") +
  scale_x_continuous(limits = c(0, 14), breaks = c(2, 4, 6, 8, 10, 12)) +
  scale_y_discrete(expand = c(0, 0.5)) +
  coord_cartesian(clip = "off") +
  geom_tile(data = dff, aes(x = 12.5, y = ID, fill = group, width = 0.5)) +
  scale_fill_manual(values = c("Differentially-expressed" = "#7294D4", "NS" = "grey")) +
  scale_color_manual(labels = paste("<span style='color:", col, "'>", unique(dff$type), "</span>"), values = col) +
  labs(x = NULL, y = NULL) +
  guides(color = guide_legend(position = "top", direction = "vertical", order = 1, title = NULL,
                              theme = theme(legend.key.height = unit(0.5, "cm"),
                                            legend.key.width = unit(1, "cm"))),
         fill = guide_legend(position = "top", direction = "vertical", order = 2, title = NULL)) +
  theme_classic() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), unit = "cm"),
        legend.text = element_markdown(),
        legend.background = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.text = element_text(color = "black"))
dev.off()

# Circular plot ("环状云雨图") to visualize gene expression distribution
dff3 <- dff %>%
  group_by(ID) %>%
  mutate(n = n(), avg = mean(value)) %>%
  arrange(-avg)

pdf("./ADNI-degcircleplot.pdf", height = 8, width = 8)
dff3 %>%
  ggplot(aes(x = ID, y = value, fill = ID)) +
  geom_half_violin(alpha = 0.9, color = "black", width = 0.7) +
  geom_half_boxplot(alpha = 0.8, color = "black", side = "r", nudge = 0.1, outlier.shape = NA, width = 0.5) +
  geom_half_point(alpha = 0.6, color = "grey", side = "r", position = position_nudge(x = 0.3), width = 0.5) +
  scale_fill_manual(values = c("#4e6d58", "#749e89", "#abccbe", "#e3cacf", "#c399a2", "#9f6e71", "#41507b",
                               "#7d87b2", "#c2cae3", "#d9700e", "#e9a00e", "#5b7314", "#c3d6ce", "#89a6bb",
                               "#efc86e", "#97c684", "#6f9969", "#aab5d5", "#808fe1", "#5c66a8", "#454a74",
                               "#b5361c", "#1c9d7c", "#31c7ba", "#369cc9", "#3a507f")) +
  ylim(0, 12) +
  coord_curvedpolar(clip = "off") +
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8, hjust = 0.5, color = "black", vjust = 0, margin = margin(l = -20, r = -30, unit = "pt")),
        axis.text.y = element_blank(),
        axis.title = element_blank())
dev.off()
