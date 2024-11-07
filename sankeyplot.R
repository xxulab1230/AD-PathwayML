#### Load necessary libraries for color palette and plotting
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(randomcoloR)

# Define custom color palettes for the Sankey plot
mycol <- brewer.pal(12, 'Paired')  # Preset 12-color palette
mycol2 <- colorRampPalette(c("#999FD0", "#93BCDC", "#ADD2D7", "#F3DAE9", "#E7BAD7"))(
  length(unique(expanded_data[expanded_data$count >= 5,]$Gene)) + 
    length(unique(expanded_data[expanded_data$count >= 5,]$classes))
)

# Filter data and reshape it to a long format for the plot
df <- expanded_data %>%
  filter(count >= 3) %>%  # Adjust filter based on dataset, e.g., rosmap:3, mayo:5, msbb:5, ADNI:4
  pivot_longer(cols = c("Gene", "classes"), 
               names_to = "stage", 
               values_to = "count")  # Reshape for Sankey plot stages

# Generate 60 distinct colors for the plot, avoiding duplicates
mycol2 <- distinctColorPalette(60)

# Create the PDF file for saving the plot
pdf("ADNI-msigdb-disease-gene_count2-sankey.pdf", height = 10, width = 10)

# Generate the Sankey plot
ggplot(data = df, aes(x = stage, y = count, group = node,
                      edge_id = edge_id, connector = connector)) +
  # Add edge information to the plot
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order = "ascending",  # Arrange by frequency in descending order
                                             v_space = "auto",
                                             width = 0.05),
                  show.legend = FALSE) +
  # Add node information to the plot
  geom_sankeynode(aes(fill = node, color = node),
                  position = position_sankey(order = "ascending",
                                             v_space = "auto",
                                             width = 0.05),
                  show.legend = FALSE) +
  scale_fill_manual(values = mycol2) +
  scale_color_manual(values = mycol2) +
  coord_cartesian(clip = "off") +  # Prevent clipping of elements beyond the plot boundary
  
  # Add labels to the left side ("from" connectors)
  geom_text(data = df %>% filter(connector == "from"),
            aes(label = node),
            stat = "sankeynode",
            position = position_sankey(width = 0.05,
                                       v_space = "auto",
                                       order = "ascending",
                                       nudge_x = -0.05),  # Left-align the text
            hjust = 1, cex = 2.8, color = "black") +
  
  # Add labels to the right side ("to" connectors)
  geom_text(data = df %>% filter(connector == "to"), aes(label = node),
            stat = "sankeynode",
            position = position_sankey(v_space = "auto",
                                       order = "ascending"),
            hjust = 0, cex = 2, color = "black") +  # Right-align the text
  
  # Customize theme to remove axes and background
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0, 3, 0, 0, unit = "cm"))

# Close the PDF device
dev.off()
