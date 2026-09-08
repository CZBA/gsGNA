# Load required libraries
library(ggplot2)
library(tidyr)
library(Cairo)  # Ensure PDF font embedding and compatibility

# ====================== 1. Build the statistical summary table stats_summary ======================
stats_summary <- data.frame(
  Algorithm = c("MCODE", "infomap", "prop", "eigen", "louvain", "walktrap", "FN", "MCL"),
  Alkalinity_Valid_Modules = c(2, 42, 4, 3, 7, 8, 4, 1),
  Aridity_Valid_Modules = c(2, 126, 3, 5, 11, 35, 6, 1),
  Cold_Valid_Modules = c(1, 128, 1, 2, 5, 2, 4, 1)
)

# Print the table as required, without row names
print(stats_summary, row.names = FALSE)

# ====================== 2. Plotting data (English column names, values synchronized) ======================
data <- data.frame(
  Algorithm = c("MCODE", "infomap", "prop", "eigen", "louvain", "walktrap", "FN", "MCL"),
  Alkalinity = c(2, 42, 4, 3, 7, 8, 4, 1),
  Aridity = c(2, 126, 3, 5, 11, 35, 6, 1),
  Cold = c(1, 128, 1, 2, 5, 2, 4, 1)
)

# Reshape the data into long format (English column names)
data_long <- pivot_longer(
  data, 
  cols = -Algorithm, 
  names_to = "Stress_Type", 
  values_to = "Valid_Module_Count"
)

# Define the research color palette (Nature-journal standard, adapted to the 3 stress types)
research_palette <- c(
  "Alkalinity" = "#0072B2",    # Alkalinity - deep sea blue
  "Aridity" = "#D55E00",       # Aridity - orange-red
  "Cold" = "#009E73"           # Cold - emerald green
)

# Draw the bar chart (all English + research color palette)
p <- ggplot(data_long, aes(x = Algorithm, y = Valid_Module_Count, fill = Stress_Type)) +
  # Bar-chart main body (research-grade width/spacing)
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  # Value labels (optimized position and style)
  geom_text(
    aes(label = Valid_Module_Count),
    position = position_dodge(width = 0.8),
    vjust = -0.3, size = 3.5,
    color = "gray30"  # Dark-gray label text for better harmony
  ) +
  # English title and labels (research-paper standard)
  labs(
    title = "Comparison of Valid Module Counts by Clustering Algorithms",
    x = "Clustering Algorithm",
    y = "Number of Valid Modules",
    fill = "Stress Factor"
  ) +
  # Research-theme optimization
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = "gray30"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "gray30"),
    axis.text.y = element_text(size = 12, color = "gray30"),
    axis.title = element_text(size = 13, face = "bold", color = "gray30"),
    legend.title = element_text(size = 12, face = "bold", color = "gray30"),
    legend.text = element_text(size = 11, color = "gray30"),
    panel.grid = element_line(color = "gray90"),  # Light-gray grid, does not dominate the subject
    panel.border = element_rect(color = "gray50") # Medium-gray border, softer
  ) +
  # Apply the research color palette
  scale_fill_manual(values = research_palette) +
  # Adjust the y-axis range to avoid labels exceeding the plot
  ylim(0, max(data_long$Valid_Module_Count) * 1.1)  

# Save as a PDF vector graphic
ggsave(
  "Valid_Module_Count_by_Algorithm.pdf", 
  plot = p, 
  width = 10, 
  height = 7,
  device = cairo_pdf,
  bg = "white"
)

# Optional: high-resolution PNG preview
# ggsave("Valid_Module_Count_by_Algorithm.png", plot = p, width = 10, height = 7, dpi = 300)
