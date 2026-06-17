library(igraph)
library(ggplot2)
library(RColorBrewer)

# 二、构建数据框
df <- data.frame(
  Pathway = c(
    rep("MAPK signaling pathway - plant", 12),
    rep("Plant hormone signal transduction", 7),
    rep("MAPK signaling pathway - plant", 14),
    rep("Plant hormone signal transduction", 9),
    rep("MAPK signaling pathway - plant", 17),
    rep("Plant hormone signal transduction", 14)
  ),
  Gene = c(
    "CATC", "BWMK1", "WRKY53", "PR1b", "MAPK6", "RBOHH", "ETR3", "Cht8", "EIN2", "ACS2", "rbohA", "PR1b",
    "PR1b", "SLR1", "MAPK6", "GH3-8", "NH1", "ORR3", "ETR3",
    "CATB", "CATC", "BWMK1", "WRKY53", "MAPK6", "rbohC", "ETR3", "RBOHH", "rbohG", "EIN2", "ETR2", "ACS2", "CaM61", "RBOHB",
    "OSRR9", "OSRR2", "MAPK6", "GH3-8", "ORR3", "ETR3", "EIN2", "ETR2", "D61",
    "CATB", "PR1b", "BWMK1", "WRKY53", "MAPK6", "rbohC", "ETR3", "RBOHH", "PR1A", "EIN2", "CaM61", "rbohG", "RBOHB", "rbohA", "Cht8", "ACS2", "ETR2",
    "ORR2", "D61", "PR1b", "OSRR1", "MAPK6", "GH3-8", "JAR1", "NH1", "ETR3", "PR1A", "EIN2", "GID1", "ORR1", "ETR2"
  ),
  stringsAsFactors = FALSE
)

# 三、构建节点与边
pathway_nodes <- unique(df$Pathway)
gene_nodes <- unique(df$Gene)
all_nodes <- c(pathway_nodes, gene_nodes)
node_type <- c(rep("Pathway", length(pathway_nodes)), rep("Gene", length(gene_nodes)))
node_df <- data.frame(id = 1:length(all_nodes), name = all_nodes, type = node_type, stringsAsFactors = FALSE)

pathway_gene_edges <- unique(df[, c("Pathway", "Gene")])
colnames(pathway_gene_edges) <- c("from", "to")
pathway_gene_edges$from_id <- match(pathway_gene_edges$from, node_df$name)
pathway_gene_edges$to_id <- match(pathway_gene_edges$to, node_df$name)
edges_list <- pathway_gene_edges[, c("from_id", "to_id")]

# 四、创建igraph对象
g <- graph_from_edgelist(as.matrix(edges_list), directed = FALSE)
V(g)$name <- node_df$name
V(g)$type <- node_df$type

# 五、可视化参数（核心：纯实心、无黑框，更换清新冷色调配色）
# 定义核心基因列表
core_mapk_genes <- c("BWMK1", "MAPK6", "WRKY53")
core_aba_genes <- c("EIN2", "ETR2", "ETR3", "JAR1", "MAPK6", "NH1", "SLR1", "GID1", "D61", "PR1b", "PR1A")

# 1. 纯实心颜色（无任何边框）
# - Pathway: 深蓝色 #1F77B4
# - 普通Gene: 薄荷绿 #2CA02C
# - MAPK核心Gene: 橙色 #FF7F0E
# - ABA核心Gene: 紫色 #9467BD
V(g)$color <- ifelse(V(g)$type == "Pathway", 
                     "#1F77B4",  # Pathway颜色
                     ifelse(V(g)$name %in% core_mapk_genes, 
                            "#FF7F0E",  # MAPK核心基因颜色
                            ifelse(V(g)$name %in% core_aba_genes, 
                                   "#9467BD",  # ABA核心基因颜色
                                   "#2CA02C")))  # 普通基因颜色

# 2. 纯实心形状（无空心/黑框）
V(g)$shape <- ifelse(V(g)$type == "Pathway", "square", "circle")

# 3. 彻底隐藏边框（边框透明）
V(g)$frame.color <- "transparent"

# 4. 节点大小（核心基因放大）
V(g)$size <- ifelse(V(g)$type == "Pathway", 
                    15,  # Pathway大小
                    ifelse(V(g)$name %in% c(core_mapk_genes, core_aba_genes), 
                           12,  # 核心基因大小
                           8))  # 普通基因大小

# 5. 边样式
E(g)$lty <- 2
E(g)$width <- 1.2
E(g)$color <- "gray50"

# 六、布局
layout <- layout_with_fr(g, niter = 1000, repulserad = vcount(g)^2.2)

# 七、绘图+保存（修改为PDF格式）
# 关键修改：png() 改为 pdf()，并调整参数适配PDF
pdf("Final_Pathway_Gene_Network_Highlighted.pdf", width = 14, height = 12)  # PDF的width/height单位是英寸
par(mar = c(1,1,2,1))
plot(
  g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  edge.arrow.mode = 0,
  main = "Pathway-Target Gene Network (Core Genes Highlighted)",
  main.cex = 1.3
)

# 图例：纯实心，无黑框（包含核心基因标注）
legend(
  "topright",
  legend = c("Pathway", "Normal Gene", "MAPK Core Gene", "ABA Core Gene"),
  fill = c("#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD"),
  pt.cex = 1.5,
  cex = 1.0,
  title = "Node Type",
  bty = "n"
)
dev.off()

# 控制台显示（保持不变）
plot(
  g,
  layout = layout,
  vertex.label = V(g)$name,
  vertex.label.cex = 0.8,
  vertex.label.color = "black",
  edge.arrow.mode = 0,
  main = "Pathway-Target Gene Network"
)

# 控制台图例（保持不变）
legend(
  "topright",
  legend = c("Pathway", "Normal Gene", "MAPK Core Gene", "ABA Core Gene"),
  fill = c("#1F77B4", "#2CA02C", "#FF7F0E", "#9467BD"),
  pt.cex = 1.5,
  cex = 1.0,
  title = "Node Type",
  bty = "n"
)

cat("标注核心基因的网络图已保存为：", file.path(getwd(), "Final_Pathway_Gene_Network_Highlighted.pdf"), "\n")