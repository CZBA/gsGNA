# 一、加载所需包（兼容所有低版本igraph，无报错）
library(igraph)
library(ggplot2)
library(RColorBrewer)

# 二、数据准备（区分已报道TF、未报道TF，解决重复节点问题）
# 1. 合并两个通路为一个核心节点（学术简洁命名）
core_pathway_name <- "MAPK + Plant hormone
signal pathway"

# 2. 读取已报道TF（你的原始TF文件）
reported_tf_list <- read.delim("核心通路筛选TF_results/Key_TFs_Annotation_Only.txt", 
                               header = FALSE, 
                               stringsAsFactors = FALSE)[,1]
reported_tf_list <- unique(reported_tf_list)  # 先对已报道TF去重

# 3. 读取未报道TF（新增：未报道的TF.txt）
unreported_tf_list <- read.delim("核心通路筛选TF_results/未报道的TF.txt", 
                                 header = FALSE, 
                                 stringsAsFactors = FALSE)[,1]
unreported_tf_list <- unique(unreported_tf_list)  # 先对未报道TF去重

# ========== 核心修正：排查并去除重复TF（已报道与未报道列表的交集） ==========
# 方法1：标记重复项，保留未报道属性（优先显示未报道，更符合科研需求）
# 提取重复的TF名称
duplicated_tf <- intersect(reported_tf_list, unreported_tf_list)
if (length(duplicated_tf) > 0) {
  cat("发现重复TF：", paste(duplicated_tf, collapse = ", "), "\n")
  cat("已自动去除已报道列表中的重复项，保留未报道属性\n")
  # 从已报道TF列表中删除重复项（避免节点名重复）
  reported_tf_list <- setdiff(reported_tf_list, duplicated_tf)
}

# 4. 合并所有TF（去重后的已报道+未报道，无重复）
all_tf_list <- c(reported_tf_list, unreported_tf_list)

# 三、构建节点表（区分核心通路、已报道TF、未报道TF，均为圆形，无重复节点）
# 先构建基础节点数据（核心通路 + 去重后的所有TF）
node_names <- c(core_pathway_name, all_tf_list)
# 定义节点类型（3类：Pathway、Reported_TF、Unreported_TF）
node_types <- vector(length = length(node_names))
node_types[node_names == core_pathway_name] <- "Pathway"
node_types[node_names %in% reported_tf_list] <- "Reported_TF"
node_types[node_names %in% unreported_tf_list] <- "Unreported_TF"
# 构建节点表（无重复名称）
nodes <- data.frame(
  name = node_names,
  type = factor(node_types, levels = c("Pathway", "Reported_TF", "Unreported_TF")),
  stringsAsFactors = FALSE
)

# 四、构建边表：所有TF（已报道+未报道）都连接核心通路（无重复边）
# 已报道TF与核心的边
edges_reported <- data.frame(
  from = core_pathway_name,
  to = reported_tf_list,
  stringsAsFactors = FALSE
)
# 未报道TF与核心的边
edges_unreported <- data.frame(
  from = core_pathway_name,
  to = unreported_tf_list,
  stringsAsFactors = FALSE
)
# 合并所有边（无重复，因为TF列表已去重）
edges <- rbind(edges_reported, edges_unreported)

# 五、创建igraph对象（无重复顶点名称，解决报错）
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# 六、非环形布局（学术首选力导向布局，节点均匀无重叠）
layout <- layout_with_fr(g)
# 可选：分层布局（适合展示层级关系）
# layout <- layout_as_tree(g, root = which(V(g)$type == "Pathway"))
layout <- layout * 1.2

# 七、可视化参数：科研配色（区分3类节点，顶刊风格）
# ========== 科研配色（适配论文/汇报，3类节点清晰区分） ==========
pathway_color <- "#0072B2"        # 核心通路：深海蓝（Nature顶刊色，不变）
reported_tf_color <- "#009E73"    # 已报道TF：橄榄绿（原TF颜色，保持一致性）
unreported_tf_color <- "#E69F00"  # 未报道TF：暖橙色（科研常用对比色，醒目不花哨）
edge_color <- "#BABABA"           # 边：中性灰（不抢焦点）

# 节点颜色赋值（按类型区分）
V(g)$color <- ifelse(V(g)$type == "Pathway", pathway_color,
                     ifelse(V(g)$type == "Reported_TF", reported_tf_color,
                            unreported_tf_color))

# 节点大小（学术规范：核心突出，TF清晰，未报道TF与已报道TF大小一致）
V(g)$size <- ifelse(V(g)$type == "Pathway", 22, 7)  # 核心22，所有TF均为7
V(g)$shape <- "circle"  # 全圆形节点，兼容低版本igraph
V(g)$frame.color <- NA  # 无边框，学术整洁样式

# 边样式（科研规范：细边+低饱和度）
E(g)$color <- edge_color
E(g)$width <- 0.7
E(g)$curved <- FALSE

# 标签样式（学术规范：字号适中，避免遮挡）
V(g)$label.cex <- ifelse(V(g)$type == "Pathway", 0.9, 0.5)
V(g)$label.color <- "black"
V(g)$label.dist <- ifelse(V(g)$type == "Pathway", 0, 0.2)

# 八、绘制并保存（PDF格式，矢量图，适合学术论文/汇报）
# 关键修改：png() 改为 pdf()，适配PDF的英寸单位（原3000x3000像素/300dpi = 10x10英寸）
pdf("Scientific_Color_TF_With_Unreported_Plot.pdf", 
    width = 10, height = 10)  # PDF单位为英寸，10x10英寸对应原300dpi的3000x3000像素
par(mar = c(0, 0, 0, 0))

# 绘制图形
plot(
  g,
  layout = layout,
  vertex.label = V(g)$name,
  main = "",
  sub = ""
)

# 添加图例（科研简洁风格，区分3类节点）
legend(
  "topright",
  legend = c("Core Pathway", "Reported TF", "Unreported TF"),
  fill = c(pathway_color, reported_tf_color, unreported_tf_color),
  pt.cex = c(2.5, 1, 1),
  cex = 0.9,
  bty = "n"
)

dev.off()

# 控制台提示
cat("区分已报道/未报道TF的可视化结果已保存：Scientific_Color_TF_With_Unreported_Plot.pdf\n")
cat("配色说明（科研规范）：\n1. 核心通路：#0072B2（深海蓝，Nature顶刊色）\n2. 已报道TF：#009E73（橄榄绿，与核心色协调）\n3. 未报道TF：#E69F00（暖橙色，醒目区分，不花哨）\n")
cat("核心特性：\n1. 错误解决：自动排查并去除重复TF，无重复顶点名称\n2. 数据区分：精准标识未报道TF，保留所有关联逻辑\n3. 样式规范：全圆形、矢量PDF、无边框，符合学术论文要求\n")
cat("PDF优势：矢量图格式，放大无模糊，可直接用于期刊投稿和高质量汇报\n")