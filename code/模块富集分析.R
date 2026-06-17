# -------------------------- 6. Infomap模块GO富集分析（筛选关键TF最多的前5模块，优先大模块） --------------------------
library(org.Osativa.eg.db)
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)  # 提供str_wrap函数
library(tidyr)
library(pheatmap)  # 热图绘制包
library(RColorBrewer) # 用于构建渐变色
library(tibble)

# 6.0 提前读取关键TF数据并统一格式
key_TFs_Alkalinity <- read_tsv("TF/P_key_TFs_Alkalinity.txt", col_names = TRUE) %>% 
  rename(GID = 1) %>% 
  pull(GID) %>% unique()

key_TFs_Aridity <- read_tsv("TF/P_key_TFs_Aridity.txt", col_names = TRUE) %>% 
  rename(GID = 1) %>% 
  pull(GID) %>% unique()

key_TFs_Cold <- read_tsv("TF/P_key_TFs_Cold.txt", col_names = TRUE) %>% 
  rename(GID = 1) %>% 
  pull(GID) %>% unique()

key_tf_list <- list(
  Alkalinity = key_TFs_Alkalinity,
  Aridity = key_TFs_Aridity,
  Cold = key_TFs_Cold
)

# 6.1 核心富集分析函数
perform_enrichment_analysis <- function(gene_list, condition_name) {
  go_mapping <- bitr(
    gene_list, 
    fromType = "GID",
    toType = "GO", 
    OrgDb = org.Osativa.eg.db
  )
  
  if (nrow(go_mapping) == 0) {
    message(sprintf("→ %s：无有效GID→GO映射", condition_name))
    return(NULL)
  }
  message(sprintf("→ %s：成功转换%d个GID到GO注释（去重后%d个GO）", 
                  condition_name, nrow(go_mapping), length(unique(go_mapping$GO))))
  
  ego <- enrichGO(
    gene          = go_mapping$GO,
    OrgDb         = org.Osativa.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    keyType       = "GO"
  )
  
  if (is.null(ego) || nrow(ego) == 0) {
    message(sprintf("→ %s：无显著富集的GO条目", condition_name))
    return(NULL)
  }
  
  return(list(GO_result = ego, GO_mapping = go_mapping))
}

# 6.2 可视化函数（修改为输出PDF）
plot_enrich_results <- function(enrich_res, save_dir, prefix) {
  # 柱状图（PDF矢量图，无需设置res，调整宽高适配GO术语）
  pdf(file.path(save_dir, paste0(prefix, "_GO_barplot.pdf")), width = 10, height = 8)
  print(barplot(enrich_res$GO_result, showCategory = 10) + 
          ggtitle(paste(prefix, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  dev.off()
  
  # 点图（PDF格式）
  pdf(file.path(save_dir, paste0(prefix, "_GO_dotplot.pdf")), width = 10, height = 8)
  print(dotplot(enrich_res$GO_result, showCategory = 10) + 
          ggtitle(paste(prefix, "GO Enrichment (BP)")) +
          theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
  dev.off()
  
  # 保存富集结果表格
  write.csv(as.data.frame(enrich_res$GO_result), 
            file.path(save_dir, paste0(prefix, "_GO_enrich.csv")), 
            row.names = FALSE)
  
  # 保存GID→GO映射表
  write.csv(enrich_res$GO_mapping, 
            file.path(save_dir, paste0(prefix, "_GID2GO.csv")), 
            row.names = FALSE)
}

# 6.3 新增：全局变量存储15个核心模块的GO数据
global_merged_go <- data.frame()

# 6.4 批量处理三个胁迫条件的Infomap模块（修改：取前5核心模块，优先大模块）
stress_names <- c("Alkalinity", "Aridity", "Cold")
method <- "infomap"
max_module_size <- 300
stats_summary <- data.frame()
top5_core_module_summary <- data.frame()  # 存储每个胁迫的前5核心模块信息

for (stress in stress_names) {
  message(sprintf("\n===== 开始处理%s胁迫下的%s模块 =====", stress, method))
  
  mod_dir <- file.path("模块聚类", sprintf("%s_module", stress), method)
  if (!dir.exists(mod_dir)) {
    message(sprintf("→ 目录不存在：%s，跳过", mod_dir))
    next
  }
  
  mod_files <- list.files(mod_dir, pattern = "module_\\d+\\.txt", full.names = TRUE)
  if (length(mod_files) == 0) {
    message("→ 该目录下无有效模块文件，跳过")
    next
  }
  
  # 创建结果保存目录
  result_root <- file.path(mod_dir, "GO_enrichment_results")
  top5_core_result_root <- file.path(mod_dir, "Top5_Core_Module_GO_enrichment")
  dir.create(result_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(top5_core_result_root, recursive = TRUE, showWarnings = FALSE)
  
  # 初始化统计计数器
  total_modules <- 0
  valid_modules <- 0
  has_enrich_modules <- 0
  no_enrich_modules <- 0
  
  # 存储每个有效模块的关键TF数量
  module_tf_count <- data.frame(
    模块ID = character(),
    模块基因数 = integer(),
    关键TF数量 = integer(),
    模块文件路径 = character(),
    stringsAsFactors = FALSE
  )
  
  # 逐个模块分析
  for (file in mod_files) {
    mod_id <- sub("module_(\\d+)\\.txt", "\\1", basename(file))
    prefix <- sprintf("%s_Infomap_Module%s", stress, mod_id)
    message(sprintf("→ 正在分析%s...", prefix))
    
    module_gids <- read.table(file, header = FALSE, stringsAsFactors = FALSE)[[1]]
    module_size <- length(module_gids)
    total_modules <- total_modules + 1
    message(sprintf("→ 模块包含%d个GID", module_size))
    
    # 过滤基因数>300的模块
    if (module_size > max_module_size) {
      message(sprintf("→ 模块%s基因数（%d）>%d，跳过", mod_id, module_size, max_module_size))
      next
    }
    valid_modules <- valid_modules + 1
    
    # 统计当前模块中的关键TF数量
    current_key_tfs <- key_tf_list[[stress]]
    tf_count <- sum(module_gids %in% current_key_tfs)
    message(sprintf("→ 模块%s包含%d个关键TF", mod_id, tf_count))
    
    # 记录模块信息
    module_tf_count <- rbind(module_tf_count, data.frame(
      模块ID = mod_id,
      模块基因数 = module_size,
      关键TF数量 = tf_count,
      模块文件路径 = file,
      stringsAsFactors = FALSE
    ))
    
    # 执行富集分析
    enrich_res <- perform_enrichment_analysis(module_gids, prefix)
    
    # 统计并保存结果
    if (!is.null(enrich_res)) {
      plot_enrich_results(enrich_res, result_root, prefix)
      message(sprintf("→ %s结果已保存到：%s", prefix, result_root))
      has_enrich_modules <- has_enrich_modules + 1
    } else {
      no_enrich_modules <- no_enrich_modules + 1
    }
  }
  
  # 筛选关键TF数量最多的前5个模块（优先大模块）
  if (nrow(module_tf_count) > 0) {
    # 排序规则：1. 关键TF数量（降序） 2. 模块基因数（降序） 3. 模块ID（升序）
    top5_core_modules <- module_tf_count %>% 
      arrange(desc(关键TF数量), desc(模块基因数), 模块ID) %>% 
      slice_head(n = min(5, nrow(module_tf_count)))  # 若模块数<5，取全部
    
    # 记录前5核心模块信息
    top5_core_module_info <- top5_core_modules %>% 
      mutate(
        胁迫条件 = stress,
        算法 = method,
        排名 = 1:nrow(.)  # 新增排名列
      ) %>% 
      select(胁迫条件, 算法, 排名, 模块ID, 模块基因数, 关键TF数量, 模块文件路径)
    
    top5_core_module_summary <- rbind(top5_core_module_summary, top5_core_module_info)
    
    # 对前5核心模块批量做富集分析并保存 + 收集GO数据到全局变量
    message(sprintf("\n→ %s胁迫的前%d核心模块（TF最多+优先大模块）：%s", 
                    stress, nrow(top5_core_modules), paste(top5_core_modules$模块ID, collapse = ",")))
    
    for (i in 1:nrow(top5_core_modules)) {
      core_mod <- top5_core_modules[i, ]
      core_mod_id <- core_mod$模块ID
      core_module_gids <- read.table(core_mod$模块文件路径, header = FALSE, stringsAsFactors = FALSE)[[1]]
      core_prefix <- sprintf("%s_Infomap_Top%d_Core_Module%s", stress, i, core_mod_id)
      message(sprintf("→ 分析Top%d核心模块%s（TF数=%d，基因数=%d）", 
                      i, core_mod_id, core_mod$关键TF数量, core_mod$模块基因数))
      
      core_enrich_res <- perform_enrichment_analysis(core_module_gids, core_prefix)
      
      if (!is.null(core_enrich_res)) {
        plot_enrich_results(core_enrich_res, top5_core_result_root, core_prefix)
        message(sprintf("→ Top%d核心模块%s富集结果已保存到：%s", i, core_mod_id, top5_core_result_root))
        
        # 收集核心模块的GO数据到全局变量（暂不处理换行，后续统一处理）
        core_go_data <- as.data.frame(core_enrich_res$GO_result) %>% 
          mutate(
            胁迫条件 = stress,
            模块编号 = core_mod_id,
            胁迫_模块 = paste(stress, paste0("mod", core_mod_id), sep = "_"),  # 热图列名
            模块关键TF数 = core_mod$关键TF数量,
            模块基因数 = core_mod$模块基因数,
            log10_qvalue = -log10(qvalue)  # 计算显著性值
          )
        global_merged_go <- rbind(global_merged_go, core_go_data)
        
      } else {
        message(sprintf("→ Top%d核心模块%s无显著富集结果", i, core_mod_id))
      }
    }
    
    # 保存模块关键TF统计表格（含排序）
    module_tf_count_sorted <- module_tf_count %>% 
      arrange(desc(关键TF数量), desc(模块基因数), 模块ID)
    write.csv(module_tf_count_sorted, 
              file.path(mod_dir, "模块关键TF数量统计（按TF数+基因数排序）.csv"), 
              row.names = FALSE, fileEncoding = "UTF-8")
    
  } else {
    message(sprintf("→ %s胁迫无有效模块，无法筛选核心模块", stress))
  }
  
  # 统计结果
  stress_stats <- data.frame(
    胁迫条件 = stress,
    算法 = method,
    总模块数 = total_modules,
    有效模块数_基因数不超过300 = valid_modules,
    有富集结果的模块数 = has_enrich_modules,
    无富集结果的模块数 = no_enrich_modules,
    stringsAsFactors = FALSE
  )
  stats_summary <- rbind(stats_summary, stress_stats)
  
  message(sprintf("===== %s胁迫的%s模块分析全部完成 =====\n", stress, method))
}

# -------------------------- 7. 绘制15模块×GO通路分组热图（红白蓝渐变色+横坐标按胁迫排序） --------------------------
message("\n===== 开始绘制15模块×核心GO通路分组热图 =====\n")

if (nrow(global_merged_go) == 0) {
  message("→ 无有效GO富集数据，无法绘制热图")
} else {
  # 7.1 筛选核心GO通路（按总显著性排序取前20）
  top_20_go <- global_merged_go %>% 
    group_by(Description) %>% 
    summarise(mean_logq = mean(log10_qvalue)) %>% 
    arrange(desc(mean_logq)) %>% 
    slice_head(n = 20) %>% 
    pull(Description)
  
  # 7.2 构建热图矩阵 + 通路名换行处理
  heatmap_data_df <- global_merged_go %>% 
    filter(Description %in% top_20_go) %>% 
    select(胁迫_模块, Description, log10_qvalue) %>% 
    distinct() %>%  # 去重：每个模块-通路只保留一个值
    mutate(Description = str_wrap(Description, width = 25)) %>%  # 通路名换行
    pivot_wider(
      names_from = 胁迫_模块,
      values_from = log10_qvalue,
      values_fill = 0  # 无富集的通路赋值为0
    )
  
  # ========== 关键修改：横坐标按胁迫排序（Alkalinity→Aridity→Cold） ==========
  # 1. 提取各胁迫对应的模块列
  alk_modules <- grep("^Alkalinity_", colnames(heatmap_data_df), value = TRUE)
  arid_modules <- grep("^Aridity_", colnames(heatmap_data_df), value = TRUE)
  cold_modules <- grep("^Cold_", colnames(heatmap_data_df), value = TRUE)
  
  # 2. 按顺序拼接列（先碱胁迫，再干旱，最后冷胁迫）
  ordered_cols <- c("Description", alk_modules, arid_modules, cold_modules)
  heatmap_data_df_ordered <- heatmap_data_df[, ordered_cols]
  
  # 3. 转换为矩阵用于绘图
  heatmap_data <- heatmap_data_df_ordered %>% 
    column_to_rownames("Description") %>% 
    as.matrix()
  
  # 7.3 构建模块分组注释（按胁迫分组，顺序与横坐标一致）
  module_groups <- data.frame(
    Stress_Group = sapply(strsplit(colnames(heatmap_data), "_"), `[`, 1),
    row.names = colnames(heatmap_data)
  )
  
  # 7.4 定义分组配色
  group_colors <- list(
    Stress_Group = c(
      Alkalinity = "#4575b4",  # 蓝色（碱胁迫）
      Aridity = "#91cf60",     # 绿色（干旱胁迫）
      Cold = "#d73027"         # 红色（冷胁迫）
    )
  )
  
  # ========== 核心修改：构建红白蓝渐变色（越富集越红） ==========
  # 颜色顺序：蓝（最低富集）→ 白（中等富集）→ 红（最高富集），生成100个渐变颜色
  heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)
  
  # 7.5 绘制分组热图（修改为输出PDF，移除dpi参数）
  pheatmap(
    mat = heatmap_data,
    annotation_col = module_groups,        # 列注释（胁迫分组）
    annotation_colors = group_colors,      # 分组配色
    color = heatmap_colors,                # 启用红白蓝渐变色（越红富集程度越高）
    border_color = "white",                # 边框白色
    treeheight_col = 0,                    # 关闭列聚类，保持手动排序
    treeheight_row = 10,                   # 行聚类（通路相似性）
    fontsize = 13,                         # 整体大字体
    fontsize_row = 11,                     # 通路名字体
    fontsize_col = 12,                     # 模块名字体
    fontsize_main = 16,                    # 标题字体
    angle_col = 45,                        # 模块标签旋转45度
    width = 20,                            # 宽度保持紧凑
    height = 14,                           # 适配大字体高度
    main = "Modules GO Pathways",
    filename = "15Modules_GO_Heatmap_OrderedX_RedWhiteBlue.pdf",  # 修改为PDF后缀
    # 移除dpi参数（PDF为矢量图，无需设置分辨率）
  )
  
  message("→ 红白蓝渐变色+横坐标按胁迫排序的热图已保存为：15Modules_GO_Heatmap_OrderedX_RedWhiteBlue.pdf")
}

# 输出最终统计汇总表
message("\n===== 三个胁迫条件的富集分析统计汇总 =====")
print(stats_summary, row.names = FALSE)

# 输出前5核心模块汇总表
message("\n===== 三个胁迫条件的前5核心模块（TF最多+优先大模块）汇总 =====")
print(top5_core_module_summary, row.names = FALSE)

# 保存统计汇总表和前5核心模块表
write.csv(stats_summary, "Infomap模块GO富集分析统计汇总.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(top5_core_module_summary, "Infomap前5核心模块（TF最多+优先大模块）汇总.csv", row.names = FALSE, fileEncoding = "UTF-8")
# 保存15模块GO数据
write.csv(global_merged_go, "15Modules_GO_Enrichment_Data.csv", row.names = FALSE, fileEncoding = "UTF-8")

message("\n→ 统计汇总表已保存为：Infomap模块GO富集分析统计汇总.csv")
message("→ 前5核心模块汇总表已保存为：Infomap前5核心模块（TF最多+优先大模块）汇总.csv")
message("→ 15模块GO富集数据已保存为：15Modules_GO_Enrichment_Data.csv")