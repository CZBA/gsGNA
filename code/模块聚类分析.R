# -------------------------- 1. 安装并加载所需包 --------------------------
# 安装基础工具包（如未安装）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!require("vroom", quietly = TRUE)) {
  install.packages("vroom")
}

# 安装并加载生物信息学相关包
bioc_packages <- c("org.Osativa.eg.db", "clusterProfiler", "enrichplot", "miRspongeR")
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg, version = "3.20")  # 匹配miRspongeR要求的版本
  }
}

# 加载所有包
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(VennDiagram)
library(igraph)
library(tidyverse)
library(miRspongeR)
library(ggplot2)
library(cowplot)
library(vroom)  
library(ReactomePA)


# -------------------------- 2. 读取数据（仅保留前三列） --------------------------
# 读取3个胁迫数据集的边数据，仅保留前三列
Z_ll_Alkalinity <- vroom("Alkalinity_significant_edges.csv", delim = ",") %>% select(1:3)
Z_ll_Aridity <- vroom("Aridity_significant_edges.csv", delim = ",") %>% select(1:3)
Z_ll_Cold <- vroom("Cold_significant_edges.csv", delim = ",") %>% select(1:3)


# -------------------------- 3. 定义核心函数（改为逐个算法处理） --------------------------
# 函数1：单个算法识别模块（每次处理一个算法）
get_single_method_modules <- function(network_data, method) {
  netModule(
    spongenetwork = network_data[, 1:2],  # 仅用前两列节点信息构建网络
    method = method,
    directed = FALSE,  # miRNA海绵网络默认无向
    modulesize = 1     # 先保留所有模块，后续筛选大小
  )
}

# 函数2：批量调用单个算法函数，逐个处理所有算法
get_modules <- function(network_data, methods) {
  module_list <- list()  # 初始化空列表存储结果
  for (method in methods) {
    cat(sprintf("→ 正在用%s算法识别模块...\n", method))
    # 逐个算法运行，添加到结果列表
    module_list[[method]] <- get_single_method_modules(network_data, method)
  }
  return(module_list)
}

# 函数3：筛选基因数>10的模块
filter_large_modules <- function(modules) {
  lapply(modules, function(alg_modules) {
    alg_modules[sapply(alg_modules, length) > 10]  # 仅保留基因数>10的模块
  })
}

# 函数4：创建"胁迫条件→聚类方法"的文件夹层级，并输出模块结果
save_modules_with_folder <- function(valid_modules, stress_name) {
  # 步骤1：创建胁迫条件的一级文件夹
  stress_folder <- sprintf("%s_module", stress_name)
  dir.create(stress_folder, recursive = TRUE, showWarnings = FALSE)
  
  # 步骤2：遍历每个聚类方法，创建二级子文件夹并输出结果
  for (method in names(valid_modules)) {
    method_folder <- file.path(stress_folder, method)
    dir.create(method_folder, recursive = TRUE, showWarnings = FALSE)
    
    method_mods <- valid_modules[[method]]
    if (length(method_mods) == 0) {
      cat(sprintf("→ 胁迫%s的%s算法无有效模块（基因数>10），跳过\n", stress_name, method))
      next
    }
    
    # 步骤3：逐个模块输出
    for (mod_idx in seq_along(method_mods)) {
      mod_genes <- method_mods[[mod_idx]]
      file_path <- file.path(method_folder, sprintf("module_%d.txt", mod_idx))
      
      write.table(
        x = data.frame(gene = mod_genes),
        file = file_path,
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
      )
    }
    cat(sprintf("→ 胁迫%s的%s算法：%d个有效模块已存入%s\n", 
                stress_name, method, length(method_mods), method_folder))
  }
  cat(sprintf("=== 胁迫%s的模块结果已全部输出到：%s\n\n", stress_name, stress_folder))
}


# -------------------------- 4. 批量处理3个胁迫数据集 --------------------------
# 定义要使用的8种聚类算法
used_methods <- c( "MCODE", "infomap", "prop", "eigen", "louvain", "walktrap","FN", "MCL")

# 4.1 处理 Alkalinity 胁迫
cat("开始处理 Alkalinity 胁迫数据...\n")
alkalinity_all_modules <- get_modules(Z_ll_Alkalinity, used_methods)  # 逐个算法运行
alkalinity_large_modules <- filter_large_modules(alkalinity_all_modules)
save_modules_with_folder(alkalinity_large_modules, stress_name = "Alkalinity")

# 4.2 处理 Aridity 胁迫
cat("开始处理 Aridity 胁迫数据...\n")
aridity_all_modules <- get_modules(Z_ll_Aridity, used_methods)  # 逐个算法运行
aridity_large_modules <- filter_large_modules(aridity_all_modules)
save_modules_with_folder(aridity_large_modules, stress_name = "Aridity")

# 4.3 处理 Cold 胁迫
cat("开始处理 Cold 胁迫数据...\n")
cold_all_modules <- get_modules(Z_ll_Cold, used_methods)  # 逐个算法运行
cold_large_modules <- filter_large_modules(cold_all_modules)
save_modules_with_folder(cold_large_modules, stress_name = "Cold")


# -------------------------- 5. 输出各胁迫的有效模块统计汇总 --------------------------
get_module_count <- function(large_modules) {
  sapply(large_modules, length)
}

stats_summary <- data.frame(
  算法 = used_methods,
  Alkalinity_有效模块数 = get_module_count(alkalinity_large_modules)[used_methods],
  Aridity_有效模块数 = get_module_count(aridity_large_modules)[used_methods],
  Cold_有效模块数 = get_module_count(cold_large_modules)[used_methods]
)

cat("=== 3个胁迫条件的有效模块数汇总（基因数>10）===\n")
print(stats_summary, row.names = FALSE)