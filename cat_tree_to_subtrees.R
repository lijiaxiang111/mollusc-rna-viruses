# -----------------------------
# ✅ 参数设置
tree_file <- "tree_no_RT.newick"        # 输入树文件
output_dir <- "subtree_iqtree"          # 输出目录
cut_height <- 2.5                        # 切割高度阈值
min_size <- 1                            # 子树最小 cluster 数
# -----------------------------

# ✅ 加载必要包（系统默认带 ape，无需额外安装）
library(ape)

# ✅ 创建输出目录
if (!dir.exists(output_dir)) dir.create(output_dir)

# ✅ 读取 Newick 树
cat("📥 读取树文件:", tree_file, "/n")
tree <- read.tree(tree_file)

# ✅ 聚类切割
hc <- as.hclust(tree)
groups <- cutree(hc, h = cut_height)

# ✅ 打印子树数量
cat("✂️ 使用高度", cut_height, "切割，共识别到", length(unique(groups)), "个子树/n")

# ✅ 按组输出
group_df <- data.frame(tip = names(groups), cluster = groups)
for (k in unique(groups)) {
     members <- group_df$tip[group_df$cluster == k]
     if (length(members) >= min_size) {
          outfile <- file.path(output_dir, sprintf("group_%04d.list", k))
          writeLines(members, outfile)
          cat(sprintf("✅ 输出: %s (%d clusters)/n", outfile, length(members)))
     }
}

cat("🎉 所有子树成员列表已生成于:", output_dir, "/n")
