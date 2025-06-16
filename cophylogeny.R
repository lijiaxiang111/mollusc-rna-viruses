# ─────────────────────────────────────────────────────────────────────────────
# 1. 加载必要的 R 包
# ─────────────────────────────────────────────────────────────────────────────
library(ape)    # 读取 Newick 树
library(paco)   # Rtapas 依赖的 PACO 函数
library(Rtapas) # max_incong() 与 tangle_gram()

# ─────────────────────────────────────────────────────────────────────────────
# 2. 设置工作目录（请根据自己的实际路径修改）
# ─────────────────────────────────────────────────────────────────────────────
setwd("cophylogeny")

# ─────────────────────────────────────────────────────────────────────────────
# 3. 读取关联矩阵 matrix.csv
#
#    要求：
#      – matrix.csv 第一列是 Host ID（会被当作行名读入），
#      – 第一行（除第一列外）是 Virus ID（自动当作列名）。
#    确保 CSV 中没有多余的空列或空白行，否则可能自动多读出空列/空行。
# ─────────────────────────────────────────────────────────────────────────────
assoc_df <- read.csv("matrix.csv",
                     row.names    = 1,             # 把第一列当作行名
                     check.names  = FALSE,         # 列名不要自动改（保留下划线、长名称等）
                     stringsAsFactors = FALSE)

# 把 data.frame 转成纯数字矩阵
assoc_mat <- as.matrix(assoc_df)

cat("关联矩阵维度 (hosts × viruses)：", dim(assoc_mat), "\n")
cat("总关联数 (1 的数量)：", sum(assoc_mat), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 4. 读取 Host 树（这里我们假定 RdRP 树是“host”，Capsid 树是“virus”）
# ─────────────────────────────────────────────────────────────────────────────
host_tree  <- read.tree("rdrp_iqtree_output.treefile")
virus_tree <- read.tree("cap_iqtree_output.treefile")

cat("Host 树 (RdRP) tip.label（前 10 个）：\n", head(host_tree$tip.label, 10), "\n\n")
cat("Virus 树 (Capsid) tip.label（前 10 个）：\n", head(virus_tree$tip.label, 10), "\n\n")

# ─────────────────────────────────────────────────────────────────────────────
# 5. 找出不匹配的 tip.label —— 这一步非常关键！
#    用 setdiff() 检查到底是哪几个名称没有对上
# ─────────────────────────────────────────────────────────────────────────────
missing_hosts  <- setdiff(host_tree$tip.label, rownames(assoc_mat))
missing_viruses <- setdiff(virus_tree$tip.label, colnames(assoc_mat))

cat("以下 Host tip 出现在 host.tree 中，却没出现在 matrix.csv 的行名里：\n")
print(missing_hosts)
cat("\n以下 Virus tip 出现在 virus.tree 中，却没出现在 matrix.csv 的列名里：\n")
print(missing_viruses)
cat("\n请根据上述缺失名单，检查并修改 matrix.csv 或 Newick 树中的标签，确保名称完全一致（包括下划线、大小写、任何前后缀都要对上）。\n\n")

# 如果发现 missing_hosts 和 missing_viruses 都为空，那么就说明行/列名已经完美对齐，可以跳过下面的修复步骤。
if (length(missing_hosts) == 0 && length(missing_viruses) == 0) {
     cat("✔️ 所有 tree tip 与 matrix.csv 名称一一对应，可以继续后面的排序与计算。\n\n")
} else {
     stop("❌ 存在不匹配的名称，请先修改 matrix.csv 或者 Newick 树文件，再重新运行本脚本。")
}

# ─────────────────────────────────────────────────────────────────────────────
# 6. 重新对齐矩阵的行列顺序：行对应 host_tree$tip.label，列对应 virus_tree$tip.label
# ─────────────────────────────────────────────────────────────────────────────
assoc_mat_ordered <- assoc_mat[host_tree$tip.label, virus_tree$tip.label]

cat("对齐后关联矩阵维度 (hosts × viruses)：", dim(assoc_mat_ordered), "\n\n")

# 再次检查行列名顺序是否正确（可选）
# stopifnot(all(rownames(assoc_mat_ordered) == host_tree$tip.label))
# stopifnot(all(colnames(assoc_mat_ordered) == virus_tree$tip.label))

# ─────────────────────────────────────────────────────────────────────────────
# 7. 运行最大不一致度算法 (max_incong)
#    – n：子样本大小 (建议 8~10)
#    – N：随机重采样次数 (建议 ≥1e4)
# ─────────────────────────────────────────────────────────────────────────────
N <- 1e4
n <- 9

LFi <- max_incong(
     HS         = assoc_mat_ordered,
     treeH      = host_tree,
     treeS      = virus_tree,
     n          = n,
     N          = N,
     method     = "paco",      # or "paraF"
     symmetric  = TRUE,        # symmetric PACo or not
     ei.correct = "sqrt.D",    # correct negative eigenvalues
     percentile = 0.99,        # top 1% for incongruence
     diff.fq    = TRUE,        # corrected freq. for incongruence
     strat      = "sequential",
     cl         = parallelly::availableCores()  # or simply cl = 1
)

# 看一下 LFi 里面的内容（list），确认算法跑通
str(LFi)

# ─────────────────────────────────────────────────────────────────────────────
# 8. 绘制 Tanglegram（缠绕图）：
#    – colscale="diverging" + colgrad = c("darkred","gray90","darkblue") 
#      表示红色代表高度不一致（host-viral 拓扑差异大），蓝色代表高度一致，灰色居中
#    – node.tag = TRUE 会在分叉处显示对应的支持度之类的标签
# ─────────────────────────────────────────────────────────────────────────────
tangle_gram(
     treeH    = host_tree,
     treeS    = virus_tree,
     HS       = assoc_mat_ordered,
     fqtab    = LFi,
     colscale = "diverging",
     colgrad  = c("darkred","gray90","darkblue"),
     node.tag = FALSE,
     pts      = FALSE,
     yspc     = 1.0
)

# ─────────────────────────────────────────────────────────────────────────────
# 9. （可选）将 tanglegram 保存为 png 文件
# png("cophylogeny_tanglegram.png", width=2000, height=2000, res=300)
# tangle_gram(host_tree, virus_tree, assoc_mat_ordered, LFi,
#             colscale="diverging", colgrad=c("darkred","gray90","darkblue"), node.tag=TRUE)
# dev.off()
# ─────────────────────────────────────────────────────────────────────────────
