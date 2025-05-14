#!/bin/bash
set -euo pipefail  # 严格错误处理

# ----------------------
# 全局变量配置
# ----------------------
# 输入输出路径（按需修改）
INPUT_CLUSTERS_DIR="clusters"         # 原始聚类FASTA目录
CLEANED_DIR="cleaned_clusters"        # 清洗后序列目录
IDMAP_DIR="idmaps"                    # ID映射表目录
MSA_EFA_DIR="msa_efa"                 # MUSCLE efa结果目录
MSA_BEST_DIR="msa_best"               # MaxCC最优比对目录
MSA_A3M_DIR="msa_a3m"                 # A3M格式目录
HMM_DIR="hmm"                         # HMM模型目录
HHR_DIR="hhr"                         # HHR比对结果目录
MATRIX_OUT="HMM_distance_matrix.csv"  # 最终距离矩阵文件

# MUSCLE参数配置
MUSCLE_ENSEMBLE_PERMS=("none" "abc" "acb" "bca")  # 排列组合
MUSCLE_PERTURB_VALUES=(0 1 2 3)                   # 扰动参数

# ----------------------
# 函数定义
# ----------------------

# 步骤1: 清洗序列并重命名ID
clean_sequences() {
    echo "=== 步骤1: 清洗序列并生成ID映射表 ==="
    mkdir -p "$CLEANED_DIR" "$IDMAP_DIR"
    
    for fasta in "$INPUT_CLUSTERS_DIR"/*.fasta; do
        cluster=$(basename "$fasta" .fasta)
        cleaned="$CLEANED_DIR/${cluster}.fasta"
        idmap="$IDMAP_DIR/${cluster}.idmap.tsv"
        
        > "$cleaned"
        echo -e "Original_ID\tNew_ID" > "$idmap"
        
        i=1
        while read line; do
            if [[ $line == ">"* ]]; then
                orig_id=$(sed 's/^>//' <<< "$line")
                new_id="${cluster}_seq${i}"
                echo ">$new_id" >> "$cleaned"
                echo -e "$orig_id\t$new_id" >> "$idmap"
                ((i++))
            else
                echo "$line" >> "$cleaned"
            fi
        done < "$fasta"
        echo "[✔] $cluster 清洗完成"
    done
}

# 步骤2: 执行MUSCLE多策略比对并提取最优
run_muscle_ensemble() {
    echo "=== 步骤2: 多策略MUSCLE比对与优化 ==="
    mkdir -p "$MSA_EFA_DIR" "$MSA_BEST_DIR"
    
    for fasta in "$CLEANED_DIR"/*.fasta; do
        cluster=$(basename "$fasta" .fasta)
        cluster_efa_dir="$MSA_EFA_DIR/$cluster"
        mkdir -p "$cluster_efa_dir"
        
        # 生成16种排列组合比对
        echo "▶ 处理 $cluster：执行16次排列比对..."
        for perm in "${MUSCLE_ENSEMBLE_PERMS[@]}"; do
            for perturb in "${MUSCLE_PERTURB_VALUES[@]}"; do
                output_afa="${cluster_efa_dir}/${perm}_${perturb}.afa"
                echo "  → $perm $perturb"
                muscle -super5 "$fasta" -perm "$perm" -perturb "$perturb" -output "$output_afa"
                # 添加标识头
                sed -i "1s/^/<${perm}_${perturb}\n/" "$output_afa"
            done
        done
        
        # 合并并提取maxcc结果
        merged_efa="${cluster_efa_dir}/merged.efa"
        cat "$cluster_efa_dir"/*.afa > "$merged_efa"
        muscle -maxcc "$merged_efa" -output "$MSA_BEST_DIR/${cluster}.afa"
        echo "[✔] $cluster 最优比对已保存"
    done
}

# 步骤3: 转换格式并构建HMM模型
build_hmm_profiles() {
    echo "=== 步骤3: 构建HMM模型 ==="
    mkdir -p "$MSA_A3M_DIR" "$HMM_DIR"
    
    # 转换AFA为A3M
    for afa in "$MSA_BEST_DIR"/*.afa; do
        base=$(basename "$afa" .afa)
        reformat.pl fas a3m "$afa" "$MSA_A3M_DIR/${base}.a3m"
        # 构建HMM
        hhmake -i "$MSA_A3M_DIR/${base}.a3m" -o "$HMM_DIR/${base}.hmm"
    done
    echo "[✔] HMM模型构建完成"
}

# 步骤4: 全对全HMM比对并生成矩阵
run_hmm_align() {
    echo "=== 步骤4: HMM全对全比对 ==="
    mkdir -p "$HHR_DIR"
    
    hmm_files=("$HMM_DIR"/*.hmm)
    total=${#hmm_files[@]}
    count=0
    
    # 遍历所有HMM对
    for ((i=0; i<total; i++)); do
        for ((j=i; j<total; j++)); do
            query="${hmm_files[$i]}"
            target="${hmm_files[$j]}"
            query_base=$(basename "$query" .hmm)
            target_base=$(basename "$target" .hmm)
            
            # 跳过重复比对（可选）
            # if [[ "$query_base" == "$target_base" ]]; then continue; fi
            
            echo "比对进度: $((++count))/$((total*(total+1)/2))"
            hhalign -i "$query" -t "$target" -o "$HHR_DIR/${query_base}_vs_${target_base}.hhr" -v 0
        done
    done
    
    # 提取得分到矩阵
    echo "Query,Target,Score" > "$MATRIX_OUT"
    for hhr in "$HHR_DIR"/*.hhr; do
        filename=$(basename "$hhr" .hhr)
        query=$(cut -d'_' -f1-2 <<< "$filename")
        target=$(cut -d'_' -f3- <<< "$filename")
        score=$(sed -n '10p' "$hhr" | cut -c58-64 | tr -d ' ')
        echo "$query,$target,$score" >> "$MATRIX_OUT"
    done
    echo "[✔] 距离矩阵已生成: $MATRIX_OUT"
}

# ----------------------
# 主流程执行
# ----------------------
main() {
    clean_sequences
    run_muscle_ensemble
    build_hmm_profiles
    run_hmm_align
    echo "✅ 全部流程完成！"
}

# 执行主函数并记录时间
time main
