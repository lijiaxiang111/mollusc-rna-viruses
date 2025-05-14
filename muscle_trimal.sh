# 输入文件
input_fasta="used.faa"

# 创建输出目录
output_dir=./RDRP_iqtree/
mkdir -p $output_dir

# 定义排列组合和扰动数
perms=("none" "abc" "acb" "bca")
perturbs=(0 1 2 3)

# 1. 运行 muscle -super5 对齐16次
echo "Running MUSCLE alignments..."
for perm in "${perms[@]}"; do
  for perturb in "${perturbs[@]}"; do
    output_afa="${output_dir}/${perm}_${perturb}.afa"
    echo "Running: muscle -super5 $input_fasta -perm $perm -perturb $perturb -output $output_afa"
    muscle -super5 "$input_fasta" -perm "$perm" -perturb "$perturb" -output "$output_afa"
  done
done

# 等待所有背景任务完成
wait

echo "All MUSCLE alignments finished."

# 2. 在每个对齐文件开头添加 <文件名> 行
echo "Modifying AFA files to add filenames as the first line..."
for perm in "${perms[@]}"; do
  for perturb in "${perturbs[@]}"; do
    output_afa="${output_dir}/${perm}_${perturb}.afa"
    sed -i "1s/^/<${perm}_${perturb}\n/" "$output_afa"
  done
done

# 3. 合并所有对齐文件
merged_file="${output_dir}/merged.efa"
echo "Merging all AFA files into $merged_file"
cat $output_dir/*.afa > $merged_file

# 4. 运行 muscle -maxcc 去除冲突
maxcc_output="${output_dir}/maxcc.afa"
echo "Running muscle -maxcc on merged file"
muscle -maxcc "$merged_file" -output "$maxcc_output"

# 5. 使用 trimal 自动修剪对齐
trimmed_output="${output_dir}/maxcc_trimal.afa"
echo "Running trimal for automated trimming"
trimal -in "$maxcc_output" -out "$trimmed_output" -gt 0.2

# 6. 运行 IQ-TREE 构建系统发育树
iqtree_output_dir="./iqtree_output"
mkdir -p $iqtree_output_dir
echo "Running IQ-TREE for phylogenetic tree construction"
iqtree -s "$trimmed_output" -bb 3000 -nt AUTO -nm 3000 -pre "$iqtree_output_dir/iqtree_output"

echo "Script finished! IQ-TREE running in background."
