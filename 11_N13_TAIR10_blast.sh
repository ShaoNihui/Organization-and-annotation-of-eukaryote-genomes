#!/bin/bash
#SBATCH --job-name=N13_TAIR10
#SBATCH -p pibu_el8
#SBATCH --output=/data/users/nshao/organization_annotation_course/gene_annotation/logs/N13_TAIR10.%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/gene_annotation/logs/N13_TAIR10.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nihui.shao@students.unibe.ch

set -e  # ← 有任何命令报错就立刻退出，不再假装成功

module load BLAST+/2.15.0-gompi-2021a

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/nshao/organization_annotation_course/gene_annotation/final"
cd "$WORKDIR"

# ✅ 正确的蛋白序列文件名（注意两个 .fasta）
QUERY="genome_N13.all.maker.proteins.fasta.renamed.fasta.filtered.fasta"

DB="${COURSEDIR}/data/TAIR10_pep_20110103_representative_gene_model"
OUT="N13_vs_TAIR10.blastp"

echo "=== Running TAIR10 BLASTP annotation for N13 ==="
echo "[0] Workdir: $(pwd)"
echo "[0] Query : $QUERY"
echo "[0] DB    : $DB"
echo

# 1. 检查输入文件是否存在
if [ ! -s "$QUERY" ]; then
    echo "ERROR: Query file not found or empty: $QUERY"
    ls -lh "$QUERY" || echo "ls cannot see this file."
    exit 1
fi

# 2. 运行 BLASTP
echo "[1] Running BLASTP..."
blastp \
  -query "$QUERY" \
  -db "$DB" \
  -num_threads 10 \
  -outfmt 6 \
  -evalue 1e-5 \
  -max_target_seqs 10 \
  -out "$OUT"

echo "[2] Sorting for best hits..."
sort -k1,1 -k12,12g "$OUT" | sort -u -k1,1 --merge > "${OUT}.besthits"

echo
echo "=== DONE! Results ==="
echo "BLAST output: $(pwd)/$OUT"
echo "Best hits:    $(pwd)/${OUT}.besthits"