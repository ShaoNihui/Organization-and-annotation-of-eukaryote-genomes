#!/bin/bash
#SBATCH --job-name=LTR_TEsorter
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/%x_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# -------- 路径（全部在 /data/users 之下） --------
WORKDIR=/data/users/nshao/organization_annotation_course
EDTA_DIR=$WORKDIR/results/EDTA_annotation
TESDIR=$WORKDIR/results/TEsorter
FIGDIR=$WORKDIR/results/figures
mkdir -p "$TESDIR" "$FIGDIR"

# TEsorter 容器和老师给的 R 脚本（只读，共享区）
TESORTER_SIF=/data/courses/assembly-annotation-course/CDS_annotation/containers/TEsorter_1.3.0.sif
R_IDENTITY=/data/courses/assembly-annotation-course/CDS_annotation/scripts/02-full_length_LTRs_identity.R

# EDTA 产生的 LTR 文件
LTR_FA=$EDTA_DIR/genome_Sha.fasta.mod.EDTA.raw/genome_Sha.fasta.mod.LTR.raw.fa
LTR_GFF=$EDTA_DIR/genome_Sha.fasta.mod.EDTA.raw/genome_Sha.fasta.mod.LTR.intact.raw.gff3

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[INFO] Node: $(hostname)"
echo "[INFO] Start: $(date)"

# -------- 基本检查 --------
if [ ! -s "$LTR_FA" ]; then
  echo "[ERROR] Missing LTR fasta: $LTR_FA"
  exit 1
fi
if [ ! -s "$LTR_GFF" ]; then
  echo "[ERROR] Missing LTR intact gff3: $LTR_GFF"
  exit 1
fi

# -------- 用节点本地盘加速 I/O --------
TMPRUN=${TMPDIR:-/tmp}/ltr_${SLURM_JOB_ID}
mkdir -p "$TMPRUN"
echo "[INFO] TMPRUN: $TMPRUN"

# 拷贝需要的文件到本地盘
cp -v "$LTR_FA"  "$TMPRUN/LTR.raw.fa"
cp -v "$LTR_GFF" "$TMPRUN/LTR.intact.raw.gff3"

# -------- 1) 跑 TEsorter 做 clade 分类 --------
# -pre LTR 会生成 LTR.liban.rexdb-plant.* 等文件
apptainer exec --bind "$WORKDIR","/data/courses","$TMPRUN" "$TESORTER_SIF" \
  TEsorter "$TMPRUN/LTR.raw.fa" -db rexdb-plant -p $SLURM_CPUS_PER_TASK -pre "$TMPRUN/LTR"

# 同步 TEsorter 输出回你的结果目录
rsync -av "$TMPRUN"/LTR.* "$TESDIR"/

# -------- 2) 运行老师的 R 脚本画 LTR identity 直方图 --------
# 脚本会从 intact.gff3 中读取 ltr_identity 信息并画图
if ! command -v Rscript >/dev/null 2>&1; then
  echo "[WARN] Rscript not found in PATH. If this step fails, try: module load R"
fi

pushd "$TMPRUN" >/dev/null
Rscript "$R_IDENTITY" "$TMPRUN/LTR.intact.raw.gff3" || echo "[WARN] R plotting step returned non-zero; check logs."
popd >/dev/null

# 拷贝 pdf/png 图到 figures 目录
rsync -av "$TMPRUN"/*.pdf "$FIGDIR"/ 2>/dev/null || true
rsync -av "$TMPRUN"/*.png "$FIGDIR"/ 2>/dev/null || true

echo "[INFO] Done: $(date)"