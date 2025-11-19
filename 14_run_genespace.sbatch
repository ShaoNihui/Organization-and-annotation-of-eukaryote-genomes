#!/bin/bash
#SBATCH --job-name=GENESPACE_N13
#SBATCH -p pibu_el8
#SBATCH --output=../logs/GENESPACE_N13.%j.out
#SBATCH --error=../logs/GENESPACE_N13.%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nihui.shao@students.unibe.ch

echo "=== [GENESPACE] job started on $(date) ==="

# 顶级工作目录（gene_annotation）
ROOT="/data/users/nshao/organization_annotation_course/gene_annotation"

# GENESPACE 目录
GENESPACE_WD="${ROOT}/genespace"

echo "Working directory = ${ROOT}"
echo "GENESPACE working directory = ${GENESPACE_WD}"

cd "${ROOT}" || { echo "Cannot cd to ROOT"; exit 1; }

# 启动 container 运行 Rscript
apptainer exec \
  --bind /data/courses/assembly-annotation-course/CDS_annotation \
  --bind /data/users/nshao/organization_annotation_course \
  --bind $SCRATCH:/temp \
  /data/courses/assembly-annotation-course/CDS_annotation/containers/genespace_latest.sif \
  Rscript ../scripts/13_genespace.R "${GENESPACE_WD}"

echo "=== [GENESPACE] finished on $(date) ==="