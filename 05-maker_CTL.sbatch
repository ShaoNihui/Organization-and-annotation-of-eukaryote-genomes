#!/bin/bash
#SBATCH --job-name=maker_CTL
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00

#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/%x_%j.err

set -euo pipefail

WORKDIR=/data/users/nshao/organization_annotation_course/gene_annotation
COURSEDIR=/data/courses/assembly-annotation-course/CDS_annotation

cd "$WORKDIR"

echo "[INFO] Running maker -CTL in $WORKDIR on node $(hostname) at $(date)"

# 注意：apptainer 只能在计算节点用，因此必须通过 sbatch 运行
apptainer exec \
  --bind /data/users/nshao/organization_annotation_course \
  "$COURSEDIR/containers/MAKER_3.01.03.sif" \
  maker -CTL

echo "[INFO] Done maker -CTL at $(date)"