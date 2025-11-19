#!/bin/bash
#SBATCH --job-name=MAKER_merge_N13
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/MAKER_merge_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/MAKER_merge_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# 课程容器目录
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
CONTAINER="$COURSEDIR/containers/MAKER_3.01.03.sif"

# 你的工作目录（有 genome_N13.maker.output）
WORKDIR="/data/users/nshao/organization_annotation_course/gene_annotation"

# 你的用户根目录（挂载整个 /data/users/nshao，保证能看到 assembly_annotation_course 等）
USERROOT="/data/users/nshao"

# RepeatMasker 目录（和跑 MAKER 时一致）
REPEATMASKER_DIR="$COURSEDIR/softwares/RepeatMasker"
export PATH="$PATH:$REPEATMASKER_DIR"

cd "$WORKDIR"

echo "[INFO] Starting MAKER merge in $WORKDIR at $(date) on $(hostname)"

module load OpenMPI/4.1.1-GCC-10.3.0
module load AUGUSTUS/3.4.0-foss-2021a

# 1) 合并所有 GFF3 为一个 genome_N13.all.gff
apptainer exec \
  --bind $SCRATCH:/TMP \
  --bind "$USERROOT" \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$CONTAINER" \
  gff3_merge -d genome_N13.maker.output/genome_N13_master_datastore_index.log

# 2) 合并所有蛋白和转录本 fasta
apptainer exec \
  --bind $SCRATCH:/TMP \
  --bind "$USERROOT" \
  --bind "$COURSEDIR" \
  --bind "$AUGUSTUS_CONFIG_PATH" \
  --bind "$REPEATMASKER_DIR" \
  "$CONTAINER" \
  fasta_merge -d genome_N13.maker.output/genome_N13_master_datastore_index.log

echo "[INFO] MAKER merge finished at $(date)"