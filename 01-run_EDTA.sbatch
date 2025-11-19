#!/bin/bash
#SBATCH --job-name=TE_annotation_EDTA
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=2-00:00
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/%x_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/%x_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

WORKDIR=/data/users/nshao/organization_annotation_course
GENOME=$WORKDIR/input/genome_N13.fasta
OUTDIR=$WORKDIR/results/EDTA_annotation
CONTAINER=/data/courses/assembly-annotation-course/CDS_annotation/containers/EDTA2.2.sif

# Source CDS (read-only)
CDS_SRC=/data/courses/assembly-annotation-course/CDS_annotation/data/TAIR10_cds_20110103_representative_gene_model_updated
# Copy to input and rename locally
CDS_LOCAL=$WORKDIR/input/TAIR10_cds.fa

mkdir -p "$OUTDIR"

# Copy CDS to local directory if not already present
if [ ! -s "$CDS_LOCAL" ]; then
  cp "$CDS_SRC" "$CDS_LOCAL"
fi

echo "[INFO] Starting EDTA at $(date)"
echo "[INFO] Genome: $GENOME"
echo "[INFO] Output: $OUTDIR"
echo "[INFO] CDS:    $CDS_LOCAL (size: $(stat -c%s "$CDS_LOCAL") bytes)"

apptainer exec --bind "$WORKDIR" "$CONTAINER" \
  EDTA.pl \
    --genome "$GENOME" \
    --species others \
    --step all \
    --sensitive 1 \
    --cds "$CDS_LOCAL" \
    --anno 1 \
    --threads $SLURM_CPUS_PER_TASK

echo "[INFO] EDTA finished at $(date)"