#!/bin/bash
#SBATCH --job-name=iprscan_N13
#SBATCH --time=1-0
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -p pibu_el8
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/08_interproscan_N13_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/08_interproscan_N13_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/nshao/organization_annotation_course"
IPRDATA="${COURSEDIR}/data/interproscan-5.70-102.0/data"

echo "[INFO] Starting InterProScan for N13 at $(date) on ${HOSTNAME}"
echo "[INFO] WORKDIR:      ${WORKDIR}/gene_annotation/final"
echo "[INFO] Protein file: ${WORKDIR}/gene_annotation/final/genome_N13.all.maker.proteins.fasta.renamed.fasta"
echo "[INFO] InterPro data: ${IPRDATA}"

cd "${WORKDIR}/gene_annotation/final" || { echo "[ERROR] Cannot cd to final dir"; exit 1; }

apptainer exec \
  --bind "${IPRDATA}:/opt/interproscan/data" \
  --bind "${WORKDIR}" \
  --bind "${COURSEDIR}" \
  --bind "${SCRATCH}:/temp" \
  "${COURSEDIR}/containers/interproscan_latest.sif" \
  /opt/interproscan/interproscan.sh \
    -appl pfam --disable-precalc -f TSV \
    --goterms --iprlookup --seqtype p \
    -T /temp \
    -cpu "${SLURM_CPUS_PER_TASK}" \
    -i genome_N13.all.maker.proteins.fasta.renamed.fasta \
    -o output.iprscan

echo "[INFO] InterProScan finished at $(date)"