#!/bin/bash
#SBATCH --time=2-0
#SBATCH --mem=32G
#SBATCH -p pibu_el8
#SBATCH --cpus-per-task=10
#SBATCH --job-name=10_N13_uniprot
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/10_N13_uniprot_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/10_N13_uniprot_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

# 路径设置
COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"
WORKDIR="/data/users/nshao/organization_annotation_course"
FINALDIR="${WORKDIR}/gene_annotation/final"
MAKERBIN="${COURSEDIR}/softwares/Maker_v3.01.03/src/bin"

# 我们的文件名（注意这个！！）
PROT="genome_N13.all.maker.proteins.fasta.renamed.fasta.filtered.fasta"
GFF="filtered.genes.renamed.gff3"

# UniProt 数据库 fasta（老师给的）
UNIPROT_FASTA="${COURSEDIR}/data/uniprot/uniprot_viridiplantae_reviewed.fa"

# BLAST 输出前缀
BLAST_OUT="N13_vs_uniprot_viridiplantae.blastp"

module load BLAST+/2.15.0-gompi-2021a

cd "${FINALDIR}"

echo "[INFO] Working dir: $(pwd)"
echo "[INFO] Protein file: ${PROT}"
echo "[INFO] GFF file:     ${GFF}"
echo "[INFO] UniProt fasta/db: ${UNIPROT_FASTA}"

# 简单安全检查
if [ ! -f "${PROT}" ]; then
    echo "[ERROR] Protein file not found: ${PROT}"
    exit 1
fi

if [ ! -f "${GFF}" ]; then
    echo "[ERROR] GFF file not found: ${GFF}"
    exit 1
fi

if [ ! -f "${UNIPROT_FASTA}" ]; then
    echo "[ERROR] UniProt fasta not found: ${UNIPROT_FASTA}"
    exit 1
fi

echo "[INFO] Start BLASTP vs UniProt reviewed viridiplantae"

blastp \
    -query "${PROT}" \
    -db "${UNIPROT_FASTA}" \
    -num_threads 10 \
    -outfmt 6 \
    -evalue 1e-5 \
    -max_target_seqs 10 \
    -out "${BLAST_OUT}"

echo "[INFO] BLASTP finished, now selecting best hits per query"

sort -k1,1 -k12,12g "${BLAST_OUT}" | sort -u -k1,1 --merge > "${BLAST_OUT}.besthits"

echo "[INFO] Copy original files for UniProt functional annotation"
cp "${PROT}" "${PROT}.Uniprot"
cp "${GFF}" "${GFF}.Uniprot.gff3"

echo "[INFO] Add UniProt functional annotation to FASTA"
${MAKERBIN}/maker_functional_fasta \
    "${UNIPROT_FASTA}" \
    "${BLAST_OUT}.besthits" \
    "${PROT}" \
    > "${PROT}.Uniprot"

echo "[INFO] Add UniProt functional annotation to GFF3"
${MAKERBIN}/maker_functional_gff \
    "${UNIPROT_FASTA}" \
    "${BLAST_OUT}.besthits" \
    "${GFF}" \
    > "${GFF}.Uniprot.gff3"

echo "[INFO] DONE: UniProt functional annotation ready."