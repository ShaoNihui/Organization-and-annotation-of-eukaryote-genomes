#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem=32G
#SBATCH -p pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --job-name=BUSCO_N13_proteins
#SBATCH --output=/data/users/nshao/organization_annotation_course/logs/BUSCO_N13_proteins_%j.out
#SBATCH --error=/data/users/nshao/organization_annotation_course/logs/BUSCO_N13_proteins_%j.err
#SBATCH --mail-user=nihui.shao@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

# ========= 基本路径设置 =========
WORKDIR="/data/users/nshao/organization_annotation_course/gene_annotation/final"
LINEAGE="brassicales_odb10"
INPUT="genome_N13.all.maker.proteins.fasta.renamed.fasta.filtered.fasta"
OUTNAME="BUSCO_N13_proteins_brassicales"

echo "==== BUSCO on N13 proteins started at $(date) ===="
echo "WORKDIR:  $WORKDIR"
echo "INPUT:    $INPUT"
echo "LINEAGE:  $LINEAGE"
echo "OUTNAME:  $OUTNAME"
echo

cd "$WORKDIR" || { echo "Cannot cd to $WORKDIR, exiting"; exit 1; }

# 检查输入文件是否存在
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file $INPUT not found in $WORKDIR"
    ls -lh
    exit 1
fi

# ========= 加载 BUSCO 模块 =========
module load BUSCO/5.4.2-foss-2021a

# ========= 运行 BUSCO =========
busco \
    -i "$INPUT" \
    -l "$LINEAGE" \
    -o "$OUTNAME" \
    -m proteins \
    --cpu 8

echo
echo "==== BUSCO finished at $(date) ===="
echo "Results should be in:"
echo "  $WORKDIR/$OUTNAME"