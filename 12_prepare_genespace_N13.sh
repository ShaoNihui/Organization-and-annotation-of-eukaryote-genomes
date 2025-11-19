#!/bin/bash
#SBATCH --job-name=prep_GENESPACE_N13
#SBATCH -p pibu_el8
#SBATCH --output=../logs/12_prep_GENESPACE_N13.%j.out
#SBATCH --error=../logs/12_prep_GENESPACE_N13.%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=nihui.shao@students.unibe.ch

set -euo pipefail

echo "=== [GENESPACE 准备：N13 + TAIR10] 开始 ==="

# 根目录（gene_annotation）
ROOT="/data/users/nshao/organization_annotation_course/gene_annotation"

# 注释结果目录（final/）
ANNOT_DIR="${ROOT}/final"

# GENESPACE 工作目录
GENESPACE_WD="${ROOT}/genespace"
BED_DIR="${GENESPACE_WD}/bed"
PEP_DIR="${GENESPACE_WD}/peptide"

mkdir -p "${BED_DIR}" "${PEP_DIR}"

cd "${ANNOT_DIR}"

########################################
# [1] 生成 N13.bed（基因坐标）
########################################
echo "[1] 生成 N13.bed（基因坐标）..."

GFF="filtered.genes.renamed.gff3"

# 只保留 gene 行
grep -P "\tgene\t" "${GFF}" > temp_genes.gff3

# GFF 是 1-based，BED 要 0-based：chr, start-1, end, geneID
# 第 9 列第一个属性为 ID=GENEID
awk 'BEGIN{OFS="\t"}{
    split($9,a,";");
    split(a[1],b,"=");
    print $1, $4-1, $5, b[2]
}' temp_genes.gff3 > "${BED_DIR}/N13.bed"

echo "  -> N13.bed: ${BED_DIR}/N13.bed"
echo "  -> 基因数: $(wc -l < ${BED_DIR}/N13.bed)"

########################################
# [2] 生成 N13.fa（每个 gene 的最长蛋白，header=基因ID）
########################################
echo "[2] 生成 N13.fa（每个 gene 的最长蛋白，header 仅基因 ID）..."

PROT_FASTA="genome_N13.all.maker.proteins.fasta.renamed.fasta.filtered.fasta"

# 2.1 计算每条蛋白序列长度
awk '
  /^>/ {
      if (id != "") {
          print id, length(seq);
      }
      id=$1;
      gsub(/^>/,"",id);
      seq="";
      next;
  }
  {
      gsub(/[ \t\r\n]/,"");
      seq = seq $0;
  }
  END {
      if (id != "") {
          print id, length(seq);
      }
  }
' "${PROT_FASTA}" > N13_protein_lengths.txt

# 2.2 对每个 gene（去掉 -R*）只保留最长 isoform
awk '
{
   id=$1; len=$2;
   gene=id;
   sub(/-R.*/, "", gene);  # 去掉 -RA/-RB 等
   if (len > maxlen[gene]) {
       maxlen[gene]=len;
       best[gene]=id;
   }
}
END{
   for (g in best) print best[g];
}
' N13_protein_lengths.txt > N13_longest_protein_ids.txt

echo "  -> 每个基因的最长 isoform 数量: $(wc -l < N13_longest_protein_ids.txt)"

# 2.3 用 awk 自己从 FASTA 里筛选这些 isoform，并把 header 换成基因 ID
awk '
  NR==FNR {
      keep[$1]=1;
      next;
  }
  /^>/{
      id=$1;
      sub(/^>/,"",id);
      if (keep[id]) {
          gene=id;
          sub(/-R.*/, "", gene);
          print ">" gene;
          printing=1;
      } else {
          printing=0;
      }
      next;
  }
  printing {print}
' N13_longest_protein_ids.txt "${PROT_FASTA}" > "${PEP_DIR}/N13.fa"

echo "  -> N13.fa: ${PEP_DIR}/N13.fa"
echo "  -> 条目数: $(grep -c '^>' ${PEP_DIR}/N13.fa)"

########################################
# [3] 拷贝 TAIR10 参考 BED 与蛋白
########################################
echo "[3] 拷贝 TAIR10 参考 BED 和蛋白..."

COURSEDIR="/data/courses/assembly-annotation-course/CDS_annotation"

cp "${COURSEDIR}/data/TAIR10.fa"  "${PEP_DIR}/TAIR10.fa"
cp "${COURSEDIR}/data/TAIR10.bed" "${BED_DIR}/TAIR10.bed"

echo "  -> TAIR10.fa:  ${PEP_DIR}/TAIR10.fa"
echo "  -> TAIR10.bed: ${BED_DIR}/TAIR10.bed"

echo "=== 完成 N13 & TAIR10 的 GENESPACE 输入准备（已修正 ID 和冒号问题，无 faSomeRecords 依赖） ==="
echo "BED 目录:    ${BED_DIR}"
echo "PEPTIDE 目录: ${PEP_DIR}"