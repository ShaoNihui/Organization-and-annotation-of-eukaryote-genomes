#!/usr/bin/env Rscript

# ---------------------------
# Full-length LTR identity histogram (N13)
# 纯 base R，无需任何 package
# ---------------------------

# -------- 正确文件路径 --------
gff_file <- "/data/users/nshao/organization_annotation_course/results/EDTA_annotation/genome_N13.fasta.mod.EDTA.intact.gff3"

cat("Reading GFF file:", gff_file, "\n")

# -------- 检查文件 --------
if (!file.exists(gff_file)) {
    stop(paste("GFF file not found:", gff_file))
}

# -------- 读取 --------
gff <- read.table(gff_file, sep="\t", header=FALSE, comment.char="#", quote="")

# 提取 attributes 列
attrs <- gff$V9

# 提取 identity=XX 信息
identity <- sapply(attrs, function(x) {
    m <- regmatches(x, regexpr("ltr_identity=[0-9.]+", x))
    if(length(m) > 0) as.numeric(sub("ltr_identity=", "", m)) else NA
})

identity <- identity[!is.na(identity)]

cat("Extracted", length(identity), "LTR identity values.\n")

# -------- 输出文件 --------
out_pdf <- "/data/users/nshao/organization_annotation_course/results/figures/N13_LTR_identity_baseR.pdf"

# -------- 画图 --------
pdf(out_pdf, width=8, height=5)
hist(identity,
     breaks=30,
     col="grey80",
     border="black",
     main="Full-length LTR-RT Identity (N13)",
     xlab="Percent identity between LTR pairs (%)"
)
dev.off()

cat("[DONE] Saved plot to:\n", out_pdf, "\n")