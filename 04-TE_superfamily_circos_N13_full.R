## 04-TE_superfamily_circos_N13_full.R
## 仅使用前若干条最长 scaffold 画 circos，避免 gap.degree 错误

cat("Working directory:", getwd(), "\n")

# -----------------------------
# 0. 路径设置
# -----------------------------
fasta  <- "input/genome_N13.fasta"
gff    <- "results/EDTA_annotation/genome_N13.fasta.mod.EDTA.TEanno.gff3"
outpdf <- "results/TE_circos_plot_N13_top10.pdf"

if (!file.exists(fasta)) stop("FASTA 不存在: ", fasta)
if (!file.exists(gff))   stop("GFF 不存在: ", gff)

# -----------------------------
# 1. 依赖包
# -----------------------------
need_pkgs <- c("data.table", "circlize")
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly=TRUE)) install.packages(p)
}
library(data.table)
library(circlize)

# -----------------------------
# 2. FASTA 读取 & 计算长度
# -----------------------------
cat("[INFO] 读取 FASTA 并计算长度...\n")

fa <- readLines(fasta)
names <- fa[grep("^>", fa)]
seq_index <- grep("^>", fa)

lens <- c()
for (i in seq_along(seq_index)) {
  start <- seq_index[i] + 1
  end <- ifelse(i == length(seq_index), length(fa), seq_index[i+1] - 1)
  len <- sum(nchar(fa[start:end]))
  lens[i] <- len
}

chroms <- gsub("^>", "", names)

ideogram <- data.table(
  chr   = chroms,
  start = 0L,
  end   = as.numeric(lens)
)

cat("[INFO] 原始 scaffold 数量：", nrow(ideogram), "\n")

## ✅ 只保留前 10 条最长的 scaffold 作为“伪染色体”
N_KEEP <- 10
ideogram <- ideogram[order(-end)][1:min(N_KEEP, .N)]
cat("[INFO] 选取前 ", nrow(ideogram), " 条最长 scaffold 用于绘图：\n")
print(ideogram)

# -----------------------------
# 3. 读取 TE GFF（兼容旧 data.table）
# -----------------------------
cat("[INFO] 读取 TE 注释...\n")

gff_lines <- readLines(gff)
gff_no_comments <- gff_lines[!grepl("^#", gff_lines)]
gff_dt <- fread(text = gff_no_comments, header = FALSE)

setnames(gff_dt, c("chr","source","feature","start","end","score","strand","phase","attr"))

# 只保留这些最长 scaffold 上的 TE
keep_chr <- ideogram$chr
gff_dt <- gff_dt[chr %in% keep_chr]

# -----------------------------
# 4. 提取 TE superfamily
# -----------------------------
extract_sf <- function(x){
  x <- as.character(x)
  sf <- sub(".*[Cc]lassification=([^;]+).*", "\\1", x)
  sf[!grepl("[Cc]lassification=", x)] <- NA
  sapply(strsplit(sf, "[/|:]"), function(v) tail(v, 1L))
}

gff_dt[, superfamily := extract_sf(attr)]

gff_dt[, class := fifelse(grepl("Gypsy", superfamily, ignore.case=TRUE),"Gypsy",
                          fifelse(grepl("Copia", superfamily, ignore.case=TRUE),"Copia",
                                  fifelse(grepl("LINE",  superfamily, ignore.case=TRUE),"LINE",
                                          fifelse(grepl("DNA",   superfamily, ignore.case=TRUE),"DNA","Other"))))]

keep <- c("Gypsy", "Copia", "DNA", "LINE")
plot_dt <- gff_dt[class %in% keep]

cat("[INFO] 有效 TE 条目：", nrow(plot_dt), "\n")

# -----------------------------
# 5. 绘制 Circos
# -----------------------------
col_map <- c(
  Gypsy = "firebrick",
  Copia = "forestgreen",
  DNA   = "steelblue",
  LINE  = "darkorange"
)

pdf(outpdf, width=8, height=8)

circos.clear()
## ⚠️ gap.degree 设得很小，避免再爆
circos.par(start.degree=90, gap.degree=0.5)
circos.genomicInitialize(ideogram, plotType=c("axis","labels"))

for(cls in keep){
  sub <- plot_dt[class == cls, .(chr, start, end)]
  if (nrow(sub) == 0) next
  circos.genomicDensity(sub, col=col_map[[cls]], track.height=0.08, window.size=100000)
}

title("TE Superfamily Density – N13 (top 10 scaffolds)")
dev.off()
circos.clear()

cat("[SUCCESS] 图已生成 → ", outpdf, "\n")