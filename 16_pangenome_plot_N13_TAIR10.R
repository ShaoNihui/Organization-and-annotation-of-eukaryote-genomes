#!/usr/bin/env Rscript

## 使用方法：
##   module load R/4.3.2-foss-2021a
##   cd /data/users/nshao/organization_annotation_course/gene_annotation
##   Rscript ../scripts/16_pangenome_plot_N13_TAIR10.R genespace

suppressPackageStartupMessages({
  # 只用 base R
})

# 关键：避免走 X11，强制使用 cairo 位图
options(bitmapType = "cairo")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("请提供 GENESPACE 工作目录，例如：Rscript 16_pangenome_plot_N13_TAIR10.R genespace")
}

wd <- args[1]
cat(">>> GENESPACE 工作目录:", wd, "\n")

# 1. 读取 summary 表
summary_file <- file.path(wd, "results", "pangenome_summary.tsv")
if (!file.exists(summary_file)) {
  stop("找不到文件: ", summary_file)
}

pg <- read.table(summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat(">>> 读入 pangenome_summary.tsv，行数:", nrow(pg), "\n")
cat(">>> 列名:", paste(colnames(pg), collapse = ", "), "\n")

# 只保留关心的两个 genome
genomes <- c("N13", "TAIR10")
pg <- pg[pg$genome %in% genomes, ]

get_counts <- function(df, g) {
  core_og    <- df$n_orthogroups[df$genome == g & df$category == "core"]
  unique_og  <- df$n_orthogroups[df$genome == g & df$category == "unique"]
  core_og[is.na(core_og)]     <- 0
  unique_og[is.na(unique_og)] <- 0
  list(core = core_og[1], unique = unique_og[1])
}

N13_counts      <- get_counts(pg, "N13")
TAIR10_counts   <- get_counts(pg, "TAIR10")

core_og         <- N13_counts$core          # 核心 orthogroups（两者共享）
unique_N13_og   <- N13_counts$unique
unique_TAIR10_og<- TAIR10_counts$unique

cat(">>> 核心 orthogroups:", core_og, "\n")
cat(">>> N13 unique orthogroups:", unique_N13_og, "\n")
cat(">>> TAIR10 unique orthogroups:", unique_TAIR10_og, "\n")

#============================#
# 2. 画 Venn 图 (orthogroups)
#============================#

venn_png <- file.path(wd, "results", "pangenome_venn_N13_TAIR10.png")
png(venn_png, width = 1400, height = 900, res = 150, type = "cairo")
par(mar = c(4, 4, 4, 4))
plot(0, 0, type = "n",
     xlim = c(0, 10), ylim = c(0, 10),
     xlab = "", ylab = "", axes = FALSE,
     main = "Pangenome orthogroups: N13 vs TAIR10")

# 两个圆的中心
x1 <- 4; y1 <- 5
x2 <- 6; y2 <- 5
r  <- 3

# 画圆
symbols(x1, y1, circles = r, inches = FALSE,
        add = TRUE, bg = adjustcolor("skyblue", alpha.f = 0.3),
        fg = "skyblue4")
symbols(x2, y2, circles = r, inches = FALSE,
        add = TRUE, bg = adjustcolor("pink", alpha.f = 0.3),
        fg = "red3")

# 标注 genome 名字
text(x1, y1 + r + 0.5, "N13",    cex = 1.4, col = "blue4")
text(x2, y2 + r + 0.5, "TAIR10", cex = 1.4, col = "red4")

# 区域数字
only_N13  <- unique_N13_og
only_TAIR <- unique_TAIR10_og
intersect <- core_og

text(x1 - 1,        y1, only_N13,  cex = 1.3)
text((x1 + x2) / 2, y1, intersect, cex = 1.3)
text(x2 + 1,        y2, only_TAIR, cex = 1.3)

legend("topleft",
       legend = c(
         paste0("N13 unique OGs: ", only_N13),
         paste0("TAIR10 unique OGs: ", only_TAIR),
         paste0("Core OGs: ", intersect)
       ),
       bty = "n", cex = 1.1)

dev.off()
cat(">>> Venn 图已输出:", venn_png, "\n")

#================================#
# 3. 画柱状图：core vs unique 数量
#================================#

bar_png <- file.path(wd, "results", "pangenome_bar_N13_TAIR10.png")
png(bar_png, width = 1400, height = 900, res = 150, type = "cairo")
par(mar = c(5, 5, 4, 2))

bar_mat <- rbind(
  core   = c(N13 = N13_counts$core,   TAIR10 = TAIR10_counts$core),
  unique = c(N13 = N13_counts$unique, TAIR10 = TAIR10_counts$unique)
)

barplot(bar_mat,
        beside = TRUE,
        col = c("steelblue3", "orange2"),
        ylim = c(0, max(bar_mat) * 1.2),
        ylab = "Number of orthogroups",
        main = "Core vs Unique orthogroups (N13 vs TAIR10)")

legend("topright",
       legend = c("Core orthogroups", "Unique orthogroups"),
       fill   = c("steelblue3", "orange2"),
       bty    = "n")

dev.off()
cat(">>> 柱状图已输出:", bar_png, "\n")
cat(">>> 完成 pangenome 可视化。\n")