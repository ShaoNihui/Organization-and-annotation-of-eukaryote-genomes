#!/usr/bin/env Rscript

library(GENESPACE)

# ----- 参数 -----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
    stop("Usage: Rscript 13_genespace.R <genespace_workdir>")
}

wd <- args[1]
message("GENESPACE working directory: ", wd)

# ----- 初始化 GENESPACE -----
gpar <- init_genespace(
    wd = wd,
    path2mcscanx = "/usr/local/bin",
    nCores = 20,
    verbose = TRUE
)

# ----- 运行 GENESPACE -----
out <- run_genespace(gpar, overwrite = TRUE)

# ----- 构建 pangenome 矩阵 -----
pangenome <- query_pangenes(
    out,
    bed = NULL,
    refGenome = "TAIR10",
    transform = TRUE,
    showArrayMem = TRUE,
    showNSOrtho = TRUE,
    maxMem2Show = Inf
)

# 保存结果
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))

message("=== GENESPACE 完成 ===")