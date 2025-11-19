#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("用法: Rscript 15_pangenome_summary_N13.R <genespace_workdir>")
}
wd <- args[1]

cat(">>> 读取 pangenome 矩阵: ", file.path(wd, "pangenome_matrix.rds"), "\n", sep = "")
pg <- readRDS(file.path(wd, "pangenome_matrix.rds"))

df <- as.data.frame(pg)
cols <- colnames(df)
cat(">>> data.frame 列名: ", paste(cols, collapse = ", "), "\n")

meta_cols <- c("pgID","interpChr","interpOrd","og","repGene",
               "genome","chr","start","end")

genome_cols <- cols[!(cols %in% meta_cols)]
cat(">>> 识别到 genome 列: ", paste(genome_cols, collapse = ", "), "\n")

if (length(genome_cols) < 2) {
  stop("基因组数量不足 (<2)")
}

n_og <- nrow(df)
n_genome <- length(genome_cols)
cat(">>> orthogroups 个数: ", n_og, "\n")
cat(">>> genome 数: ", n_genome, " -> ", paste(genome_cols, collapse = ", "), "\n")

presence_fun <- function(x) {
  if (is.null(x)) return(FALSE)
  if (length(x) == 0) return(FALSE)

  if (is.numeric(x)) return(any(x > 0, na.rm = TRUE))
  if (is.logical(x)) return(any(x, na.rm = TRUE))

  x <- as.character(x)
  if (all(is.na(x))) return(FALSE)
  any(x != "" & !is.na(x))
}

presence_mat <- matrix(FALSE, nrow = n_og, ncol = n_genome,
                       dimnames = list(NULL, genome_cols))

for (j in seq_along(genome_cols)) {
  col <- genome_cols[j]
  presence_mat[, j] <- vapply(df[[col]],
                              presence_fun,
                              logical(1))
}

n_present <- rowSums(presence_mat)

core_idx      <- which(n_present == n_genome)
accessory_idx <- which(n_present > 1 & n_present < n_genome)
unique_idx    <- which(n_present == 1)

cat(">>> 全局 core orthogroups: ", length(core_idx), "\n")
cat(">>> 全局 accessory orthogroups: ", length(accessory_idx), "\n")
cat(">>> 全局 unique orthogroups: ", length(unique_idx), "\n")

count_genes_in_set <- function(row_idx, genome_col) {
  if (length(row_idx) == 0) return(0L)
  vals <- df[[genome_col]][row_idx]

  s <- 0L
  for (v in vals) {
    if (is.null(v) || length(v) == 0) next
    if (is.numeric(v)) {
      s <- s + sum(v, na.rm = TRUE)
    } else {
      v <- as.character(v)
      v <- v[!is.na(v) & v != ""]
      s <- s + length(v)
    }
  }
  return(s)
}

res <- list()

for (g in genome_cols) {
  core_g      <- core_idx[presence_mat[core_idx, g]]
  acc_g       <- accessory_idx[presence_mat[accessory_idx, g]]
  uniq_g      <- unique_idx[presence_mat[unique_idx, g]]

  res[[length(res)+1]] <- data.frame(genome=g, category="core",
                                     n_orthogroups=length(core_g),
                                     n_genes=count_genes_in_set(core_g, g))

  res[[length(res)+1]] <- data.frame(genome=g, category="accessory",
                                     n_orthogroups=length(acc_g),
                                     n_genes=count_genes_in_set(acc_g, g))

  res[[length(res)+1]] <- data.frame(genome=g, category="unique",
                                     n_orthogroups=length(uniq_g),
                                     n_genes=count_genes_in_set(uniq_g, g))
}

res_df <- do.call(rbind, res)

outdir <- file.path(wd, "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(outdir, "pangenome_summary.tsv")
write.table(res_df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)

cat(">>> 结果已写入: ", outfile, "\n")
cat(">>> 预览:\n")
print(res_df)