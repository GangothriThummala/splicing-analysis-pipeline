#!/usr/bin/env Rscript
# edgeR diffSpliceDGE differential exon usage, one run per row of contrasts.csv.
# usage: edger_deu.R counts.txt samplesheet.csv contrasts.csv group_col pair_col outdir
suppressPackageStartupMessages(library(edgeR))

a <- commandArgs(TRUE)
counts_f <- a[1]; ss_f <- a[2]; ct_f <- a[3]; grp <- a[4]; pair <- a[5]; out <- a[6]
dir.create(out, showWarnings = FALSE, recursive = TRUE)

ss <- read.csv(ss_f, stringsAsFactors = FALSE)
ct <- read.csv(ct_f, stringsAsFactors = FALSE)
stopifnot(all(c("b1", "b2", "paired") %in% names(ct)))

cnt <- read.delim(counts_f, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
fixed <- 1:6
names(cnt)[-fixed] <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(names(cnt)[-fixed]))
anno <- data.frame(GeneID = cnt$Geneid,
                   ExonID = paste(cnt$Geneid, cnt$Chr, cnt$Start, cnt$End, sep = ":"),
                   stringsAsFactors = FALSE)
mat <- as.matrix(cnt[, -fixed]); rownames(mat) <- anno$ExonID

for (i in seq_len(nrow(ct))) {
  g1 <- ct$b1[i]; g2 <- ct$b2[i]
  paired <- tolower(as.character(ct$paired[i])) %in% c("yes", "true", "1")
  nm <- paste0(g1, "_vs_", g2)

  si <- ss[ss[[grp]] %in% c(g1, g2) & ss$sample %in% colnames(mat), ]
  if (paired) {
    both <- intersect(si[[pair]][si[[grp]] == g1], si[[pair]][si[[grp]] == g2])
    si <- si[si[[pair]] %in% both, ]
  }
  condition <- factor(si[[grp]], levels = c(g1, g2))
  animal <- factor(si[[pair]])
  design <- if (paired) model.matrix(~ animal + condition) else model.matrix(~ condition)

  y <- DGEList(counts = mat[, si$sample], genes = anno)
  y <- y[filterByExpr(y, design), , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  ds <- diffSpliceDGE(fit, coef = ncol(design), geneid = "GeneID", exonid = "ExonID")

  write.table(topSpliceDGE(ds, test = "Simes", number = Inf), file.path(out, paste0(nm, ".edger.genes_simes.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(topSpliceDGE(ds, test = "exon", number = Inf), file.path(out, paste0(nm, ".edger.exons.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
