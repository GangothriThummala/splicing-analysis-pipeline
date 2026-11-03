#!/usr/bin/env Rscript
# DEXSeq differential exon usage, one run per row of contrasts.csv (columns b1,b2,paired).
# usage: dexseq_deu.R counts.txt flattened.gtf samplesheet.csv contrasts.csv group_col pair_col threads outdir
# Fold change direction: b2 relative to b1 (b1 = reference, e.g. D0).
# This is the OPPOSITE sign convention to rMATS (b1 - b2).
suppressPackageStartupMessages({ library(DEXSeq); library(BiocParallel) })

a <- commandArgs(TRUE)
counts_f <- a[1]; flat_f <- a[2]; ss_f <- a[3]; ct_f <- a[4]
grp <- a[5]; pair <- a[6]; nt <- as.integer(a[7]); out <- a[8]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
bp <- MulticoreParam(nt)

ss <- read.csv(ss_f, stringsAsFactors = FALSE)
ct <- read.csv(ct_f, stringsAsFactors = FALSE)
stopifnot(all(c("b1", "b2", "paired") %in% names(ct)))

# featureCounts table: Geneid Chr Start End Strand Length <one column per BAM>
cnt <- read.delim(counts_f, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
fixed <- 1:6
names(cnt)[-fixed] <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(names(cnt)[-fixed]))

# Exon-bin numbers and transcripts come from the flattened GTF; match rows by gene + start + end
gtf <- read.delim(flat_f, header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE)
gtf <- gtf[gtf$V3 == "exonic_part", ]
g_gene <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
g_num  <- sub('.*exonic_part_number "([^"]+)".*', "\\1", gtf$V9)
g_tx   <- sub('.*transcripts "([^"]+)".*', "\\1", gtf$V9)
m <- match(paste(cnt$Geneid, cnt$Start, cnt$End), paste(g_gene, gtf$V4, gtf$V5))
stopifnot(!anyNA(m))
cnt$featureID <- g_num[m]
tx <- setNames(g_tx[m], paste(cnt$Geneid, cnt$featureID))

for (i in seq_len(nrow(ct))) {
  g1 <- ct$b1[i]; g2 <- ct$b2[i]
  paired <- tolower(as.character(ct$paired[i])) %in% c("yes", "true", "1")
  nm <- paste0(g1, "_vs_", g2)

  si <- ss[ss[[grp]] %in% c(g1, g2) & ss$sample %in% names(cnt), ]
  if (paired) {                                  # keep animals present in both groups
    both <- intersect(si[[pair]][si[[grp]] == g1], si[[pair]][si[[grp]] == g2])
    si <- si[si[[pair]] %in% both, ]
  }
  sdata <- data.frame(row.names = si$sample,
                      condition = factor(si[[grp]], levels = c(g1, g2)),
                      animal    = factor(si[[pair]]))
  message(nm, ": ", nrow(sdata), " samples, paired = ", paired)

  cm <- as.matrix(cnt[, rownames(sdata), drop = FALSE])
  storage.mode(cm) <- "integer"
  keep_bins <- rowSums(cm) >= 10                 # drop near-empty bins (speed + power)
  full    <- if (paired) ~ sample + exon + animal:exon + condition:exon else ~ sample + exon + condition:exon
  reduced <- if (paired) ~ sample + exon + animal:exon else ~ sample + exon

  dxd <- DEXSeqDataSet(countData = cm[keep_bins, , drop = FALSE], sampleData = sdata, design = full,
                       featureID = cnt$featureID[keep_bins], groupID = cnt$Geneid[keep_bins],
                       featureRanges = GRanges(cnt$Chr[keep_bins],
                                               IRanges(cnt$Start[keep_bins], cnt$End[keep_bins]),
                                               strand = cnt$Strand[keep_bins]))
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd, formula = full, BPPARAM = bp)
  dxd <- testForDEU(dxd, fullModel = full, reducedModel = reduced, BPPARAM = bp)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = bp)

  res <- DEXSeqResults(dxd)
  pgq <- perGeneQValue(res)
  keep <- setdiff(colnames(res), c("genomicData", "countData", "transcripts"))
  df <- as.data.frame(res[, keep])
  df$transcripts <- tx[paste(df$groupID, df$featureID)]

  write.table(df, file.path(out, paste0(nm, ".dexseq.exon_bins.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(gene = names(pgq), gene_qval = pgq), file.path(out, paste0(nm, ".dexseq.genes.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
