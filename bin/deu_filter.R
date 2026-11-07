#!/usr/bin/env Rscript
# Post-processing filter for DEU results (like rmats_filter.py).
# Keeps genes with gene-level q below --max-q that have at least one exon bin with
# padj below --max-bin-padj and |usage change| of at least --min-usage-change.
# usage change = (mean over animals of the bin's share of its gene's counts at b2)
#                minus (the same at b1), i.e. later minus earlier. Paired contrasts use only
#                animals present at both timepoints.
# usage: deu_filter.R --deu-dir DIR --samplesheet CSV --contrasts CSV [--group-column timepoint]
#        [--pair-column animal] [--min-usage-change 0.1] [--max-q 0.05] [--max-bin-padj 0.05] [--out deu_filtered]
# DIR is the --outdir/deu folder written by deu_only.nf.
args <- commandArgs(TRUE)
opt <- list("deu-dir"=NA, "samplesheet"=NA, "contrasts"=NA, "group-column"="timepoint", "pair-column"="animal",
            "min-usage-change"="0.1", "max-q"="0.05", "max-bin-padj"="0.05", "out"="deu_filtered")
i <- 1
while (i <= length(args)) {
  k <- sub("^--", "", args[i])
  if (!k %in% names(opt) || i == length(args)) stop("bad option: ", args[i])
  opt[[k]] <- args[i + 1]; i <- i + 2
}
for (k in c("deu-dir", "samplesheet", "contrasts")) if (is.na(opt[[k]])) stop("missing --", k)
mu <- as.numeric(opt[["min-usage-change"]]); maxq <- as.numeric(opt[["max-q"]]); maxp <- as.numeric(opt[["max-bin-padj"]])
deu <- opt[["deu-dir"]]; grp <- opt[["group-column"]]; pair <- opt[["pair-column"]]
dir.create(opt[["out"]], showWarnings = FALSE, recursive = TRUE)

ss <- read.csv(opt[["samplesheet"]], stringsAsFactors = FALSE)
ct <- read.csv(opt[["contrasts"]], stringsAsFactors = FALSE)
stopifnot(all(c("b1", "b2", "paired") %in% names(ct)))

cnt <- read.delim(file.path(deu, "counts", "deu_exonic_part_counts.txt"), comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
names(cnt)[-(1:6)] <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(names(cnt)[-(1:6)]))
gtf <- read.delim(file.path(deu, "annotation", "flattened.dexseq.gtf"), header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE)
gtf <- gtf[gtf$V3 == "exonic_part", ]
gg <- sub('.*gene_id "([^"]+)".*', "\\1", gtf$V9)
gn <- sub('.*exonic_part_number "([^"]+)".*', "\\1", gtf$V9)
m <- match(paste(cnt$Geneid, cnt$Start, cnt$End), paste(gg, gtf$V4, gtf$V5)); stopifnot(!anyNA(m))
cnt$featureID <- gn[m]
smp <- setdiff(names(cnt), c("Geneid", "Chr", "Start", "End", "Strand", "Length", "featureID"))
mat <- as.matrix(cnt[, smp])
prop <- mat / rowsum(mat, cnt$Geneid)[cnt$Geneid, , drop = FALSE]   # bin share of its gene's counts, per sample
key <- paste(cnt$Geneid, cnt$featureID)

summ <- list()
for (i in seq_len(nrow(ct))) {
  g1 <- ct$b1[i]; g2 <- ct$b2[i]; nm <- paste0(g1, "_vs_", g2)
  paired <- tolower(as.character(ct$paired[i])) %in% c("yes", "true", "1")
  si <- ss[ss[[grp]] %in% c(g1, g2) & ss$sample %in% smp, ]
  if (paired) {
    both <- intersect(si[[pair]][si[[grp]] == g1], si[[pair]][si[[grp]] == g2])
    si <- si[si[[pair]] %in% both, ]
  }
  s1 <- si$sample[si[[grp]] == g1]; s2 <- si$sample[si[[grp]] == g2]
  delta <- rowMeans(prop[, s2, drop = FALSE], na.rm = TRUE) - rowMeans(prop[, s1, drop = FALSE], na.rm = TRUE)

  g <- read.delim(file.path(deu, "dexseq", "dexseq", paste0(nm, ".dexseq.genes.tsv")), stringsAsFactors = FALSE)
  b <- read.delim(file.path(deu, "dexseq", "dexseq", paste0(nm, ".dexseq.exon_bins.tsv")), colClasses = c(featureID = "character"), stringsAsFactors = FALSE)
  b$usage_change <- delta[match(paste(b$groupID, b$featureID), key)]
  genes_ok <- g$gene[!is.na(g$gene_qval) & g$gene_qval < maxq]
  hit <- b[!is.na(b$padj) & b$padj < maxp & !is.na(b$usage_change) & abs(b$usage_change) >= mu & b$groupID %in% genes_ok, ]
  hit <- hit[order(-abs(hit$usage_change)), ]
  top <- hit[!duplicated(hit$groupID), ]
  gene_tab <- data.frame(gene = top$groupID, gene_qval = g$gene_qval[match(top$groupID, g$gene)],
                         n_bins = as.vector(table(hit$groupID)[top$groupID]), top_bin = top$featureID,
                         top_usage_change = round(top$usage_change, 4), top_bin_padj = top$padj)
  lfc <- grep("^log2fold", names(hit), value = TRUE)[1]
  cols <- intersect(c("groupID", "featureID", "exonBaseMean", "padj", lfc, "usage_change", "transcripts"), names(hit))
  write.table(gene_tab, file.path(opt[["out"]], paste0(nm, ".deu_filtered_genes.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(hit[, cols], file.path(opt[["out"]], paste0(nm, ".deu_filtered_bins.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  summ[[nm]] <- data.frame(contrast = nm, samples = nrow(si), genes_tested = nrow(g), genes_q = length(genes_ok),
                           genes_pass = nrow(gene_tab), genes_pass_split_ids = length(unique(unlist(strsplit(gene_tab$gene, "[+]")))),
                           min_usage_change = mu, max_q = maxq, max_bin_padj = maxp)
}
write.table(do.call(rbind, summ), file.path(opt[["out"]], "deu_filter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
