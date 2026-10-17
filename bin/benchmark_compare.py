#!/usr/bin/env python3
"""Summarise an alignment benchmark (bin/benchmark_align.sh output).

Per arm and sample: wall time, peak memory, uniquely mapped %. Against the
reference arm (the first CPU arm in timing.tsv, else the first arm): uniquely mapped reads, splice junctions (from
SJ.out.tab, junctions with >= 5 unique reads; Jaccard index and the fraction of
reference junctions found), and the Pearson correlation of log gene counts
(ReadsPerGene column 2, unstranded).

Usage: benchmark_compare.py <benchmark_outdir>
"""
import csv
import math
import os
import sys


def log_final(path):
    out = {}
    for line in open(path):
        if "|" in line:
            k, v = line.split("|", 1)
            out[k.strip()] = v.strip()
    return out


def junctions(path, min_reads=5):
    return {tuple(f[:4]) for f in (l.split("\t") for l in open(path)) if int(f[6]) >= min_reads}


def gene_counts(path):
    return {f[0]: float(f[1]) for f in (l.split("\t") for l in open(path)) if not f[0].startswith("N_")}


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sab = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    sa = math.sqrt(sum((x - ma) ** 2 for x in a))
    sb = math.sqrt(sum((y - mb) ** 2 for y in b))
    return sab / (sa * sb)


def main(out):
    rows = list(csv.DictReader(open(os.path.join(out, "timing.tsv")), delimiter="\t"))
    arms = list(dict.fromkeys(r["arm"] for r in rows))
    ref = next((a for a in arms if a.startswith("cpu")), arms[0])
    arms = [ref] + [a for a in arms if a != ref]
    files = lambda arm, s, suffix: os.path.join(out, arm, s, f"{s}.{suffix}")

    print(f"Reference arm: {ref}\n")
    print("arm\tsample\twall_min\tpeak_GB\tunique_%\tspeedup_vs_ref")
    wall = {(r["arm"], r["sample"]): float(r["wall_s"]) for r in rows}
    for r in rows:
        lf = log_final(files(r["arm"], r["sample"], "Log.final.out"))
        sp = wall.get((ref, r["sample"]), float("nan")) / float(r["wall_s"])
        print(f"{r['arm']}\t{r['sample']}\t{float(r['wall_s']) / 60:.1f}\t"
              f"{int(r['max_rss_kb']) / 1e6:.1f}\t{lf.get('Uniquely mapped reads %', 'NA')}\t{sp:.2f}")

    print("\narm\tsample\tunique_reads_vs_ref\tjunction_jaccard\tref_junctions_found\tgene_count_r(log)")
    for arm in arms[1:]:
        for s in dict.fromkeys(r["sample"] for r in rows if r["arm"] == arm):
            if (ref, s) not in wall:
                continue
            a = log_final(files(arm, s, "Log.final.out"))
            b = log_final(files(ref, s, "Log.final.out"))
            ua, ub = int(a["Uniquely mapped reads number"]), int(b["Uniquely mapped reads number"])
            ja, jb = junctions(files(arm, s, "SJ.out.tab")), junctions(files(ref, s, "SJ.out.tab"))
            ga, gb = gene_counts(files(arm, s, "ReadsPerGene.out.tab")), gene_counts(files(ref, s, "ReadsPerGene.out.tab"))
            genes = [g for g in gb if g in ga and (ga[g] + gb[g]) > 0]
            r = pearson([math.log1p(ga[g]) for g in genes], [math.log1p(gb[g]) for g in genes])
            print(f"{arm}\t{s}\t{ua / ub:.4f}\t{len(ja & jb) / len(ja | jb):.4f}\t"
                  f"{len(ja & jb) / len(jb):.4f}\t{r:.5f}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
