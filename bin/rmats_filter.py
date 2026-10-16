#!/usr/bin/env python3
"""Filter one rMATS contrast and recompute FDR on the events that remain.

rMATS tests every event, including events with a handful of reads, and applies
Benjamini-Hochberg over all of them. Events that cannot reach significance only
dilute the correction, and low-count events give unstable PSI estimates.

The filters use all samples of the contrast together and ignore the group
labels (independent filtering, Bourgon et al. 2010, PNAS), so they do not bias
the test; FDR is then recomputed over the events that pass:
  - p-value is not NA;
  - mean junction reads (IJC + SJC) per sample >= --min-reads;
  - mean PSI over all samples within [--min-psi, 1 - --min-psi]
    (events that are almost always included or always skipped carry no signal).
An event is significant when FDR_filtered < --fdr and |IncLevelDifference| >= --min-dpsi.
IncLevelDifference is group 1 - group 2 (b1 - b2 in the contrast table).

Usage: rmats_filter.py <contrast_dir> <out_dir> [--min-reads 10] [--min-psi 0.05]
                       [--fdr 0.05] [--min-dpsi 0.1]
Writes <out_dir>/<type>.MATS.<JC|JCEC>.filtered.txt (passing events, with MeanReads,
MeanPSI and FDR_filtered) and <out_dir>/summary.tsv.
"""
import argparse
import csv
import os

TYPES = ["SE", "A5SS", "A3SS", "MXE", "RI"]


def numbers(s):
    return [float(v) for v in s.split(",") if v not in ("NA", "")]


def bh(pvalues):
    """Benjamini-Hochberg adjusted p-values, in input order."""
    m = len(pvalues)
    order = sorted(range(m), key=lambda i: pvalues[i])
    q, running = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * m / rank)
        q[i] = running
    return q


def filter_file(path, out_path, a):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields, rows = reader.fieldnames, list(reader)
    kept = []
    for r in rows:
        if r["PValue"] == "NA":
            continue
        reads = [sum(x) for x in zip(*(numbers(r[c]) for c in
                 ("IJC_SAMPLE_1", "SJC_SAMPLE_1")))] + \
                [sum(x) for x in zip(*(numbers(r[c]) for c in
                 ("IJC_SAMPLE_2", "SJC_SAMPLE_2")))]
        psi = numbers(r["IncLevel1"]) + numbers(r["IncLevel2"])
        if not reads or not psi:
            continue
        mean_reads, mean_psi = sum(reads) / len(reads), sum(psi) / len(psi)
        if mean_reads < a.min_reads or not a.min_psi <= mean_psi <= 1 - a.min_psi:
            continue
        r["MeanReads"], r["MeanPSI"] = f"{mean_reads:.1f}", f"{mean_psi:.4f}"
        kept.append(r)

    for r, q in zip(kept, bh([float(r["PValue"]) for r in kept])):
        r["FDR_filtered"] = f"{q:.6g}"
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["MeanReads", "MeanPSI", "FDR_filtered"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept)

    sig = [r for r in kept if float(r["FDR_filtered"]) < a.fdr
           and abs(float(r["IncLevelDifference"])) >= a.min_dpsi]
    up = sum(float(r["IncLevelDifference"]) > 0 for r in sig)
    tested = sum(r["PValue"] != "NA" for r in rows)
    return [len(rows), tested, len(kept), len(sig), up, len(sig) - up]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("contrast_dir")
    p.add_argument("out_dir")
    p.add_argument("--min-reads", type=float, default=10)
    p.add_argument("--min-psi", type=float, default=0.05)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min-dpsi", type=float, default=0.1)
    a = p.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    with open(os.path.join(a.out_dir, "summary.tsv"), "w") as out:
        out.write(f"# min_reads={a.min_reads:g} min_psi={a.min_psi:g} fdr={a.fdr:g} "
                  f"min_dpsi={a.min_dpsi:g}; up = higher inclusion in group 1\n")
        out.write("event_type\tcounts\tevents\ttested\tpass_filter\tsignificant\tup\tdown\n")
        for t in TYPES:
            for counts in ("JC", "JCEC"):
                src = os.path.join(a.contrast_dir, f"{t}.MATS.{counts}.txt")
                dst = os.path.join(a.out_dir, f"{t}.MATS.{counts}.filtered.txt")
                n = filter_file(src, dst, a)
                out.write("\t".join([t, counts] + list(map(str, n))) + "\n")


if __name__ == "__main__":
    main()
