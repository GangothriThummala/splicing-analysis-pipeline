#!/usr/bin/env python3
"""Collect rMATS prep read outcomes (*read_outcomes_by_bam.txt) into one table.

Each file is a BAM path, then "OUTCOME: count" lines (USED, NOT_NH_1, CLIPPED,
...). Fails if any read was rejected as NOT_EXPECTED_READ_LENGTH: with the raw
read length and --variable-read-length that count must be 0, so a non-zero
value means --readLength is wrong.

Usage: rmats_read_outcomes.py <out.tsv> <*read_outcomes_by_bam.txt>...
"""
import sys


def main(out_tsv, files):
    rows, keys = [], []
    for path in sorted(files):
        with open(path) as fh:
            bam = fh.readline().strip()
            counts = {}
            for line in fh:
                k, v = line.split(":")
                counts[k.strip()] = int(v)
                if k.strip() not in keys:
                    keys.append(k.strip())
        sample = bam.split("/")[-1].split(".Aligned")[0]
        rows.append((sample, counts))

    with open(out_tsv, "w") as out:
        out.write("\t".join(["sample"] + keys + ["used_frac"]) + "\n")
        for s, c in rows:
            used = c["USED"] / c["TOTAL_FOR_BAM"] if c.get("TOTAL_FOR_BAM") else 0
            out.write("\t".join([s] + [str(c.get(k, 0)) for k in keys] + [f"{used:.4f}"]) + "\n")

    bad = [f"{s}: {c['NOT_EXPECTED_READ_LENGTH']}" for s, c in rows
           if c.get("NOT_EXPECTED_READ_LENGTH", 0) > 0]
    if bad:
        sys.exit("Reads rejected for length (check --readLength):\n  " + "\n  ".join(bad))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
