#!/usr/bin/env python3
"""Raw read length for rMATS --readLength, from FastQC reports on the raw reads.

FastQC reports "Sequence length" as one number (101) or a range (35-101). The
longest read of every file must be the same: rMATS counts a read only if its
length matches --readLength (or, with --variable-read-length, is not longer),
so a wrong value silently drops reads. Fails if the files disagree.

Usage: read_length.py <out.tsv> <*_fastqc.zip>...
Writes a per-file table and prints the read length on stdout.
"""
import os
import sys
import zipfile


def sequence_length(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        name = next(n for n in z.namelist() if n.endswith("/fastqc_data.txt"))
        for line in z.read(name).decode().splitlines():
            if line.startswith("Sequence length\t"):
                return line.split("\t")[1]
    sys.exit(f"No 'Sequence length' in {zip_path}")


def main(out_tsv, zips):
    rows = []
    for z in sorted(zips):
        rng = sequence_length(z)
        rows.append((os.path.basename(z).removesuffix("_fastqc.zip"), rng, int(rng.split("-")[-1])))

    with open(out_tsv, "w") as out:
        out.write("file\tsequence_length\tmax_length\n")
        for f, rng, mx in rows:
            out.write(f"{f}\t{rng}\t{mx}\n")

    lengths = {mx for *_, mx in rows}
    if len(lengths) != 1:
        bad = [f"{f} {rng}" for f, rng, _ in rows]
        sys.exit("Raw read lengths differ between files:\n  " + "\n  ".join(bad))
    print(lengths.pop())


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:])
