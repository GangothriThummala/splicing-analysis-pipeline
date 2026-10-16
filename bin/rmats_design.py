#!/usr/bin/env python3
"""Turn the samplesheet and contrast table into rMATS inputs.

rMATS post runs once over every sample, so all contrasts share one event set
(same event IDs). Each contrast then selects its samples by 0-based index into
that post run (rMATS_P/prepare_stat_inputs.py).

Paired contrasts: PAIRADISE pairs sample i of group 1 with sample i of group 2,
by position only. So both groups are restricted to animals present in both and
put in the same animal order. Animals present in only one group are dropped
and listed.

Usage: rmats_design.py <samplesheet.csv> <contrasts.csv> <group_column> <pair_column>
Writes post_b1.txt (BAM names, samplesheet order) and contrasts.tsv.
"""
import csv
import sys

BAM = "{}.Aligned.sortedByCoord.out.bam"   # name written by bin/star_align.sh
YES, NO = {"yes", "true", "1"}, {"no", "false", "0"}


def main(sheet, contrasts, group_col, pair_col):
    samples = list(csv.DictReader(open(sheet)))
    for col in (group_col, pair_col):
        if col not in samples[0]:
            sys.exit(f"Column '{col}' not in {sheet}")
    groups = {s[group_col] for s in samples}

    with open("post_b1.txt", "w") as fh:
        fh.write(",".join(BAM.format(s["sample"]) for s in samples))

    out = open("contrasts.tsv", "w")
    out.write("contrast\tpaired\tn1\tn2\tg1_indices\tg2_indices\tg1_samples\tg2_samples\tdropped\n")
    for c in csv.DictReader(open(contrasts)):
        b1, b2, paired = c["b1"], c["b2"], c["paired"].strip().lower()
        for g in (b1, b2):
            if g not in groups:
                sys.exit(f"Contrast {b1} vs {b2}: no samples with {group_col} = '{g}'")
        if paired not in YES | NO:
            sys.exit(f"Contrast {b1} vs {b2}: paired must be yes or no, got '{c['paired']}'")
        i1 = [i for i, s in enumerate(samples) if s[group_col] == b1]
        i2 = [i for i, s in enumerate(samples) if s[group_col] == b2]
        dropped = []

        if paired in YES:
            by_pair = []
            for g, idx in ((b1, i1), (b2, i2)):
                m = {}
                for i in idx:
                    p = samples[i][pair_col]
                    if p in m:
                        sys.exit(f"{g}: two samples for {pair_col} {p} "
                                 f"({samples[m[p]]['sample']}, {samples[i]['sample']})")
                    m[p] = i
                by_pair.append(m)
            common = sorted(set(by_pair[0]) & set(by_pair[1]))
            dropped = sorted(set(by_pair[0]) ^ set(by_pair[1]))
            if len(common) < 2:
                sys.exit(f"Contrast {b1} vs {b2}: {len(common)} pairs, need at least 2")
            i1 = [by_pair[0][p] for p in common]
            i2 = [by_pair[1][p] for p in common]

        out.write("\t".join([
            f"{b1}_vs_{b2}", "yes" if paired in YES else "no", str(len(i1)), str(len(i2)),
            ",".join(map(str, i1)), ",".join(map(str, i2)),
            ",".join(samples[i]["sample"] for i in i1),
            ",".join(samples[i]["sample"] for i in i2),
            ",".join(dropped) or "-"]) + "\n")
    out.close()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
