// rMATS-turbo 4.4.0 in three stages, so every contrast shares one event set:
//   RMATS_PREP  per BAM: junction and exon read counts (.rmats)
//   RMATS_POST  once, over all BAMs: event set (fromGTF.*) and counts per sample
//   RMATS_STAT  per contrast: picks its samples out of the post output
//               (prepare_stat_inputs.py) and runs the statistical model.
// Settings shared by prep and post. --allow-clipping keeps soft-clipped reads
// (5.5% more usable reads on s001); NH != 1 reads are always dropped by rMATS.
def rmats_args(read_length, lib_type) {
    "-t paired --readLength ${read_length} --variable-read-length --libType ${lib_type} --allow-clipping"
}

// Samples to compare, by position in RMATS_POST's BAM list (bin/rmats_design.py).
process RMATS_DESIGN {
    publishDir "${params.outdir}/rmats", mode: 'copy'

    input:
    path samplesheet
    path contrasts

    output:
    path 'post_b1.txt',   emit: post_b1
    path 'contrasts.tsv', emit: contrasts

    script:
    """
    rmats_design.py ${samplesheet} ${contrasts} ${params.group_column} ${params.pair_column}
    """
}

process RMATS_PREP {
    tag "${meta.sample}"

    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    val read_length
    val lib_type

    output:
    path "${meta.sample}_*.rmats",                    emit: rmats
    path "${meta.sample}_*read_outcomes_by_bam.txt", emit: outcomes

    script:
    // The BAM is listed by file name only. rMATS post matches prep output to
    // BAMs by this exact string, and the file name is the same in every task.
    """
    printf '%s' ${bam} > b1.txt
    rmats.py --task prep --b1 b1.txt --gtf ${gtf} --od od --tmp tmp \\
        ${rmats_args(read_length, lib_type)} --nthread 1 > prep.log 2>&1 \\
        || { cat prep.log >&2; exit 1; }
    # prep names its files by timestamp; prefix the sample so they cannot collide.
    for f in tmp/*.rmats tmp/*read_outcomes_by_bam.txt; do mv "\$f" "${meta.sample}_\$(basename "\$f")"; done
    """
}

process RMATS_POST {
    publishDir "${params.outdir}/rmats", mode: 'copy', pattern: '{post,read_outcomes.tsv}'

    input:
    path rmats_files, stageAs: 'tmp/*'
    path outcomes
    path post_b1
    path gtf
    val read_length
    val lib_type

    output:
    path 'post',              emit: post
    path 'read_outcomes.tsv'

    script:
    // --task inte first: every BAM in the list has exactly one prep file,
    // and every prep file belongs to a listed BAM.
    """
    rmats_read_outcomes.py read_outcomes.tsv ${outcomes}
    rmats.py --task inte --b1 ${post_b1} --gtf ${gtf} --od post --tmp tmp ${rmats_args(read_length, lib_type)}
    rmats.py --task post --b1 ${post_b1} --gtf ${gtf} --od post --tmp tmp \\
        ${rmats_args(read_length, lib_type)} --nthread ${task.cpus} --statoff
    rm -rf post/tmp
    """
}

process RMATS_STAT {
    tag "${contrast}"
    publishDir "${params.outdir}/rmats/contrasts", mode: 'copy'

    input:
    tuple val(contrast), val(paired), val(g1), val(g2)
    path post

    output:
    tuple val(contrast), path("${contrast}")

    script:
    def model = paired == 'yes' ? '--paired-stats' : ''
    // rMATS exits 0 even when the paired model (PAIRADISE) fails; it only
    // prints "error in paired model". Check for that, and for p-values.
    """
    P=\$(dirname "\$(readlink -f "\$(which rmats.py)")")/rMATS_P/prepare_stat_inputs.py
    python "\$P" --new-output-dir ${contrast} --old-output-dir ${post} \\
        --group-1-indices ${g1} --group-2-indices ${g2}
    rmats.py --task stat --od ${contrast} --tmp stat_tmp --nthread ${task.cpus} ${model} \\
        > ${contrast}/stat.log 2>&1 || { cat ${contrast}/stat.log >&2; exit 1; }
    if grep -q 'error in paired model' ${contrast}/stat.log; then
        cat ${contrast}/stat.log >&2; exit 1
    fi
    awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "PValue") c = i; next }
                 \$c != "NA" { n++ } END { exit n > 0 ? 0 : 1 }' ${contrast}/SE.MATS.JC.txt \\
        || { echo "No p-values in ${contrast}/SE.MATS.JC.txt" >&2; exit 1; }
    rm -rf ${contrast}/tmp
    """
}

// Label-blind coverage and PSI filters, then FDR recomputed on the events that
// pass (bin/rmats_filter.py). The rMATS tables above are left as they are.
process RMATS_FILTER {
    tag "${contrast}"
    publishDir "${params.outdir}/rmats/filtered", mode: 'copy'

    input:
    tuple val(contrast), path(contrast_dir, stageAs: 'in/*')

    output:
    path "${contrast}"

    script:
    """
    rmats_filter.py ${contrast_dir} ${contrast} --min-reads ${params.filter_min_reads} \\
        --min-psi ${params.filter_min_psi} --fdr ${params.filter_fdr} --min-dpsi ${params.filter_min_dpsi}
    """
}
