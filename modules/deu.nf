// Differential Exon Usage (DEU) branch, modelled on nf-core/rnasplice (red path):
//   flatten GTF -> featureCounts on exonic parts -> DEXSeq + edgeR diffSpliceDGE (per contrast)
// No conda directive: like modules/rmats.nf, tools come from the activated environment.

process DEU_ANNOTATION {
    tag 'flatten_gtf'
    publishDir "${params.outdir}/deu/annotation", mode: 'copy'

    input:
    path gtf

    output:
    path 'flattened.dexseq.gtf', emit: flat_gtf

    script:
    """
    DEXSEQ_PY=\$(Rscript -e 'cat(system.file("python_scripts", package="DEXSeq"))')
    python \$DEXSEQ_PY/dexseq_prepare_annotation.py ${gtf} flattened.dexseq.gtf
    """
}

process DEU_COUNT {
    tag 'featureCounts_exonic_parts'
    cpus 8
    publishDir "${params.outdir}/deu/counts", mode: 'copy'

    input:
    path bams
    path flat_gtf
    val  lib_type

    output:
    path 'deu_exonic_part_counts.txt',         emit: counts
    path 'deu_exonic_part_counts.txt.summary', emit: summary

    script:
    // rMATS libType -> featureCounts -s (0 unstranded, 1 forward, 2 reverse)
    def strand = lib_type == 'fr-firststrand' ? 2 : (lib_type == 'fr-secondstrand' ? 1 : 0)
    """
    featureCounts -T ${task.cpus} \\
        -p --countReadPairs -B -C \\
        -s ${strand} \\
        -f -O \\
        -F GTF -t exonic_part -g gene_id \\
        -a ${flat_gtf} \\
        -o deu_exonic_part_counts.txt \\
        ${bams}
    """
}

process DEU_DEXSEQ {
    tag 'DEXSeq'
    cpus 8
    publishDir "${params.outdir}/deu/dexseq", mode: 'copy'

    input:
    path counts
    path flat_gtf
    path samplesheet
    path contrasts

    output:
    path 'dexseq/*', emit: results

    script:
    """
    dexseq_deu.R ${counts} ${flat_gtf} ${samplesheet} ${contrasts} \\
        ${params.group_column} ${params.pair_column} ${task.cpus} dexseq
    """
}

process DEU_EDGER {
    tag 'edgeR_diffSpliceDGE'
    publishDir "${params.outdir}/deu/edger", mode: 'copy'

    input:
    path counts
    path samplesheet
    path contrasts

    output:
    path 'edger/*', emit: results

    script:
    """
    edger_deu.R ${counts} ${samplesheet} ${contrasts} \\
        ${params.group_column} ${params.pair_column} edger
    """
}
