// Test driver for the DEU branch: starts from existing BAMs (--aligned_dir); no alignment, no rMATS.
// Optional: --samples a,b,c to use a subset; --lib_type (default fr-unstranded).
include { DEU_ANNOTATION; DEU_COUNT; DEU_DEXSEQ; DEU_EDGER } from './modules/deu'

workflow {
    def wanted   = params.samples ? params.samples.toString().tokenize(',') : null
    def lib_type = params.lib_type ?: 'fr-unstranded'

    bams = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .filter { row -> wanted == null || row.sample in wanted }
        .map { row -> file("${params.aligned_dir}/align/star/${row.sample}/${row.sample}.Aligned.sortedByCoord.out.bam", checkIfExists: true) }
        .collect()

    DEU_ANNOTATION(file(params.gtf, checkIfExists: true))
    DEU_COUNT(bams, DEU_ANNOTATION.out.flat_gtf, lib_type)
    DEU_DEXSEQ(DEU_COUNT.out.counts, DEU_ANNOTATION.out.flat_gtf,
               file(params.samplesheet), file(params.contrasts))
    DEU_EDGER(DEU_COUNT.out.counts, file(params.samplesheet), file(params.contrasts))
}
