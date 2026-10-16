#!/usr/bin/env nextflow
// Splicing analysis pipeline:
//   FastQC (raw) -> fastp + STAR -> strandedness and read-length checks -> MultiQC
//   -> rMATS prep (per BAM) -> post (all BAMs, one event set) -> stat (per contrast)
//   -> filter (coverage and PSI, FDR recomputed)
//
//   nextflow run main.nf -params-file datasets/<virus>/params.yaml \
//       --raw_dir <fastq_dir> --star_index <index_dir> --gtf <annotation.gtf> --outdir <outdir>
//
// --aligned_dir <dir> starts from a finished alignment instead (the --outdir of an
// earlier run, or of bin/align_all.sh): BAMs, gene counts, FastQC and fastp
// reports are read from there, and --raw_dir/--star_index are not needed.

include { FASTQC       } from './modules/fastqc'
include { TRIM_ALIGN   } from './modules/trim_align'
include { STRANDEDNESS } from './modules/strandedness'
include { MULTIQC      } from './modules/multiqc'
include { READ_LENGTH  } from './modules/read_length'
include { RMATS_DESIGN; RMATS_PREP; RMATS_POST; RMATS_STAT; RMATS_FILTER } from './modules/rmats'

workflow {
    def required = ['dataset', 'samplesheet', 'contrasts', 'group_column', 'pair_column', 'gtf', 'outdir']
    if (!params.aligned_dir) required += ['adapter', 'raw_dir', 'star_index']
    required.each { p ->
        if (!params[p]) error "Missing --${p} (see README, 'Running a dataset')"
    }

    // One item per samplesheet row: [meta, row]. meta is the row without the
    // file names (sample, animal, species, timepoint, ...), used later to group
    // and pair samples for rMATS.
    rows = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // Paired-end only for now; single-end samples need their own run.
            if (!row.fastq_2) error "${row.sample}: no fastq_2; single-end is not supported yet"
            tuple(row.findAll { k, _v -> !(k in ['fastq_1', 'fastq_2']) }, row)
        }

    // Either align now, or take the outputs of an earlier alignment. Both give
    // the same channels; everything after this block is shared.
    if (!params.aligned_dir) {
        samples = rows.map { meta, row ->
            tuple(meta,
                  file("${params.raw_dir}/${row.fastq_1}", checkIfExists: true),
                  file("${params.raw_dir}/${row.fastq_2}", checkIfExists: true))
        }
        FASTQC(samples)
        TRIM_ALIGN(samples, params.adapter, file(params.star_index, type: 'dir', checkIfExists: true))
        bams        = TRIM_ALIGN.out.bam
        gene_counts = TRIM_ALIGN.out.gene_counts
        star_logs   = TRIM_ALIGN.out.star_log
        fastp_json  = TRIM_ALIGN.out.fastp_json
        fastqc_zips = FASTQC.out.reports.flatten().filter { f -> f.name.endsWith('.zip') }
    } else {
        def a = "${params.aligned_dir}/align"
        bams = rows.map { meta, _row ->
            def s = meta.sample
            tuple(meta,
                  file("${a}/star/${s}/${s}.Aligned.sortedByCoord.out.bam", checkIfExists: true),
                  file("${a}/star/${s}/${s}.Aligned.sortedByCoord.out.bam.bai", checkIfExists: true))
        }
        gene_counts = rows.map { meta, _row ->
            file("${a}/star/${meta.sample}/${meta.sample}.ReadsPerGene.out.tab", checkIfExists: true) }
        star_logs   = rows.map { meta, _row ->
            file("${a}/star/${meta.sample}/${meta.sample}.Log.final.out", checkIfExists: true) }
        fastp_json  = rows.map { meta, _row ->
            file("${a}/fastp/${meta.sample}.fastp.json", checkIfExists: true) }
        fastqc_zips = rows.flatMap { _meta, row ->
            [row.fastq_1, row.fastq_2].collect { fq ->
                file("${a}/fastqc_raw/${fq.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')}_fastqc.zip", checkIfExists: true)
            }
        }
    }

    // Dataset-wide checks; each fails the run if samples disagree.
    STRANDEDNESS(gene_counts.collect())
    READ_LENGTH(fastqc_zips.collect())
    lib_type    = STRANDEDNESS.out.lib_type.map { f -> f.text.trim() }
    read_length = READ_LENGTH.out.read_length.map { f -> f.text.trim() }

    MULTIQC(fastqc_zips.mix(fastp_json, star_logs, gene_counts).collect())

    // rMATS
    gtf = file(params.gtf, checkIfExists: true)
    RMATS_DESIGN(file(params.samplesheet), file(params.contrasts, checkIfExists: true))
    RMATS_PREP(bams, gtf, read_length, lib_type)
    RMATS_POST(RMATS_PREP.out.rmats.collect(), RMATS_PREP.out.outcomes.collect(),
               RMATS_DESIGN.out.post_b1, gtf, read_length, lib_type)
    contrasts = RMATS_DESIGN.out.contrasts
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.contrast, r.paired, r.g1_indices, r.g2_indices) }
    RMATS_STAT(contrasts, RMATS_POST.out.post)
    RMATS_FILTER(RMATS_STAT.out)
}
