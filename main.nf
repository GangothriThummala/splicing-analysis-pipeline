#!/usr/bin/env nextflow
// Splicing analysis pipeline:
//   FastQC (raw) -> fastp + STAR -> strandedness and read-length checks -> MultiQC
//   -> rMATS prep (per BAM) -> post (all BAMs, one event set) -> stat (per contrast)
//
//   nextflow run main.nf -params-file datasets/<virus>/params.yaml \
//       --raw_dir <fastq_dir> --star_index <index_dir> --gtf <annotation.gtf> --outdir <outdir>

include { FASTQC       } from './modules/fastqc'
include { TRIM_ALIGN   } from './modules/trim_align'
include { STRANDEDNESS } from './modules/strandedness'
include { MULTIQC      } from './modules/multiqc'
include { READ_LENGTH  } from './modules/read_length'
include { RMATS_DESIGN; RMATS_PREP; RMATS_POST; RMATS_STAT } from './modules/rmats'

workflow {
    ['dataset', 'samplesheet', 'contrasts', 'group_column', 'pair_column', 'adapter',
     'raw_dir', 'star_index', 'gtf', 'outdir'].each { p ->
        if (!params[p]) error "Missing --${p} (see README, 'Running a dataset')"
    }

    // One item per sample: [meta, R1, R2]. meta is the samplesheet row without
    // the file names (sample, animal, species, timepoint, ...), used later to
    // group and pair samples for rMATS.
    samples = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // Paired-end only for now; single-end samples need their own run.
            if (!row.fastq_2) error "${row.sample}: no fastq_2; single-end is not supported yet"
            def meta = row.findAll { k, _v -> !(k in ['fastq_1', 'fastq_2']) }
            tuple(meta,
                  file("${params.raw_dir}/${row.fastq_1}", checkIfExists: true),
                  file("${params.raw_dir}/${row.fastq_2}", checkIfExists: true))
        }

    FASTQC(samples)
    TRIM_ALIGN(samples, params.adapter, file(params.star_index, type: 'dir', checkIfExists: true))

    // Dataset-wide checks; each fails the run if samples disagree.
    STRANDEDNESS(TRIM_ALIGN.out.gene_counts.collect())
    READ_LENGTH(FASTQC.out.reports.flatten().filter { f -> f.name.endsWith('.zip') }.collect())
    lib_type    = STRANDEDNESS.out.lib_type.map { f -> f.text.trim() }
    read_length = READ_LENGTH.out.read_length.map { f -> f.text.trim() }

    MULTIQC(
        FASTQC.out.reports
            .mix(TRIM_ALIGN.out.fastp_json, TRIM_ALIGN.out.star_log, TRIM_ALIGN.out.gene_counts)
            .collect()
    )

    // rMATS
    gtf = file(params.gtf, checkIfExists: true)
    RMATS_DESIGN(file(params.samplesheet), file(params.contrasts, checkIfExists: true))
    RMATS_PREP(TRIM_ALIGN.out.bam, gtf, read_length, lib_type)
    RMATS_POST(RMATS_PREP.out.rmats.collect(), RMATS_PREP.out.outcomes.collect(),
               RMATS_DESIGN.out.post_b1, gtf, read_length, lib_type)
    contrasts = RMATS_DESIGN.out.contrasts
        .splitCsv(header: true, sep: '\t')
        .map { r -> tuple(r.contrast, r.paired, r.g1_indices, r.g2_indices) }
    RMATS_STAT(contrasts, RMATS_POST.out.post)
}
