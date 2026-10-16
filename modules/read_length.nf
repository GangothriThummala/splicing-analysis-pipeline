// Raw read length for rMATS --readLength, from the FastQC reports.
// Fails if files disagree (bin/read_length.py).
process READ_LENGTH {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path fastqc_zips

    output:
    path 'read_length.tsv'
    path 'read_length.txt', emit: read_length

    script:
    """
    read_length.py read_length.tsv ${fastqc_zips} > read_length.txt
    """
}
