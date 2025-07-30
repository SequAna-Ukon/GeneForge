process TRNASCAN_SE {
    tag "$meta.id"
    label 'process_medium'
    conda "bioconda::trnascan-se"

    publishDir "${params.outdir}/tRNA_scan", mode: 'copy', pattern: "${meta.id}_trna_annotation.gff"
    publishDir "${params.outdir}/tRNA_scan", mode: 'copy', pattern: "${meta.id}_highconf.tbl"

    input:
    tuple val(meta), path(fasta)
    path script_dir

    output:
    tuple val(meta), path("${meta.id}_highconf.tbl"), emit: highconf
    tuple val(meta), path("${meta.id}_trnascan-se.out"), emit: out
    tuple val(meta), path("${meta.id}_trnascan-se.tbl"), emit: tbl
    tuple val(meta), path("${meta.id}_trnascan-se.log"), emit: log
    tuple val(meta), path("${meta.id}_eukconf"), emit: eukconf
    tuple val(meta), path("${meta.id}_trna_annotation.gff"), emit: gff
    path "versions.yml", emit: versions
    path "${meta.id}_error.log", emit: error

    script:
    def prefix = "${meta.id}"
    def cpus = task.cpus
    """
    # Check write permissions to /tmp/
    touch /tmp/test_write 2>> ${prefix}_error.log && rm /tmp/test_write 2>> ${prefix}_error.log || { echo "Cannot write to /tmp/" > ${prefix}_error.log; exit 1; }

    # Validate FASTA file
    if ! grep -q '^>' ${fasta} || ! grep -q '[ACGTN]' ${fasta}; then
        echo "Invalid FASTA format: missing headers or sequences" > ${prefix}_error.log
        exit 1
    fi

    # Log tRNAscan-SE start
    echo "Starting tRNAscan-SE at \$(date)" >> ${prefix}_error.log

    tRNAscan-SE \\
        -E \\
        -I \\
        -H \\
        --detail \\
        --thread ${cpus} \\
        -o ${prefix}_trnascan-se.out \\
        -f ${prefix}_trnascan-se.tbl \\
        -m ${prefix}_trnascan-se.log \\
        ${fasta} 2>> ${prefix}_error.log

    # Check if tRNAscan-SE produced output
    if [ ! -s ${prefix}_trnascan-se.log ]; then
        echo "tRNAscan-SE log is empty" >> ${prefix}_error.log
        exit 1
    fi

    # Check if tRNAscan-SE output exists
    if [ ! -s ${prefix}_trnascan-se.out ]; then
        echo "tRNAscan-SE output file is empty" >> ${prefix}_error.log
        exit 1
    fi

    EukHighConfidenceFilter \\
        -i ${prefix}_trnascan-se.out \\
        -s ${prefix}_trnascan-se.tbl \\
        -o ${prefix}_eukconf \\
        -p ${prefix}_filt 2>> ${prefix}_error.log

    # Filter high-confidence tRNAs
    perl ${script_dir}/filter_highconf_tRNAs.pl \\
        ${prefix}_eukconf/${prefix}_filt.out \\
        ${prefix}_highconf.tbl 2>> ${prefix}_error.log

    # Convert to GFF3
    perl ${script_dir}/convert_tRNAScanSE_to_gff3.pl --input=${prefix}_highconf.tbl > ${prefix}_trna_annotation.gff 2>> ${prefix}_error.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: \$(tRNAscan-SE -h | grep tRNAscan-SE | awk '{print \$2}')
        perl: \$(perl --version | grep -oP 'v\\d+\\.\\d+\\.\\d+' | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.id}"
    """
    touch ${prefix}_trnascan-se.out
    touch ${prefix}_trnascan-se.tbl
    touch ${prefix}_trnascan-se.log
    mkdir ${prefix}_eukconf
    touch ${prefix}_eukconf/${prefix}_filt.out
    touch ${prefix}_highconf.tbl
    touch ${prefix}_trna_annotation.gff
    touch ${prefix}_error.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tRNAscan-SE: 2.0.12
        perl: 5.32.1
    END_VERSIONS
    """
}
