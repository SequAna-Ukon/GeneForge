process TRNASCAN_SE {
    tag "$meta.id"
    label 'process_medium'
    container 'docker://abdoallahsharaf/geneforge-trnascan:2.1'

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
    path "${meta.id}_error.log", emit: error
    

    script:
    def prefix = "${meta.id}"
    def cpus = task.cpus
    """
    #!/bin/bash
    set -euo pipefail

    # Initialize empty required output files
    touch ${prefix}_highconf.tbl ${prefix}_trnascan-se.out ${prefix}_trnascan-se.tbl ${prefix}_trnascan-se.log ${prefix}_trna_annotation.gff
    mkdir -p ${prefix}_eukconf

    export PATH=/opt/conda/bin:/opt/conda/envs/trnascan/bin:/usr/local/bin:/usr/bin:/bin:\$PATH


    if ! grep -q '^>' ${fasta} || ! grep -q '[ACGTN]' ${fasta}; then
        echo "Invalid FASTA format: missing headers or sequences" >> ${prefix}_error.log
        exit 1
    fi

    rm -f ${prefix}_trnascan-se.out ${prefix}_trnascan-se.tbl ${prefix}_trnascan-se.log
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

    if [ ! -s ${prefix}_trnascan-se.log ]; then
        echo "tRNAscan-SE log is empty" >> ${prefix}_error.log
        exit 1
    fi

    if [ ! -s ${prefix}_trnascan-se.out ]; then
        echo "tRNAscan-SE output file is empty" >> ${prefix}_error.log
        exit 1
    fi

    EukHighConfidenceFilter \\
        -i ${prefix}_trnascan-se.out \\
        -s ${prefix}_trnascan-se.tbl \\
        -o ${prefix}_eukconf \\
        -p ${prefix}_filt 2>> ${prefix}_error.log

    perl ${script_dir}/filter_highconf_tRNAs.pl \\
        ${prefix}_eukconf/${prefix}_filt.out \\
        ${prefix}_highconf.tbl 2>> ${prefix}_error.log

    perl ${script_dir}/convert_tRNAScanSE_to_gff3.pl --input=${prefix}_highconf.tbl > ${prefix}_trna_annotation.gff 2>> ${prefix}_error.log
    
    """
}
