process BRAKER_POST {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*_braker.gff3'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*.braker.prot.fasta'
    publishDir "${params.outdir}/braker", mode: 'copy', pattern: '*_busco_braker.txt'

    conda "bioconda::agat=1.4.0 bioconda::gffread=0.12.7 bioconda::busco=5.4.7"

    input:
    tuple val(meta),  // Metadata is first element
          path(braker_gtf),
          path(braker_dir),
          path(genome_unmasked),
          path(trna_gff),
          path(script_dir)

    output:
    tuple val(meta), path("${meta.id}_braker.gff3"), emit: gff3
    tuple val(meta), path("${meta.id}.braker.prot.fasta"), emit: proteins
    tuple val(meta), path("${meta.id}_busco_braker.txt"), emit: busco_summary

    script:
    def prefix = meta.id
    def busco_db = meta.busco_db
    

    return """
    #!/bin/bash
    set -euo pipefail
    
    if [[ -f "${braker_dir}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff" ]]; then
        echo ">>> Adding UTRs..."
        python3.8 ${script_dir}/stringtie2utr.py \\
            -g ${braker_gtf} \\
            -s ${braker_dir}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff \\
            -o braker_with_utrs.gtf
        cat braker_with_utrs.gtf | gtf2gff.pl --gff3 -o braker.gff3    
    else
        echo ">>> no RNASeq..."
        cat ${braker_gtf} | gtf2gff.pl --gff3 -o braker.gff3
    fi
  
    echo ">>> Merging tRNAs..."
    agat_sp_merge_annotations.pl --gff braker.gff3 --gff ${trna_gff} --out merged.gff

    echo ">>> Post-filtering..."
    agat_sp_filter_by_ORF_size.pl -g merged.gff -s 50 -o ${prefix}_filtered.gff
    agat_sp_fix_overlaping_genes.pl -f ${prefix}_filtered_sup50.gff -o ${prefix}_braker.gff3

    echo ">>> Extracting proteins..."
    gffread ${prefix}_braker.gff3 -g ${genome_unmasked} -y ${prefix}.braker.prot.fasta

    echo ">>> Running BUSCO..."
    busco -i ${prefix}.braker.prot.fasta -o ${prefix}_busco -m proteins -l ${busco_db} -c ${task.cpus}
    
    # Move final outputs
    mv "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_braker.txt
    """
}
