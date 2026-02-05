process MERGE_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_high'
    
    // Conda with required tools (same as BRAKER_POST + busco)
    conda "bioconda::agat=1.4.0 bioconda::gffread=0.12.7 bioconda::busco=5.4.7"

    // Publish only the final key files, named by sample ID
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*_GeneForge.gff3'
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*.GeneForge.prot.fasta'
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*_busco_geneforge.txt'

    input:
    tuple val(meta),
          path(best_choice_file, stageAs: 'best_choice_file.txt'),
          path(fa_gff),
          path(fa_proteins),
          path(br_gff),
          path(br_proteins),
          path(trna_gff_file),
          path(genome)

    output:
    tuple val(meta), path("${meta.id}_GeneForge.gff3")         , emit: final_gff
    tuple val(meta), path("${meta.id}.GeneForge.prot.fasta")   , emit: final_proteins
    tuple val(meta), path("${meta.id}_busco_geneforge.txt")    , emit: busco_summary
    path "${meta.id}_merge_annotations.log"                    , emit: log

    script:
    def prefix     = meta.id
    def busco_db   = meta.busco_db
    
    """
    #!/bin/bash
    set -euo pipefail

    echo "Starting GeneForge merging for ${prefix}" > ${prefix}_merge_annotations.log
    
    # Read the best tool FROM BASH, not from Groovy
    # The file is staged as 'best_choice_file.txt' 
    best_tool=\$(cat "best_choice_file.txt" | tr -d '[:space:]')
    echo "Selected best annotation (via BUSCO): \${best_tool}" >> ${prefix}_merge_annotations.log

    # Determine reference and alternative GFF - USE BASH VARIABLE \$best_tool
    if [[ "\${best_tool}" == "braker" ]]; then
        REF_GFF=${br_gff}
        ALT_GFF=${fa_gff}
        echo "Using BRAKER as reference" >> ${prefix}_merge_annotations.log
    elif [[ "\${best_tool}" == "funannotate" ]]; then
        REF_GFF=${fa_gff}
        ALT_GFF=${br_gff}
        echo "Using Funannotate as reference" >> ${prefix}_merge_annotations.log
    else
        echo "ERROR: Invalid content in best_annotation.txt: \${best_tool}" >> ${prefix}_merge_annotations.log
        exit 1
    fi

    # 1. Complement reference with alternative
    agat_sp_complement_annotations.pl --ref \$REF_GFF --add \$ALT_GFF -o complemented.gff3 \
        >> ${prefix}_merge_annotations.log 2>&1

    # 2. Fix duplicated/overlapping features
    agat_sp_fix_features_locations_duplicated.pl --gff complemented.gff3 -o cleaned.gff \
        >> ${prefix}_merge_annotations.log 2>&1

    # 3. Remove tRNAs from protein-coding annotation
    grep -v \$'\ttRNA\t' cleaned.gff > cleaned_no_trna.gff

    # 4. Merge independent tRNA annotation
    agat_sp_merge_annotations.pl --gff cleaned_no_trna.gff --gff ${trna_gff_file} --out with_trna.gff \
        >> ${prefix}_merge_annotations.log 2>&1

    # 5. Remove any residual funannotate tRNAs (safety)
    grep -v "funannotate"\$'\tRNA\t' with_trna.gff | grep -v "funannotate"\$'\ttRNA\t' > ${prefix}_GeneForge.gff3

    # 6. Extract protein sequences
    gffread ${prefix}_GeneForge.gff3 -g ${genome} -y ${prefix}.GeneForge.prot.fasta \
        >> ${prefix}_merge_annotations.log 2>&1

    # 7. Run BUSCO on final proteins
    echo "Running BUSCO on final GeneForge proteins..." >> ${prefix}_merge_annotations.log
    busco -i ${prefix}.GeneForge.prot.fasta \
          -o ${prefix}_busco \
          -m proteins \
          -l ${busco_db} \
          -c ${task.cpus} \
          >> ${prefix}_merge_annotations.log 2>&1

    # Extract and rename the short summary
    if [[ -f "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ]]; then
        mv "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_geneforge.txt
    elif [[ -f "${prefix}_busco/short_summary.generic.${busco_db}.${prefix}_busco.txt" ]]; then
        mv "${prefix}_busco/short_summary.generic.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_geneforge.txt
    else
        echo "WARNING: BUSCO summary not found!" >> ${prefix}_merge_annotations.log
        touch ${prefix}_busco_geneforge.txt  # empty placeholder
    fi

    echo "GeneForge pipeline completed for ${prefix}" >> ${prefix}_merge_annotations.log
    """
}
