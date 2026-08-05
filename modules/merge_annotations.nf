process MERGE_ANNOTATIONS {
    tag "${meta.id}"
    label 'process_high'
    
    container 'docker://abdoallahsharaf/geneforge-braker3:2.1'

    // Publish final evaluated files
    publishDir "${params.outdir}/busco_comparison/${meta.id}", mode: 'copy', pattern: '*_busco_comparison.log'
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*_GeneForge.gff3'
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*.GeneForge.prot.fasta'
    publishDir "${params.outdir}/geneforge", mode: 'copy', pattern: '*_busco_geneforge.txt'

    input:
    tuple val(meta), 
          path(fa_busco), 
          path(br_busco), 
          path(fa_proteins), 
          path(fa_gff), 
          path(br_proteins), 
          path(br_gff),
          path(trna_gff_file),
          path(genome)
    path script_dir

    output:
    tuple val(meta), path("${meta.id}_GeneForge.gff3")         , emit: final_gff
    tuple val(meta), path("${meta.id}.GeneForge.prot.fasta")   , emit: final_proteins
    tuple val(meta), path("${meta.id}_busco_geneforge.txt")    , emit: busco_summary
    path "${meta.id}_busco_comparison.log"                     , emit: comparison_log
    path "${meta.id}_merge_annotations.log"                    , emit: merge_log

    script:
    def prefix   = meta.id
    def busco_db = meta.busco_db
    
    """
    #!/bin/bash
    set -euo pipefail

    # Point to your new environment paths
    export PATH=/opt/conda/envs/agat_busco/bin:/usr/local/bin:/usr/bin:/bin:$PATH

    # ===============================================================
    # PHASE 1: BUSCO COMPARISON
    # ===============================================================
    echo "Starting BUSCO comparison for ${prefix}" > ${prefix}_busco_comparison.log

    if [ ! -f "${script_dir}/compare_busco.py" ]; then
        echo "ERROR: compare_busco.py not found in ${script_dir}" >> ${prefix}_busco_comparison.log
        exit 1
    fi

    # Log file validation
    for file in ${fa_busco} ${br_busco} ${fa_proteins} ${fa_gff} ${br_proteins} ${br_gff}; do
        if [ -s "\$file" ]; then
            echo "Valid: \$file" >> ${prefix}_busco_comparison.log
        else
            echo "Empty/Missing: \$file" >> ${prefix}_busco_comparison.log
        fi
    done

    echo "Running compare_busco.py..." >> ${prefix}_busco_comparison.log
    python3 ${script_dir}/compare_busco.py \\
        --funannotate_busco \$(realpath ${fa_busco}) \\
        --braker_busco \$(realpath ${br_busco}) \\
        --funannotate_proteins \$(realpath ${fa_proteins}) \\
        --funannotate_gff \$(realpath ${fa_gff}) \\
        --braker_proteins \$(realpath ${br_proteins}) \\
        --braker_gff \$(realpath ${br_gff}) \\
        >> ${prefix}_busco_comparison.log 2>&1

    if [ ! -f busco_comparison.txt ]; then
        echo "ERROR: compare_busco.py failed to create busco_comparison.txt" >> ${prefix}_busco_comparison.log
        exit 1
    fi
    mv busco_comparison.txt ${prefix}_busco_comparison.txt

    if ! grep -q "Selected tool:" ${prefix}_busco_comparison.log; then
        echo "ERROR: Python script didn't select a tool!" >> ${prefix}_busco_comparison.log
        exit 1
    fi
    
    best_tool=\$(grep "Selected tool:" ${prefix}_busco_comparison.log | awk '{print \$NF}' | tr -d '[:space:]')
    echo "Python selected: \${best_tool}" >> ${prefix}_busco_comparison.log

    # ===============================================================
    # PHASE 2: STRUCTURAL ANNOTATION MERGING
    # ===============================================================
    echo "Starting GeneForge merging for ${prefix}" > ${prefix}_merge_annotations.log
    echo "Selected best annotation (via BUSCO): \${best_tool}" >> ${prefix}_merge_annotations.log
    
    # Wipe out pre-existing report files 
    rm -f cleaned_report.txt complemented_report.txt

    if [[ "\${best_tool}" == "braker" ]]; then
        REF_GFF=${br_gff}
        ALT_GFF=${fa_gff}
        echo "Using BRAKER as reference" >> ${prefix}_merge_annotations.log
    elif [[ "\${best_tool}" == "funannotate" ]]; then
        REF_GFF=${fa_gff}
        ALT_GFF=${br_gff}
        echo "Using Funannotate as reference" >> ${prefix}_merge_annotations.log
    else
        echo "ERROR: Invalid tool configuration resolved: \${best_tool}" >> ${prefix}_merge_annotations.log
        exit 1
    fi

    # 1. Complement reference with alternative
    agat_sp_complement_annotations.pl --ref \$REF_GFF --add \$ALT_GFF -o complemented.gff3 \\
        >> ${prefix}_merge_annotations.log 2>&1

    # 2. Fix duplicated/overlapping features
    agat_sp_fix_features_locations_duplicated.pl --gff complemented.gff3 -o cleaned.gff \\
        >> ${prefix}_merge_annotations.log 2>&1

    # 3. Remove tRNAs from protein-coding annotation
    grep -v \$'\ttRNA\t' cleaned.gff > cleaned_no_trna.gff

    # 4. Merge independent tRNA annotation
    agat_sp_merge_annotations.pl --gff cleaned_no_trna.gff --gff ${trna_gff_file} --out with_trna.gff \\
        >> ${prefix}_merge_annotations.log 2>&1

    # 5. Remove residual funannotate tRNAs and fix ncRNA designations
    grep -v "funannotate"\$'\tRNA\t' with_trna.gff | grep -v "funannotate"\$'\ttRNA\t' > ${prefix}_GeneForge.gff3
    sed -i 's/\tRNA\t/\tncRNA\t/g' ${prefix}_GeneForge.gff3
    sed -i 's/;anticodon=[^;]*//g' ${prefix}_GeneForge.gff3

    # 6. Extract protein sequences
    gffread ${prefix}_GeneForge.gff3 -g ${genome} -y ${prefix}.GeneForge.prot.fasta \\
        >> ${prefix}_merge_annotations.log 2>&1

    # 6b. Sanitize translated proteins
    
    n_dots=\$(grep -v '^>' ${prefix}.GeneForge.prot.fasta | grep -o '\\.' | wc -l || true)
    if [ "\$n_dots" -gt 0 ]; then
        echo "WARNING: \$n_dots untranslatable dot residues found; replacing with X" >> ${prefix}_merge_annotations.log
        sed -i '/^>/!s/\\./X/g' ${prefix}.GeneForge.prot.fasta
    fi


    # 7. Run BUSCO on final proteins
    echo "Running BUSCO on final GeneForge proteins..." >> ${prefix}_merge_annotations.log
    /opt/conda/envs/busco6_env/bin/python3 /opt/conda/envs/busco6_env/bin/busco \\
          -i ${prefix}.GeneForge.prot.fasta \\
          -o ${prefix}_busco \\
          -m proteins \\
          -l ${busco_db} \\
          -c ${task.cpus} \\
          >> ${prefix}_merge_annotations.log 2>&1

    if [[ -f "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ]]; then
        mv "${prefix}_busco/short_summary.specific.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_geneforge.txt
    elif [[ -f "${prefix}_busco/short_summary.generic.${busco_db}.${prefix}_busco.txt" ]]; then
        mv "${prefix}_busco/short_summary.generic.${busco_db}.${prefix}_busco.txt" ${prefix}_busco_geneforge.txt
    else
        echo "WARNING: BUSCO summary not found!" >> ${prefix}_merge_annotations.log
        touch ${prefix}_busco_geneforge.txt
    fi

    echo "GeneForge pipeline completed for ${prefix}" >> ${prefix}_merge_annotations.log
    """
}

