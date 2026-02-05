process COMPARE_BUSCO {
    tag "${meta.id}"
    label 'process_low'
    
    publishDir "${params.outdir}/busco_comparison/${meta.id}", mode: 'copy', pattern: '*_busco_comparison.log'

    input:
    tuple val(meta), path(fa_busco), path(br_busco), path(fa_proteins), path(fa_gff), path(br_proteins), path(br_gff)
    path script_dir

    output:
    tuple val(meta), path("${meta.id}_busco_comparison.log")      , emit: log
    tuple val(meta), path("${meta.id}_busco_comparison.txt")      , emit: comparison_table
    tuple val(meta), path("${meta.id}_best_annotation.txt")       , emit: best_choice

    script:
    def sample_id = meta.id
    """
    #!/bin/bash
    set -euo pipefail
    
    echo "Starting BUSCO comparison for ${sample_id}" > ${sample_id}_busco_comparison.log

    if [ ! -f "${script_dir}/compare_busco.py" ]; then
        echo "ERROR: compare_busco.py not found in ${script_dir}" >> ${sample_id}_busco_comparison.log
        exit 1
    fi

    # Log file validation
    for file in ${fa_busco} ${br_busco} ${fa_proteins} ${fa_gff} ${br_proteins} ${br_gff}; do
        if [ -s "\$file" ]; then
            echo "Valid: \$file" >> ${sample_id}_busco_comparison.log
        else
            echo "Empty/Missing: \$file" >> ${sample_id}_busco_comparison.log
        fi
    done

    # Run the comparison script - CAPTURE BOTH STDOUT AND STDERR
    echo "Running compare_busco.py..." >> ${sample_id}_busco_comparison.log
    ${script_dir}/compare_busco.py \
        --funannotate_busco \$(realpath ${fa_busco}) \
        --braker_busco \$(realpath ${br_busco}) \
        --funannotate_proteins \$(realpath ${fa_proteins}) \
        --funannotate_gff \$(realpath ${fa_gff}) \
        --braker_proteins \$(realpath ${br_proteins}) \
        --braker_gff \$(realpath ${br_gff}) \
        >> ${sample_id}_busco_comparison.log 2>&1

    # Check if comparison was created
    if [ ! -f busco_comparison.txt ]; then
        echo "ERROR: compare_busco.py failed to create busco_comparison.txt" >> ${sample_id}_busco_comparison.log
        exit 1
    fi

    mv busco_comparison.txt ${sample_id}_busco_comparison.txt
    echo "Created comparison file: ${sample_id}_busco_comparison.txt" >> ${sample_id}_busco_comparison.log

    # Parse the results - get what the python script actually says
    echo "Parsing results from python output..." >> ${sample_id}_busco_comparison.log
    
    # Get the selected tool from python output
    if ! grep -q "Selected tool:" ${sample_id}_busco_comparison.log; then
        echo "ERROR: Python script didn't select a tool!" >> ${sample_id}_busco_comparison.log
        exit 1
    fi
    
    best_tool=\$(grep "Selected tool:" ${sample_id}_busco_comparison.log | awk '{print \$NF}' | tr -d '[:space:]')
    echo "Python selected: \${best_tool}" >> ${sample_id}_busco_comparison.log
    
    # Get the actual scores
    if [ "\${best_tool}" = "braker" ]; then
        if grep -q "BRAKER BUSCO Complete:" ${sample_id}_busco_comparison.log; then
            best_score=\$(grep "BRAKER BUSCO Complete:" ${sample_id}_busco_comparison.log | awk '{print \$NF}' | sed 's/%//')
        else
            echo "WARNING: BRAKER score not found" >> ${sample_id}_busco_comparison.log
            best_score="NA"
        fi
    else
        if grep -q "Funannotate BUSCO Complete:" ${sample_id}_busco_comparison.log; then
            best_score=\$(grep "Funannotate BUSCO Complete:" ${sample_id}_busco_comparison.log | awk '{print \$NF}' | sed 's/%//')
        else
            echo "WARNING: Funannotate score not found" >> ${sample_id}_busco_comparison.log
            best_score="NA"
        fi
    fi

    # Write the result
    echo "\${best_tool}" > ${sample_id}_best_annotation.txt
    echo "Selected best annotation: \${best_tool} (Complete BUSCO: \${best_score}%)" >> ${sample_id}_busco_comparison.log
    
    # Show both scores for transparency
    echo "" >> ${sample_id}_busco_comparison.log
    echo "=== ACTUAL BUSCO SCORES ===" >> ${sample_id}_busco_comparison.log
    if grep -q "Funannotate BUSCO Complete:" ${sample_id}_busco_comparison.log; then
        funa_score=\$(grep "Funannotate BUSCO Complete:" ${sample_id}_busco_comparison.log | awk '{print \$NF}')
        echo "Funannotate: \${funa_score}" >> ${sample_id}_busco_comparison.log
    fi
    if grep -q "BRAKER BUSCO Complete:" ${sample_id}_busco_comparison.log; then
        braker_score=\$(grep "BRAKER BUSCO Complete:" ${sample_id}_busco_comparison.log | awk '{print \$NF}')
        echo "BRAKER: \${braker_score}" >> ${sample_id}_busco_comparison.log
    fi
    
    echo "BUSCO comparison completed for ${sample_id}" >> ${sample_id}_busco_comparison.log
    """
}
