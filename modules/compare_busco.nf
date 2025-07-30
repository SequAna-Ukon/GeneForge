process COMPARE_BUSCO {
    tag "${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}/busco_comparison", mode: 'copy'

    input:
    tuple val(meta), path(fa_busco), path(br_busco), path(fa_proteins), path(fa_gff), path(br_proteins), path(br_gff)
    path script_dir

    output:
    path "busco_comparison.log", emit: log

    script:
    """
    # Verify script existence
    if [ ! -f "${script_dir}/compare_busco.py" ]; then
        echo "ERROR: compare_busco.py not found in ${script_dir}" > busco_comparison.log
        exit 1
    fi

    # Debug: Check if input files are valid
    for file in ${fa_busco} ${br_busco} ${fa_proteins} ${fa_gff} ${br_proteins} ${br_gff}; do
        if [ -s "\$file" ]; then
            echo "Valid file: \$file" >> busco_comparison.log
        else
            echo "Invalid or empty file: \$file" >> busco_comparison.log
        fi
    done

    # Run comparison with absolute paths
    ${script_dir}/compare_busco.py \
        --funannotate_busco \$(realpath ${fa_busco}) \
        --braker_busco \$(realpath ${br_busco}) \
        --funannotate_proteins \$(realpath ${fa_proteins}) \
        --funannotate_gff \$(realpath ${fa_gff}) \
        --braker_proteins \$(realpath ${br_proteins}) \
        --braker_gff \$(realpath ${br_gff}) \
        2>> busco_comparison.log

    # Check if busco_comparison.txt was created
    if [ ! -f "busco_comparison.txt" ]; then
        echo "ERROR: busco_comparison.txt not created" >> busco_comparison.log
        exit 1
    fi
    """
}
