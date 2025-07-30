process FUNCTIONAL_ANNOTATION {
    tag "${meta.id}"
    label 'process_high'

    conda 'bioconda::funannotate=1.8.17 bioconda::eggnog-mapper=2.1.9'

    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    tuple val(meta), path(protein_fasta), path(gff_file), path(genome_unmasked), path(funanno_db), path(eggnog_db), path(genemark_key), path(genemark_tar), path(phobius_tarball), path(signalp_tarball)

    output:
    tuple val(meta), path("${meta.id}_functional_annotation/annotate_results/*"), emit: dir, optional: true

    script:
    def prefix = meta.id
    def species = meta.species
    def busco_db_fun = meta.busco_db_fun
    def funanno_db = meta.funanno_DB
    def use_setup = (!funanno_db || !file(funanno_db).exists())
    """
    #!/bin/bash
    set -euo pipefail

    # Setup EggNOG database first to define eggnog_db_dir
    eggnog_db_dir=""
    if [ -n "${eggnog_db}" ] && [ -d "${eggnog_db}" ]; then
        echo "Using provided EggNOG database: ${eggnog_db}" > ${prefix}_error.log
        eggnog_db_dir=\$(realpath "${eggnog_db}")
    else
        echo "EggNOG database directory not found at ${eggnog_db}. Downloading EggNOG database..." >> ${prefix}_error.log
        mkdir -p eggnog_db
        if ! download_eggnog_data.py -M -y --data_dir ./eggnog_db 2>> ${prefix}_error.log; then
             echo "ERROR: Failed to download EggNOG database" >> ${prefix}_error.log
       
             exit 1
        fi
        eggnog_db_dir=\$(realpath ./eggnog_db)
    fi


    # Validate inputs
    if [ ! -s "${protein_fasta}" ] || [ ! -s "${gff_file}" ]; then
        echo "ERROR: Protein FASTA or GFF file is missing or empty. Skipping functional annotation." >> ${prefix}_error.log
        exit 1
    fi

    # Setup GeneMark
    mkdir -p gmes
    if [ -s "${genemark_key}" ]; then
        export HOME=\$(pwd)
        gunzip -c "${genemark_key}" > "\$HOME/.gm_key"
    else
        echo "ERROR: GeneMark key missing: ${genemark_key}" >> ${prefix}_error.log
        exit 1
    fi

    if [ -s "${genemark_tar}" ]; then
        tar -xzf ${genemark_tar} -C gmes
        subdir=\$(find gmes -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -n "\$subdir" ]; then
            mv "\$subdir"/* gmes/
            rmdir "\$subdir"
        fi
        chmod +x gmes/gmes_petap.pl
        (cd gmes && perl change_path_in_perl_scripts.pl \$(which perl))
    else
        echo "ERROR: GeneMark tar missing: ${genemark_tar}" >> ${prefix}_error.log
        exit 1
    fi
    export GENEMARK_PATH=\$(realpath gmes)
    export PATH=\$GENEMARK_PATH:\$PATH

    # Setup Phobius
    if [ -s "${phobius_tarball}" ]; then
        tar -zxf ${phobius_tarball};
    else
        echo "WARNING: Phobius tarball missing, skipping Phobius annotation" >> ${prefix}_error.log
        touch phobius.results.txt
    fi

    # SignalP6 setup
    if [ -s "${signalp_tarball}" ]; then
        echo "Installing SignalP from tarball: ${signalp_tarball}" >> ${prefix}_error.log
        mkdir -p signalp6
        if ! tar -xzf "${signalp_tarball}" -C signalp6; then
            echo "ERROR: Failed to extract ${signalp_tarball}" >> ${prefix}_error.log
            exit 1
        fi
        # Create and activate virtual environment
        python3 -m venv signalp_venv
        source signalp_venv/bin/activate
        if ! pip install signalp6/signalp6_fast/signalp-6-package/ >> ${prefix}_error.log 2>&1; then
            echo "ERROR: Failed to install signalp module via pip" >> ${prefix}_error.log
            exit 1
        fi
        if ! pip install 'numpy<2' >> ${prefix}_error.log 2>&1; then
            echo "ERROR: Failed to install numpy<2 via pip" >> ${prefix}_error.log
            exit 1
        fi
        SP_DIR=\$(python3 -c 'import signalp; import os; print(os.path.dirname(signalp.__file__))' 2>> ${prefix}_error.log || echo "")
        if [ -n "\$SP_DIR" ]; then
            echo "SignalP module found at: \$SP_DIR" >> ${prefix}_error.log
            if ! cp -r signalp6/signalp6_fast/signalp-6-package/models/* "\$SP_DIR/model_weights/" >> ${prefix}_error.log 2>&1; then
                echo "WARNING: Failed to copy SignalP model weights to \$SP_DIR/model_weights/" >> ${prefix}_error.log
            fi
        else
            echo "ERROR: SignalP module not found after installation" >> ${prefix}_error.log
            exit 1
        fi
        # Run SignalP
        if ! signalp6 --output_dir ${prefix}_signalp -org euk --mode fast -format txt -fasta "${protein_fasta}" --write_procs ${task.cpus} 2>> ${prefix}_error.log; then
            echo "WARNING: SignalP6 execution failed" >> ${prefix}_error.log
        fi
        deactivate
    else
        echo "WARNING: SignalP tarball missing, skipping SignalP annotation" >> ${prefix}_error.log
    fi

    # funannotate_db setup
    if ${use_setup}; then
        mkdir -p ./funannotate_db
        funannotate setup --install all -b ${busco_db_fun} -f --database ./funannotate_db
        export FUNANNOTATE_DB=\$(realpath ./funannotate_db)
    else
        export FUNANNOTATE_DB="${funanno_db}"
    fi

    # Fix permissions for aux_scripts
    SITE_PACKAGES=\$(python3 -c 'import site; print(site.getsitepackages()[0])' 2>/dev/null || echo "")
    echo "SITE_PACKAGES: \$SITE_PACKAGES" >> ${prefix}_error.log
    if [ -n "\$SITE_PACKAGES" ]; then
        AUX_SCRIPTS_DIR="\$SITE_PACKAGES/funannotate/aux_scripts"
        echo "AUX_SCRIPTS_DIR: \$AUX_SCRIPTS_DIR" >> ${prefix}_error.log
        if [ -d "\$AUX_SCRIPTS_DIR" ]; then
            chmod -R +x "\$AUX_SCRIPTS_DIR" 2>> ${prefix}_error.log || {
                echo "WARNING: Failed to set executable permissions for \$AUX_SCRIPTS_DIR" >> ${prefix}_error.log
            }
        else
            echo "WARNING: Directory \$AUX_SCRIPTS_DIR not found" >> ${prefix}_error.log
        fi
    else
        echo "WARNING: Could not determine site-packages directory for Python" >> ${prefix}_error.log
    fi


    # Run Phobius
    if [ -d "phobius" ]; then
        phobius/phobius.pl -short "${protein_fasta}" > phobius.results.txt 2>> ${prefix}_error.log
    else
        echo "WARNING: Phobius directory not found, skipping Phobius annotation" >> ${prefix}_error.log
    fi

    # Run InterProScan
    funannotate iprscan -i "${protein_fasta}" -m docker -c ${task.cpus} -o ${prefix}_iprscan.xml 2>> ${prefix}_error.log

    # Run EggNOG
    emapper.py --cpu ${task.cpus} -m mmseqs --data_dir "\$eggnog_db_dir" -i "${protein_fasta}" -o ${prefix}_eggnog 2>> ${prefix}_error.log

    # Run Funannotate annotate
    funannotate annotate \\
        --gff "${gff_file}" \\
        --fasta "${genome_unmasked}" \\
        --species "${species}" \\
        --busco_db ${busco_db_fun} \\
        --eggnog ${prefix}_eggnog.emapper.annotations \\
        --iprscan ${prefix}_iprscan.xml \\
        --phobius phobius.results.txt \\
        --signalp ${prefix}_signalp/prediction_results.txt \\
        --cpus ${task.cpus} \\
        -o ${prefix}_functional_annotation 2>> ${prefix}_error.log
    """
}
