process FUNCTIONAL_ANNOTATION {
    tag "$meta.id"
    label 'process_high'

    conda 'bioconda::funannotate=1.8.17 bioconda::eggnog-mapper=2.1.9'
    publishDir "${params.outdir}/functional_annotation", mode: 'copy'

    input:
    tuple val(meta), path(proteins), path(gff), path(genome), path(funannotate_db), path(eggnog_db), path(gm_key), path(gmes_tar), path(phobius_tar), path(signalp_tar)

    output:
    path("${meta.id}_functional_annotation"), emit: annotation_dir
    path("${meta.id}_error.log"), emit: log

    script:
    """
    #!/bin/bash
    set -euo pipefail
    
    # -------------------------------------------------------------------------
    # 1. EggNOG database path resolution & installation
    # -------------------------------------------------------------------------
    # Your original exact logic—resolves the Nextflow input symlink to the real path
    export EGGNOG_DB_DIR="\$(realpath ${eggnog_db})"
    
    if [[ ! -d "\$EGGNOG_DB_DIR" ]]; then
        mkdir -p "\$EGGNOG_DB_DIR"
    fi

    # Foolproof file-level check inside the real path
    if [[ ! -f "\$EGGNOG_DB_DIR/eggnog.db" || ! -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
        (
            flock -x 201
            if [[ ! -f "\$EGGNOG_DB_DIR/eggnog.db" || ! -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
                echo ">>> Database missing. Downloading directly into: \$EGGNOG_DB_DIR" >> ${meta.id}_error.log
                
                BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
                
                # Download straight into the real resolved path using -P
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog.db.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog_proteins.dmnd.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog.taxa.tar.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/mmseqs.tar.gz"
                echo ">>> Extracting archives directly within cache path..." >> ${meta.id}_error.log
                
                # Decompress inside the real destination folder
                gunzip -f "\$EGGNOG_DB_DIR/eggnog.db.gz"
                gunzip -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd.gz"
                
                tar -xzf "\$EGGNOG_DB_DIR/eggnog.taxa.tar.gz" -C "\$EGGNOG_DB_DIR"
                rm -f "\$EGGNOG_DB_DIR/eggnog.taxa.tar.gz"

                tar -xzf "\$EGGNOG_DB_DIR/mmseqs.tar.gz" -C "\$EGGNOG_DB_DIR"
                rm -f "\$EGGNOG_DB_DIR/mmseqs.tar.gz"

                echo ">>> Manual EggNOG database installation complete." >> ${meta.id}_error.log
            else
                echo ">>> Reusing existing EggNOG database (populated by parallel process)" >> ${meta.id}_error.log
            fi
        ) 201>"\$EGGNOG_DB_DIR/.download.lock"
    else
        echo ">>> Verified existing EggNOG database files at \$EGGNOG_DB_DIR" >> ${meta.id}_error.log
    fi
    
    # -------------------------------------------------------------------------
    # 2. Validate inputs
    # -------------------------------------------------------------------------
    if [ ! -s "${proteins}" ] || [ ! -s "${gff}" ]; then
        echo "ERROR: Protein FASTA or GFF file is missing or empty. Skipping functional annotation." >> ${meta.id}_error.log
        exit 1
    fi
    
    # -------------------------------------------------------------------------
    # 3. GeneMark setup
    # -------------------------------------------------------------------------
    mkdir -p gmes
    if [ -s "${gm_key}" ] && [ "${gm_key.name}" != "NO_FILE" ]; then
        export HOME=\$(pwd)
        gunzip -c "${gm_key}" > "\$HOME/.gm_key"
    else
        echo "ERROR: GeneMark key missing." >> ${meta.id}_error.log
        exit 1
    fi
    
    if [ -s "${gmes_tar}" ] && [ "${gmes_tar.name}" != "NO_FILE" ]; then
        tar -xzf "${gmes_tar}" -C gmes
        subdir=\$(find gmes -mindepth 1 -maxdepth 1 -type d | head -n1)
        if [ -n "\$subdir" ]; then
            mv "\$subdir"/* gmes/
            rmdir "\$subdir"
        fi
        chmod +x gmes/gmes_petap.pl
        (cd gmes && perl change_path_in_perl_scripts.pl \$(which perl))
    else
        echo "ERROR: GeneMark tarball missing." >> ${meta.id}_error.log
        exit 1
    fi
    export GENEMARK_PATH=\$(realpath gmes)
    export PATH=\$GENEMARK_PATH:\$PATH
    
    # -------------------------------------------------------------------------
    # 4. Phobius setup
    # -------------------------------------------------------------------------
    if [ -s "${phobius_tar}" ] && [ "${phobius_tar.name}" != "NO_PHOBIUS_TARBALL.empty" ]; then
        tar -zxf "${phobius_tar}"
    else
        echo "WARNING: Phobius tarball missing, skipping Phobius annotation" >> ${meta.id}_error.log
        touch phobius.results.txt
    fi
    
    # -------------------------------------------------------------------------
    # 5. SignalP6 setup
    # -------------------------------------------------------------------------
    if [ -s "${signalp_tar}" ] && [ "${signalp_tar.name}" != "NO_SIGNALP_TARBALL.empty" ]; then
        echo "Installing SignalP from tarball..." >> ${meta.id}_error.log
        mkdir -p signalp6
        tar -xzf "${signalp_tar}" -C signalp6
        python3 -m venv signalp_venv
        source signalp_venv/bin/activate
        pip install signalp6/signalp6_fast/signalp-6-package/ >> ${meta.id}_error.log 2>&1
        pip install 'numpy<2' >> ${meta.id}_error.log 2>&1
        
        SP_DIR=\$(python3 -c 'import signalp; import os; print(os.path.dirname(signalp.__file__))' 2>> ${meta.id}_error.log || echo "")
        if [ -n "\$SP_DIR" ]; then
            cp -r signalp6/signalp6_fast/signalp-6-package/models/* "\$SP_DIR/model_weights/" >> ${meta.id}_error.log 2>&1
        fi
        
        if ! signalp6 --output_dir ${meta.id}_signalp -org euk --mode fast -format txt -fasta "${proteins}" --write_procs ${task.cpus} 2>> ${meta.id}_error.log; then
            echo "WARNING: SignalP6 execution failed" >> ${meta.id}_error.log
        fi
        deactivate
    else
        echo "WARNING: SignalP tarball missing, skipping SignalP annotation" >> ${meta.id}_error.log
    fi
    
    # -------------------------------------------------------------------------
    # 6. Funannotate database setup
    # -------------------------------------------------------------------------
    export FUNANNOTATE_DB="\$(realpath ${funannotate_db})"
    
    if [[ ! -d "\$FUNANNOTATE_DB" ]]; then
        if [[ "\$FUNANNOTATE_DB" == *"/work/funannotate_db" ]]; then
            mkdir -p "\$FUNANNOTATE_DB"
        else
            echo "ERROR: Provided DB path \$FUNANNOTATE_DB does not exist." >> ${meta.id}_error.log
            echo "       Please run: funannotate setup --install all -b ${meta.busco_db_fun ?: 'metazoa'} --wget -f --database \$FUNANNOTATE_DB" >> ${meta.id}_error.log
            exit 1
        fi
    fi

    if [[ -f "\$FUNANNOTATE_DB/funannotate-db-info.txt" ]]; then
        echo ">>> Using existing funannotate database at \$FUNANNOTATE_DB" >> ${meta.id}_error.log
    else
        (
            flock -x 200
            if [[ ! -f "\$FUNANNOTATE_DB/funannotate-db-info.txt" ]]; then
                echo ">>> Installing funannotate database to \$FUNANNOTATE_DB" >> ${meta.id}_error.log
                funannotate setup --install all -b ${meta.busco_db_fun ?: 'metazoa'} --wget -f --database "\$FUNANNOTATE_DB" >> ${meta.id}_error.log 2>&1
            else
                echo ">>> Reusing existing funannotate database at \$FUNANNOTATE_DB" >> ${meta.id}_error.log
            fi
        ) 200>"\$FUNANNOTATE_DB/.install.lock"
    fi
    
    # -------------------------------------------------------------------------
    # 7. Fix permissions for aux_scripts
    # -------------------------------------------------------------------------
    SITE_PACKAGES=\$(python3 -c 'import site; print(site.getsitepackages()[0])' 2>/dev/null || echo "")
    if [ -n "\$SITE_PACKAGES" ]; then
        AUX_SCRIPTS_DIR="\$SITE_PACKAGES/funannotate/aux_scripts"
        if [ -d "\$AUX_SCRIPTS_DIR" ]; then
            chmod -R +x "\$AUX_SCRIPTS_DIR" 2>> ${meta.id}_error.log || true
        fi
    fi
    
    # -------------------------------------------------------------------------
    # 8. Phobius run execution
    # -------------------------------------------------------------------------
    if [ -d "phobius" ]; then
        phobius/phobius.pl -short "${proteins}" > phobius.results.txt 2>> ${meta.id}_error.log
    fi
    
    # -------------------------------------------------------------------------
    # 9. InterProScan
    # -------------------------------------------------------------------------
    funannotate iprscan -i "${proteins}" -m docker -c ${task.cpus} -o ${meta.id}_iprscan.xml 2>> ${meta.id}_error.log
    
    # -------------------------------------------------------------------------
    # 10. EggNOG Map
    # -------------------------------------------------------------------------
    emapper.py --cpu ${task.cpus} -m mmseqs --data_dir "\$EGGNOG_DB_DIR" -i "${proteins}" -o ${meta.id}_eggnog 2>> ${meta.id}_error.log
    
    # -------------------------------------------------------------------------
    # 11. Funannotate final synthesis
    # -------------------------------------------------------------------------
    # Build arguments safely into a bash array
    ANNOTATE_ARGS=(
        --gff "${gff}"
        --fasta "${genome}"
        --species "${meta.species}"
        --busco_db "${meta.busco_db_fun ?: 'metazoa'}"
        --eggnog "${meta.id}_eggnog.emapper.annotations"
        --iprscan "${meta.id}_iprscan.xml"
        --cpus "${task.cpus}"
        -o "${meta.id}_functional_annotation"
    )

    if [ -f "${meta.id}_signalp/prediction_results.txt" ]; then
        ANNOTATE_ARGS+=("--signalp" "${meta.id}_signalp/prediction_results.txt")
    fi

    if [ -s "phobius.results.txt" ]; then
        ANNOTATE_ARGS+=("--phobius" "phobius.results.txt")
    fi

    funannotate annotate "\${ANNOTATE_ARGS[@]}" 2>> ${meta.id}_error.log
    """
}
