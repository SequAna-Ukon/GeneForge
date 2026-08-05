process FUNCTIONAL_ANNOTATION {
    tag "$meta.id"
    label 'process_high'

    container 'docker://abdoallahsharaf/geneforge-funannotate-func:2.1'
    publishDir "${params.outdir}/functional_annotation/${meta.id}", mode: 'copy', saveAs: { filename -> file(filename).name }

    input:
    tuple val(meta), path(proteins), path(gff), path(genome), val(funannotate_db), val(eggnog_db), val(interproscan_db), path(gm_key), path(gmes_tar), path(phobius_tar), path(signalp_tar)

    output:
    tuple val(meta), path("${meta.id}_functional_annotation/annotate_results/*"), emit: annotation_results

    script:
    """
    #!/bin/bash
    set -euo pipefail
        
    export PATH=/opt/conda/envs/funannotate_func/bin:/usr/local/bin:/usr/bin:/bin:\$PATH
    
    # -------------------------------------------------------------------------
    # 0. Global Environmental Variable Fixes
    # -------------------------------------------------------------------------
    export AUGUSTUS_CONFIG_PATH=/opt/conda/envs/funannotate_func/config

    # -------------------------------------------------------------------------
    # 1. EggNOG database setup
    # -------------------------------------------------------------------------
    export EGGNOG_DB_DIR="\$(realpath ${eggnog_db})"
    
    if [[ ! -d "\$EGGNOG_DB_DIR" ]]; then
        mkdir -p "\$EGGNOG_DB_DIR"
    fi

    if [[ ! -f "\$EGGNOG_DB_DIR/eggnog.db" || ! -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
        (
            flock -x 201
            if [[ ! -f "\$EGGNOG_DB_DIR/eggnog.db" || ! -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd" ]]; then
                echo ">>> Database missing. Downloading directly into: \$EGGNOG_DB_DIR" >> ${meta.id}_error.log
                
                BASE_URL="http://eggnog5.embl.de/download/emapperdb-5.0.2"
                
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog.db.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog_proteins.dmnd.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/eggnog.taxa.tar.gz"
                wget -q -P "\$EGGNOG_DB_DIR" "\$BASE_URL/mmseqs.tar.gz"
                
                echo ">>> Extracting EggNOG archives..." >> ${meta.id}_error.log
                
                gunzip -f "\$EGGNOG_DB_DIR/eggnog.db.gz"
                gunzip -f "\$EGGNOG_DB_DIR/eggnog_proteins.dmnd.gz"
                
                tar -xzf "\$EGGNOG_DB_DIR/eggnog.taxa.tar.gz" -C "\$EGGNOG_DB_DIR"
                rm -f "\$EGGNOG_DB_DIR/eggnog.taxa.tar.gz"

                tar -xzf "\$EGGNOG_DB_DIR/mmseqs.tar.gz" -C "\$EGGNOG_DB_DIR"
                rm -f "\$EGGNOG_DB_DIR/mmseqs.tar.gz"

                echo ">>> EggNOG database installation complete." >> ${meta.id}_error.log
            else
                echo ">>> Reusing existing EggNOG database (populated by parallel process)" >> ${meta.id}_error.log
            fi
        ) 201>"\$EGGNOG_DB_DIR/.download.lock"
    else
        echo ">>> Verified existing EggNOG database at \$EGGNOG_DB_DIR" >> ${meta.id}_error.log
    fi

    # -------------------------------------------------------------------------
    # 2. InterProScan setup
    # FIXED: Use full 64-bit distribution (not data-only alt/ tarball)
    #        which provides interproscan.sh + all binaries alongside data/
    # -------------------------------------------------------------------------
    export IPRSCAN_DB_DIR="\$(realpath ${interproscan_db})"

    if [[ ! -d "\$IPRSCAN_DB_DIR" ]]; then
        mkdir -p "\$IPRSCAN_DB_DIR"
    fi

    if [[ ! -f "\$IPRSCAN_DB_DIR/interproscan-5.67-99.0/interproscan.sh" ]]; then
        (
            flock -x 202
            if [[ ! -f "\$IPRSCAN_DB_DIR/interproscan-5.67-99.0/interproscan.sh" ]]; then
                echo ">>> InterProScan missing. Downloading full 64-bit distribution..." >> ${meta.id}_error.log
                
                wget -q -P "\$IPRSCAN_DB_DIR" \
                    "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/interproscan-5.67-99.0-64-bit.tar.gz"
                
                echo ">>> Extracting InterProScan..." >> ${meta.id}_error.log
                tar -xzf "\$IPRSCAN_DB_DIR/interproscan-5.67-99.0-64-bit.tar.gz" -C "\$IPRSCAN_DB_DIR"
                rm -f "\$IPRSCAN_DB_DIR/interproscan-5.67-99.0-64-bit.tar.gz"
                
                echo ">>> InterProScan installation complete." >> ${meta.id}_error.log
            else
                echo ">>> Reusing existing InterProScan installation (populated by parallel process)" >> ${meta.id}_error.log
            fi
        ) 202>"\$IPRSCAN_DB_DIR/.ipr_download.lock"
    else
        echo ">>> Verified existing InterProScan at \$IPRSCAN_DB_DIR" >> ${meta.id}_error.log
    fi
    
    # -------------------------------------------------------------------------
    # 3. Validate inputs
    # -------------------------------------------------------------------------
    if [ ! -s "${proteins}" ] || [ ! -s "${gff}" ]; then
        echo "ERROR: Protein FASTA or GFF file is missing or empty." >> ${meta.id}_error.log
        exit 1
    fi
    
    # -------------------------------------------------------------------------
    # 4. GeneMark setup
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
    # 5. Phobius setup
    # -------------------------------------------------------------------------
    if [ -s "${phobius_tar}" ] && [ "${phobius_tar.name}" != "NO_PHOBIUS_TARBALL.empty" ]; then
        tar -zxf "${phobius_tar}"
    else
        echo "WARNING: Phobius tarball missing, skipping Phobius annotation" >> ${meta.id}_error.log
        touch phobius.results.txt
    fi
    
    # -------------------------------------------------------------------------
    # 6. SignalP6 setup
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
    # 7. Funannotate database setup
    # -------------------------------------------------------------------------
    export FUNANNOTATE_DB="\$(realpath ${funannotate_db})"
    
    if [[ ! -d "\$FUNANNOTATE_DB" ]]; then
        if [[ "\$FUNANNOTATE_DB" == *"/work/funannotate_db" ]]; then
            mkdir -p "\$FUNANNOTATE_DB"
        else
            echo "ERROR: Provided DB path \$FUNANNOTATE_DB does not exist." >> ${meta.id}_error.log
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
    # 8. Phobius execution
    # -------------------------------------------------------------------------
    if [ -d "phobius" ]; then
        REAL_PROTEINS=\$(realpath "${proteins}")
        phobius/phobius.pl -short "\$REAL_PROTEINS" > phobius.results.txt 2>> ${meta.id}_error.log
    fi
    
    # -------------------------------------------------------------------------
    # 9. InterProScan
    # -------------------------------------------------------------------------
    IPRSCAN_SH="\$IPRSCAN_DB_DIR/interproscan-5.67-99.0/interproscan.sh"

    if [[ ! -f "\$IPRSCAN_SH" ]]; then
        echo "ERROR: interproscan.sh not found at \$IPRSCAN_SH" >> ${meta.id}_error.log
        exit 1
    fi
    chmod +x "\$IPRSCAN_SH"

    # InterProScan rejects sequences containing stop-codon asterisks
    tr -d '*' < "${proteins}" > proteins_for_ipr.fa


    "\$IPRSCAN_SH" \
        -i proteins_for_ipr.fa \
        -b ${meta.id}_iprscan \
        -f XML \
        -goterms \
        -pa \
        -dp \
        --cpu ${task.cpus} \
        2>> ${meta.id}_error.log

    if [[ ! -s "${meta.id}_iprscan.xml" ]]; then
        echo "ERROR: InterProScan produced no output XML" >> ${meta.id}_error.log
        exit 1
    fi
    
    # -------------------------------------------------------------------------
    # 9.5. Pre-parse InterProScan XML
    # -------------------------------------------------------------------------
    mkdir -p "${meta.id}_functional_annotation/annotate_misc"
    cp "${meta.id}_iprscan.xml" "${meta.id}_functional_annotation/annotate_misc/iprscan.xml"

    python3 "${projectDir}/scripts/iprscan2annotations_streaming.py" \
    "${meta.id}_functional_annotation/annotate_misc/iprscan.xml" \
    "${meta.id}_functional_annotation/annotate_misc/annotations.iprscan.txt" \
    2>> ${meta.id}_error.log

    if [[ ! -s "${meta.id}_functional_annotation/annotate_misc/annotations.iprscan.txt" ]]; then
        echo "WARNING: streaming InterProScan parse produced no output, funannotate will attempt its own parsing" >> ${meta.id}_error.log
    fi



    # -------------------------------------------------------------------------
    # 10. EggNOG Map
    # -------------------------------------------------------------------------
    emapper.py --cpu ${task.cpus} -m mmseqs --data_dir "\$EGGNOG_DB_DIR" -i "${proteins}" -o ${meta.id}_eggnog 2>> ${meta.id}_error.log
    
    # -------------------------------------------------------------------------
    # 11. Funannotate final synthesis
    # -------------------------------------------------------------------------
    ANNOTATE_ARGS=(
        --gff "${gff}"
        --fasta "${genome}"
        --species "${meta.species}"
        --busco_db "${meta.busco_db_fun ?: 'metazoa'}"
        --eggnog "${meta.id}_eggnog.emapper.annotations"
        --iprscan "${meta.id}_functional_annotation/annotate_misc/iprscan.xml"
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
