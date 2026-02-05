#!/usr/bin/env python3
import argparse
import re
import os
import sys
import logging

def setup_logging():
    """Set up logging to file and console."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('busco_comparison.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def parse_busco_summary(file_path):
    """Parse BUSCO summary file and extract completeness percentage."""
    if not file_path or not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        logging.warning(f"BUSCO summary file {file_path} is missing or empty")
        return 0.0
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        match = re.search(r'C:(\d+\.\d+)%', content)
        if match:
            score = float(match.group(1))
            logging.info(f"Parsed BUSCO score from {file_path}: {score}%")
            return score
        logging.warning(f"No BUSCO completeness score found in {file_path}")
        return 0.0
    except Exception as e:
        logging.error(f"Error parsing {file_path}: {e}")
        return 0.0

def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Compare BUSCO summary files and select the best tool.")
    parser.add_argument('--funannotate_busco', help="Path to Funannotate BUSCO summary file")
    parser.add_argument('--braker_busco', help="Path to BRAKER BUSCO summary file")
    parser.add_argument('--funannotate_proteins', help="Path to Funannotate proteins FASTA")
    parser.add_argument('--funannotate_gff', help="Path to Funannotate GFF3")
    parser.add_argument('--braker_proteins', help="Path to BRAKER proteins FASTA")
    parser.add_argument('--braker_gff', help="Path to BRAKER GFF3")
    args = parser.parse_args()

    # Parse BUSCO scores
    fa_score = parse_busco_summary(args.funannotate_busco)
    br_score = parse_busco_summary(args.braker_busco)

    logging.info(f"Funannotate BUSCO Complete: {fa_score}%")
    logging.info(f"BRAKER BUSCO Complete: {br_score}%")

    # Check file existence and validity
    fa_proteins_valid = args.funannotate_proteins and os.path.exists(args.funannotate_proteins) and os.path.getsize(args.funannotate_proteins) > 0
    fa_gff_valid = args.funannotate_gff and os.path.exists(args.funannotate_gff) and os.path.getsize(args.funannotate_gff) > 0
    br_proteins_valid = args.braker_proteins and os.path.exists(args.braker_proteins) and os.path.getsize(args.braker_proteins) > 0
    br_gff_valid = args.braker_gff and os.path.exists(args.braker_gff) and os.path.getsize(args.braker_gff) > 0

    # Select the best tool
    if fa_score > 0 and fa_proteins_valid and fa_gff_valid and (fa_score >= br_score or br_score == 0):
        selected_tool = "funannotate"
        selected_proteins = args.funannotate_proteins
        selected_gff = args.funannotate_gff
    elif br_score > 0 and br_proteins_valid and br_gff_valid:
        selected_tool = "braker"
        selected_proteins = args.braker_proteins
        selected_gff = args.braker_gff
    else:
        logging.warning("No valid BUSCO summaries or output files provided. Defaulting to BRAKER dummy outputs.")
        selected_tool = "braker"
        selected_proteins = "NO_BR_PROTEINS.fa"
        selected_gff = "NO_BR_GFF.gff3"

    # Ensure dummy files exist if selected
    if selected_proteins == "NO_BR_PROTEINS.fa" and not os.path.exists(selected_proteins):
        with open(selected_proteins, 'w') as f:
            f.write("")
    if selected_gff == "NO_BR_GFF.gff3" and not os.path.exists(selected_gff):
        with open(selected_gff, 'w') as f:
            f.write("")

    logging.info(f"Selected tool: {selected_tool}")
    logging.info(f"Selected proteins: {selected_proteins}")
    logging.info(f"Selected GFF: {selected_gff}")

    # Write to file for Nextflow
    with open("busco_comparison.txt", "w") as f:
        f.write(f"tool={selected_tool}\n")
        f.write(f"proteins={selected_proteins}\n")
        f.write(f"gff={selected_gff}\n")

if __name__ == "__main__":
    main()
