#!/usr/bin/env python3

"""
Script Name: stringtie2utr.py
Description: This script decorates a gtf file with genes generated by AUGUSTUS,
             BRAKER, or TSEBRA, which UTR features from a stringtie gtf file.
Author: Katharina J. Hoff
Email: katharina.hoff@uni-greifswald.de
Date: September 29th 2023

Copyright (C) The year of copyright, Katharina J. Hoff, University of Greifswald

This program is free software; you can redistribute it and/or modify
it under the terms of the Artistic License.
"""

import argparse
import re
from intervaltree import IntervalTree, Interval


def read_gtf(gtf_file):
    """
    Reads a GTF file and extracts gene and non-gene features.

    Args:
    - gtf_file (str): Path to the GTF file.

    Returns:
    tuple: (non_gene_dict, gene_dict, tx_to_gene_dict, tx_dict) where:
    - non_gene_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.
    """
    non_gene_dict = {}
    gene_dict = {}
    tx_to_gene_dict = {}
    tx_dict = {}
    transcript_id_pattern = re.compile(r'transcript_id "([^"]+)"')
    gene_id_pattern = re.compile(r'gene_id "([^"]+)"')

    with open(gtf_file, 'r') as f:
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                feature_type = fields[2]
                last_field = fields[-1]

                if feature_type == "gene":
                    gene_dict[last_field] = line.strip()  # Use the entire last field as gene_id
                elif feature_type == "transcript":
                    tx_dict[last_field] = line.strip()  # Use the entire last field as transcript_id
                else:
                    transcript_id_match = transcript_id_pattern.search(last_field)
                    gene_id_match = gene_id_pattern.search(last_field)
                    gene_id = None
                    transcript_id = None
                    if transcript_id_match:
                        transcript_id = transcript_id_match.group(1)
                    if gene_id_match:
                        gene_id = gene_id_match.group(1)
                    if transcript_id:  # Ensure we have a transcript ID
                        non_gene_dict.setdefault(transcript_id, []).append(line.strip())
                        if transcript_id not in tx_to_gene_dict:
                            tx_to_gene_dict[transcript_id] = gene_id

    return non_gene_dict, gene_dict, tx_to_gene_dict, tx_dict




def add_intron_features(gtf_dict):
    """
    Adds intron feature lines based on the exon feature lines.
    
    Args:
    - gtf_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Updated GTF dictionary with intron feature lines added.
    """
    
    for transcript_id, entries in gtf_dict.items():
        sorted_exons = sorted([entry for entry in entries if entry.split('\t')[2] == 'exon'], key=lambda x: int(x.split('\t')[3]))
        
        new_entries = []
        for i in range(len(sorted_exons) - 1):
            curr_exon = sorted_exons[i]
            next_exon = sorted_exons[i + 1]
            
            curr_end = int(curr_exon.split('\t')[4])
            next_start = int(next_exon.split('\t')[3])
            
            if next_start - curr_end > 1:
                # Construct intron feature line based on exon line, but replace the feature name with "intron"
                intron_line = curr_exon.split('\t')
                intron_line[2] = 'intron'
                intron_line[3] = str(curr_end + 1)
                intron_line[4] = str(next_start - 1)
                
                # Remove exon_number attribute
                intron_line[8] = re.sub(r'exon_number "\d+";', '', intron_line[8]).strip()

                # Add the intron line to new entries
                new_entries.append('\t'.join(intron_line))

        # Extend the original entries with new intron feature lines
        gtf_dict[transcript_id].extend(new_entries)

    return gtf_dict



def create_introns_hash(non_gene_dict):
    """
    Creates a dictionary with intron strings as keys and transcript IDs as values.
    
    Args:
    - non_gene_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Dictionary with intron strings as keys and transcript IDs as values.
    """
    
    introns_hash = {}
    
    for transcript_id, entries in non_gene_dict.items():
        for entry in entries:
            # Check if the feature is an intron
            if entry.split('\t')[2] == 'intron':
                seqname = entry.split('\t')[0]
                start = entry.split('\t')[3]
                end = entry.split('\t')[4]
                strand = entry.split('\t')[6]
                intron_key = f"{seqname}_{start}_{end}_{strand}"

                # Store the transcript ID associated with the intron in the hash
                introns_hash[intron_key] = transcript_id

    return introns_hash


def find_matching_transcripts(intron_hash1, intron_hash2):
    """
    Find matching transcript IDs based on intron patterns.
    This only works for multi-exon genes.
    
    Args:
    - intron_hash1 (dict): Dictionary with intron strings as keys and transcript IDs from the first dataset as values.
    - intron_hash2 (dict): Dictionary with intron strings as keys and transcript IDs from the second dataset as values.

    Returns:
    dict: Dictionary with transcript IDs from intron_hash1 as keys and lists of matching transcript IDs from intron_hash2 as values.
    """
    
    # Reverse the hashes for easy lookup of intron patterns for each transcript
    reverse_hash1 = {}
    for intron, transcript in intron_hash1.items():
        if transcript not in reverse_hash1:
            reverse_hash1[transcript] = []
        reverse_hash1[transcript].append(intron)

    reverse_hash2 = {}
    for intron, transcript in intron_hash2.items():
        if transcript not in reverse_hash2:
            reverse_hash2[transcript] = []
        reverse_hash2[transcript].append(intron)

    matching_transcripts = {}

    for transcript1, introns1 in reverse_hash1.items():
        matches = set()
        
        for transcript2, introns2 in reverse_hash2.items():
            if all(intron in introns2 for intron in introns1):
                matches.add(transcript2)

        if matches:
            matching_transcripts[transcript1] = list(matches)

    return matching_transcripts


def find_single_exon_genes(non_gene_dict):
    """
    Find single exon genes in non_gene_dict.
    
    Args:
    - non_gene_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: List of transcript IDs for single exon genes.
    """
    
    single_exon_genes = {}
    
    for transcript_id, entries in non_gene_dict.items():
        # Check if the transcript has only one exon
        if len([entry for entry in entries if entry.split('\t')[2] == 'exon']) == 1:
            single_exon_genes[transcript_id] = True

    return single_exon_genes

def tx_len(non_gene_dict):
    """
    Calculate the length of each transcript.

    Args:
    - non_gene_dict (dict): Dictionary with transcript_id as key and a list of 
    GTF entries as values.

    Returns:
    dict: Dictionary with transcript IDs as keys and their lengths as values.
    """
    
    transcript_lengths = {}

    for transcript_id, entries in non_gene_dict.items():
        total_length = 0
        
        for entry in entries:
            fields = entry.split('\t')
            start = int(fields[3])
            end = int(fields[4])
            
            total_length += (end - start + 1)  # +1 because both start and end are inclusive

        transcript_lengths[transcript_id] = total_length

    return transcript_lengths


def select_longest_tx(matched_transcripts, tx_lens_stringtie):
    """
    For each key in matched_transcripts, select the longest transcript based on lengths provided in tx_lens_stringtie.

    Args:
    - matched_transcripts (dict): Dictionary with keys as transcript IDs and values as lists of matching transcript IDs.
    - tx_lens_stringtie (dict): Dictionary with transcript IDs as keys and their lengths as values.

    Returns:
    dict: Dictionary with keys from matched_transcripts and values as the IDs of the longest transcripts.
    """
    
    selected_transcripts = {}

    for original_tx, matched_tx_list in matched_transcripts.items():
        # Find the longest transcript in the matched_tx_list
        longest_tx = max(matched_tx_list, key=lambda tx: tx_lens_stringtie.get(tx, 0))

        selected_transcripts[original_tx] = longest_tx

    return selected_transcripts


def overlap(exon_start, exon_end, cds_start, cds_end):
    """Check if exon overlaps with CDS."""
    return exon_start <= cds_end and exon_end >= cds_start


def merge_features(tsebra_gtf, stringtie_gtf, selected_transcripts):
    for tsebra_tx, stringtie_tx in selected_transcripts.items():
        # Retrieve features for current transcripts
        tsebra_features = tsebra_gtf[tsebra_tx]
        stringtie_exons = [f for f in stringtie_gtf[stringtie_tx] if f.split('\t')[2] == 'exon']
        tsebra_cds_features = [f for f in tsebra_features if f.split('\t')[2] == 'CDS']
        
        for exon in stringtie_exons:
            exon_parts = exon.split('\t')
            exon_start, exon_end = int(exon_parts[3]), int(exon_parts[4])
            
            # Flag to determine if current exon should be merged
            merge_exon = True
            
            for cds in tsebra_cds_features:
                cds_parts = cds.split('\t')
                cds_start, cds_end = int(cds_parts[3]), int(cds_parts[4])
                
                if overlap(exon_start, exon_end, cds_start, cds_end):
                    # Exon overlaps with CDS. Check the length condition.
                    if exon_end - exon_start + 1 < cds_end - cds_start + 1:
                        merge_exon = False
                        break
            
            # Add exon to tsebra_features if it passed the checks
            if merge_exon:
                tsebra_features.append(exon)
        
        # Update the tsebra_gtf dictionary with new features
        tsebra_gtf[tsebra_tx] = tsebra_features

    return tsebra_gtf


def compute_utr_features(tsebra_gtf):
    """
    Compute the UTR features for each transcript in tsebra_gtf based on strand information.

    Args:
    - tsebra_gtf (dict): Dictionary with transcript IDs as keys and lists of GTF feature lines as values.

    Returns:
    dict: Updated tsebra_gtf dictionary with UTR features added.
    """
    for transcript_id, features in tsebra_gtf.items():
        # Sort features by start position
        features.sort(key=lambda x: int(x.split('\t')[3]))
        
        utr_features = []
        for feature in features:
            fields = feature.split('\t')
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # If feature is an exon, check if there are overlapping CDS features
            if feature_type == "exon":
                overlapping_cds = [f for f in features if f.split('\t')[2] == "CDS" and int(f.split('\t')[3]) <= end and int(f.split('\t')[4]) >= start]
                
                if overlapping_cds:
                    cds_start = int(overlapping_cds[0].split('\t')[3])
                    cds_end = int(overlapping_cds[0].split('\t')[4])

                    # Check for UTR based on strand
                    if strand == "+":
                        if start < cds_start:
                            utr5 = "\t".join(fields[:2] + ["five_prime_UTR"] + [fields[3]] + [str(cds_start - 1)] + fields[5:])
                            utr_features.append(utr5)

                        if end > cds_end:
                            utr3 = "\t".join(fields[:2] + ["three_prime_UTR"] + [fields[3]] + [str(cds_end + 1)] + fields[5:])
                            utr_features.append(utr3)
                    
                    elif strand == "-":
                        if end > cds_end:
                            utr5 = "\t".join(fields[:2] + ["five_prime_UTR"] + [str(cds_end + 1)] + fields[4:])
                            utr_features.append(utr5)

                        if start < cds_start:
                            utr3 = "\t".join(fields[:2] + ["three_prime_UTR"] + [fields[3]] + [str(cds_start - 1)] + fields[5:])
                            utr_features.append(utr3)
                
        # Add the computed UTR features to the list of features for this transcript
        features.extend(utr_features)
        features.sort(key=lambda x: int(x.split('\t')[3]))

        # Modify StringTie exons to UTR features if needed
        seen_cds = False
        for idx, feature in enumerate(features):
            fields = feature.split('\t')
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            featuretype = fields[2]
            if featuretype == "CDS":
                seen_cds = True
            if strand == "+" and not seen_cds:
                utr_type = "five_prime_UTR"
            elif strand == "-" and not seen_cds:
                utr_type = "three_prime_UTR"
            elif strand == "+" and seen_cds:
                utr_type = "three_prime_UTR"
            else:
                utr_type = "five_prime_UTR"
            
            if "StringTie" in fields[1]:  # Assuming this condition identifies StringTie exons
                # Check if it has the same start as any of its neighbors
                same_start = False
                for offset in [-1, 1]:
                    if idx + offset >= 0 and idx + offset < len(features):
                        neighbor = features[idx + offset].split('\t')
                        if int(neighbor[3]) == start:
                            same_start = True
                            break
                
                if same_start:
                    # Remove the StringTie feature
                    features.pop(idx)
                else:
                    fields[2] = utr_type
                    features[idx] = "\t".join(fields)

    return tsebra_gtf


def fix_feature_coordinates(gtf_dict, gene_dict, tx_to_gene_dict, tx_dict):
    """
    The list of features in gtf_dict has been expanded compared to the original version. We need to identify the left most and right most coordinate
     for each transcript. Then, update the tx_dict lines with these new coordinates. In the end, loop over all tx_dict entries and identify the left most and right
      most coordinate per gene, then update the gene lines in gene_dict.
    
    Args:
    - gtf_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.

    Returns:
    gene_dict, tx_dict with updated coordinates
    """
    # 1. Update transcript coordinates based on feature coordinates in gtf_dict
    for tx_id, features in gtf_dict.items():
        starts = []
        ends = []
        for feature in features:
            fields = feature.split('\t')
            starts.append(int(fields[3]))
            ends.append(int(fields[4]))
        
        # Find the min and max coordinates for the transcript
        min_coord = min(starts)
        max_coord = max(ends)
        
        tx_fields = tx_dict[tx_id].split('\t')
        tx_fields[3] = str(min_coord)
        tx_fields[4] = str(max_coord)
        tx_dict[tx_id] = '\t'.join(tx_fields)
    
    # 2. Update gene coordinates based on updated transcript coordinates in tx_dict
    for gene_id in gene_dict:
        associated_transcripts = [tx for tx, gid in tx_to_gene_dict.items() if gid == gene_id]
        
        starts = []
        ends = []
        for tx_id in associated_transcripts:
            tx_fields = tx_dict[tx_id].split('\t')
            starts.append(int(tx_fields[3]))
            ends.append(int(tx_fields[4]))
        
        # Find the min and max coordinates for the gene
        min_coord = min(starts)
        max_coord = max(ends)
        
        gene_fields = gene_dict[gene_id].split('\t')
        gene_fields[3] = str(min_coord)
        gene_fields[4] = str(max_coord)
        gene_dict[gene_id] = '\t'.join(gene_fields)
    
    return gene_dict, tx_dict


def print_gtf(filename, gtf_dict, gene_dict, tx_to_gene_dict, tx_dict):
    """
    Print GTF lines based on gene_dict and gtf_dict.
    
    Args:
    - filename (str): Path to the output file.
    - gtf_dict (dict): Dictionary with transcript IDs as keys and lists of non-gene feature lines as values.
    - gene_dict (dict): Dictionary with gene IDs extracted from the last column of gene feature lines as keys and the corresponding gene line as value.
    - tx_to_gene_dict (dict): Dictionary with transcript IDs as keys and the corresponding gene ID as value, extracted from the last column of transcript feature lines.
    - tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value. The entire last column is used as the transcript ID.

    Returns:
    None: Prints the GTF lines to stdout.
    """
    printed_gene = {}
    try:
        with open(filename, 'w') as f:
            # Iterate over gene_dict entries
            for tx_id, tx_line in tx_dict.items():
                gene_id = tx_to_gene_dict[tx_id]
                if not gene_id in printed_gene:
                    f.write(gene_dict[tx_to_gene_dict[tx_id]] + "\n")
                    printed_gene[gene_id] = True
                f.write(tx_line + "\n")
                sorted_features = sorted(gtf_dict.get(tx_id, []), key=lambda x: int(x.split('\t')[3]))
                for feature in sorted_features:
                    # If the feature is a UTR line, remove the exon_number
                    if "UTR" in feature:
                        # split line into fields
                        fields = feature.split("\t")
                        # build new gtf line
                        f.write(fields[0] + "\tstringtie2utr\t" + "\t".join(fields[2:8]) + "\ttranscript_id \"" + tx_id + "\"; gene_id \"" + tx_to_gene_dict[tx_id] + "\";" + "\n")
                    else:
                        f.write(feature + "\n")
    except IOError:
        print("Could not write to file: " + filename)
        exit(1)

def build_tree(data):
    """Build an interval tree from data. We will use that to quickly find the overlapping single exon genes/transcript pairs."""
    tree = IntervalTree()
    for item in data:
        start, end, value = item
        tree.addi(start, end, value)
    return tree

def construct_transcript_tree(tx_dict):
    """Contructs a dictionary that has sequence name as first key, then strand as second key, and transcript interval tree as value.
    Internally, before calling build_tree(data), a data structure for each sequence and strand must be constructed. It looks like this:
    transcripts = [(50, 550, "tx1"), (580, 950, "tx2")]
    
    Args: 
    - tx_dict (dict): Dictionary with transcript IDs as keys and lists of transcript feature lines as values.
    """
    seq_strand_to_transcripts = {}
    # Parse gtf_dict to segregate transcript data by sequence name and strand
    for transcript_id, transcript_data in tx_dict.items():
        fields = transcript_data.split('\t')
        seq_name = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        # Use setdefault to initialize a nested dictionary
        seq_strand_dict = seq_strand_to_transcripts.setdefault(seq_name, {})
        # Use setdefault to initialize a list for the current strand
        data_list = seq_strand_dict.setdefault(strand, [])
        data_list.append((start, end, transcript_id))
    
    # Build interval trees for each seq_name-strand combination
    for seq_name, strand_data in seq_strand_to_transcripts.items():
        for strand, data_list in strand_data.items():
            seq_strand_to_transcripts[seq_name][strand] = build_tree(data_list)

    return seq_strand_to_transcripts

def construct_gene_tree(tx_dict, single_exons):
    """Contructs a dictionary that has sequence name as first key, then strand as second key, and transcript interval tree as value.
    Internally, before calling build_tree(data), a data structure for each sequence and strand must be constructed. It looks like this:
    transcripts = [(50, 550, "tx1"), (580, 950, "tx2")]
    
    Args: 
    - tx_dict (dict): Dictionary with transcript IDs as keys and lists of transcript feature lines as values.
    - single_exons (dict): Dictionary with single-exon transcript IDs. We only want to build the tree for single-exon transcripts.
    """
    seq_strand_to_genes = {}
    # Parse tx_dict to segregate gene data by sequence name and strand, but only for single-exon transcripts
    for transcript_id, line in tx_dict.items():
        # Check if the transcript_id is in single_exons
        if transcript_id in single_exons:
            fields = line.split('\t')
            seq_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            # Use setdefault to initialize a nested dictionary
            seq_strand_dict = seq_strand_to_genes.setdefault(seq_name, {})
            # Use setdefault to initialize a list for the current strand
            data_list = seq_strand_dict.setdefault(strand, [])
            data_list.append((start, end, transcript_id))
    
    # Build interval trees for each seq_name-strand combination
    for seq_name, strand_data in seq_strand_to_genes.items():
        for strand, data_list in strand_data.items():
            seq_strand_to_genes[seq_name][strand] = build_tree(data_list)

    return seq_strand_to_genes

def find_overlapping_transcripts(gene_tree, transcript_tree):
    """Find transcripts that overlap with genes."""
    gene_to_transcripts = {}

    # loop over the sequences in outer dict:
    for seq_name, strand_data in gene_tree.items():
        # loop over the strands in the inner dict:
        for strand, gene_tree in strand_data.items():
            # loop over the genes in the gene_tree:
            for gene_interval in gene_tree:
                if seq_name in transcript_tree:
                    if strand in transcript_tree[seq_name]:
                        overlapping_transcripts = transcript_tree[seq_name][strand].overlap(gene_interval.begin, gene_interval.end)
                        if overlapping_transcripts:
                            gene_id = gene_interval.data
                            gene_to_transcripts[gene_id] = [tx.data for tx in overlapping_transcripts]
                    
    return gene_to_transcripts

def extract_introns_from_dict(tx_dict, tx_id):
    """Extract introns directly from stringtie_tx_dict."""
    introns = []
    for feature_line in tx_dict[tx_id]:
        feature_type, start, end = feature_line.split("\t")[2:5]
        if feature_type == "intron":
            introns.append((int(start), int(end)))
    return introns

def check_overlap_compatibility(gene_to_transcripts, tsebra_tx_dict, stringtie_tx_dict):
    """Check if the overlapping transcripts are compatible with the gene models.
    In particular, we want to make sure that no stringtie transcript contains an intron within
    the CDS feature range of the braker transcripts. Only verfy this for those gene/transcript list matches that
    are in gene_to_transcripts. If an intron breaks the gene/transcript association, the transcript should be removed from the list in gene_to_transcripts.
    If the list is empty, also the key should be removed.

    Args:
    - gene_to_transcipts (dict): has transcript IDs from braker as keys and lists of transcript IDs from stringtie as values.
    - tsebra_tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value.
    - stringtie_tx_dict (dict): Dictionary with transcript IDs as keys and the corresponding transcript line as value.

    Returns:
    dict: Updated gene_to_transcripts dictionary.
    """
    for tsebra_tx_id, stringtie_tx_ids in list(gene_to_transcripts.items()):  # Using list() to allow dictionary changes during iteration
        tsebra_tx = tsebra_tx_dict[tsebra_tx_id]
        for line in tsebra_tx:
            if re.search(r'\tCDS\t', line):
                fields = line.split('\t')
                start = int(fields[3])
                end = int(fields[4])
                break # there's only one CDS because it's a single-exon transcript
        # Assuming the CDS range is at the end of the transcript line in the format start-end
        tsebra_cds_start, tsebra_cds_end = [start,end]
        
        # Store compatible transcripts
        compatible_tx = []
        
        for st_tx_id in stringtie_tx_ids:
            tx_id = re.search(r'gene_id "[^"]+"; transcript_id "([^"]+)";', st_tx_id).group(1)
            introns = extract_introns_from_dict(stringtie_tx_dict, tx_id)
            
            # Check if any intron overlaps with the CDS range
            overlap = any(intron_start <= tsebra_cds_end and intron_end >= tsebra_cds_start for intron_start, intron_end in introns)
            
            if not overlap:
                compatible_tx.append(tx_id)
                
        # Update the dictionary
        if compatible_tx:
            gene_to_transcripts[tsebra_tx_id] = compatible_tx
        else:
            del gene_to_transcripts[tsebra_tx_id]
    
    return gene_to_transcripts
    
    
def main():
    parser = argparse.ArgumentParser(description="Updates GaiusAugustus gene models with UTRs from a StringTie assembly.")
    
    # Mandatory input arguments
    parser.add_argument("-g", "--genes", required=True, help="File with gene gene models in GTF format, typically the output of Augustus, BRAKER, or TSEBRA.")
    parser.add_argument("-s", "--stringtie", required=True, help="File with StringTie transcript models in GFF format")
    parser.add_argument("-o", "--output", required=True, help="Output file name, file is in GTF format.")

    args = parser.parse_args()

    tsebra_file = args.genes 
    stringtie_file = args.stringtie

    # read the tsebra_file with gene models
    tsebra_none_gene_dict, tsebra_gene_line_dict, tsebra_tx_to_gene_dict, tsebra_tx_dict = read_gtf(tsebra_file)
    # read the stringtie_file
    stringtie_none_gene_dict, stringtie_gene_line_dict, stringtie_tx_to_gene_dict, stringtie_tx_dict = read_gtf(stringtie_file)
    # add intron features to the stringtie_none_gene_dict
    stringtie_none_gene_dict = add_intron_features(stringtie_none_gene_dict)

    # in order to match genes and transcripts efficiently, we will take 2 different approaches. For multi-exon genes,
    # we create intron hashes for matching transcripts. For single-exon genes, we construct interval trees.
    # first the multi-exon gene matching
    tsebra_introns_hash = create_introns_hash(tsebra_none_gene_dict)
    stringtie_introns_hash = create_introns_hash(stringtie_none_gene_dict)
    matched_transcripts = find_matching_transcripts(tsebra_introns_hash, stringtie_introns_hash)

    # find single exon genes, these now do not have UTRs yet
    single_exon_genes = find_single_exon_genes(tsebra_none_gene_dict)
    
    # construct transcript trees
    seq_strand_to_transcripts = construct_transcript_tree(stringtie_tx_dict)
    seq_strand_to_genes = construct_gene_tree(tsebra_tx_dict, single_exon_genes)
    overlaps = find_overlapping_transcripts(seq_strand_to_genes, seq_strand_to_transcripts)
    # the overlaps might include matches where stringtie introns break cds, remove these
    overlaps = check_overlap_compatibility(overlaps, tsebra_none_gene_dict, stringtie_none_gene_dict)
    # merge the single exon matches on top of the multi exon matches
    matched_transcripts.update(overlaps)

    # We may now have several alternative RNA-Seq inferred transcripts that match a single predicted transcript
    # we will select the longest of these. For this calculate the length of each transcript
    tx_lens_stringtie = tx_len(stringtie_none_gene_dict)
    final_matching_tx = select_longest_tx(matched_transcripts, tx_lens_stringtie)

    # merge features from stringtie_gtf into tsebra_gtf based on the selected transcripts, this is a simple concatenation,
    # UTR features are not inferred, yet
    tsebra_gtf = merge_features(tsebra_none_gene_dict, stringtie_none_gene_dict, final_matching_tx)
    # compute UTR features, remove StringTie exon features
    tsebra_gtf = compute_utr_features(tsebra_gtf)

    # expanding gene/transcripts by UTRs shifts the coordinates of these features, we need to update the coordinates in the gene and transcript lines
    tsebra_gene_line_dict, tsebra_tx_dict = fix_feature_coordinates(tsebra_gtf, tsebra_gene_line_dict, tsebra_tx_to_gene_dict, tsebra_tx_dict)

    # print the updated tsebra_gtf
    print_gtf(args.output, tsebra_gtf, tsebra_gene_line_dict, tsebra_tx_to_gene_dict, tsebra_tx_dict)

if __name__ == "__main__":
    main()
