#!/usr/bin/env python3
"""
Repeat Structure Classifier
Classifies the structure of reads containing repeat expansions by comparing
sliding windows to known unit sequences using minimap2.

Author: Freddie
Date: 2026-02-04
"""

import argparse
import subprocess
import tempfile
import os
from pathlib import Path
from collections import defaultdict


def parse_fasta(fasta_file):
    """
    Parse FASTA file and return dictionary of {header: sequence}
    """
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:].split()[0]  # Get first word after >
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


def sanitize_name(name):
    """
    Sanitize a name to be safe for use in filenames
    Replace problematic characters with underscores
    """
    # Replace forward slashes, backslashes, and other problematic characters
    return name.replace('/', '_').replace('\\', '_').replace(' ', '_')


def write_fasta(header, sequence, filename):
    """
    Write a single sequence to a FASTA file
    """
    with open(filename, 'w') as f:
        f.write(f">{header}\n")
        # Write sequence in 80-character lines
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + '\n')


def run_minimap2(query_fasta, ref_fasta, preset='map-hifi'):
    """
    Run minimap2 and return PAF output as list of lines
    """
    cmd = [
        'minimap2',
        '-c',  # Output in PAF format
        '-x', preset,
        ref_fasta,
        query_fasta
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return result.stdout.strip().split('\n') if result.stdout.strip() else []
    except subprocess.CalledProcessError as e:
        print(f"Warning: minimap2 failed with error: {e}")
        return []


def parse_paf_line(paf_line):
    """
    Parse a single PAF line and return relevant information
    Returns: dict with query_name, query_len, query_start, query_end, strand,
             target_name, target_len, target_start, target_end, matches, aln_len, mapq
    """
    if not paf_line:
        return None
    
    fields = paf_line.split('\t')
    if len(fields) < 12:
        return None
    
    return {
        'query_name': fields[0],
        'query_len': int(fields[1]),
        'query_start': int(fields[2]),
        'query_end': int(fields[3]),
        'strand': fields[4],
        'target_name': fields[5],
        'target_len': int(fields[6]),
        'target_start': int(fields[7]),
        'target_end': int(fields[8]),
        'matches': int(fields[9]),
        'aln_len': int(fields[10]),
        'mapq': int(fields[11])
    }


def calculate_identity(matches, aln_len):
    """
    Calculate percent identity from PAF alignment
    """
    if aln_len == 0:
        return 0.0
    return (matches / aln_len) * 100


def is_uppercase_assignment(target_start, target_len, threshold=0.3):
    """
    Determine if assignment should be uppercase based on alignment position
    Returns True (uppercase) if alignment includes first 30% of reference
    """
    # If alignment starts within first 30% of reference, it's uppercase
    return (target_start / target_len) < threshold


def classify_window(window_seq, window_name, unit_seqs, unit_names, 
                   temp_dir, identity_threshold=40.0, uppercase_threshold=0.3,
                   preset='map-hifi'):
    """
    Classify a single window against all unit sequences
    
    Returns:
        - assigned_unit: The assigned unit (uppercase or lowercase) or 'U'
        - scores_dict: Dictionary with alignment info for each unit
        - strand: The strand of the best alignment ('+', '-', or None)
    """
    # Write window to temporary FASTA (sanitize name for filesystem)
    safe_window_name = sanitize_name(window_name)
    window_fasta = os.path.join(temp_dir, f"{safe_window_name}.fa")
    write_fasta(window_name, window_seq, window_fasta)
    
    # Store results for each unit
    scores_dict = {}
    best_identity = 0.0
    best_unit = None
    best_target_start = None
    best_target_len = None
    best_strand = None
    
    # Run minimap2 against each unit separately
    for unit_name in unit_names:
        unit_fasta = os.path.join(temp_dir, f"unit_{unit_name}.fa")
        write_fasta(unit_name, unit_seqs[unit_name], unit_fasta)
        
        paf_output = run_minimap2(window_fasta, unit_fasta, preset)
        
        if paf_output and paf_output[0]:
            # Parse the best alignment (first line)
            aln = parse_paf_line(paf_output[0])
            if aln:
                identity = calculate_identity(aln['matches'], aln['aln_len'])
                scores_dict[unit_name] = {
                    'identity': identity,
                    'target_start': aln['target_start'],
                    'target_end': aln['target_end'],
                    'target_len': aln['target_len'],
                    'strand': aln['strand']
                }
                
                # Track best alignment
                if identity > best_identity:
                    best_identity = identity
                    best_unit = unit_name
                    best_target_start = aln['target_start']
                    best_target_len = aln['target_len']
                    best_strand = aln['strand']
            else:
                scores_dict[unit_name] = {
                    'identity': 0.0,
                    'target_start': 'NA',
                    'target_end': 'NA',
                    'target_len': len(unit_seqs[unit_name]),
                    'strand': 'NA'
                }
        else:
            # No alignment found
            scores_dict[unit_name] = {
                'identity': 0.0,
                'target_start': 'NA',
                'target_end': 'NA',
                'target_len': len(unit_seqs[unit_name]),
                'strand': 'NA'
            }
    
    # Determine final assignment
    if best_identity < identity_threshold or best_unit is None:
        assigned_unit = 'U'
    else:
        # Check if should be uppercase
        if is_uppercase_assignment(best_target_start, best_target_len, uppercase_threshold):
            assigned_unit = best_unit.upper()
        else:
            assigned_unit = best_unit.lower()
    
    return assigned_unit, scores_dict, best_strand


def process_read(read_name, read_seq, unit_seqs, unit_names, temp_dir,
                window_size=1000, identity_threshold=40.0, 
                uppercase_threshold=0.3, preset='map-hifi'):
    """
    Process a single read by sliding window analysis
    
    Returns:
        - structure: String like "A-a-a-B-C-U-D"
        - detailed_scores: List of dictionaries with scores for each window
        - read_strand: Overall strand of the read ('+', '-', or 'mixed')
    """
    structure = []
    detailed_scores = []
    strand_votes = []  # Track strand for each assigned window
    
    read_len = len(read_seq)
    window_num = 0
    
    # Slide window across read
    for start in range(0, read_len, window_size):
        window_num += 1
        end = min(start + window_size, read_len)
        window_seq = read_seq[start:end]
        window_name = f"{read_name}_window{window_num}"
        
        # Classify this window
        assigned_unit, scores, strand = classify_window(
            window_seq, window_name, unit_seqs, unit_names, temp_dir,
            identity_threshold, uppercase_threshold, preset
        )
        
        structure.append(assigned_unit)
        
        # Track strand for non-U assignments
        if assigned_unit != 'U' and strand is not None:
            strand_votes.append(strand)
        
        # Store detailed scores
        score_entry = {
            'read_name': read_name,
            'window_num': window_num,
            'window_start': start,
            'window_end': end,
            'assigned_unit': assigned_unit,
            'strand': strand if strand else 'NA'
        }
        
        # Add scores for each unit
        for unit_name in unit_names:
            score_entry[f'{unit_name}_identity'] = scores[unit_name]['identity']
            score_entry[f'{unit_name}_ref_start'] = scores[unit_name]['target_start']
            score_entry[f'{unit_name}_ref_end'] = scores[unit_name]['target_end']
            score_entry[f'{unit_name}_strand'] = scores[unit_name]['strand']
        
        detailed_scores.append(score_entry)
    
    # Determine overall read strand by majority vote
    if not strand_votes:
        read_strand = 'unknown'
    else:
        plus_count = strand_votes.count('+')
        minus_count = strand_votes.count('-')
        
        if plus_count > minus_count * 2:  # Strongly forward
            read_strand = '+'
        elif minus_count > plus_count * 2:  # Strongly reverse
            read_strand = '-'
        else:
            read_strand = 'mixed'
    
    structure_string = '-'.join(structure)
    return structure_string, detailed_scores, read_strand


def interpret_structure(structure_string, strand='+'):
    """
    Convert structure string like 'b-C-C-C-c-c-c-c-c-c-C-C-C-c-c-c-c-c'
    into simplified format like 'B-C×2'
    
    If strand is '-', the structure is reversed before interpretation to give
    the true 5'->3' structure of the read.
    
    Logic:
    1. If reverse strand: reverse the structure string
    2. Group consecutive windows of same unit
    3. Count tandem copies by looking for lowercase→uppercase transitions
    4. Simplify consecutive Us to single U
    5. Output format: Unit×count (if count > 1)
    """
    # If reverse strand, reverse the structure to get true 5'->3' orientation
    if strand == '-':
        elements = structure_string.split('-')
        elements.reverse()
        structure_string = '-'.join(elements)
    
    elements = structure_string.split('-')
    
    result = []
    i = 0
    
    while i < len(elements):
        current_unit = elements[i].upper()
        
        # Special handling for U (undetermined)
        if current_unit == 'U':
            # Skip all consecutive Us
            j = i
            while j < len(elements) and elements[j].upper() == 'U':
                j += 1
            result.append('U')
            i = j
            continue
        
        # Count tandem copies for this unit
        tandem_count = 1  # Start with 1 (first occurrence)
        j = i
        prev_was_lower = False
        
        while j < len(elements) and elements[j].upper() == current_unit:
            # Count transition from lowercase to uppercase as new tandem
            if elements[j].isupper() and prev_was_lower:
                tandem_count += 1
            
            prev_was_lower = elements[j].islower()
            j += 1
        
        # Format output
        if tandem_count > 1:
            result.append(f"{current_unit}×{tandem_count}")
        else:
            result.append(current_unit)
        
        i = j
    
    return '-'.join(result)


def has_internal_u(structure_list):
    """
    Returns True if the read has U windows worth investigating.
    A single trailing U is ignored — this is typically a short leftover
    window at the end of the read with too few bases to align.
    A run of 2+ consecutive trailing Us is kept, as this likely represents
    a genuine undetermined region rather than a truncation artefact.
    """
    trimmed = list(structure_list)
    # Only strip a single trailing U, not a whole run
    if trimmed and trimmed[-1] == 'U' and (len(trimmed) < 2 or trimmed[-2] != 'U'):
        trimmed.pop()
    return 'U' in trimmed


def get_u_regions(read_name, read_seq, structure_list, window_size):
    """
    Given a read structure list and its sequence, return:
      masked_seq : full read where internal-U-window bases are UPPERCASE,
                   all other bases lowercase
      u_regions  : list of dicts, one per contiguous internal-U stretch:
                     {region_idx, start_bp, end_bp, seq, window_indices}
    A single trailing U is excluded. A run of 2+ trailing Us is kept.
    """
    trimmed = list(structure_list)
    # Only strip a single trailing U, not a whole run
    if trimmed and trimmed[-1] == 'U' and (len(trimmed) < 2 or trimmed[-2] != 'U'):
        trimmed.pop()
    internal_u_indices = {i for i, unit in enumerate(trimmed) if unit == 'U'}

    # Build masked sequence
    read_len = len(read_seq)
    masked_chars = []
    for i, unit in enumerate(structure_list):
        start = i * window_size
        end = min(start + window_size, read_len)
        window_seq = read_seq[start:end]
        if i in internal_u_indices:
            masked_chars.append(window_seq.upper())
        else:
            masked_chars.append(window_seq.lower())
    masked_seq = ''.join(masked_chars)

    # Group contiguous internal-U windows into regions
    u_regions = []
    if not internal_u_indices:
        return masked_seq, u_regions

    sorted_indices = sorted(internal_u_indices)
    region_start_idx = sorted_indices[0]
    prev_idx = sorted_indices[0]

    for idx in sorted_indices[1:]:
        if idx == prev_idx + 1:
            prev_idx = idx
        else:
            bp_start = region_start_idx * window_size
            bp_end = min((prev_idx + 1) * window_size, read_len)
            u_regions.append({
                'region_idx': len(u_regions) + 1,
                'start_bp': bp_start,
                'end_bp': bp_end,
                'seq': read_seq[bp_start:bp_end],
                'window_indices': list(range(region_start_idx, prev_idx + 1))
            })
            region_start_idx = idx
            prev_idx = idx

    bp_start = region_start_idx * window_size
    bp_end = min((prev_idx + 1) * window_size, read_len)
    u_regions.append({
        'region_idx': len(u_regions) + 1,
        'start_bp': bp_start,
        'end_bp': bp_end,
        'seq': read_seq[bp_start:bp_end],
        'window_indices': list(range(region_start_idx, prev_idx + 1))
    })

    return masked_seq, u_regions


def run_blat_gfclient(query_fasta, output_psl, host, port, genome_dir,
                      min_score=20, min_identity=70):
    """
    Run gfClient against a running gfServer and return parsed PSL hits.

    Parameters
    ----------
    query_fasta : path to FASTA file with sequences to query
    output_psl  : path to write PSL output
    host        : hostname where gfServer is running
    port        : port number (int or str)
    genome_dir  : directory on the server containing the .2bit file
    min_score   : minimum BLAT score to report
    min_identity: minimum % identity to report

    Returns list of dicts sorted by score descending:
        {q_name, chrom, strand, t_start, t_end, matches, mismatches,
         percent_identity, q_size, score, block_count}
    """
    cmd = [
        'gfClient',
        host,
        str(port),
        genome_dir,
        query_fasta,
        output_psl,
        '-t=dna',
        '-q=dna',
        f'-minScore={min_score}',
        f'-minIdentity={min_identity}',
        '-nohead'
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        # gfClient writes "gfClient error exit" to stderr on failure but may
        # still write a valid (empty) PSL — check for hard errors
        if 'error' in result.stderr.lower() and 'mustOpen' in result.stderr:
            print(f"  gfClient error: {result.stderr.strip()}")
            return []
    except FileNotFoundError:
        print("  Error: gfClient not found in PATH")
        return []

    # Parse PSL output
    hits = []
    if not os.path.exists(output_psl):
        return hits

    with open(output_psl) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('psLayout') \
                    or line.startswith('match') or line.startswith('---'):
                continue
            fields = line.split('\t')
            if len(fields) < 21:
                continue
            try:
                matches    = int(fields[0])
                mismatches = int(fields[1])
                aln_len    = matches + mismatches
                pct_id     = (matches / aln_len * 100) if aln_len > 0 else 0.0
                score      = matches - mismatches
                t_start = int(fields[15])
                t_end   = int(fields[16])
                hits.append({
                    'q_name':           fields[9],
                    'chrom':            fields[13],
                    'strand':           fields[8],
                    't_start':          t_start,
                    't_end':            t_end,
                    'span':             t_end - t_start,
                    'matches':          matches,
                    'mismatches':       mismatches,
                    'percent_identity': round(pct_id, 2),
                    'q_size':           int(fields[10]),
                    'score':            score,
                    'block_count':      int(fields[17])
                })
            except (ValueError, IndexError):
                continue

    hits.sort(key=lambda x: -x['score'])
    return hits



def main():
    parser = argparse.ArgumentParser(
        description='Classify repeat structure in long reads using minimap2',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                       help='Input FASTA file with unassigned reads')
    parser.add_argument('-u', '--units', required=True,
                       help='FASTA file with unit sequences')
    parser.add_argument('-o', '--output', required=True,
                       help='Output prefix for results files')
    parser.add_argument('-w', '--window-size', type=int, default=1000,
                       help='Window size in base pairs')
    parser.add_argument('-t', '--identity-threshold', type=float, default=40.0,
                       help='Minimum percent identity for assignment (below this = U)')
    parser.add_argument('-c', '--uppercase-threshold', type=float, default=0.3,
                       help='Fraction of unit length for uppercase assignment (0.3 = first 30%%)')
    parser.add_argument('-p', '--preset', default='map-hifi',
                       choices=['map-ont', 'map-hifi', 'map-pb', 'asm20', 'sr'],
                       help='Minimap2 preset')

    # ── BLAT / U-region options ────────────────────────────────────────────
    blat_group = parser.add_argument_group('BLAT / U-region analysis')
    blat_group.add_argument('--blat', action='store_true',
                       help='Run local gfClient BLAT on internal U regions')
    blat_group.add_argument('--blat-host', default=None,
                       help='Hostname where gfServer is running (required with --blat)')
    blat_group.add_argument('--blat-port', type=int, default=7777,
                       help='gfServer port (default: 7777)')
    blat_group.add_argument('--blat-genome-dir', default=None,
                       help='Directory on server containing the .2bit file (required with --blat)')
    blat_group.add_argument('--blat-min-identity', type=float, default=70.0,
                       help='Min %% identity for BLAT hits (default: 70.0)')
    blat_group.add_argument('--blat-min-score', type=int, default=20,
                       help='Min score for BLAT hits (default: 20)')
    blat_group.add_argument('--blat-top-hits', type=int, default=3,
                       help='Number of top BLAT hits to report per U region (default: 3)')

    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} not found")
        return 1
    if not os.path.exists(args.units):
        print(f"Error: Units file {args.units} not found")
        return 1
    if args.blat:
        if not args.blat_host:
            print("Error: --blat requires --blat-host (hostname where gfServer is running)")
            print("  Start gfServer with: start_blat_server.lsf, then check blat_server_host.txt")
            return 1
        if not args.blat_genome_dir:
            print("Error: --blat requires --blat-genome-dir (directory containing .2bit file)")
            return 1
    
    # Parse input files
    print("Reading input sequences...")
    reads = parse_fasta(args.input)
    units = parse_fasta(args.units)
    
    print(f"Found {len(reads)} reads to process")
    print(f"Found {len(units)} unit sequences: {', '.join(units.keys())}")
    
    unit_names = sorted(units.keys())
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Define output file paths
        structure_file = f"{args.output}_structures.txt"
        scores_file    = f"{args.output}_detailed_scores.txt"
        summary_file   = f"{args.output}_structure_summary.txt"
        u_context_file = f"{args.output}_u_context.fasta"
        u_regions_file = f"{args.output}_u_regions.fasta"
        blat_hits_file = f"{args.output}_u_blat_hits.tsv"

        structure_counts    = defaultdict(int)
        strand_distribution = defaultdict(int)
        reads_with_internal_u = 0

        # Collect structures in memory so we can annotate with BLAT hits later
        # {read_name: (read_strand, structure, simplified)}
        structures_cache = {}

        with open(scores_file, 'w') as df, \
             open(u_context_file, 'w') as ucf, \
             open(u_regions_file, 'w') as urf:

            header_parts = ['read_name', 'window_num', 'window_start', 'window_end', 'strand']
            for unit_name in unit_names:
                header_parts.extend([
                    f'{unit_name}_identity',
                    f'{unit_name}_ref_start',
                    f'{unit_name}_ref_end',
                    f'{unit_name}_strand'
                ])
            header_parts.append('assigned_unit')
            df.write('\t'.join(header_parts) + '\n')

            # Process each read
            for idx, (read_name, read_seq) in enumerate(reads.items(), 1):
                print(f"Processing read {idx}/{len(reads)}: {read_name}")

                structure, detailed, read_strand = process_read(
                    read_name, read_seq, units, unit_names, temp_dir,
                    args.window_size, args.identity_threshold,
                    args.uppercase_threshold, args.preset
                )

                simplified = interpret_structure(structure, read_strand)
                structures_cache[read_name] = (read_strand, structure, simplified)

                structure_counts[simplified] += 1
                strand_distribution[read_strand] += 1

                # Write detailed scores
                for score_entry in detailed:
                    row_parts = [
                        score_entry['read_name'],
                        str(score_entry['window_num']),
                        str(score_entry['window_start']),
                        str(score_entry['window_end']),
                        str(score_entry['strand'])
                    ]
                    for unit_name in unit_names:
                        row_parts.extend([
                            f"{score_entry[f'{unit_name}_identity']:.2f}",
                            str(score_entry[f'{unit_name}_ref_start']),
                            str(score_entry[f'{unit_name}_ref_end']),
                            str(score_entry[f'{unit_name}_strand'])
                        ])
                    row_parts.append(score_entry['assigned_unit'])
                    df.write('\t'.join(row_parts) + '\n')

                # ── U-region analysis ──────────────────────────────────────
                structure_list = structure.split('-')

                if has_internal_u(structure_list):
                    reads_with_internal_u += 1
                    masked_seq, u_regions = get_u_regions(
                        read_name, read_seq, structure_list, args.window_size
                    )

                    # Context FASTA: full read with U windows in uppercase
                    ucf.write(f">{read_name} strand={read_strand} "
                              f"structure={simplified}\n")
                    for i in range(0, len(masked_seq), 80):
                        ucf.write(masked_seq[i:i+80] + '\n')

                    # U-regions FASTA: one entry per contiguous U stretch
                    for region in u_regions:
                        region_header = (
                            f"{read_name}_Uregion{region['region_idx']}"
                            f"_bp{region['start_bp']}-{region['end_bp']}"
                            f"_windows{region['window_indices'][0]+1}"
                            f"-{region['window_indices'][-1]+1}"
                        )
                        urf.write(f">{region_header}\n")
                        seq = region['seq']
                        for i in range(0, len(seq), 80):
                            urf.write(seq[i:i+80] + '\n')

        # ── BLAT all U regions ─────────────────────────────────────────────
        # {read_name: "chr13:16682663-16689663;chr5:1000-2000"} — top hit per region
        blat_hit_by_read = {}

        if args.blat:
            with open(u_regions_file) as f:
                n_regions = sum(1 for l in f if l.startswith('>'))

            if n_regions == 0:
                print("\nNo internal U regions found — nothing to BLAT.")
            else:
                print(f"\nRunning local BLAT on {n_regions} U region(s) "
                      f"via gfClient ({args.blat_host}:{args.blat_port})...")
                psl_out = f"{args.output}_u_regions.psl"
                hits = run_blat_gfclient(
                    u_regions_file, psl_out,
                    host=args.blat_host,
                    port=args.blat_port,
                    genome_dir=args.blat_genome_dir,
                    min_score=args.blat_min_score,
                    min_identity=args.blat_min_identity
                )

                # Group hits by query name
                from collections import OrderedDict
                hits_by_query = OrderedDict()
                for hit in hits:
                    hits_by_query.setdefault(hit['q_name'], []).append(hit)

                # Build per-read top-hit string for structures.txt
                # q_name format: {read_name}_Uregion{n}_bp{start}-{end}_windows{a}-{b}
                for q_name, q_hits in hits_by_query.items():
                    # Extract original read name (everything before _Uregion)
                    read_name = q_name.rsplit('_Uregion', 1)[0]
                    top = q_hits[0]
                    loc = f"{top['chrom']}:{top['t_start']}-{top['t_end']}"
                    if read_name not in blat_hit_by_read:
                        blat_hit_by_read[read_name] = []
                    blat_hit_by_read[read_name].append(loc)

                # Write detailed BLAT hits file
                all_query_names = set()
                with open(u_regions_file) as f:
                    for line in f:
                        if line.startswith('>'):
                            all_query_names.add(line[1:].split()[0])

                with open(blat_hits_file, 'w') as bhf:
                    bhf.write("region_id\trank\tchrom\tstrand\t"
                              "t_start\tt_end\tmatches\tmismatches\t"
                              "percent_identity\tq_size\tspan\tscore\tblock_count\n")
                    for q_name, q_hits in hits_by_query.items():
                        for rank, hit in enumerate(q_hits[:args.blat_top_hits], 1):
                            bhf.write(
                                f"{hit['q_name']}\t{rank}\t"
                                f"{hit['chrom']}\t{hit['strand']}\t"
                                f"{hit['t_start']}\t{hit['t_end']}\t"
                                f"{hit['matches']}\t{hit['mismatches']}\t"
                                f"{hit['percent_identity']}\t"
                                f"{hit['q_size']}\t{hit['span']}\t"
                                f"{hit['score']}\t{hit['block_count']}\n"
                            )
                    for q_name in sorted(all_query_names - set(hits_by_query.keys())):
                        bhf.write(f"{q_name}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n")

                print(f"  PSL output:  {psl_out}")
                print(f"  BLAT hits:   {blat_hits_file} "
                      f"({len(hits_by_query)} regions with hits)")

        # ── Write structures file (with BLAT column if available) ──────────
        with open(structure_file, 'w') as sf:
            if args.blat:
                sf.write("read_name\tstrand\tfull_structure\t"
                         "simplified_structure\tu_blat_top_hits\n")
            else:
                sf.write("read_name\tstrand\tfull_structure\t"
                         "simplified_structure\n")
            for read_name, (read_strand, structure, simplified) in structures_cache.items():
                if args.blat:
                    hits_str = ';'.join(blat_hit_by_read.get(read_name, ['.']))
                    sf.write(f"{read_name}\t{read_strand}\t{structure}\t"
                             f"{simplified}\t{hits_str}\n")
                else:
                    sf.write(f"{read_name}\t{read_strand}\t{structure}\t"
                             f"{simplified}\n")

        # Generate structure summary
        total_reads = len(reads)
        with open(summary_file, 'w') as sumf:
            sumf.write("# Structure Summary\n")
            sumf.write(f"# Total reads processed: {total_reads}\n")
            sumf.write(f"# Window size: {args.window_size} bp\n")
            sumf.write(f"# Identity threshold: {args.identity_threshold}%\n")
            sumf.write(f"# Reads with internal U regions: {reads_with_internal_u}\n")
            sumf.write("#\n")
            sumf.write("# Strand Distribution:\n")
            for strand in sorted(strand_distribution.keys()):
                count = strand_distribution[strand]
                pct = (count / total_reads) * 100
                sumf.write(f"#   {strand}: {count} reads ({pct:.1f}%)\n")
            sumf.write("#\n")
            sumf.write("simplified_structure\tcount\tpercentage\n")
            sorted_structures = sorted(
                structure_counts.items(),
                key=lambda x: (-x[1], x[0])
            )
            for structure, count in sorted_structures:
                percentage = (count / total_reads) * 100
                sumf.write(f"{structure}\t{count}\t{percentage:.2f}\n")

    print(f"\nComplete! Results written to:")
    print(f"  Structures:        {structure_file}")
    print(f"  Detailed scores:   {scores_file}")
    print(f"  Structure summary: {summary_file}")
    print(f"  U context FASTA:   {u_context_file}  ({reads_with_internal_u} reads)")
    print(f"  U regions FASTA:   {u_regions_file}")
    if args.blat:
        print(f"  BLAT hits:         {blat_hits_file}")
        print(f"  Raw PSL:           {args.output}_u_regions.psl")
    print(f"\nStructure Summary:")
    print(f"  Total reads:           {total_reads}")
    print(f"  Reads with internal U: {reads_with_internal_u}")
    print(f"  Unique structures:     {len(structure_counts)}")
    if sorted_structures:
        top_structure, top_count = sorted_structures[0]
        top_pct = (top_count / total_reads) * 100
        print(f"  Most common: {top_structure} ({top_count} reads, {top_pct:.1f}%)")

    return 0


if __name__ == '__main__':
    exit(main())
