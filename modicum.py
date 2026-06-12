import argparse
import sys
import os
import gzip
import csv
import re
from collections import defaultdict, Counter
from Bio import SeqIO

try:
    import pyhmmer
except ImportError:
    sys.exit("Error: pyhmmer not found. Install: pip install pyhmmer")

try:
    import matplotlib.pyplot as plt
    from sankeyflow import Sankey
    SANKEY_AVAILABLE = True
except ImportError:
    SANKEY_AVAILABLE = False

def load_hmm_mapping(tsv_path):
    """
    Loads a mapping file with support for optional columns:
    Format: HMM_ID \t Gene \t Pathway \t Type (Marker/Accessory) [\t Priority] [\t HMM_Variant] [\t New_Pathway]
    """
    hmm_map = defaultdict(list)
    pathway_accessories = defaultdict(set)
    child_to_parent = {}
    
    if not os.path.exists(tsv_path):
        sys.exit(f"[!] Error: Mapping file not found: {tsv_path}")
        
    with open(tsv_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        
        # Dynamically determine indices based on headers, falling back to positional defaults
        idx_hmm, idx_gene, idx_pathway, idx_type = 0, 1, 2, 3
        idx_priority, idx_variant, idx_new_pathway = -1, -1, -1
        
        if header:
            for idx, col in enumerate(header):
                c = col.strip().lower()
                if c in ['hmm_id', 'hmm']: idx_hmm = idx
                elif c == 'gene': idx_gene = idx
                elif c == 'pathway': idx_pathway = idx
                elif c == 'type': idx_type = idx
                elif c == 'priority': idx_priority = idx
                elif c == 'hmm_variant': idx_variant = idx
                elif c == 'new_pathway': idx_new_pathway = idx
                
        for row in reader:
            if not row or row[0].startswith('#'): continue
            if len(row) > max(idx_hmm, idx_gene, idx_pathway, idx_type):
                hmm_id = row[idx_hmm].strip()
                gene = row[idx_gene].strip()
                pathway = row[idx_pathway].strip()
                gene_type = row[idx_type].strip()
                
                # Priority configuration
                priority = 1
                if idx_priority != -1 and idx_priority < len(row):
                    priority_val = row[idx_priority].strip()
                    if priority_val != "":
                        try:
                            priority = int(priority_val)
                        except ValueError:
                            priority = 1
                
                # Variant and New Pathway rule parsing
                variant_rules = []
                if idx_variant != -1 and idx_new_pathway != -1 and idx_variant < len(row) and idx_new_pathway < len(row):
                    var_str = row[idx_variant].strip()
                    np_str = row[idx_new_pathway].strip()
                    if var_str and np_str:
                        var_parts = [v.strip() for v in var_str.split(';')]
                        np_parts = [n.strip() for n in np_str.split(';')]
                        
                        for v_part, n_part in zip(var_parts, np_parts):
                            if v_part and n_part:
                                conds = v_part.split('&')
                                valid_rule = True
                                parsed_conds = []
                                for cond in conds:
                                    cond = cond.strip()
                                    m = re.match(r"^(\d+)(.+)$", cond)
                                    if m:
                                        p_num = int(m.group(1))
                                        res_options = set(m.group(2).split('/'))
                                        parsed_conds.append((p_num, res_options))
                                    else:
                                        valid_rule = False
                                if valid_rule:
                                    # store (conditions, target_pathway, specificity_score)
                                    variant_rules.append((parsed_conds, n_part, len(parsed_conds)))
                                    if n_part.lower() != 'canonical':
                                        child_to_parent[n_part] = pathway
                
                hmm_map[hmm_id].append((gene, pathway, gene_type, priority, variant_rules))
                
                if gene_type.lower() in ['accessory', 'acc']:
                    pathway_accessories[pathway].add(gene)
                elif pathway not in pathway_accessories:
                    pathway_accessories[pathway] = set()
                
    total_acc_per_pathway = {p: len(genes) for p, genes in pathway_accessories.items()}
    for child, parent in child_to_parent.items():
        if parent in total_acc_per_pathway:
            total_acc_per_pathway[child] = total_acc_per_pathway[parent]
    print(f"[*] Loaded mappings for {len(hmm_map)} unique HMM profiles from {tsv_path}")
    return hmm_map, total_acc_per_pathway, child_to_parent

def count_sequences(fasta_path):
    """Counts the number of sequences (lines starting with '>') in a FASTA/FAA file in a memory-efficient way."""
    is_gzip = False
    try:
        with open(fasta_path, 'rb') as f:
            magic = f.read(2)
            is_gzip = (magic == b'\x1f\x8b')
    except Exception:
        pass
    
    open_func = gzip.open if is_gzip else open
    count = 0
    with open_func(fasta_path, 'rb') as f:
        found = False
        last_char = b''
        while True:
            chunk = f.read(65536)
            if not chunk:
                break
            idx = chunk.find(b'>')
            if idx != -1:
                count = 1
                found = True
                count += chunk[idx+1:].count(b'\n>')
                last_char = chunk[-1:]
                break
        
        if found:
            while True:
                chunk = f.read(65536)
                if not chunk:
                    break
                if last_char == b'\n' and chunk.startswith(b'>'):
                    count += 1
                count += chunk.count(b'\n>')
                last_char = chunk[-1:]
    return count

def parse_fasta_into_genomes(fasta_path, all_hmm_hits):
    """Groups sequences by parent genome/phage ID based on FASTA headers, storing only sequences with HMM hits in memory."""
    phage_groups = defaultdict(list)
    all_phages = []
    seen_phages = set()
    total_seqs = 0
    
    is_gzip = False
    try:
        with open(fasta_path, 'rb') as f:
            magic = f.read(2)
            is_gzip = (magic == b'\x1f\x8b')
    except Exception:
        pass
        
    open_func = gzip.open if is_gzip else open
    with open_func(fasta_path, "rt", encoding="utf-8", errors="ignore") as f_in:
        for seq in SeqIO.parse(f_in, "fasta"):
            total_seqs += 1
            parts = seq.id.rsplit('_', 1)
            if len(parts) > 1 and parts[1].isdigit():
                phage_id = parts[0]
            else:
                phage_id = os.path.basename(fasta_path).rsplit('.', 1)[0]
                
            if phage_id not in seen_phages:
                seen_phages.add(phage_id)
                all_phages.append(phage_id)
                
            if seq.id in all_hmm_hits:
                phage_groups[phage_id].append(seq)
    
    print(f"[*] Parsed {total_seqs} sequences across {len(all_phages)} genomes.")
    return phage_groups, all_phages

def run_pyhmmer_scan(fasta_path, hmm_db_path, min_coverage=0.8):
    """Runs pyhmmer and extracts aligned residues per HMM match state coordinate position."""
    hits_dict = defaultdict(list)
    
    total_seqs = count_sequences(fasta_path)
    if total_seqs == 0:
        return hits_dict
        
    with pyhmmer.plan7.HMMFile(hmm_db_path) as hmm_file:
        hmms = list(hmm_file)
        
    with pyhmmer.easel.SequenceFile(fasta_path, digital=True) as seq_file:
        while True:
            seq_block = seq_file.read_block(sequences=100000)
            if not seq_block:
                break
                
            for hmm, top_hits in zip(hmms, pyhmmer.hmmsearch(hmms, seq_block, Z=total_seqs, domZ=None)):
                raw_id = hmm.accession if hmm.accession else hmm.name
                hmm_id_base = raw_id.split('.')[0]
                hmm_length = hmm.M 
                
                for hit in top_hits:
                    if hit.included: 
                        covered_positions = set()
                        hmm_residues = {}
                        
                        for domain in hit.domains:
                            if domain.included:
                                covered_positions.update(range(domain.alignment.hmm_from, domain.alignment.hmm_to + 1))
                                
                                # Map HMM positions directly to aligned target residues
                                position = domain.alignment.hmm_from
                                for hmm_letter, amino_acid in zip(domain.alignment.hmm_sequence, domain.alignment.target_sequence):
                                    if hmm_letter != ".":
                                        hmm_residues[position] = amino_acid
                                        position += 1
                        
                        coverage = len(covered_positions) / hmm_length
                        
                        if coverage >= min_coverage:
                            seq_id = hit.name
                            hits_dict[seq_id].append((hmm_id_base, hit.evalue, hit.score, coverage, hmm_residues))
                            
    return hits_dict

def analyze_single_phage(phage_id, proteins, all_hmm_hits, active_hmm_map, total_acc_per_pathway, child_to_parent):
    """Analyzes hits, updates pathways based on residue variants, and resolves priority contentions."""
    hit_pathways = defaultdict(lambda: defaultdict(set))
    gene_to_type = {}
    pathway_marker_priorities = defaultdict(int)
    annotated_seqs = []
    gene_hmm_seqs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    
    for seq in proteins:
        seq_id = seq.id
        if seq_id in all_hmm_hits:
            hit_annotations = set()
            
            # Sort HMM hits by bitscore in descending order (highest score first)
            sorted_hits = sorted(all_hmm_hits[seq_id], key=lambda x: x[2], reverse=True)
            selected_hit = None
            
            for hit in sorted_hits:
                hit_hmm_id, evalue, bitscore, coverage, hmm_residues = hit
                if hit_hmm_id in active_hmm_map:
                    # Check if this HMM hit has gaps in its variant positions
                    has_gap = False
                    for gene, pathway, gene_type, priority, variant_rules in active_hmm_map[hit_hmm_id]:
                        if variant_rules:
                            for parsed_conds, _, _ in variant_rules:
                                for pos, _ in parsed_conds:
                                    if pos not in hmm_residues or hmm_residues[pos] == '-':
                                        has_gap = True
                                        break
                                if has_gap:
                                    break
                        if has_gap:
                            break
                    
                    if not has_gap:
                        selected_hit = hit
                        break
            
            # Fallback to the top hit in the map if all hits had gaps
            if selected_hit is None:
                for hit in sorted_hits:
                    if hit[0] in active_hmm_map:
                        selected_hit = hit
                        break
            
            if selected_hit:
                hit_hmm_id, evalue, bitscore, coverage, hmm_residues = selected_hit
                for gene, pathway, gene_type, priority, variant_rules in active_hmm_map[hit_hmm_id]:
                    
                    # Evaluate variant assignments
                    assigned_pathway = pathway
                    if variant_rules:
                        matching_rules = []
                        for parsed_conds, new_pathway, spec_score in variant_rules:
                            match_all = True
                            for pos, allowed in parsed_conds:
                                if pos in hmm_residues and hmm_residues[pos] in allowed:
                                    continue
                                else:
                                    match_all = False
                                    break
                            if match_all:
                                matching_rules.append((spec_score, new_pathway))
                        
                        if matching_rules:
                            # Prioritize the rule with the highest specificity score.
                            # If specificity score is equal, prioritize 'canonical' to correctly skip wild-types.
                            matching_rules.sort(key=lambda x: (x[0], x[1].strip().lower() == 'canonical'), reverse=True)
                            assigned_pathway = matching_rules[0][1]
                        
                    if assigned_pathway.strip().lower() == "canonical":
                        continue
                        
                    hit_pathways[assigned_pathway][gene].add(seq_id)
                    gene_hmm_seqs[assigned_pathway][gene][hit_hmm_id].add(seq_id)
                    gene_to_type[gene] = gene_type
                    
                    # Monitor marker priorities per pathway
                    if gene_type.lower() in ['marker', 'primary']:
                        if assigned_pathway not in pathway_marker_priorities:
                            pathway_marker_priorities[assigned_pathway] = priority
                        else:
                            if priority > pathway_marker_priorities[assigned_pathway]:
                                pathway_marker_priorities[assigned_pathway] = priority
                    
                    metric_string = f"E={evalue:.1e}|S={bitscore:.1f}|C={coverage:.2f}"
                    hit_annotations.add(f"{assigned_pathway}|{gene}|{gene_type}|{metric_string}")
            
            if hit_annotations:
                seq.description = f"Hits: {', '.join(sorted(list(hit_annotations)))}"
                annotated_seqs.append(seq)

    # Copy accessory genes from parent to children if children have hits
    for child, parent in child_to_parent.items():
        if child in hit_pathways and parent in hit_pathways:
            for gene, seq_ids in hit_pathways[parent].items():
                if gene_to_type.get(gene, '').lower() in ['accessory', 'acc']:
                    hit_pathways[child][gene].update(seq_ids)
                    for hmm, seqs in gene_hmm_seqs[parent][gene].items():
                        gene_hmm_seqs[child][gene][hmm].update(seqs)

    best_prediction = "Unknown"
    status = "None"
    evidence_strings = []
    gene_details_output = {}
    hmm_names_output = {}
    pathway_statuses = {} 
    
    if hit_pathways:
        predictions = []
        statuses = []
        
        # Compile pathways triggering traditional marker classifications
        strong_pathways = []
        for pathway, genes_dict in hit_pathways.items():
            markers_hit = sorted([g for g in genes_dict.keys() if gene_to_type[g].lower() in ['marker', 'primary']])
            if markers_hit:
                prio = pathway_marker_priorities.get(pathway, 1)
                strong_pathways.append((pathway, prio, markers_hit))
        
        # Apply priority filters
        final_strong_pathways = []
        promoted_to_putative_from_priority_0 = []
        
        if strong_pathways:
            max_prio = max(p for _, p, _ in strong_pathways)
            if max_prio > 0:
                for b_pathway, p_val, m_hit in strong_pathways:
                    if p_val == max_prio:
                        final_strong_pathways.append((b_pathway, m_hit))
            else:
                # If only priority 0 pathways exist, downgrade designation to Putative
                for b_pathway, p_val, m_hit in strong_pathways:
                    promoted_to_putative_from_priority_0.append((b_pathway, m_hit))
        
        pathway_results = {}
        for pathway, genes_dict in hit_pathways.items():
            markers_hit = sorted([g for g in genes_dict.keys() if gene_to_type[g].lower() in ['marker', 'primary']])
            accessories_hit = sorted([g for g in genes_dict.keys() if gene_to_type[g].lower() in ['accessory', 'acc']])
            total_acc = total_acc_per_pathway.get(pathway, 0)
            
            gene_counts = {g: len(seq_ids) for g, seq_ids in genes_dict.items()}
            gene_details_output[pathway] = [f"{g}({c})" for g, c in sorted(gene_counts.items())]
            
            hmm_parts = []
            for g, c in sorted(gene_counts.items()):
                hmm_seqs_dict = gene_hmm_seqs[pathway][g]
                hmm_subparts = [f"{hmm}({len(seqs)})" for hmm, seqs in sorted(hmm_seqs_dict.items())]
                hmm_parts.append(",".join(hmm_subparts))
            hmm_names_output[pathway] = hmm_parts
            
            is_strong = any(pathway == sp[0] for sp in final_strong_pathways)
            is_prio_0_putative = any(pathway == pp[0] for pp in promoted_to_putative_from_priority_0)
            
            p_status = "None"
            p_pred = None
            p_evidence = None
            
            if is_strong:
                p_status = "Strong"
                p_pred = pathway
                p_evidence = f"Primary: {','.join(markers_hit)} | Acc: {len(accessories_hit)}/{total_acc}"
            elif is_prio_0_putative:
                p_status = "Putative"
                p_pred = f"{pathway}?"
                p_evidence = f"Primary: {','.join(markers_hit)} | Acc: {len(accessories_hit)}/{total_acc}"
            else:
                if not markers_hit:
                    if total_acc > 1 and len(accessories_hit) > (total_acc / 2.0):
                        p_status = "Putative"
                        p_pred = f"{pathway}?"
                        p_evidence = f"Primary: None | Acc: {len(accessories_hit)}/{total_acc}"
            
            if p_status != "None":
                pathway_results[pathway] = {
                    'status': p_status,
                    'prediction': p_pred,
                    'evidence': p_evidence
                }
                
        # Suppress parent pathway prediction if a child pathway is predicted
        for child, parent in child_to_parent.items():
            if child in pathway_results:
                if parent in pathway_results:
                    del pathway_results[parent]
                    
        # Now rebuild predictions, statuses, evidence_strings, pathway_statuses
        predictions = []
        statuses = []
        evidence_strings = []
        
        for pathway in sorted(pathway_results.keys()):
            p_info = pathway_results[pathway]
            predictions.append(p_info['prediction'])
            statuses.append(p_info['status'])
            evidence_strings.append(p_info['evidence'])
            pathway_statuses[pathway] = p_info['status']

        if "Strong" in statuses:
            status = "Strong"
        elif "Putative" in statuses:
            status = "Putative"
        elif not statuses:
            status = "None"
        else:
            status = "Unknown"

        best_prediction = " / ".join(predictions) if predictions else "Unknown"

    return {
        'PhageID': phage_id,
        'Prediction': best_prediction,
        'Status': status,
        'Evidence': " || ".join(evidence_strings),
        'Gene_Details': gene_details_output,
        'HMM_Names': hmm_names_output,
        'Pathway_Statuses': pathway_statuses,
        'Annotated_Seqs': annotated_seqs
    }

def build_sankey_plot(final_results, output_file, child_to_parent):
    if not SANKEY_AVAILABLE:
        print("[!] Matplotlib or SankeyFlow not found. Skipping chart generation.")
        return
        
    raw_flow_counts = Counter()
    for res in final_results:
        if res['Status'] == "None": 
            continue
        for pathway, p_status in res.get('Pathway_Statuses', {}).items():
            if pathway in child_to_parent:
                parent = child_to_parent[pathway]
                raw_flow_counts[(p_status, parent)] += 1
                raw_flow_counts[(parent, pathway)] += 1
            else:
                raw_flow_counts[(p_status, pathway)] += 1

    if not raw_flow_counts:
        print("[!] No confident or putative predictions to chart.")
        return

    flows = []
    for (src, des), count in raw_flow_counts.items():
        flows.append((src, des, count))

    plt.style.use('ggplot')

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica', 'Liberation Sans', 'sans-serif']
    plt.rcParams['text.color'] = '#2D3748'

    fig, ax = plt.subplots(figsize=(10, 6), dpi=600)
    fig.patch.set_facecolor('#F7FAFC')
    ax.set_facecolor('#F7FAFC')

    from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_hex
    
    all_nodes = set()
    for src, des, count in flows:
        all_nodes.add(src)
        all_nodes.add(des)
    num_nodes = max(len(all_nodes), 1)
    
    base_colors = [
        '#1b4965', '#e2b13c', '#457b9d', '#e9c46a', 
        '#1d3557', '#f2b963', '#62b6cb', '#fcc647', 
        '#91b3d5', '#a8dadc', '#bee9e8'
    ]
    
    if num_nodes <= len(base_colors):
        colors_list = base_colors[:num_nodes]
    else:
        cmap_temp = LinearSegmentedColormap.from_list('blue_yellow_custom', base_colors, N=num_nodes)
        colors_list = [to_hex(cmap_temp(i / (num_nodes - 1))) for i in range(num_nodes)]
        
    custom_cmap = ListedColormap(colors_list)

    s = Sankey(
        flows=flows,
        align_y='tree clamp', # Tree clamp alignment groups node structures nicely
        flow_color_mode='source', # Blends ribbons automatically based on originating node
        flow_color_mode_alpha=0.45, # Soft flow opacity for clean layering
        cmap=custom_cmap,
        node_width=0.025,
        node_opts={
            'edgecolor': 'white', # White border separates node from flows
            'linewidth': 1.5,
            'label_opts': {'fontsize': 12, 'weight': 'normal', 'color': '#4A5568'}
        }
    )

    # Adjust node label positions to prevent overlapping with flow ribbons
    for level, node_level in enumerate(s.nodes):
        for node in node_level:
            node.label_format = '{label} ({value:,.0f})'
            if level == 0:
                node.label_pos = 'left'
                node.label_pad_x = 0.015
            else:
                node.label_pos = 'right'
                node.label_pad_x = 0.015

    s.draw(ax=ax)

    unknown_count = sum(1 for res in final_results if res.get('Prediction') == "Unknown")
    ax.set_title("MODICUM pathway predictions", 
                 fontsize=15, pad=15, fontweight='normal', color='#1A202C', loc='left')
    ax.text(0.5, -0.05, f"Genomes with no MODICUM classification: {unknown_count}", 
            transform=ax.transAxes, ha='center', fontsize=14, color='#4A5568', fontweight='normal')

    # Remove axis spines and ticks while keeping the background color panel
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.savefig(output_file, format='svg', bbox_inches='tight', facecolor=fig.get_facecolor())
    plt.close()
    print(f"[*] Sankey plot written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MODICUM: Dynamic HMM Pathway Predictor (designed originally for DNA modifications)")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA of phage proteins")
    parser.add_argument("-d", "--database", required=True, help="HMM database (.hmm file)")
    parser.add_argument("-m", "--map", required=True, help="4-to-7 column TSV mapping file")
    parser.add_argument("-o", "--output", default="modicum_results", help="Output prefix")
    parser.add_argument("-c", "--coverage", type=float, default=0.8, help="Minimum HMM coverage threshold (0.0 to 1.0). Default: 0.8")
    args = parser.parse_args()

    active_hmm_map, total_acc_per_pathway, child_to_parent = load_hmm_mapping(args.map)

    print(f"[*] Scanning {args.input} against {args.database} with >= {args.coverage*100}% coverage threshold...")
    all_hmm_hits = run_pyhmmer_scan(args.input, args.database, min_coverage=args.coverage)

    phage_groups, all_phages = parse_fasta_into_genomes(args.input, all_hmm_hits)
    
    final_results = []
    all_hit_sequences = []
    
    for phage_id in all_phages:
        proteins = phage_groups[phage_id]
        result = analyze_single_phage(phage_id, proteins, all_hmm_hits, active_hmm_map, total_acc_per_pathway, child_to_parent)
        final_results.append(result)
        all_hit_sequences.extend(result['Annotated_Seqs'])

    tsv_file = f"{args.output}.tsv"
    with open(tsv_file, 'w', newline='') as tsvfile:
        fieldnames = ['PhageID', 'Prediction', 'Status', 'Evidence', 'Gene_Details', 'HMM_Names']
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        
        for row in final_results:
            row_copy = row.copy()
            formatted_details = []
            for path, gene_counts in row['Gene_Details'].items():
                formatted_details.extend(gene_counts)
            
            seen_details = set()
            deduped_details = []
            for item in formatted_details:
                if item not in seen_details:
                    seen_details.add(item)
                    deduped_details.append(item)
            row_copy['Gene_Details'] = "; ".join(deduped_details)
            
            formatted_hmms = []
            for path, hmm_strings in row.get('HMM_Names', {}).items():
                formatted_hmms.extend(hmm_strings)
            
            seen_hmms = set()
            deduped_hmms = []
            for item in formatted_hmms:
                if item not in seen_hmms:
                    seen_hmms.add(item)
                    deduped_hmms.append(item)
            row_copy['HMM_Names'] = "; ".join(deduped_hmms)
            
            writer.writerow(row_copy)
            
    print(f"[*] Results written to: {tsv_file}")

    fasta_file = f"{args.output}_annotated_hits.fasta"
    if all_hit_sequences:
        SeqIO.write(all_hit_sequences, fasta_file, "fasta")
        print(f"[*] Annotated sequences written to: {fasta_file}")
    else:
        print("[!] No hits found, skipping FASTA generation.")

    if SANKEY_AVAILABLE and all_hit_sequences:
        build_sankey_plot(final_results, f"{args.output}_sankey.svg", child_to_parent)