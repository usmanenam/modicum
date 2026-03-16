import argparse
import sys
import os
import csv
from collections import defaultdict, Counter
from Bio import SeqIO

try:
    import pyhmmer
except ImportError:
    sys.exit("Error: pyhmmer not found. Install: pip install pyhmmer")

try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

def load_hmm_mapping(tsv_path):
    """
    Loads a 4-column TSV mapping file.
    Format: HMM_ID \t Gene \t Pathway \t Type (Marker/Accessory)
    """
    hmm_map = defaultdict(list)
    pathway_accessories = defaultdict(set)
    
    if not os.path.exists(tsv_path):
        sys.exit(f"[!] Error: Mapping file not found: {tsv_path}")
        
    with open(tsv_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None) # Skip header
        
        for row in reader:
            if not row or row[0].startswith('#'): continue
            if len(row) >= 4:
                hmm_id, gene, pathway, gene_type = row[0].strip(), row[1].strip(), row[2].strip(), row[3].strip()
                hmm_map[hmm_id].append((gene, pathway, gene_type))
                
                # Keep track of unique accessory genes for threshold logic
                if gene_type.lower() in ['accessory', 'acc']:
                    pathway_accessories[pathway].add(gene)
                elif pathway not in pathway_accessories:
                    pathway_accessories[pathway] = set()
                
    total_acc_per_pathway = {p: len(genes) for p, genes in pathway_accessories.items()}
    
    print(f"[*] Loaded mappings for {len(hmm_map)} unique HMM profiles from {tsv_path}")
    return hmm_map, total_acc_per_pathway

def parse_fasta_into_genomes(fasta_path):
    """Groups sequences by parent genome/phage ID based on FASTA headers."""
    phage_groups = defaultdict(list)
    for seq in SeqIO.parse(fasta_path, "fasta"):
        parts = seq.id.rsplit('_', 1)
        if len(parts) > 1 and parts[1].isdigit():
            phage_id = parts[0]
        else:
            phage_id = os.path.basename(fasta_path).rsplit('.', 1)[0]
        phage_groups[phage_id].append(seq)
    
    print(f"[*] Parsed {sum(len(v) for v in phage_groups.values())} sequences across {len(phage_groups)} genomes.")
    return phage_groups

def run_pyhmmer_scan(fasta_path, hmm_db_path, min_coverage=0.8):
    """Runs pyhmmer and returns a dict of hits per sequence ID, applying a coverage threshold."""
    hits_dict = defaultdict(list)
    
    with pyhmmer.easel.SequenceFile(fasta_path, digital=True) as seq_file:
        sequences = list(seq_file)
        
    with pyhmmer.plan7.HMMFile(hmm_db_path) as hmm_file:
        hmms = list(hmm_file)
        
    for hmm, top_hits in zip(hmms, pyhmmer.hmmsearch(hmms, sequences)):
        raw_id = hmm.accession.decode() if hmm.accession else hmm.name.decode()
        hmm_id_base = raw_id.split('.')[0]
        hmm_length = hmm.M 
        
        for hit in top_hits:
            if hit.included: 
                covered_positions = set()
                for domain in hit.domains:
                    if domain.included:
                        covered_positions.update(range(domain.alignment.hmm_from, domain.alignment.hmm_to + 1))
                
                coverage = len(covered_positions) / hmm_length
                
                if coverage >= min_coverage:
                    seq_id = hit.name.decode()
                    hits_dict[seq_id].append((hmm_id_base, hit.evalue, hit.score, coverage))
                    
    return hits_dict

def analyze_single_phage(phage_id, proteins, all_hmm_hits, active_hmm_map, total_acc_per_pathway):
    """Analyzes hits and deduplicates redundant HMM matches to the same gene."""
    hit_pathways = defaultdict(lambda: defaultdict(set))
    gene_to_type = {}
    annotated_seqs = []
    
    for seq in proteins:
        seq_id = seq.id
        if seq_id in all_hmm_hits:
            hit_annotations = set()
            
            for hit_hmm_id, evalue, bitscore, coverage in all_hmm_hits[seq_id]:
                if hit_hmm_id in active_hmm_map:
                    for gene, pathway, gene_type in active_hmm_map[hit_hmm_id]:
                        hit_pathways[pathway][gene].add(seq_id)
                        gene_to_type[gene] = gene_type
                        
                        metric_string = f"E={evalue:.1e}|S={bitscore:.1f}|C={coverage:.2f}"
                        hit_annotations.add(f"{pathway}|{gene}|{gene_type}|{metric_string}")
            
            if hit_annotations:
                seq.description = f"Hits: {', '.join(sorted(list(hit_annotations)))}"
                annotated_seqs.append(seq)

    best_prediction = "Unknown"
    status = "None"
    evidence_strings = []
    gene_details_output = {}
    pathway_statuses = {} # Track individual pathway statuses for Sankey plot
    
    if hit_pathways:
        predictions = []
        statuses = []
        
        for pathway, genes_dict in hit_pathways.items():
            markers_hit = sorted([g for g in genes_dict.keys() if gene_to_type[g].lower() in ['marker', 'primary']])
            accessories_hit = sorted([g for g in genes_dict.keys() if gene_to_type[g].lower() in ['accessory', 'acc']])
            
            total_acc = total_acc_per_pathway.get(pathway, 0)
            
            # ALWAYS log the genes hit into Gene_Details, even if they fail the threshold
            gene_counts = {g: len(seq_ids) for g, seq_ids in genes_dict.items()}
            gene_details_output[pathway] = [f"{g}({c})" for g, c in sorted(gene_counts.items())]
            
            if markers_hit:
                # Strong hit
                predictions.append(pathway)
                statuses.append("Strong")
                pathway_statuses[pathway] = "Strong"
                evidence_strings.append(f"Primary: {','.join(markers_hit)} | Acc: {len(accessories_hit)}/{total_acc}")
                
            else:
                # Putative check
                if total_acc > 1 and len(accessories_hit) > (total_acc / 2.0):
                    predictions.append(f"{pathway}?")
                    statuses.append("Putative")
                    pathway_statuses[pathway] = "Putative"
                    evidence_strings.append(f"Primary: None | Acc: {len(accessories_hit)}/{total_acc}")
                else:
                    # failed threshold: Do nothing for prediction/evidence/Sankey
                    pass

        if "Strong" in statuses:
            status = "Strong"
        elif "Putative" in statuses:
            status = "Putative"
        elif not statuses:
            status = "None" # Occurs if everything failed the threshold
        else:
            status = "Unknown"

        best_prediction = " / ".join(predictions) if predictions else "Unknown"

    return {
        'PhageID': phage_id,
        'Prediction': best_prediction,
        'Status': status,
        'Evidence': " || ".join(evidence_strings),
        'Gene_Details': gene_details_output,
        'Pathway_Statuses': pathway_statuses,
        'Annotated_Seqs': annotated_seqs
    }

def build_sankey_plot(final_results, output_file):
    if not PLOTLY_AVAILABLE: return
    
    labels = []
    source = []
    target = []
    value = []

    def get_or_add_label(name):
        if name not in labels:
            labels.append(name)
        return labels.index(name)

    # Count aggregated flows from Prediction Status -> Specific Pathway
    flow_counts = Counter()
    
    for res in final_results:
        if res['Status'] == "None": continue
        
        for pathway, p_status in res.get('Pathway_Statuses', {}).items():
            flow_counts[(p_status, pathway)] += 1

    if not flow_counts:
        print("[!] No confident or putative predictions to chart.")
        return

        # Build Sankey Links
    for (status_name, pathway_name), count in flow_counts.items():
        s_idx = get_or_add_label(status_name)
        t_idx = get_or_add_label(pathway_name)
        source.append(s_idx)
        target.append(t_idx)
        value.append(count)

    # Calculate totals for each node to display next to the label
    node_totals = {}
    for s, t, v in zip(source, target, value):
        node_totals[s] = node_totals.get(s, 0) + v
        node_totals[t] = node_totals.get(t, 0) + v
        
    updated_labels = [f"{name} ({node_totals.get(i, 0)})" for i, name in enumerate(labels)]

    fig = go.Figure(data=[go.Sankey(
        node = dict(pad=15, thickness=20, line=dict(color="black", width=0.5), label=updated_labels),
        link = dict(source=source, target=target, value=value)
    )])
    fig.update_layout(title_text="Modicum: Prediction Status to Pathway Breakdown", font_size=12)
     # Calculate unknowns
    unknown_count = sum(1 for res in final_results if res.get('Prediction') == "Unknown")
    
    # Generate HTML string and inject custom text at the bottom
    html_content = fig.to_html(full_html=True)
    custom_footer = f"<div style='text-align: center; font-family: Arial; margin-top: 20px; font-size: 16px;'>Total Unclassified: {unknown_count}</div>"
    html_content = html_content.replace("</body>", f"{custom_footer}\n</body>")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Modicum: Dynamic DNA Modification Pathway Predictor")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA of phage proteins")
    parser.add_argument("-d", "--database", required=True, help="HMM database (.hmm file)")
    parser.add_argument("-m", "--map", required=True, help="4-column TSV mapping file (HMM_ID, Gene, Pathway, Type)")
    parser.add_argument("-o", "--output", default="modicum_results", help="Output prefix")
    parser.add_argument("-c", "--coverage", type=float, default=0.8, help="Minimum HMM coverage threshold (0.0 to 1.0). Default: 0.8")
    args = parser.parse_args()

    active_hmm_map, total_acc_per_pathway = load_hmm_mapping(args.map)

    print(f"[*] Scanning {args.input} against {args.database} with >= {args.coverage*100}% coverage threshold...")
    all_hmm_hits = run_pyhmmer_scan(args.input, args.database, min_coverage=args.coverage)

    print("[*] Analyzing pathways...")
    phage_groups = parse_fasta_into_genomes(args.input)
    
    final_results = []
    all_hit_sequences = []
    
    for phage_id, proteins in phage_groups.items():
        result = analyze_single_phage(phage_id, proteins, all_hmm_hits, active_hmm_map, total_acc_per_pathway)
        final_results.append(result)
        all_hit_sequences.extend(result['Annotated_Seqs'])

    csv_file = f"{args.output}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = ['PhageID', 'Prediction', 'Status', 'Evidence', 'Gene_Details']
        # extrasaction='ignore' ensures dictionary keys like 'Pathway_Statuses' are safely omitted from the CSV
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        
        for row in final_results:
            row_copy = row.copy()
            
            formatted_details = []
            for path, gene_counts in row['Gene_Details'].items():
                formatted_details.extend(gene_counts)
            row_copy['Gene_Details'] = "; ".join(formatted_details)
            
            writer.writerow(row_copy)
            
    print(f"[*] Results written to: {csv_file}")

    fasta_file = f"{args.output}_annotated_hits.fasta"
    if all_hit_sequences:
        SeqIO.write(all_hit_sequences, fasta_file, "fasta")
        print(f"[*] Annotated sequences written to: {fasta_file}")
    else:
        print("[!] No hits found, skipping FASTA generation.")

    if PLOTLY_AVAILABLE and all_hit_sequences:
        build_sankey_plot(final_results, f"{args.output}_sankey.html")