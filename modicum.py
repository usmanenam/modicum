import argparse
import sys
import os
import csv
from collections import defaultdict, Counter
from Bio import SeqIO
import matplotlib.pyplot as plt

try:
    import pyhmmer
except ImportError:
    sys.exit("Error: pyhmmer not found. Install: pip install pyhmmer")

try:
    import plotly.graph_objects as go
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# ==========================================
# 1. CONFIGURATION & MAPPINGS
# ==========================================

HMM_MAP = {
    # --- Z-BASE ---
    "NF038379": ("PurZ", "Z-Base"),
    "PF00709":  ("PurA/PurZ_Generic", "Z-Base"),
    "PF04447":  ("MazG-like", "Z-Base"),
    "PF18909":  ("MazG-like-Nt", "Z-Base"),
    "PF12917":  ("DatZ", "Z-Base"),
    # "NF038380": ("DpoZ", "Z-Base"),

    # --- 7-dG ---
    "PF23859":   ("DpdA (TGT)", "7-dG"),
    "NF041059":  ("DpdA (TGT)", "7-dG"),
    "TIGR00364": ("QueC", "7-dG"),
    "NF008317":  ("QueC", "7-dG"),
    "PF06508":   ("QueC", "7-dG"),
    "PF01242":   ("QueD", "7-dG"),
    "TIGR03367": ("QueD", "7-dG"),
    "PF01639":   ("QueE", "7-dG"),
    "TIGR03365": ("QueE", "7-dG"),
    "TIGR04322": ("QueE", "7-dG"),
    "NF006824":  ("FolE", "7-dG"),
    "PF02649":   ("FolE", "7-dG"),
    "PF14489":   ("QueF", "7-dG"),
    "TIGR03138": ("QueF", "7-dG"),
    "TIGR03139": ("QueF", "7-dG"),
    "PF13522":   ("Gat", "7-dG"),
    "NF040592":  ("arcS", "7-dG"),

    # --- 5(hm)Y / Glucosylation ---
    "PF00303":  ("thyA", "5(hm)Y"), # Careful: thyA is canonical too. Needs logic.
    "NF002499": ("thyA", "5(hm)Y"),
    "PF09198":  ("Beta-GT", "5(hm)Y"),
    "PF11440":  ("Alpha-GT", "5(hm)Y"),
    "PF00201":  ("GT-family", "5(hm)Y"),
    "PF00535":  ("GT-family-2", "5(hm)Y"),
    
    # --- dU --- #
    "PF22769":  ("dcd", "dU"),
    "PF00383":  ("dCMP_cyt_deam_1", "dU"),
    "PF18880":  ("Uracil-glycosylase-inhibitor", "dU"),

    # --- dI --- #
    # https://journals.asm.org/doi/10.1128/jvi.01111-19#T1
    "TIGR00551": ("dGMP_reductase", "dI"), # potential dGMP reductase but is otherwise nadB; L-aspartate oxidase. L-aspartate oxidase is the B protein, NadB, of the quinolinate synthetase complex.
    "TIGR01078": ("dAMP_deaminase", "dI"), # potential dAMP deaminase but is otherwise arcA; arginine deiminase. Arginine deiminase is the first enzyme of the arginine deiminase pathway
    "TIGR00576": ("dGTPase", "dI"), # potential dGTAPase but is actually a dUTPase
    # --- momylation --- #
    "PF25680":  ("Mom", "Momylation"),

    # --- Other --- #
    "PF12851":  ("TET/JBP", "Base J/5hmU"),
}

# Pathway Logic
# Note: Added empty set for accessory to avoid key errors if missing
PATHWAY_DEFINITIONS = {
    "dZ": {
        "primary": {"PurZ", "PurA/PurZ_Generic"}, 
        "accessory": {"MazG-like", "MazG-like-Nt", "DatZ"}
    },
    "7-dG": {
        "primary": {"DpdA (TGT)"}, 
        "accessory": {"QueC", "QueD", "QueE", "FolE", "QueF", "Gat", "arcS"}
    },
    "5(hm)Y": {
        "primary": {"thyA", "dUMP_Hmase", "TET/JBP"}, # Assuming you'll add these HMMs
        "accessory": {"Alpha-GT", "Beta-GT", "GT-family", "GT-family-2"}
    },
    "dU": {
        "primary": {"Uracil-glycosylase-inhibitor", "dcd"}, 
        "accessory": {"dCMP_cyt_deam_1"}
    },
    "dI": {
        "primary": {"dGMP_reductase", "dAMP_deaminase"}, 
        "accessory": {"dGTPase"}
    },
    "Mom": {
        "primary": {"Mom"}, 
        "accessory": set()
    },
    # "Hypermod-T": {
    #     "primary": {"5-HMUDK"}, 
    #     "accessory": set()
    # }
}

# ==========================================
# 2. CORE PROCESSING
# ==========================================

def get_phage_id(header):
    """
    Extracts PhageID from header formats like:
    >PhageID_cds001 annotation...
    >PhageID_prot1...
    """
    # Split by _cds or first space
    simple_id = header.split()[0] # Remove description
    if "_cds" in simple_id:
        return simple_id.split("_cds")[0]
    elif "_prot" in simple_id:
        return simple_id.split("_prot")[0]
    else:
        # Fallback: Assume the whole ID is the phage ID if no delimiter found
        # Or remove the last underscore suffix (common in Prodigal)
        return simple_id.rsplit('_', 1)[0]

def group_proteins_by_phage(fasta_path):
    """Parses FASTA and groups SeqRecords by PhageID."""
    phage_groups = defaultdict(list)
    for record in SeqIO.parse(fasta_path, "fasta"):
        phage_id = get_phage_id(record.id)
        phage_groups[phage_id].append(record)
    return phage_groups

def run_hmm_scan(fasta_path, hmm_path, evalue=1e-5):
    """
    Scans the ENTIRE combined fasta once.
    Returns dict: { protein_id: [list of accession_hits] }
    """
    results = defaultdict(list)
    print("[*] Loading HMM database...")
    
    with pyhmmer.easel.SequenceFile(fasta_path, digital=True) as seq_file:
        proteins = list(seq_file)
    
    print(f"[*] Scanning {len(proteins)} proteins against HMMs...")
    
    with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
        for hmm in hmm_file:
            # Robust Accession Extraction
            if hmm.accession:
                h_acc = hmm.accession.decode().split('.')[0]
            else:
                h_acc = hmm.name.decode()
            
            # Run Search
            for hits in pyhmmer.hmmer.hmmsearch(hmm, proteins, E=evalue):
                for hit in hits:
                    if hit.included:
                        pid = hit.name.decode()
                        results[pid].append(h_acc)
    return results

def analyze_single_phage(phage_id, protein_records, all_hmm_hits):
    """
    Analyzes one specific phage using its subset of proteins and their hits.
    Returns analysis results AND a list of annotated SeqRecords for hits.
    """
    found_genes_map = defaultdict(list)
    found_gene_names = set()
    annotated_hits = [] # List to store SeqRecords of interest

    # 1. Gather Evidence (HMM)
    for record in protein_records:
        pid = record.id
        protein_hit_names = set()
        
        # HMM Evidence
        if pid in all_hmm_hits:
            for acc in all_hmm_hits[pid]:
                if acc in HMM_MAP:
                    gene_name, _ = HMM_MAP[acc]
                    found_genes_map[gene_name].append(pid)
                    found_gene_names.add(gene_name)
                    protein_hit_names.add(gene_name)
            
        # If this protein is relevant, save it for FASTA output
        if protein_hit_names:
            # Create a copy of the record to annotate
            annotated_rec = record[:] 
            # Append the hit info to the description
            hit_str = "+".join(sorted(protein_hit_names))
            annotated_rec.description = f"{record.description} [Modicum_Hit: {hit_str}]"
            # Ensure ID is preserved correctly
            annotated_rec.id = record.id 
            annotated_hits.append(annotated_rec)

    # 2. Evaluate Pathways
    best_pathway = "Unknown"
    status = "None"
    evidence_str = ""

    # Priority Check (First strong match wins)
    for p_name, rules in PATHWAY_DEFINITIONS.items():
        primary_hits = rules["primary"].intersection(found_gene_names)
        accessory_hits = rules["accessory"].intersection(found_gene_names)
        
        if primary_hits:
            best_pathway = p_name
            status = "Strong"
            evidence_str = f"Primary: {','.join(primary_hits)} | Acc: {len(accessory_hits)}"
            break # Stop after first strong hit
        
        elif len(accessory_hits) >= 2 and status == "None":
            best_pathway = p_name
            status = "Partial/Incomplete"
            evidence_str = f"Missing Primary | Acc: {','.join(accessory_hits)}"

    return {
        "PhageID": phage_id,
        "Prediction": best_pathway,
        "Status": status,
        "Evidence": evidence_str,
        "Gene_Details": dict(found_genes_map),
        "Annotated_Seqs": annotated_hits
    }

# ==========================================
# 3. VISUALIZATION
# ==========================================

def generate_sankey_plot(results_list, output_base):
    """Creates a Sankey diagram mapping Status -> Prediction."""
    if not PLOTLY_AVAILABLE:
        print("[!] Plotly not found. Skipping Sankey plot. Install: pip install plotly")
        return

    # Count flows: (Status) -> (Prediction)
    flows = defaultdict(int)
    for r in results_list:
        source = r["Status"] if r["Status"] else "Unclassified"
        target = r["Prediction"] if r["Prediction"] else "Unknown"
        flows[(source, target)] += 1

    # Create node lists
    all_sources = list(set(k[0] for k in flows.keys()))
    all_targets = list(set(k[1] for k in flows.keys()))
    all_nodes = all_sources + all_targets
    
    node_map = {name: i for i, name in enumerate(all_nodes)}

    # Create link lists
    sources = []
    targets = []
    values = []

    for (src, tgt), count in flows.items():
        sources.append(node_map[src])
        targets.append(node_map[tgt])
        values.append(count)

    # Define colors
    node_colors = ["#81B1C2" if n in all_sources else "#C78484" for n in all_nodes]

    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values
        ))])

    fig.update_layout(title_text="Phage DNA Modification Flows (Status -> Type)", font_size=10)
    
    out_file = f"{output_base}_sankey.html"
    fig.write_html(out_file)
    print(f"[*] Sankey plot saved to: {out_file}")

def generate_pie_chart(results_list, output_base):
    """Creates a pie chart of the modification distribution."""
    # Count predictions
    counts = Counter([r["Prediction"] for r in results_list])
    
    # Data for plotting
    labels = [f"{k} ({v})" for k, v in counts.items()]
    sizes = list(counts.values())
    
    # Visual Tweaks - Palette inspired by phage genome map (Pastels)
    # Pink, Light Blue, Light Green, Yellow/Orange, Purple, Grey
    colors = [
        "#C78484", # Pastel Pink (Salmon-like)
        "#81B1C2", # Pastel Blue
        "#9EBE7D", # Pastel Green
        "#CF923B", # Pastel Orange
        '#B39EB5', # Pastel Purple
        '#FFF68F', # Pastel Yellow
        "#7AB1FA", # Blue Light Grey (for Unknown/Others)
        '#FF6961', # Soft Red
    ]
    
    plt.figure(figsize=(10, 7))
    patches, texts = plt.pie(sizes, colors=colors[:len(sizes)], startangle=90, radius=1.2)
    plt.legend(patches, labels, loc="best")
    plt.axis('equal')
    plt.title(f"Phage DNA Modification Distribution (N={len(results_list)})")
    
    out_file = f"{output_base}_summary.png"
    plt.savefig(out_file)
    print(f"[*] Pie chart saved to: {out_file}")

# ==========================================
# 4. MAIN
# ==========================================

def main():
    parser = argparse.ArgumentParser(description="MODICUM: Predicting Phage DNA modifications")
    parser.add_argument("-i", "--input", required=True, help="Combined Protein FASTA")
    parser.add_argument("-d", "--database", required=True, help="Custom HMM Database (.hmm)")
    parser.add_argument("-o", "--output", default="modicum_results", help="Output basename (no extension)")
    parser.add_argument("--sankey", action="store_true", help="Generate a Sankey plot instead of a Pie chart")

    args = parser.parse_args()

    if not os.path.exists(args.input): sys.exit("Input file error.")
    if not os.path.exists(args.database): sys.exit("Database file error.")

    # 1. Parse Input & Group
    print(f"[*] Grouping proteins by Phage ID...")
    phage_groups = group_proteins_by_phage(args.input)
    print(f"    > Found {len(phage_groups)} unique phages.")

    # 2. Run Scan
    all_hmm_hits = run_hmm_scan(args.input, args.database)

    # 3. Analyze Each Phage
    final_results = []
    all_hit_sequences = [] # Collect all annotated proteins here
    
    print("[*] Analyzing pathways...")
    for phage_id, proteins in phage_groups.items():
        result = analyze_single_phage(phage_id, proteins, all_hmm_hits)
        final_results.append(result)
        all_hit_sequences.extend(result['Annotated_Seqs'])

    # 4. Write CSV
    csv_file = f"{args.output}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = ['PhageID', 'Prediction', 'Status', 'Evidence', 'Gene_Details']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in final_results:
            # Flatten gene details for cleaner CSV
            row_copy = row.copy()
            del row_copy['Annotated_Seqs'] # Don't write the seq objects to CSV
            row_copy['Gene_Details'] = "; ".join([f"{k}({len(v)})" for k,v in row['Gene_Details'].items()])
            writer.writerow(row_copy)
    print(f"[*] Results written to: {csv_file}")

    # 5. Write Annotated FASTA
    fasta_file = f"{args.output}_annotated_hits.fasta"
    if all_hit_sequences:
        SeqIO.write(all_hit_sequences, fasta_file, "fasta")
        print(f"[*] Annotated sequences written to: {fasta_file}")
    else:
        print("[!] No hits found, skipping FASTA generation.")

    # 6. Generate Chart
    if args.sankey:
        generate_sankey_plot(final_results, args.output)
    else:
        generate_pie_chart(final_results, args.output)

    print("\nDone! Modicum analysis complete.")

if __name__ == "__main__":
    main()