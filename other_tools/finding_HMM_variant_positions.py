import pyhmmer
import csv
import re
import os

def update_hmm_mapping_file(ref_fasta_path, hmm_database_path, input_tsv, output_tsv, ref_variants_str, new_pathways_str):
    """
    Aligns a reference protein to HMMs to extract variant nodes,
    and updates an existing Modicum mapping file with the specific HMM variants.
    """
    alphabet = pyhmmer.easel.Alphabet.amino()
    
    # Automatically extract target residues from the variant string
    target_residues = list(set(map(int, re.findall(r'\d+', ref_variants_str))))
    if not target_residues:
        print("[!] Error: Could not find any residue positions (numbers) in ref_variants_str.")
        return

    # Load the Reference sequence
    if not os.path.exists(ref_fasta_path):
        print(f"[!] Error: Reference FASTA not found at {ref_fasta_path}")
        return
        
    with pyhmmer.easel.SequenceFile(ref_fasta_path, digital=True, alphabet=alphabet) as seq_file:
        ref_seqs = seq_file.read_block()
        
    if len(ref_seqs) == 0:
        print("[!] Error: No sequence found in the provided Reference FASTA.")
        return

    pattern = re.compile(r'(?<!\d)(' + '|'.join(map(str, target_residues)) + r')(?!\d)')
    updates_dict = {}

    # Search Reference against EVERY HMM and build our dictionary of updates
    print("[*] Scanning HMMs and mapping active site nodes...")
    with pyhmmer.plan7.HMMFile(hmm_database_path) as hmm_file:
        for hmm in hmm_file:
            
            hmm_name = hmm.name.decode() if isinstance(hmm.name, bytes) else (hmm.name or "")
            hmm_acc = hmm.accession.decode() if hmm.accession and isinstance(hmm.accession, bytes) else (hmm.accession or "")
            
            possible_keys = set()
            if hmm_name:
                possible_keys.add(hmm_name)
                possible_keys.add(hmm_name.split('.')[0]) 
            if hmm_acc:
                possible_keys.add(hmm_acc)
                possible_keys.add(hmm_acc.split('.')[0]) 
                
            pipeline = pyhmmer.plan7.Pipeline(alphabet)
            hits = pipeline.search_hmm(hmm, ref_seqs)
            
            if not hits:
                continue
                
            ref_hit = hits[0]
            mapped_nodes = {res: "Not_Found" for res in target_residues}
            
            for domain in ref_hit.domains:
                if not domain.included:
                    continue
                    
                ali = domain.alignment
                hmm_pos = ali.hmm_from
                target_pos = ali.target_from
                
                for hmm_char, target_char in zip(ali.hmm_sequence, ali.target_sequence):
                    if target_pos in target_residues and target_char != '-' and hmm_char != '.':
                        mapped_nodes[target_pos] = hmm_pos
                    
                    if target_char != '-':
                        target_pos += 1
                    if hmm_char != '.':
                        hmm_pos += 1
            
            if "Not_Found" in mapped_nodes.values():
                for k in possible_keys:
                    updates_dict[k] = ("", "")
            else:
                def replace_with_mapped_node(match):
                    res_num = int(match.group(1))
                    return str(mapped_nodes[res_num])
                    
                final_variants = pattern.sub(replace_with_mapped_node, ref_variants_str)
                for k in possible_keys:
                    updates_dict[k] = (final_variants, new_pathways_str)

    # Update the existing Mapping File
    if not os.path.exists(input_tsv):
        print(f"[!] Error: Input TSV mapping file not found at {input_tsv}")
        return

    print(f"[*] Updating mapping file: {input_tsv}")
    updated_rows = []
    
    with open(input_tsv, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        header = [h.strip() for h in header]
        
        if 'HMM_Variant' not in header:
            header.append('HMM_Variant')
        if 'New_Pathway' not in header:
            header.append('New_Pathway')
            
        var_idx = header.index('HMM_Variant')
        path_idx = header.index('New_Pathway')
        
        hmm_id_idx = header.index('HMM_ID') if 'HMM_ID' in header else 0 
        
        processed_count = 0
        
        for row in reader:
            if not row or row[0].startswith('#'):
                updated_rows.append(row)
                continue
                
            while len(row) < len(header):
                row.append("")
                
            hmm_id = row[hmm_id_idx].strip()
            
            if hmm_id in updates_dict:
                variant_val, pathway_val = updates_dict[hmm_id]
                row[var_idx] = variant_val
                row[path_idx] = pathway_val
                if variant_val != "":
                    processed_count += 1
            else:
                pass
                
            updated_rows.append(row)

    with open(output_tsv, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        writer.writerows(updated_rows)
        
    print(f"[*] Successfully mapped and populated variants for {processed_count} HMMs.")
    print(f"[*] Updated mapping file saved to: {output_tsv}")


if __name__ == "__main__":
    
    # This is an example of what I did to make a map of the relevant thyA variants

    REF_PROT_FASTA = "ecoli_thya.faa"  
    HMM_DB = "../latest_hmm_db/pfam_ncbi.hmm"   
    
    INPUT_MAPPING_TSV = "../latest_hmm_db/pfam_ncbi_v1.2_map.tsv"   
    OUTPUT_MAPPING_TSV = "../latest_hmm_db/pfam_ncbi_v1.3_map.tsv"
    
    REF_VARIANTS = "177N&169D;177N&169N;177N&143L/N/I;177N&143D/E;177N&143D/E&169D;177D&169D;177D&169N;177D&169D&143L/N/I;177D&143D/E;177D&169D&143D/E"
    NEW_PATHWAYS = "5(h)mdU;5hdU;canonical;5h(m)dU;5hmdU;5(h)mdC;5hdC;5mdC;5h(m)dC;5hmdC"
    
    update_hmm_mapping_file(
        ref_fasta_path=REF_PROT_FASTA, 
        hmm_database_path=HMM_DB, 
        input_tsv=INPUT_MAPPING_TSV,
        output_tsv=OUTPUT_MAPPING_TSV,
        ref_variants_str=REF_VARIANTS, 
        new_pathways_str=NEW_PATHWAYS
    )