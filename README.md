# MODICUM: Predicting DNA modifications in bacteriophage genomes by base biosynthetic pathway analysis

MODICUM predicts non-canonical DNA modifications in bacteriophage genomes by searching any given genome for tell-tale genes found in bacteriophage DNA modification pathways. The genes involved in each pathway have been curated from literature and will be updated (probably not-so-regularly). The tool is only as good as my reading of the literature as of 6/7/2026.
This tool can also technically support automated detection of any curated pathways of your choice so long as you make a mapping file. More on that below.

<img width="2523" height="1024" alt="Screenshot 2026-03-16 152108" src="https://github.com/user-attachments/assets/2f35b628-7f17-4def-85d4-4f01bf62f814" />


## Supported Pathways

*   **dZ (2-aminoadenine)**: PurZ pathway.
*   **7-dG (7-deazaguanine)**: DpdA/Que cluster.
*   **5(hm)dY**: Thymidylate synthase-like genes; Hydroxymethylated, hydroxylated or methylated pyrimidines + glycosylation. With the latest version it will now also predict whether the modification is 5hmdU, 5hmdC, 5mdC or 5hdC by inspecting key catalytic residues.
*   **dU (Deoxyuracil)**: UGI / dCMP deaminase / Pyrophosphotases (This is perhaps the weakest and hardest pathway prediction. There aren't any good marker genes or accessory genes that are unique to dU phages and so prediction is difficult. Lots of false positives). This will almost always be given putative status
*   **dI (Deoxyinosine)**: dGMP_reductase / dAMP_deaminase (not tested; this is perhaps the second weakest prediction 1) because the pathway has not been mechanistically proven and 2) because there are a fair number of false positives). This will almost always be given putative status
*   **Momylation**: Mu-like Mom.

## Databases

I have compiled and curated two different different databases containing HMMs for genes involved in DNA modifications. You can use either or both. Modicum runs quite fast thanks to pyhmmer so running both is not a bad idea but I would start with the EnVhog database.

### EnVhog

hmm file: ```/latest_hmm_db/envhogs_dnamods_v4.hmm```

mapping file: ```/latest_hmm_db/envhogs_dnamods_v4.3_map.tsv```

The EnVhog database comprises clusters of families of viral proteins compiled from environmental sequencing data and filtered for cellular sequences. See paper here: Pérez-Bucio, R.; Enault, F.; Galiez, C. EnVhogDB: an extended view of the viral protein families on Earth through a vast collection of HMM profiles. Peer Community Journal, Volume 5 (2025), article no. e100. https://doi.org/10.24072/pcjournal.627

I took representative DNA modification genes and extracted HMM profiles of the most similar clusters from this database with coverage of at least 0.8 and hhblits probability of at least 95%. 

### PFAM/TIGRFAM

hmm file: ```/latest_hmm_db/pfam_ncbi.hmm```

mapping file: ```/latest_hmm_db/pfam_ncbi_v1.3_map.tsv```

These are general PFAM and TIGRFAM HMMs that are commonly found in bacteriophage DNA modification pathways. In some cases, these offer greater sensitivity than the other database (above) because they are models that have been made across the tree of life rather than just for bacteriophages.

Using the EnVhog database for Modicum will often result in fewer false positive hits at the cost of sensitivity (i.e. you may be unable to pick up more divergent versions of genes). However, there is a greater diversity of HMMs that are covered in the EnVhog database.

## Installation

### 1. Dependencies

We recommend using Mamba (or Conda) to manage dependencies. Create and activate the modicum environment:

```bash
mamba create -n modicum -c conda-forge -c bioconda python biopython matplotlib pyhmmer hmmer
mamba activate modicum
pip install sankeyflow
```

These are the versions of the dependencies that Modicum has currently been tested with:
```
python 3.13.13
biopython 1.87
matplotlib 3.10.9
pyhmmer 0.12.0
hmmer 3.4
sankeyflow 0.4.1
```

### 2. Script and Databases

Download the ```modicum.py``` python script and the ```/latest_hmm_db/``` folder using either git or manual download


## Usage

### Input Format

Input must be a Protein FASTA (`.faa`). If this is a multi-phage fasta, then headers must use an underscore separator to group proteins by phage.

**Format:** `>UniquePhageName_GeneNumber`

**Example:**
```fasta
>PhageAlpha_cds001 putative terminase
MKYL...
>PhageAlpha_cds002 capsid protein
MSTR...
```

These protein fasta files can be generated via Pharokka, Prodigal or Prodigal-gv amongst other tools. You can put DNA fasta file containing multiple phage genomes into these programs and they will spit out a protein fasta of the format above. 

### Running the Tool

```bash
python modicum.py -i phages.faa -d latest_hmm_db/pfam_ncbi.hmm -m latest_hmm_db/pfam_ncbi_map.tsv -o modicum_results -c 0.8
```

### Arguments

*   `-i, --input`: Combined protein FASTA file (Required).
*   `-d, --database`: HMMER3 Pressed `db.hmm` file (Required).
*   `-m, --map`: TSV mapping file [HMM_ID, Gene, Pathway, Marker/Accessory] (Required).
*   `-o, --output`: Output basename (Default: `modicum_results`).
*   `-c, --coverage`: Minimum HMM coverage threshold (0.0 to 1.0; Default: 0.8)

Note regarding coverage: The lower the value, the greater the sensitivity at the cost of an increase in false positives. I have found bona fide hits at low coverage ~0.2-0.3 but also false positives. The best thing to do to validate a hit if you are using a low coverage value (or otherwise) is to AlphaFold the sequence and run it through foldseek. 

## Outputs

1.  **Results CSV (`modicum_results.csv`)**
    *   **PhageID**: Identifier extracted from FASTA header.
    *   **Prediction**: Predicted pathway (e.g., dZ, 7-dG, '?' added if putative).
    *   **Status**: Confidence (Strong: Primary marker(s) found; Putative: > half of accessory and no marker genes found; Unknown: <= half of accessory and no marker genes found).
    *   **Evidence**: Summary of hits.
    *   **Gene_Details**: All hits listed regardless of if the prediction threshold was reached.
    *   **HMM_Names**: The accession numbers or names of the HMMs that were hit in the same order as **Gene_Details**

2.  **Annotated FASTA (`modicum_results_annotated_hits.fasta`)**
    *   Protein FASTA sequences of identified modification genes. Headers appended with each HMM that was hit `Pathway|HMM name|Marker/Accessory|E-val|Bit-Score|Coverage, ...`.

3.  **Visualization**
    *   **Sankey Plot (`modicum_results_sankey.svg`)**: Interactive flow diagram (Status -> Prediction).

## Logic & Classification

*   **Strong**: Primary marker detected (e.g., PurZ, DpdA). From my reading of the literature, these genes are either always present in phages studied to contain a particular DNA modification or are sufficient to result in a particular DNA modification.
*   **Partial**: Primary marker missing, but greater than half of the accessory genes for this pathway mentioned in the map file present. Could suggest incomplete genome or divergent pathway.
*   **Unknown**: No significant modification hits found or no primary and less than half of accessory genes found. If you know for a fact that your genome is modified then this is potentially an interesting result containing a novel DNA modification biosynthesis pathway.

## Making your own database

### HMM profiles

As mentioned earlier, this tool can technically be use to predict any pathway that you have curated as long as you provide a HMMER3 `hmmpress` database containing HMMs of genes involved in your pathways of interest. All you have to do is make a mapping file (as described below). 

I will use this space to describe how I constructed/curated the EnVhog and PFAM/NCBI DNA modification databases soon. 

### Map file

Note: The mapping file for DNA modifications is already provided but you may also generate your own.

This is a TSV file containing the mapping between HMMs, gene names, pathways and whether they are an accessory or marker gene. It should be of the following format:

```
HMM_ID   Gene   Pathway   Type   
NF038379   PurZ   dZ   Marker
PF00709   PurA   dZ   Marker
PF04447   MazG-like   dZ   Accessory
PF18909   MazG-like-Nt   dZ   Accessory
PF12917   DatZ   dZ   Accessory
PF23859   DpdA   7-dG   Marker
NF041059   DpdA   7-dG   Marker
TIGR00364   QueC   7-dG   Accessory
.
.
.
```
   - The TSV file must start with a header as the first line is always skipped
   - The `HMM_ID` field should contain the specific alphanumeric value associated with the accession field of the HMMs.
   - The `Gene` field contains the gene name you associate with that HMM
   - The `Pathway` field contains the pathway that gene belongs to
   - The `Marker/Accessory` field should contain either Marker or Accessory depending on how unique and sufficient that gene is to the pathway

#### Some advanced options
The following functionalities are totally optional and are not required for the map file.
1. **Priority**
  - This is applied to marker genes so that if two markers from different pathways are hit in the same genome, then the higher priority one is presented in the output.
  - All marker genes by default are assigned a priority of 1
  - To assign a priority just add a **Priority** column to the map file and fill in the appropriate priority for each marker gene (again, blank means priority of 1)
  - A priority of 0 is a special status. It will always be lower priority than any other markers and if a genome containing a marker gene with this priority is found, and no other marker genes are found, then the output pathway prediction will always be putative. You can assign this if you think your curated evidence for the pathway is weak as I have done for dU and dI.
2. **HMM_Variants**
    - Sometimes it is necessary to distinguish different HMM-protein hits by key/characteristic residues. For example, in Thymidylate Synthase-like genes, the residue at position 177 (relative to E. coli ThyA) determines whether the substrate is a dUMP (E. coli (canonical), SPOI phage...) or a dCMP (T4 phage, RB69 phage...).
    - A standard TS-like dCMP Hmase HMM will also hit a canonical ThyA very well and so to filter these out it can be useful to look at the identity of the key residues in the alignment of the HMM and protein hit.
    - To allow for this, I have added functionality that does this within the modicum.py code if it finds the optional **HMM_Variant** AND **New_Pathway** columns in the map file.
    - Here you can specify which HMM residue variants you are looking for, and what the final prediction should be. Here is an example map entry:
    - ```
      HMM_ID	Gene	Pathway	Type	Priority	HMM_Variant	New_Pathway
      1HJf8	TS_hdC	5(hm)dY	Marker		185N&177D;185N&177N;185N&146L/N/I;185N&146D/E;185N&146D/E&177D;185D&177D;185D&177N;185D&177D&146L/N/I;185D&146D/E;185D&177D&146D/E	5(h)mdU;5hdU;canonical;5h(m)dU;5hmdU;5(h)mdC;5hdC;5mdC;5h(m)dC;5hmdC
      ```
   - Notice that the "parent" pathway prediction is 5(hm)dY. So if the HMM 1HJf8 hits a similar protein, that may even be canonical ThyA, the prediction originally would have been 5(hm)dY. However, additional logic in the **HMM_Variant** and **New_Pathway** field will result in more specific, "children" predictions.
   - For example, if the protein contains an **N** at position 185 AND a **D** at position 177 relative to the HMM, the new prediction will be 5(h)mdU.
   - Parentheticals in the pathway are treated as "may or may not contain". So the product may or may not contain a hydroxy group but it will contain a methyl group.
   - The code will always prioritize a prediction that is more specific than less specific if it can. I.e. if the protein satisfies the 185N & 146D/E & 177D logic for 5hmdU, as well as the 185N & 177D logic for 5(h)mdU, the final predicition will be 5hmdU.
   - The "child" prediction will also inherit all of the accesory gene predictions of the "parent" predicition.
   - Here's how to generate the mapping file with the HMM_Variant and New_Pathway column:
      - I have included a helper tool in the ```./other_tools/``` folder called ```finding_HMM_variant_positions.py```.
      - To run this tool you need:
         - A reference protein sequence (this is the protein that you base all of your variants on as informed by your specific literature search. For me, it was the E. coli ThyA protein, which has multiple studied residue variants that confer different functionality to the protein)
         - Your compiled HMM database file
         - Your basic mapping file that you have already created
         - ```REF_VARIANTS``` Replace with your string of the specific variants. ```&``` implies both variants have to be present, ```\``` implies the identity of the residue can be either, and each variant or combination of variants should be punctuated by a ```;```
         - ```NEW_PATHWAYS``` Replace with your string of new "children" predictions. This should be in the same order as the ```REF_VARIANTS``` string. Again, parentheticals here ```()``` have a special role where they imply the optionality of the element within the parenthicals. In the code logic, 5hdC will always be preferred over 5h(m)dC if the logic for both is satisfied.
         - Use "canonical" for any prediction that you do not want to be output in the final modicum predictions. For example, canonical ThyAs are not associated with DNA modifications and so their discovery in a genome is not revealing. Hence, if one of my TS-like HMMs hits a TS-like protein and it finds 185N & 146L/N/I in the protein, then it will toss this gene altogether and not include it in the final prediction. 

## Citation

Please cite the repository if used (publication in works). Also cite PFAM/NCBI, the makers of EnVhog, HMMER3, HH-suite and the papers which I used to curate the database (adding soon)

The phage drawings above are made by BioRender and I think I'm supposed to say that.

*This code hasn’t been heavily benchmarked yet or tested. There are probably a lot of bugs and lots of things left incomplete. Also note that a lot of this code was written with the help of Google Gemini. I will add citations for the pathways that this information was drawn from soon.*
