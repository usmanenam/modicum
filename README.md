# MODICUM: Predicting DNA modifications in bacteriophage genomes by base biosynthetic pathway analysis

MODICUM predicts non-canonical DNA modifications in bacteriophage genomes by searching any given genome for tell-tale genes found in bacteriophage DNA modification pathways. The genes involved in each pathway have been curated from literature and will be updated (probably not-so-regularly). The tool is only as good as my reading of the literature as of 2/17/2026.

## Supported Pathways

*   **dZ (2-aminoadenine)**: PurZ pathway.
*   **7-dG (7-deazaguanine)**: DpdA/Que cluster.
*   **5(hm)Y**: Hydroxymethylated, hydroxylated or methylated pyrimidines + glycosylation.
*   **dU (Deoxyuracil)**: UGI / dcd (not tested).
*   **dI (Deoxyinosine)**: dGMP_reductase / dAMP_deaminase (not tested).
*   **Momylation**: Mu-like Mom (not tested).
*   **Base J**: JBP/TET pathway (not tested).

## Installation

### 1. Dependencies

We recommend using Mamba (or Conda) to manage dependencies. Create and activate the modicum environment:

```bash
mamba create -n modicum -c conda-forge -c bioconda python biopython matplotlib plotly pyhmmer hmmer
mamba activate modicum
```

### 2. Prepare Database (`custom_db.hmm`)

**Note:** This has already been compiled for you and is included in the folder `latest_hmm_db`. The following is what I did to compile this database.

Create `custom_db.hmm` using `hmmfetch` (part of the HMMER suite installed above) and the source databases (`Pfam-A.hmm` and `hmm_PGAP.hmm`). Pfam-A was downloaded from Interpro and hmm_PGAP from NCBI.

#### Step A: Create Key Files

Create `pfam_list.txt` with these accessions:

```text
PF04447
PF18909
PF00709
PF12917
PF23859
PF06508
PF01242
PF01639
PF02649
PF13522
PF00585
PF11440
PF09198
PF00201
PF00535
PF12851
PF00303
PF18880
PF22769
PF00383
PF25680
```

Create `ncbi_list.txt` with these accessions:

```text
NF038380.1
NF038379.2
NF041059.1
TIGR00364.1
NF008317.0
TIGR03367.1
TIGR03365.1
TIGR04322.1
TIGR03138.1
TIGR03139.1
NF006824.0
NF040592.1
NF002499.0
TIGR00551.1
TIGR01078.1
TIGR00576.1
```

#### Step B: Fetch and Compile

```bash
# 1. Fetch Pfam profiles (auto-detect versions)
grep "^ACC" Pfam-A.hmm | grep -F -f pfam_list.txt | awk '{print $2}' > pfam_versions.txt
hmmfetch -f Pfam-A.hmm pfam_versions.txt > part1.hmm

# 2. Fetch NCBI/TIGRFAM profiles
hmmfetch -f hmm_PGAP.hmm ncbi_list.txt > part2.hmm

# 3. Combine and Press
cat part1.hmm part2.hmm > custom_db.hmm
hmmpress custom_db.hmm
```

## Usage

### Input Format

Input must be a Protein FASTA (`.faa`). If this is a multi-phage fasta, then headers must use an underscore separator (e.g., `_cds` or `_prot`) to group proteins by phage.

**Format:** `>UniquePhageName_GeneNumber`

**Example:**
```fasta
>PhageAlpha_cds001 putative terminase
MKYL...
>PhageAlpha_cds002 capsid protein
MSTR...
```

### Running the Tool

**Basic Run (CSV + Pie Chart):**
```bash
python modicum.py -i phages.faa -d latest_hmm_db/pfam_ncbi_021726.hmm -o my_results
```

**Sankey Plot Run:**
```bash
python modicum.py -i phages.faa -d latest_hmm_db/pfam_ncbi_021726.hmm -o my_results --sankey
```

### Arguments

*   `-i, --input`: Combined protein FASTA file (Required).
*   `-d, --database`: Pressed `custom_db.hmm` file (Required).
*   `-o, --output`: Output basename (Default: `modicum_results`).
*   `--sankey`: Generates an interactive HTML Sankey diagram instead of a static Pie chart.

## Outputs

1.  **Results CSV (`my_results.csv`)**
    *   **PhageID**: Identifier from FASTA header.
    *   **Prediction**: Predicted pathway (e.g., dZ, 7-dG).
    *   **Status**: Confidence (Strong = Primary marker found; Partial = Accessory only).
    *   **Evidence**: Summary of hits.
    *   **Gene_Details**: Specific genes detected.

2.  **Annotated FASTA (`my_results_annotated_hits.fasta`)**
    *   Sequences of identified modification genes. Headers appended with `[Modicum_Hit: GeneName]`.

3.  **Visualization**
    *   **Summary Chart (`.png`)**: Distribution pie chart.
    *   **Sankey Plot (`.html`)**: Interactive flow diagram (Status -> Prediction).

## Logic & Classification

*   **Strong**: Primary marker detected (e.g., PurZ, DpdA). From my reading of the literature, these genes are always present in phages studied to contain a particular DNA modification.
*   **Partial/Incomplete**: Primary marker missing, but multiple accessory genes present. Suggests variants or fragmentation.
*   **Unknown**: No significant modification hits found. If you know for a fact that your genome is modified then this is potentially an interesting result containing a novel DNA modification biosynthesis pathway.

## Citation

Please cite the repository if used (publication in works).

*This code hasn’t been heavily benchmarked yet or tested. There are probably a lot of bugs and lots of things left incomplete. Also note that a lot of this code was written with the help of Google Gemini. I will add citations for the pathways that this information was drawn from soon.*