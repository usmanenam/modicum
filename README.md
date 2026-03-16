# MODICUM: Predicting DNA modifications in bacteriophage genomes by base biosynthetic pathway analysis

MODICUM predicts non-canonical DNA modifications in bacteriophage genomes by searching any given genome for tell-tale genes found in bacteriophage DNA modification pathways. The genes involved in each pathway have been curated from literature and will be updated (probably not-so-regularly). The tool is only as good as my reading of the literature as of 2/17/2026.
This tool can also technically support automated detection of any curated pathways of your choice so long as you make a mapping file. More on that below.

## Supported Pathways

*   **dZ (2-aminoadenine)**: PurZ pathway.
*   **7-dG (7-deazaguanine)**: DpdA/Que cluster.
*   **5(hm)Y**: Thymidylate synthase-like genes; Hydroxymethylated, hydroxylated or methylated pyrimidines + glycosylation.
*   **dU (Deoxyuracil)**: UGI / dCMP deaminase / Pyrophosphotases (This is perhaps the weakest and hardest pathway prediction. There aren't any good marker genes or accessory genes that are unique to dU phages and so prediction is difficult. Lots of false positives).
*   **dI (Deoxyinosine)**: dGMP_reductase / dAMP_deaminase (not tested; this is perhaps the second weakest prediction 1) because the pathway has not been mechanistically proven and 2) because there are a fair number of false positives).
*   **Momylation**: Mu-like Mom.

## Databases

I have compiled and curated two different different databases containing HMMs for genes involved in DNA modifications. You can use either or both. Modicum runs quite fast thanks to pyhmmer so running both is not a bad idea.

### PFAM/TIGRFAM

hmm file: ```/latest_hmm_db/pfam_ncbi.hmm```

mapping file: ```/latest_hmm_db/pfam_ncbi_map.tsv```

These are general PFAM and TIGRFAM HMMs that are commonly found in bacteriophage DNA modification pathways. In some cases, these offer greater sensitivity than the other database (below) because they are models that have been made across the tree of life rather than just for bacteriophages.

### EnVhog

hmm file: ```/latest_hmm_db/envhogs_dnamods_v4.hmm```

mapping file: ```/latest_hmm_db/envhogs_dnamods_v4_map.tsv```

The EnVhog database comprises clusters of families of viral proteins compiled from environmental sequencing data and filtered for cellular sequences. See paper here: Pérez-Bucio, R.; Enault, F.; Galiez, C. EnVhogDB: an extended view of the viral protein families on Earth through a vast collection of HMM profiles. Peer Community Journal, Volume 5 (2025), article no. e100. https://doi.org/10.24072/pcjournal.627

I took representative DNA modification genes and extracted HMM profiles of the most similar clusters from this database with coverage of at least 0.8 and hhblits probability of at least 95%. 

Using the EnVhog database for Modicum will often result in fewer false positive hits at the cost of sensitivity (i.e. you may be unable to pick up more divergent versions of genes)

## Installation

### 1. Dependencies

We recommend using Mamba (or Conda) to manage dependencies. Create and activate the modicum environment:

```bash
mamba create -n modicum -c conda-forge -c bioconda python biopython matplotlib plotly pyhmmer hmmer
mamba activate modicum
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

2.  **Annotated FASTA (`modicum_results_annotated_hits.fasta`)**
    *   Protein FASTA sequences of identified modification genes. Headers appended with each HMM that was hit `Pathway|HMM name|Marker/Accessory|E-val|Bit-Score|Coverage, ...`.

3.  **Visualization**
    *   **Sankey Plot (`modicum_results.html`)**: Interactive flow diagram (Status -> Prediction).

## Logic & Classification

*   **Strong**: Primary marker detected (e.g., PurZ, DpdA). From my reading of the literature, these genes are either always present in phages studied to contain a particular DNA modification or are sufficient to result in a particular DNA modification.
*   **Partial/Incomplete**: Primary marker missing, but greater than half of the accessory genes for this pathway mentioned in the map file present. Could suggest incomplete genome or divergent pathway.
*   **Unknown**: No significant modification hits found or no primary and less than half of accessory genes found. If you know for a fact that your genome is modified then this is potentially an interesting result containing a novel DNA modification biosynthesis pathway.

## Making your own database

### HMM profiles

As mentioned earlier, this tool can technically be use to predict any pathway that you have curated as long as you provide a HMMER3 `hmmpress` database containing HMMs of genes involved in your pathways of interest. All you have to do is make a mapping file (as described below). 

I will use this space to describe how I constructed/curated the EnVhog and PFAM/NCBI DNA modification databases soon. 

### Map file

Note: The mapping file for DNA modifications is already provided but you may also generate your own.

This is a TSV file containing the mapping between HMMs, gene names, pathways and whether they are an accessory or marker gene. It should be of the following format:

```
HMM_ID   Gene   Pathway   Marker/Accessory
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
The TSV file must start with a header as the first line is always skipped
The `HMM_ID` field should contain the specific alphanumeric value associated with the accession field of the HMMs.
The `Gene` field contains the gene name you associate with that HMM
The `Pathway` field contains the pathway that gene belongs to
The `Marker/Accessory` field should contain either Marker or Accessory depending on how unique and sufficient that gene is to the pathway

## Citation

Please cite the repository if used (publication in works). Also cite PFAM/NCBI, the makers of EnVhog, HMMER3, HH-suite and the papers which I used to curate the database (adding soon)

*This code hasn’t been heavily benchmarked yet or tested. There are probably a lot of bugs and lots of things left incomplete. Also note that a lot of this code was written with the help of Google Gemini. I will add citations for the pathways that this information was drawn from soon.*
