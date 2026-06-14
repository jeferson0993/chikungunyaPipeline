# Pipeline technical reference

This document describes each step of the `ChikungunyaPipeline` class in detail.

## 1. `fetch_data_ncbi()`

```
Query: Chikungunya virus[Organism] AND Brazil[Geo Location] AND 2015:2024[Date - Collection]
Database: nucleotide (GenBank)
Retrieval: esearch → efetch (rettype="gb", retmode="text")
```

**Logic:**
- Searches NCBI with a hardcoded query targeting Brazilian CHIKV isolates from 2015–2024.
- Downloads GenBank records and parses metadata: **accession**, **collection year**, **country**, **sequence length**.
- Filters out records where length ≤ 100 bp or collection date is missing.
- Sorts chronologically and writes a multi-FASTA file (`raw_sequences.fasta`).
- Returns `None` if no records match, skipping all downstream steps.

## 2. `align_sequences()`

```
Command: mafft --auto --quiet <input> > aligned.fasta
```

**Logic:**
- Calls MAFFT as a subprocess via `os.system()`.
- `--auto` lets MAFFT auto-select the alignment strategy.
- `--quiet` suppresses verbose output.
- Reads the aligned output back into a Biopython `MultipleSeqAlignment` object.
- Requires MAFFT to be installed as a system binary (done via `apt-get` in the Dockerfile).

## 3. `analyze_molecular_evolution()`

**Logic:**
- Strips gaps (`-`) from each aligned sequence.
- Trims to the longest in-frame length (multiple of 3).
- Uses the **first sequence** as the reference.
- For every other sequence, walks codon-by-codon (3 nt windows):
  - Skips codons containing `N` (ambiguous base).
  - Compares translation (Biopython `Seq.translate()`):
    - **Different amino acid** → counts as a non-synonymous mutation (N).
    - **Same amino acid** → counts as a synonymous mutation (S).
- Computes dN/dS = N/S (0 if S = 0).
- Writes `mutation_data.csv` with columns: `id, date, synonymous, nonsynonymous, total_mutations, dn_ds`.

**Limitation:** This is a simplified per-sequence pairwise dN/dS against a single reference, not a full phylogenetic-aware method (e.g. PAML, HyPhy). It is useful for relative comparisons across samples.

## 4. `generate_plots()`

### Temporal regression (`plot_temporal.png`)
- Scatter plot of `total_mutations` vs `collection year` for each sample.
- Overlays a linear regression line (via `seaborn.regplot`).
- Title: "Molecular Clock: Mutations vs Time".

### Phylogenetic tree (`plot_tree.png`, `tree.xml`)
- Computes a distance matrix from the alignment using the `identity` model (proportion of differing sites).
- Builds a Neighbor-Joining tree (Biopython `DistanceTreeConstructor`).
- Roots at midpoint.
- Saves tree as both **PNG** (matplotlib rendering) and **phyloXML** (`tree.xml`) for external viewers (e.g. Archaeopteryx, FigTree).
- Errors are caught and logged but do not halt the pipeline.

## Output files summary

| File | Step | Format | Contents |
|------|------|--------|----------|
| `raw_sequences.fasta` | 1 | FASTA | Downloaded GenBank sequences |
| `aligned.fasta` | 2 | FASTA | MAFFT-aligned sequences |
| `mutation_data.csv` | 3 | CSV | Per-sample dN/dS metrics |
| `plot_temporal.png` | 4 | PNG 300dpi | Molecular clock scatter plot |
| `plot_tree.png` | 4 | PNG 300dpi | NJ phylogenetic tree |
| `tree.xml` | 4 | phyloXML | Machine-readable tree |

## Environment variables

| Variable | Default | Purpose |
|----------|---------|---------|
| `NCBI_EMAIL` | `jeferson0993@gmail.com` | Required by NCBI Entrez API |
| `PYTHONDONTWRITEBYTECODE` | `1` | Prevents `.pyc` files in the container |
| `PYTHONUNBUFFERED` | `1` | Flushes stdout/stderr immediately |

## Docker details

- **Base image:** `python:3.9-slim` (Debian-based).
- **System packages:** `mafft`, `build-essential`.
- **Working directory:** `/app`.
- **Data volume:** `./data:/app/data` — all output files.
- **Container lifecycle:** Runs once and exits (`restart: "no"`).
