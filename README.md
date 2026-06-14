# chikungunyaPipeline

**Dockerized genomic pipeline for evolutionary analysis of Chikungunya virus.**

Queries NCBI for CHIKV sequences from Brazil (2015–2024), aligns them with MAFFT, estimates dN/dS ratios, and generates temporal regression plots and a phylogenetic tree.

## Prerequisites

- Docker and Docker Compose

## Quick start

```bash
docker compose up --build
```

Output files land in `./data/` on the host (mounted at `/app/data` in the container).

## Pipeline steps

| Step | Method | Description |
|------|--------|-------------|
| 1. Fetch | NCBI Entrez (`esearch` + `efetch`) | Searches GenBank for *Chikungunya virus* sequences from Brazil (2015–2024), filters by length >100 bp and records with collection dates |
| 2. Align | MAFFT (`--auto --quiet`) | Multiple sequence alignment, output as FASTA |
| 3. Evol. Analysis | Custom codon translation | Picks first sequence as reference, counts synonymous (S) vs non-synonymous (N) codon changes per sample, computes dN/dS ratio |
| 4. Plot | matplotlib + seaborn | Temporal regression (mutations vs collection date) + Neighbor-Joining tree (rooted at midpoint) |

## Output files

All written to `./data/`:

| File | Description |
|------|-------------|
| `raw_sequences.fasta` | Raw downloaded sequences (FASTA) |
| `aligned.fasta` | MAFFT-aligned sequences |
| `mutation_data.csv` | Per-sample dN/dS statistics |
| `plot_temporal.png` | Molecular clock scatter plot |
| `plot_tree.png` | Phylogenetic tree (NJ, rooted) |
| `tree.xml` | Tree in phyloXML format |

## Configuration

Set the `NCBI_EMAIL` environment variable (default in `docker-compose.yml` is `jeferson0993@gmail.com`). NCBI requires a valid email address for API access.

## Customization

Edit the query in `src/main.py:38` to change search terms or `max_records` (default 30 in `run()`, 50 in `fetch_data_ncbi()`).

## Tech stack

- **Python 3.9** (slim Docker image)
- **Biopython** — NCBI Entrez API, sequence I/O, alignment parsing, phylogenetics
- **MAFFT** — multiple sequence alignment (system binary)
- **pandas / numpy** — data manipulation
- **matplotlib / seaborn** — plotting
- **scipy** — implicit statistical dependencies
- **reportlab** — unused import dependency

## License

BSD 3-Clause — see [LICENSE](LICENSE).
