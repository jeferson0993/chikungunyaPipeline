# AGENTS.md — chikungunyaPipeline

**Pipeline genômico Dockerizado para análise evolutiva do Vírus Chikungunya.**

## Quick start

```bash
docker compose up --build
```

Data output lands in `./data/` (mounted at `/app/data` in container).

## Architecture

Single-module pipeline at `src/main.py`. Class `ChikungunyaPipeline` runs four sequential steps:

1. `fetch_data_ncbi()` — queries NCBI Entrez for Chikungunya virus sequences
2. `align_sequences()` — calls **MAFFT** (system binary, installed via `apt` in Docker)
3. `analyze_molecular_evolution()` — crude dN/dS via codon translation
4. `generate_plots()` — temporal regression + NJ phylogenetic tree

Entrypoint: `python src/main.py` (set as `CMD` in `Dockerfile`).

## Key requirements

- `NCBI_EMAIL` env var (default: `jeferson0993@gmail.com`, set in `docker-compose.yml`)
- MAFFT must be installed — only available inside the Docker container
- Python 3.9, dependencies in `requirements.txt` (biopython, pandas, matplotlib, seaborn, scipy, reportlab, numpy)
- `matplotlib.use('Agg')` is hardcoded for headless environments

## No tests / no CI

No test framework, no linting, no type checking configured. The `.gitignore` includes `.ruff_cache/` but no Ruff config exists.

## Constraints

- Single container, no restart policy (`restart: "no"`)
- Volume mount `./data:/app/data` — output files persist there
- Container runs once and exits — not a long-running service
