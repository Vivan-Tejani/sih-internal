# Eukaryon: Deep-Sea eDNA Taxonomy and Biodiversity Pipeline

Smart India Hackathon 2026 | Problem Statement DJS_26_SW_15
Team: The Fellowship of the Ring (DJSW 115) | Theme: Blue Economy | Category: Software

A database-query-free pipeline that identifies eukaryotic marine organisms directly from raw eDNA reads. The pipeline clusters sequences first and assigns labels second, so novel and database-absent species are discovered rather than discarded.

## Overview

Deep-sea biodiversity is critical to study but poorly understood. Most organisms in this environment remain unidentified, and existing classification tools depend on reference databases in which deep-sea and Indian Ocean taxa are barely represented. Forcing every read into a known label mislabels or discards the organisms most worth discovering.

Eukaryon replaces the standard classify-first approach with a cluster-first workflow:

1. Group similar DNA sequences together before querying any reference database.
2. Check those groups for internal consistency.
3. Identify known matches and flag potentially undiscovered organisms separately.

This allows novel or under-represented organisms to form their own cluster instead of being forced into an incorrect label, while every accepted match still carries a calibrated confidence score.

## Key Features

| Feature | Description |
|---|---|
| Database-free identification | Detects species missing from reference databases instead of marking them unassigned |
| Compute-efficient clustering | Identifies and annotates only cluster representatives, not every read |
| Calibrated confidence matching | Every match carries a confidence score calibrated from real reference-similarity data, not a fixed cutoff |
| Neighborhood cross-check | Each top match is verified against its 4 next-nearest neighbors at genus level; disagreements are downgraded |
| Novel organism detection | Separates confident, uncertain, and candidate-novel organisms into distinct categories |
| Automated biodiversity metrics | Computes Shannon, Simpson, and Chao1 diversity indices automatically |

## Architecture

```
Preprocessing (Stage 0)             Embedding Backbone (Stage 1)
Raw paired-end FASTQ                 DNABERT-S (pretrained)
  -> DADA2: filterAndTrim,             -> 768-dim embeddings
     learnErrors, dada,                -> Pure PyTorch attention fallback
     mergePairs
  -> ASV Table and FASTA

Offline Reference Build (Stage 2)   Clustering Architecture (Stage 3)
PR2 Database                         L2-Normalize -> UMAP (dimensionality
  -> embedded once via DNABERT-S        reduction) -> HDBSCAN (density-based
  -> FAISS IVF index                    clustering) -> bootstrap
                                         reproducibility scoring

Species Annotation (Stage 4)
Representative selection (1 per cluster)
  -> FAISS top-5 nearest-neighbor search
  -> k-NN genus-level agreement check (cross-checked against 4 nearest neighbors)
  -> BLAST cross-validation
  -> Biodiversity metrics, final reports, dashboard
```

Matching a sequence without an index requires scanning all 240,000+ references individually, comparable to checking every shelf in a library to find one book. Eukaryon groups the reference database by similarity in a one-time offline build, so each search only checks the relevant section. This reduces lookup time from 365.98 seconds (BLAST) to 2.66 seconds.

## Technology Stack
 
**Libraries and Frameworks**
- fastp / DADA2: read trimming, quality filtering, and denoising into ASVs
- DNABERT-S: pretrained model that converts DNA sequences into numeric embedding vectors
- HDBSCAN: density-based unsupervised clustering of embeddings
- UMAP: dimensionality reduction prior to clustering
- FAISS: approximate nearest-neighbor search against the reference index
- Python dataclasses: structured configuration and pipeline stage storage
**Algorithms**
- HDBSCAN density-based clustering, with no fixed cluster count; "no confident group" is a valid outcome
- FAISS nearest-neighbor search (IVF/HNSW index) against the pre-embedded reference database
- Empirically calibrated confidence thresholding, derived from the observed reference similarity distribution
- k-NN neighborhood cross-check, verifying each top match against its 4 next-nearest neighbors before acceptance
- UMAP manifold approximation for high-dimensional embeddings

## Repository Structure

```
.
├── src/                                    Frontend (React, Vite, TypeScript): dashboard
├── public/
├── PIPELINE_ARCHITECTURE.md
├── context.md
├── index.html
├── package.json / package-lock.json
├── tailwind.config.js
├── postcss.config.js
├── vite.config.js
└── model-codes/                            Pipeline and model branch
    ├── DADA2_PREPROCESSING.R                Stage 0: read trimming, denoising, ASVs
    ├── preprocessing.py
    ├── Embedding and clustering.py          Stage 1 and 3: DNABERT-S embedding, UMAP, HDBSCAN
    ├── Annotion and labelling.py            Stage 4: annotation and taxa labeling
    ├── reference_index(Library analogy).py  Stage 2: offline FAISS reference index build
    ├── comparison with blast.py             Benchmarking against BLAST
    ├── biodiversity_summary.json            Shannon, Simpson, and Chao1 outputs
    ├── candidate_novel_taxa.csv             Flagged novel-candidate organisms
    └── final_report.csv                     Final per-sample taxonomy report
```

The pipeline code resides on the `model-codes` branch. The dashboard (React, Vite, TypeScript frontend) resides on `main`.

## Getting Started

### Prerequisites
- Python 3.10 or later
- R (for DADA2 preprocessing)
- Node.js 18 or later (for the dashboard)
- GPU recommended for the one-time offline reference index build; CPU is sufficient for per-sample inference

### Pipeline (model-codes branch)

```bash
git clone https://github.com/Vivan-Tejani/sih-internal.git
git checkout model-codes
cd sih-internal

# Install Python dependencies
pip install -r requirements.txt

# Stage 0: Preprocessing (DADA2)
Rscript DADA2_PREPROCESSING.R --input <raw_fastq_dir> --output asv_table.fasta

# Stage 1 and 3: Embedding and clustering
python "Embedding and clustering.py" --input asv_table.fasta --output embeddings.npy

# Stage 2: Build offline reference index (one-time)
python "reference_index(Library analogy).py" --db PR2_database.fasta --output reference.index

# Stage 4: Annotation, novel-taxa flagging, biodiversity metrics
python "Annotion and labelling.py" --embeddings embeddings.npy --index reference.index --output final_report.csv
```

### Dashboard (main branch)

```bash
git clone https://github.com/Vivan-Tejani/sih-internal.git
cd sih-internal
npm install
npm run dev
```

Live dashboard: [sih-internal-one.vercel.app](https://sih-internal-one.vercel.app)

## Results (Real Sample Run)

Tested on 503,014 raw eDNA reads (300 bp, 18S rRNA gene) from NCBI SRA run SRR17870216.

| Metric | BLAST | Eukaryon |
|---|---|---|
| Coverage | 87.8% (125 of 1,024 clusters returned zero hit) | 100%; rescued 79 of those clusters with a real similarity signal |
| Clusters with zero hit | Forces every read to a label; novel organisms get mislabeled | 1,024 clusters formed by similarity before any label is assigned |
| Divergent taxa rescued | Misapplied cutoffs flagged 93% of clusters as novel | Thresholds calibrated on 240K reference sequences' real similarity distribution |
| Neighbor-validated confidence | No safeguard; 896 of 1,024 confident calls taken at face value | k-NN cross-check downgraded 46 of 80 confident calls lacking support |
| Lookup time | 365.98 seconds | 2.66 seconds |

Summary: 1,024 groups identified from 157,333 reads. 100% cluster coverage. 214 novel candidates flagged and retained, not discarded. 46 of 80 confident calls verified or downgraded on cross-check.
