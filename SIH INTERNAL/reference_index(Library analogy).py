%%writefile build_reference_index.py
#!/usr/bin/env python3
"""
build_reference_index.py  (STEP 3 PREREQUISITE -- run ONCE, before step3)

ONE-TIME, OFFLINE precompute step. Run this ONCE before Step 3 -- its output
(reference.faiss + reference_metadata.json) is what Step 3 uses for species
annotation.

Takes a reference database FASTA (e.g. SILVA, PR2) with accompanying taxonomy
metadata, embeds every reference sequence with the SAME embedding model used
in the main pipeline (DNABERT-S by default, or k-mer fallback), and builds a
FAISS approximate-nearest-neighbor index (IVF or HNSW) from those embeddings.

This only needs to be run ONCE per reference database version. The resulting
index + metadata files are then reused by pipeline.py for every sample run --
the reference database is never re-embedded per sample.

Usage:
    python build_reference_index.py \
        --ref-fasta silva_138_ssu_nr99.fasta \
        --ref-taxonomy silva_138_taxonomy.tsv \
        --outdir ref_index/ \
        --index-type ivf --nlist 100

    # then, in pipeline.py:
    python pipeline.py ... --ref-index ref_index/reference.faiss \
                            --ref-metadata ref_index/reference_metadata.json
"""

import argparse
import json
import logging
import random
import sys
from pathlib import Path

import numpy as np
import faiss
from Bio import SeqIO


# ===============================================================================
#  EDIT THESE PATHS -- sab kuch /kaggle/input/ se load hoga, sirf output
#  /kaggle/working/ mein jayega (jaisa aap chahte ho)
# ===============================================================================

# Jis folder mein aapki s2.py (patched DNABERT-S wali) hai -- input dataset ka path.
# Ye sys.path mein add hoga taaki "import s2" kaam kare (input dirs Python ke
# module search path mein by default nahi hoti).
S2_MODULE_DIR = "/kaggle/input/datasets/divyashkigf/s2data"
REF_FASTA = "/kaggle/input/datasets/krishs23/sih-database/SIH-databse/pr2_version_5.1.1_SSU_mothur.fasta"
REF_TAXONOMY = "/kaggle/input/datasets/krishs23/sih-database/SIH-databse/pr2_version_5.1.1_SSU_mothur.tax"
OUTPUT_DIR = "/kaggle/working/ref_index"

MAX_REF_SEQUENCES = 240000   # cap so a large DB doesn't need full embedding
INDEX_TYPE = "ivf"
NLIST = 100

# LOCAL patched DNABERT-S model folder -- SAME path jo s2.py use karti hai.
# Isse config.json fold mila. Zaroori hai warna default HF Hub id use hoga,
# poora fresh download hoga, aur unpatched Triton trans_b crash phir aayega.
DNABERT_MODEL_PATH = "/kaggle/input/datasets/anshikabarai/anshika/dnabert-s"
# ===============================================================================

sys.path.insert(0, S2_MODULE_DIR)
try:
    from s2 import dnabert_embedding, kmer_embedding
except ImportError as exc:
    raise ImportError(
        f"'s2' module nahi mila S2_MODULE_DIR='{S2_MODULE_DIR}' mein. "
        f"Is script ke top mein S2_MODULE_DIR ko us dataset/folder ke exact "
        f"path se update karo jahan aapki s2.py hai (Kaggle 'Add Input' se "
        f"jo path mila wahi -- e.g. /kaggle/input/<dataset-name>)."
    ) from exc


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build a FAISS reference index from a taxonomy database")
    p.add_argument("--ref-fasta", default=REF_FASTA, type=Path,
                    help="Reference database FASTA (e.g. PR2/SILVA SSU sequences)")
    p.add_argument("--ref-taxonomy", default=REF_TAXONOMY, type=Path,
                    help="TSV mapping sequence ID -> taxonomy lineage string (id<TAB>lineage)")
    p.add_argument("--outdir", default=OUTPUT_DIR, type=Path,
                    help="Directory to write the FAISS index + metadata")
    p.add_argument("--max-ref-sequences", type=int, default=MAX_REF_SEQUENCES,
                    help="Randomly sample at most this many reference sequences "
                         "(keeps download/embedding cost bounded on large DBs). "
                         "Pass 0 to use every sequence.")
    p.add_argument("--sample-seed", type=int, default=42,
                    help="Random seed for --max-ref-sequences sampling (reproducibility)")
    p.add_argument("--no-dnabert", action="store_true",
                    help="Use k-mer fallback embedding instead of DNABERT-S "
                         "(must match the flag used later in pipeline.py)")
    p.add_argument("--dnabert-model", default=DNABERT_MODEL_PATH,
                    help="Local model path (or HF Hub id if you want a fresh download)")
    p.add_argument("--kmer-size", type=int, default=6)
    p.add_argument("--index-type", choices=["ivf", "hnsw"], default=INDEX_TYPE,
                    help="FAISS index type: 'ivf' (IndexIVFFlat) or 'hnsw' (IndexHNSWFlat)")
    p.add_argument("--nlist", type=int, default=NLIST,
                    help="Number of IVF clusters (only used for --index-type ivf; "
                         "rule of thumb: ~sqrt(num_reference_sequences))")
    p.add_argument("--hnsw-m", type=int, default=32,
                    help="HNSW graph connectivity (only used for --index-type hnsw)")
    return p.parse_args()


def load_taxonomy(path: Path) -> dict:
    tax_map = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                tax_map[parts[0]] = parts[1]
    return tax_map


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                         datefmt="%H:%M:%S")
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    logging.info("Loading reference FASTA: %s", args.ref_fasta)
    records = list(SeqIO.parse(str(args.ref_fasta), "fasta"))
    logging.info("Reference FASTA has %d total sequences", len(records))

    # ---- Cap the reference set so a large DB doesn't force a full download/embed ----
    if args.max_ref_sequences and args.max_ref_sequences > 0 and len(records) > args.max_ref_sequences:
        logging.info("Sampling %d sequences out of %d (--max-ref-sequences cap, seed=%d)",
                     args.max_ref_sequences, len(records), args.sample_seed)
        random.seed(args.sample_seed)
        records = random.sample(records, args.max_ref_sequences)

    seq_ids = [r.id for r in records]
    sequences = [str(r.seq) for r in records]
    logging.info("Using %d reference sequences for this index", len(sequences))

    tax_map = load_taxonomy(args.ref_taxonomy)

    logging.info("Generating embeddings for reference sequences (one-time cost)...")
    if args.no_dnabert:
        embeddings = kmer_embedding(sequences, k=args.kmer_size)
    else:
        try:
            embeddings = dnabert_embedding(sequences, args.dnabert_model)
        except Exception as exc:
            import traceback
            logging.warning("DNABERT-S failed: %s: %r", type(exc).__name__, exc)
            traceback.print_exc()
            logging.warning("Falling back to k-mer embedding.")
            embeddings = kmer_embedding(sequences, k=args.kmer_size)

    embeddings = np.ascontiguousarray(embeddings.astype(np.float32))
    dim = embeddings.shape[1]

    # normalize so inner-product search == cosine similarity
    faiss.normalize_L2(embeddings)

    logging.info("Building FAISS %s index (dim=%d, n=%d)...", args.index_type, dim, len(embeddings))
    if args.index_type == "ivf":
        # IVF training needs nlist <= number of training vectors (roughly nlist*39
        # points minimum is FAISS's own recommendation) -- cap automatically so a
        # small reference set never crashes the index build
        safe_nlist = min(args.nlist, max(1, len(embeddings) // 39))
        if safe_nlist < args.nlist:
            logging.warning("Reducing --nlist from %d to %d (too few reference "
                            "sequences for the requested nlist)", args.nlist, safe_nlist)
        quantizer = faiss.IndexFlatIP(dim)
        index = faiss.IndexIVFFlat(quantizer, dim, safe_nlist, faiss.METRIC_INNER_PRODUCT)
        index.train(embeddings)
        index.add(embeddings)
        index.nprobe = max(1, args.nlist // 10)
    else:  # hnsw
        index = faiss.IndexHNSWFlat(dim, args.hnsw_m, faiss.METRIC_INNER_PRODUCT)
        index.hnsw.efConstruction = 200
        index.add(embeddings)

    index_path = args.outdir / "reference.faiss"
    faiss.write_index(index, str(index_path))

    metadata = [
        {"seq_id": sid, "lineage": tax_map.get(sid, "Unknown lineage")}
        for sid in seq_ids
    ]
    with open(args.outdir / "reference_metadata.json", "w") as fh:
        json.dump(metadata, fh, indent=2)

    config = {
        "index_type": args.index_type,
        "embedding_model": "kmer" if args.no_dnabert else args.dnabert_model,
        "kmer_size": args.kmer_size if args.no_dnabert else None,
        "dim": dim,
        "n_sequences": len(sequences),
        "nlist": args.nlist if args.index_type == "ivf" else None,
        "hnsw_m": args.hnsw_m if args.index_type == "hnsw" else None,
    }
    with open(args.outdir / "reference_config.json", "w") as fh:
        json.dump(config, fh, indent=2)

    logging.info("Reference index built: %s", index_path)
    logging.info("Metadata: %s", args.outdir / "reference_metadata.json")
    logging.info("IMPORTANT: use the SAME --dnabert-model / --kmer-size (or --no-dnabert) "
                 "flags in pipeline.py, otherwise query and reference embeddings won't "
                 "live in the same vector space.")


if __name__ == "__main__":
    main()