%%writefile /kaggle/working/step3_annotate_species.py 
#!/usr/bin/env python3
"""
===============================================================================
STEP 3 of 3 -- Species Annotation + Biodiversity Report (Pure Python)
===============================================================================
Kaam: Step 2 ke clusters ko ek NAAM (species label) deta hai -- reference
      database (SILVA/PR2) ke saath embedding-space similarity search
      (FAISS) karke. Agar match confident nahi hai to cluster ko
      "Candidate novel taxon" mark kar deta hai -- force-fit NAHI karta.
      Optional: BLAST se cross-check bhi kar sakte ho (secondary check).
      Phir abundance aur biodiversity indices (Shannon/Simpson/Chao1) bhi
      calculate karta hai.

INPUT  : Step 2 ka output -- asv_with_clusters.csv, cluster_representatives.fasta,
                              representative_embeddings.npy
         + precomputed reference index (build_reference_index.py se, ek baar banega)
OUTPUT : final_report.csv          -- final species labels + confidence + abundance
         candidate_novel_taxa.csv  -- sirf un clusters ki list jo "novel" flag hue
         biodiversity_summary.json -- Shannon, Simpson, richness, Chao1

Kaggle mein test karne ka tarika:
    python step3_annotate_species.py \
        --clusters-csv /kaggle/working/step2_output/asv_with_clusters.csv \
        --rep-embeddings /kaggle/working/step2_output/representative_embeddings.npy \
        --rep-ids /kaggle/working/step2_output/representative_ids.csv \
        --ref-index /kaggle/working/ref_index/reference.faiss \
        --ref-metadata /kaggle/working/ref_index/reference_metadata.json \
        --outdir /kaggle/working/step3_output
===============================================================================
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


# ===============================================================================
#  EDIT THESE PATHS -- notebook cell mein bina CLI args ke chalane ke liye
# ===============================================================================
CLUSTERS_CSV = "/kaggle/input/datasets/mistridaivya/step2-output-v2/step2_output_v2/asv_with_clusters_v2.csv"
REP_EMBEDDINGS = "/kaggle/input/datasets/mistridaivya/step2-output-v2/step2_output_v2/representative_embeddings_v2.npy"
REP_IDS = "/kaggle/input/datasets/mistridaivya/step2-output-v2/step2_output_v2/representative_ids_v2.csv"
OUTPUT_DIR = "/kaggle/working/step3_output_v2"

REF_INDEX = "/kaggle/input/datasets/divyashkigf/pr2ref/pr2ref/reference.faiss"
REF_METADATA = "/kaggle/input/datasets/divyashkigf/pr2ref/pr2ref/reference_metadata.json"

CONFIDENT_THRESHOLD = 0.80
UNCERTAIN_THRESHOLD = 0.65

   # ispar se kam = "Candidate novel taxon"
# NOTE: Ye thresholds DNABERT-S cosine-similarity ke actual observed distribution
# se calibrate kiye hain (poore 240k PR2 DB pe run karne ke baad: max similarity
# ~85%, 95th percentile ~80%, median ~72%). BLAST-style % identity intuition
# (95%/80%) yahan directly apply nahi hoti -- embedding-space similarity ka
# natural range alag hota hai, chahe species well-known ho.
INFER_ECOLOGICAL_ROLES = False
# ===============================================================================


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="STEP 3: Annotate clusters + biodiversity report")

    p.add_argument("--clusters-csv", default=CLUSTERS_CSV, type=Path,
                    help="Step 2 ka asv_with_clusters.csv")
    p.add_argument("--rep-embeddings", default=REP_EMBEDDINGS, type=Path,
                    help="Step 2 ka representative_embeddings.npy")
    p.add_argument("--rep-ids", default=REP_IDS, type=Path,
                    help="Step 2 ka representative_ids.csv")
    p.add_argument("--outdir", default=OUTPUT_DIR, type=Path)

    # PRIMARY annotation: FAISS reference index
    p.add_argument("--ref-index", type=Path, default=REF_INDEX,
                    help="build_reference_index.py se bana FAISS index (.faiss file)")
    p.add_argument("--ref-metadata", type=Path, default=REF_METADATA,
                    help="build_reference_index.py se bana reference_metadata.json")
    p.add_argument("--faiss-nprobe", type=int, default=10)

    # SECONDARY, optional: BLAST cross-check
    p.add_argument("--enable-blast-crossvalidation", action="store_true")
    p.add_argument("--blast-db", type=Path, default=None)
    p.add_argument("--blastn-path", default="blastn")
    p.add_argument("--threads", type=int, default=4)

    # thresholds
    p.add_argument("--confident-threshold", type=float, default=CONFIDENT_THRESHOLD,
                    help="Isse zyada similarity = 'Confident match'")
    p.add_argument("--uncertain-threshold", type=float, default=UNCERTAIN_THRESHOLD,
                    help="Isse kam similarity = 'Candidate novel taxon'")

    p.add_argument("--infer-ecological-roles", action="store_true", default=INFER_ECOLOGICAL_ROLES)

    return p.parse_args()




# -----------------------------------------------------------------------------
# PRIMARY: embedding-similarity annotation (FAISS)
# -----------------------------------------------------------------------------

def load_faiss_index(ref_index: Optional[Path], ref_metadata: Optional[Path], nprobe: int):
    if ref_index is None or ref_metadata is None:
        return None, None
    import faiss
    index = faiss.read_index(str(ref_index))
    with open(ref_metadata) as fh:
        metadata = json.load(fh)
    if hasattr(index, "nprobe"):
        index.nprobe = nprobe
    return index, metadata

# === CROSSOVER FEATURE (priority 3, CROSSOVER_FEATURES.md): genus-level agreement helper ===
# Extracts the genus field (second-to-last position) from a PR2-style lineage
# string. Used to compare k-NN neighbor agreement at genus level instead of
# full lineage string, since full-string comparison over-penalizes cases where
# the top hit is an unresolved "_sp." placeholder (different placeholder tags
# within the same genus were being counted as disagreement).
def get_genus(lineage: str) -> str:
    parts = lineage.strip(";").split(";")
    return parts[-2] if len(parts) >= 2 else lineage


def annotate_via_embeddings(cluster_ids: list[str], rep_embeddings: np.ndarray,
                             ref_index, ref_metadata,
                             confident_thresh: float, uncertain_thresh: float) -> pd.DataFrame:
    if ref_index is None:
        out = pd.DataFrame({"Cluster_ID": cluster_ids})
        out["Best_Match"] = "N/A (--ref-index nahi diya gaya)"
        out["Similarity_Pct"] = 0.0
        out["Status"] = "Candidate novel taxon"
        return out

    import faiss
    query = np.ascontiguousarray(rep_embeddings.astype(np.float32))
    faiss.normalize_L2(query)

    # === CROSSOVER FEATURE (priority 3, CROSSOVER_FEATURES.md): k-NN cross-check ===
    # Instead of trusting only the single closest match (k=1), pull the top-5
    # nearest neighbors. This lets us sanity-check the top-1 hit against its
    # neighborhood before accepting it at face value.
    K = 5
    similarities, indices = ref_index.search(query, K)

    rows = []
    for i, cluster_id in enumerate(cluster_ids):
        top1_sim = float(similarities[i][0])
        top1_idx = int(indices[i][0])
        sim_pct = max(0.0, top1_sim) * 100.0
        lineage = (ref_metadata[top1_idx]["lineage"]
                   if 0 <= top1_idx < len(ref_metadata) else "Unknown lineage")

        # === CROSSOVER FEATURE: check agreement across top-5 neighbors ===
        # Count how many of the top-5 neighbors share the same lineage as top-1.
        # If top-1 looks confident but the rest of its neighborhood disagrees,
        # that's a signal the top-1 hit may be a spurious/coincidental match
        # rather than a reliable identification.
        neighbor_lineages = [
            ref_metadata[int(indices[i][k])]["lineage"]
            if 0 <= int(indices[i][k]) < len(ref_metadata) else "Unknown lineage"
            for k in range(K)
        ]
        # Genus-level comparison instead of full-string match (see get_genus above)
        top1_genus = get_genus(lineage)
        neighbor_genera = [get_genus(lin) for lin in neighbor_lineages]
        agreement_count = sum(1 for g in neighbor_genera if g == top1_genus)
        knn_agreement_pct = round((agreement_count / K) * 100.0, 1)

        if sim_pct >= confident_thresh * 100:
            status = "Confident match"
        elif sim_pct >= uncertain_thresh * 100:
            status = "Low confidence / possible divergent"
        else:
            status = "Candidate novel taxon"

        # === CROSSOVER FEATURE: downgrade confident calls with low neighbor agreement ===
        # A "Confident match" whose neighbors mostly disagree is untrustworthy
        # even though its raw similarity score looks good — downgrade it.
        if status == "Confident match" and agreement_count <= 1:
            status = "Low confidence / possible divergent (neighbor disagreement)"

        rows.append({"Cluster_ID": cluster_id, "Best_Match": lineage,
                      "Similarity_Pct": round(sim_pct, 2), "Status": status,
                      "kNN_Agreement_Pct": knn_agreement_pct})

    return pd.DataFrame(rows)

# -----------------------------------------------------------------------------
# SECONDARY, optional: BLAST cross-check
# -----------------------------------------------------------------------------

def run_blast_crosscheck(representatives_fasta: Path, cluster_ids: list[str],
                          blast_db: Path, blastn_path: str, threads: int,
                          outdir: Path, confident_thresh: float, uncertain_thresh: float) -> pd.DataFrame:
    import subprocess
    blast_out = outdir / "blast_results.tsv"
    cmd = [blastn_path, "-query", str(representatives_fasta), "-db", str(blast_db),
           "-out", str(blast_out),
           "-outfmt", "6 qseqid sseqid pident length evalue stitle",
           "-max_target_seqs", "1", "-num_threads", str(threads)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"BLAST fail hua:\n{result.stderr}")

    columns = ["Cluster_ID", "Subject_ID", "BLAST_Similarity_Pct", "Align_Length", "Evalue", "BLAST_Best_Match"]
    if blast_out.stat().st_size == 0:
        hits = pd.DataFrame(columns=columns)
    else:
        hits = pd.read_csv(blast_out, sep="\t", names=columns)
    hits = hits.sort_values("BLAST_Similarity_Pct", ascending=False).drop_duplicates("Cluster_ID")

    merged = pd.DataFrame({"Cluster_ID": cluster_ids}).merge(hits, on="Cluster_ID", how="left")
    merged["BLAST_Similarity_Pct"] = merged["BLAST_Similarity_Pct"].fillna(0.0)
    merged["BLAST_Best_Match"] = merged["BLAST_Best_Match"].fillna("No match found")

    def classify(pct):
        if pct >= confident_thresh * 100:
            return "Confident match"
        elif pct >= uncertain_thresh * 100:
            return "Low confidence / possible divergent"
        return "Candidate novel taxon"

    merged["BLAST_Status"] = merged["BLAST_Similarity_Pct"].apply(classify)
    return merged[["Cluster_ID", "BLAST_Best_Match", "BLAST_Similarity_Pct", "BLAST_Status"]]


# -----------------------------------------------------------------------------
# Abundance + biodiversity indices
# -----------------------------------------------------------------------------

def compute_abundance(asv_df: pd.DataFrame) -> pd.DataFrame:
    agg = asv_df.groupby("Cluster_ID", as_index=False)["Count"].sum()
    agg = agg.rename(columns={"Count": "Total_Reads"})
    total = agg["Total_Reads"].sum()
    agg["Pct_of_Sample"] = (agg["Total_Reads"] / total * 100).round(2) if total else 0.0
    return agg


def compute_diversity_metrics(abundance: pd.DataFrame) -> dict:
    counted = abundance[abundance["Cluster_ID"] != "Unclustered"].copy()
    total = counted["Total_Reads"].sum()
    if total == 0 or len(counted) == 0:
        return {"species_richness": 0, "shannon_index": 0.0,
                "simpson_index": 0.0, "chao1_estimate": 0.0, "total_reads_used": 0}

    p = counted["Total_Reads"] / total
    shannon = float(-(p * np.log(p)).sum())
    simpson = float(1 - (p ** 2).sum())
    richness = int(len(counted))
    f1 = int((counted["Total_Reads"] == 1).sum())
    f2 = int((counted["Total_Reads"] == 2).sum())
    chao1 = richness + (f1 * (f1 - 1)) / (2 * (f2 + 1)) if richness > 0 else 0.0

    return {"species_richness": richness, "shannon_index": round(shannon, 4),
            "simpson_index": round(simpson, 4), "chao1_estimate": round(float(chao1), 2),
            "total_reads_used": int(total)}


ECOLOGICAL_ROLE_KEYWORDS = {
    "Protist": ["protist", "alveolata", "stramenopile", "rhizaria", "excavata"],
    "Cnidarian": ["cnidaria", "anthozoa", "hydrozoa", "scyphozoa", "jellyfish"],
    "Metazoan (other)": ["metazoa", "annelida", "mollusca", "arthropoda", "porifera"],
    "Fungi": ["fungi", "ascomycota", "basidiomycota"],
}


def infer_ecological_role(lineage: str) -> str:
    lineage_lower = str(lineage).lower()
    for role, keywords in ECOLOGICAL_ROLE_KEYWORDS.items():
        if any(kw in lineage_lower for kw in keywords):
            return role
    return "Unclassified / Unknown role"


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                         datefmt="%H:%M:%S")
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("========================================")
    print("STEP 3: Species Annotation + Biodiversity Report")
    print("========================================")

    asv_df = pd.read_csv(args.clusters_csv)
    rep_ids_df = pd.read_csv(args.rep_ids)
    rep_embeddings = np.load(args.rep_embeddings)
    cluster_ids = rep_ids_df["Cluster_ID"].tolist()

    logging.info("Loaded %d ASVs across %d clusters", len(asv_df), len(cluster_ids))

    # ---- PRIMARY: embedding-similarity annotation ----
    logging.info("PRIMARY annotation: FAISS embedding-similarity search ...")
    ref_index_path = args.ref_index if (args.ref_index and Path(args.ref_index).exists()) else None
    ref_metadata_path = args.ref_metadata if (args.ref_metadata and Path(args.ref_metadata).exists()) else None
    if args.ref_index and ref_index_path is None:
        logging.warning("--ref-index path '%s' nahi mila -- annotation skip hogi, "
                        "sab clusters 'Candidate novel taxon' mark honge.", args.ref_index)
    ref_index, ref_metadata = load_faiss_index(ref_index_path, ref_metadata_path, args.faiss_nprobe)
    annotation = annotate_via_embeddings(
        cluster_ids, rep_embeddings, ref_index, ref_metadata,
        args.confident_threshold, args.uncertain_threshold
    )

    # ---- SECONDARY, optional: BLAST cross-check ----
    blast_crosscheck = None
    if args.enable_blast_crossvalidation:
        if args.blast_db is None:
            logging.warning("--enable-blast-crossvalidation diya but --blast-db nahi -- skip.")
        else:
            logging.info("SECONDARY cross-check: BLAST ...")
            rep_fasta = args.rep_ids.parent / "cluster_representatives.fasta"
            blast_crosscheck = run_blast_crosscheck(
                rep_fasta, cluster_ids, args.blast_db, args.blastn_path,
                args.threads, args.outdir, args.confident_threshold, args.uncertain_threshold
            )

    # ---- Abundance + biodiversity ----
    logging.info("Abundance + biodiversity indices calculate kar rahe hain ...")
    abundance = compute_abundance(asv_df)
    diversity = compute_diversity_metrics(abundance)

    with open(args.outdir / "biodiversity_summary.json", "w") as fh:
        json.dump(diversity, fh, indent=2)
    print(f"\nBiodiversity summary: {diversity}\n")

    # ---- Final report ----
    n_per_cluster = (
        asv_df.groupby("Cluster_ID", as_index=False)["ASV_ID"]
              .count().rename(columns={"ASV_ID": "Num_ASVs"})
    )
    report = annotation.merge(abundance, on="Cluster_ID", how="left")
    report = report.merge(n_per_cluster, on="Cluster_ID", how="left")
    if blast_crosscheck is not None:
        report = report.merge(blast_crosscheck, on="Cluster_ID", how="left")
    if args.infer_ecological_roles:
        report["Ecological_Role"] = report["Best_Match"].apply(infer_ecological_role)

    report = report.sort_values("Total_Reads", ascending=False)

    report_path = args.outdir / "final_report.csv"
    report.to_csv(report_path, index=False)

    novel = report[report["Status"] == "Candidate novel taxon"]
    novel_path = args.outdir / "candidate_novel_taxa.csv"
    novel.to_csv(novel_path, index=False)

    print(f"\n== STEP 3 DONE ==")
    print(f"Final report      : {report_path}")
    print(f"Novel taxa flagged: {len(novel)} -> {novel_path}")
    print(">>> Poora pipeline complete! final_report.csv hi aapka main deliverable hai <<<")


if __name__ == "__main__":
    main()