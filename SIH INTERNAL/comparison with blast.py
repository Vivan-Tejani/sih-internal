%%writefile /kaggle/working/compare_annotation_methods.py
#!/usr/bin/env python3
"""
===============================================================================
COMPARISON SCRIPT: BLAST vs FAISS-Top1 vs FAISS+kNN-Crosscheck
===============================================================================
Kaam: SAME cluster representatives + SAME PR2 reference database pe teeno
      annotation methods chalata hai, aur common metrics pe compare karta hai
      -- PPT ke "our approach is better than BLAST" slide ke liye.

Teeno methods:
  1. BLAST (alignment-based, purana/traditional approach)
  2. FAISS Top-1 (embedding-similarity, sirf nearest neighbor -- humara v1)
  3. FAISS + kNN cross-check (embedding-similarity + neighbor agreement -- humara v2/latest)

METRICS measured (bina ground-truth ke, genuinely proxy metrics):
  - Time taken (wall-clock, per method, total + per-cluster average)
  - Coverage % (kitne clusters ko "kuch bhi" match mila vs "No hit")
  - Divergent-zone handling (50-70% similarity zone mein kitne clusters
    BLAST "no hit" bolta hai jabki FAISS meaningful closest-match deta hai
    -- ye PS-15 ka core claim hai)
  - kNN self-consistency rate (sirf method 3 ke liye -- confident matches
    mein se kitne % genuinely neighbor-agreement se bhi confirm hue)

INPUT  : cluster representatives (FASTA + embeddings, Step 2 ka output)
         PR2 FAISS index (build_reference_index.py se)
         PR2 BLAST-formatted DB (makeblastdb se, alag se banana padega)
OUTPUT : comparison_report.csv  -- per-cluster comparison, teeno methods
         comparison_summary.json -- aggregate metrics for the PPT
===============================================================================
"""

import json
import subprocess
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


# ===============================================================================
#  EDIT THESE PATHS
# ===============================================================================
REP_FASTA = "/kaggle/working/step2_output_v2/cluster_representatives_v2.fasta"
REP_EMBEDDINGS = "/kaggle/working/step2_output_v2/representative_embeddings_v2.npy"
REP_IDS = "/kaggle/working/step2_output_v2/representative_ids_v2.csv"

REF_INDEX = "/kaggle/input/datasets/divyashkigf/pr2ref/pr2ref/reference.faiss"
REF_METADATA = "/kaggle/input/datasets/divyashkigf/pr2ref/pr2ref/reference_metadata.json"

BLAST_DB = "/kaggle/working/pr2_blast_db/pr2"
BLASTN_PATH = "blastn"

CONFIDENT_THRESHOLD = 0.80
UNCERTAIN_THRESHOLD = 0.65
K_NEIGHBORS = 5

OUTPUT_DIR = "/kaggle/working/comparison_output"
# ===============================================================================

Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)


def classify(pct: float, confident: float, uncertain: float) -> str:
    if pct >= confident * 100:
        return "Confident match"
    elif pct >= uncertain * 100:
        return "Low confidence"
    return "No usable match"


def get_genus(lineage: str) -> str:
    parts = lineage.strip(";").split(";")
    return parts[-2] if len(parts) >= 2 else lineage


# -----------------------------------------------------------------------------
# METHOD 1: BLAST (traditional alignment-based)
# -----------------------------------------------------------------------------

def run_method_blast(rep_fasta: Path, cluster_ids: list[str]) -> tuple[pd.DataFrame, float]:
    print("\n--- Method 1: BLAST (alignment-based) ---")
    if BLAST_DB is None:
        print("BLAST_DB not configured -- skipping BLAST method.")
        out = pd.DataFrame({"Cluster_ID": cluster_ids})
        out["BLAST_Match"] = "SKIPPED (no BLAST_DB configured)"
        out["BLAST_Similarity_Pct"] = np.nan
        out["BLAST_Status"] = "N/A"
        return out, 0.0

    n_queries = len(cluster_ids)
    blast_out = Path(OUTPUT_DIR) / "blast_results.tsv"
    # start with a clean file so the line-count watcher below reads only
    # this run's output, not leftovers from a previous run
    if blast_out.exists():
        blast_out.unlink()

    cmd = [BLASTN_PATH, "-query", str(rep_fasta), "-db", str(BLAST_DB),
           "-out", str(blast_out),
           "-outfmt", "6 qseqid sseqid pident length evalue stitle",
           "-max_target_seqs", "1", "-num_threads", "4"]

    # BLAST writes its tabular output incrementally, one line per hit (roughly
    # one line per query since we asked for max_target_seqs=1). We poll the
    # output file's line count in a background thread to drive a live tqdm
    # bar -- this is an approximation (not every query is guaranteed exactly
    # one line, e.g. "no hit" queries write zero lines), but it gives a much
    # better sense of "is this stuck or progressing" than total silence.
    stop_event = threading.Event()
    pbar = tqdm(total=n_queries, desc="BLAST (est. queries completed)", unit="qry")

    def watch_progress():
        last_count = 0
        while not stop_event.is_set():
            if blast_out.exists():
                try:
                    with open(blast_out) as f:
                        count = sum(1 for _ in f)
                except (OSError, UnicodeDecodeError):
                    count = last_count
                if count > last_count:
                    pbar.update(min(count, n_queries) - last_count)
                    last_count = count
            time.sleep(2)

    watcher = threading.Thread(target=watch_progress, daemon=True)

    start = time.perf_counter()
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    watcher.start()
    stdout, stderr = process.communicate()   # blocks here until BLAST actually finishes
    elapsed = time.perf_counter() - start

    stop_event.set()
    watcher.join(timeout=3)
    pbar.n = n_queries   # snap to 100% once truly done, regardless of the estimate
    pbar.refresh()
    pbar.close()

    if process.returncode != 0:
        print(f"BLAST failed: {stderr}")
        out = pd.DataFrame({"Cluster_ID": cluster_ids})
        out["BLAST_Match"] = "ERROR"
        out["BLAST_Similarity_Pct"] = np.nan
        out["BLAST_Status"] = "N/A"
        return out, elapsed

    columns = ["Cluster_ID", "Subject_ID", "BLAST_Similarity_Pct", "Align_Length", "Evalue", "BLAST_Match"]
    if blast_out.stat().st_size == 0:
        hits = pd.DataFrame(columns=columns)
    else:
        hits = pd.read_csv(blast_out, sep="\t", names=columns)
    hits = hits.sort_values("BLAST_Similarity_Pct", ascending=False).drop_duplicates("Cluster_ID")

    merged = pd.DataFrame({"Cluster_ID": cluster_ids}).merge(hits, on="Cluster_ID", how="left")
    # BLAST's real behavior: if no significant alignment found, there's simply
    # NO ROW for that query -- that IS the "no hit" case, not a low score.
    merged["BLAST_Similarity_Pct"] = merged["BLAST_Similarity_Pct"]
    merged["BLAST_Status"] = merged["BLAST_Similarity_Pct"].apply(
        lambda p: "No hit" if pd.isna(p) else classify(p, CONFIDENT_THRESHOLD, UNCERTAIN_THRESHOLD)
    )
    merged["BLAST_Match"] = merged["BLAST_Match"].fillna("No hit")

    print(f"BLAST took {elapsed:.2f}s for {n_queries} clusters")
    return merged[["Cluster_ID", "BLAST_Match", "BLAST_Similarity_Pct", "BLAST_Status"]], elapsed


# -----------------------------------------------------------------------------
# METHOD 2: FAISS Top-1 (embedding-similarity, single nearest neighbor)
# -----------------------------------------------------------------------------

def run_method_faiss_top1(rep_embeddings: np.ndarray, cluster_ids: list[str],
                           index, metadata) -> tuple[pd.DataFrame, float]:
    print("\n--- Method 2: FAISS Top-1 (embedding similarity) ---")
    import faiss

    query = np.ascontiguousarray(rep_embeddings.astype(np.float32))
    faiss.normalize_L2(query)

    n = len(cluster_ids)
    batch_size = max(1, n // 20)   # ~20 progress updates regardless of n
    all_sims, all_idxs = [], []

    start = time.perf_counter()
    for i in tqdm(range(0, n, batch_size), desc="FAISS Top-1 search", unit="batch"):
        sims, idxs = index.search(query[i:i + batch_size], 1)
        all_sims.append(sims)
        all_idxs.append(idxs)
    elapsed = time.perf_counter() - start

    similarities = np.vstack(all_sims)
    indices = np.vstack(all_idxs)

    rows = []
    for i, cid in enumerate(cluster_ids):
        sim_pct = max(0.0, float(similarities[i][0])) * 100.0
        ref_idx = int(indices[i][0])
        lineage = metadata[ref_idx]["lineage"] if 0 <= ref_idx < len(metadata) else "Unknown"
        status = classify(sim_pct, CONFIDENT_THRESHOLD, UNCERTAIN_THRESHOLD)
        rows.append({"Cluster_ID": cid, "FAISS1_Match": lineage,
                     "FAISS1_Similarity_Pct": round(sim_pct, 2), "FAISS1_Status": status})

    print(f"FAISS Top-1 took {elapsed:.4f}s for {n} clusters")
    return pd.DataFrame(rows), elapsed


# -----------------------------------------------------------------------------
# METHOD 3: FAISS + kNN cross-check (our latest/current approach)
# -----------------------------------------------------------------------------

def run_method_faiss_knn(rep_embeddings: np.ndarray, cluster_ids: list[str],
                          index, metadata) -> tuple[pd.DataFrame, float]:
    print("\n--- Method 3: FAISS + kNN cross-check (current approach) ---")
    import faiss

    query = np.ascontiguousarray(rep_embeddings.astype(np.float32))
    faiss.normalize_L2(query)

    n = len(cluster_ids)
    batch_size = max(1, n // 20)
    all_sims, all_idxs = [], []

    start = time.perf_counter()
    for i in tqdm(range(0, n, batch_size), desc="FAISS+kNN search", unit="batch"):
        sims, idxs = index.search(query[i:i + batch_size], K_NEIGHBORS)
        all_sims.append(sims)
        all_idxs.append(idxs)
    elapsed = time.perf_counter() - start

    similarities = np.vstack(all_sims)
    indices = np.vstack(all_idxs)

    rows = []
    for i, cid in enumerate(tqdm(cluster_ids, desc="Computing kNN agreement", unit="cluster")):
        top1_sim = float(similarities[i][0])
        top1_idx = int(indices[i][0])
        sim_pct = max(0.0, top1_sim) * 100.0
        lineage = metadata[top1_idx]["lineage"] if 0 <= top1_idx < len(metadata) else "Unknown"

        neighbor_lineages = [
            metadata[int(indices[i][k])]["lineage"] if 0 <= int(indices[i][k]) < len(metadata) else "Unknown"
            for k in range(K_NEIGHBORS)
        ]
        top1_genus = get_genus(lineage)
        agreement = sum(1 for lin in neighbor_lineages if get_genus(lin) == top1_genus)
        knn_pct = round((agreement / K_NEIGHBORS) * 100.0, 1)

        status = classify(sim_pct, CONFIDENT_THRESHOLD, UNCERTAIN_THRESHOLD)
        if status == "Confident match" and agreement <= 1:
            status = "Low confidence (neighbor disagreement)"

        rows.append({"Cluster_ID": cid, "FAISS2_Match": lineage,
                     "FAISS2_Similarity_Pct": round(sim_pct, 2), "FAISS2_Status": status,
                     "FAISS2_kNN_Agreement_Pct": knn_pct})

    print(f"FAISS + kNN took {elapsed:.4f}s for {n} clusters")
    return pd.DataFrame(rows), elapsed


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    print("Loading cluster representatives...")
    rep_ids_df = pd.read_csv(REP_IDS)
    cluster_ids = rep_ids_df["Cluster_ID"].tolist()
    rep_embeddings = np.load(REP_EMBEDDINGS)
    print(f"{len(cluster_ids)} cluster representatives loaded")

    import faiss
    index = faiss.read_index(REF_INDEX)
    with open(REF_METADATA) as fh:
        metadata = json.load(fh)

    blast_df, blast_time = run_method_blast(Path(REP_FASTA), cluster_ids)
    faiss1_df, faiss1_time = run_method_faiss_top1(rep_embeddings, cluster_ids, index, metadata)
    faiss2_df, faiss2_time = run_method_faiss_knn(rep_embeddings, cluster_ids, index, metadata)

    # merge everything into one comparison table
    comparison = blast_df.merge(faiss1_df, on="Cluster_ID").merge(faiss2_df, on="Cluster_ID")
    comparison_path = Path(OUTPUT_DIR) / "comparison_report.csv"
    comparison.to_csv(comparison_path, index=False)

    # ---- aggregate metrics for the PPT ----
    n = len(cluster_ids)

    def coverage_pct(status_col, no_hit_labels):
        return round(100 * (~comparison[status_col].isin(no_hit_labels)).sum() / n, 1)

    blast_coverage = coverage_pct("BLAST_Status", ["No hit", "N/A"])
    faiss1_coverage = coverage_pct("FAISS1_Status", [])   # FAISS always returns something
    faiss2_coverage = coverage_pct("FAISS2_Status", [])

    # divergent zone: clusters where BLAST says "No hit" but FAISS still gives
    # a similarity-based closest match in the 50-70% range (a real, usable signal)
    if BLAST_DB is not None:
        blast_no_hit = comparison["BLAST_Status"].isin(["No hit"])
        faiss_has_signal = comparison["FAISS1_Similarity_Pct"].between(50, 70)
        divergent_rescued = int((blast_no_hit & faiss_has_signal).sum())
    else:
        divergent_rescued = None

    knn_confident = (comparison["FAISS2_Status"] == "Confident match").sum()
    knn_downgraded = comparison["FAISS2_Status"].str.contains("disagreement", na=False).sum()

    summary = {
        "n_clusters_tested": n,
        "time_seconds": {
            "BLAST": round(blast_time, 3),
            "FAISS_Top1": round(faiss1_time, 4),
            "FAISS_kNN_crosscheck": round(faiss2_time, 4),
        },
        "speedup_vs_blast": {
            "FAISS_Top1": round(blast_time / faiss1_time, 1) if faiss1_time > 0 and BLAST_DB else "N/A (BLAST skipped)",
            "FAISS_kNN_crosscheck": round(blast_time / faiss2_time, 1) if faiss2_time > 0 and BLAST_DB else "N/A (BLAST skipped)",
        },
        "coverage_pct": {
            "BLAST": blast_coverage if BLAST_DB else "N/A (BLAST skipped)",
            "FAISS_Top1": faiss1_coverage,
            "FAISS_kNN_crosscheck": faiss2_coverage,
        },
        "divergent_sequences_rescued_by_FAISS": divergent_rescued,
        "kNN_confident_matches": int(knn_confident),
        "kNN_downgraded_due_to_neighbor_disagreement": int(knn_downgraded),
    }

    with open(Path(OUTPUT_DIR) / "comparison_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    print(json.dumps(summary, indent=2))
    print(f"\nFull comparison table: {comparison_path}")


if __name__ == "__main__":
    main()




    
       
    