#!/usr/bin/env python3
"""
===============================================================================
STEP 2 of 3 -- Embedding + Clustering (Pure Python, Deep Learning + HDBSCAN)
===============================================================================
Kaam: Step 1 ke ASV sequences ko numbers (embeddings) mein convert karta hai
      (DNABERT-S pretrained model se, ya k-mer fallback), phir HDBSCAN se
      similar sequences ko groups (clusters) mein baantha hai.
      KOI reference database yahan use NAHI hoti -- sirf ASVs aapas mein
      compare hote hain.

INPUT  : Step 1 ka output -- asv_sequences.fasta + asv_counts.csv
OUTPUT : asv_with_clusters.csv -- har ASV ke saath uska Cluster_ID
         embeddings.npy         -- (bonus) raw embedding vectors, Step 3 reuse karega
         cluster_representatives.fasta -- har cluster ka ek representative sequence

Kaggle mein test karne ka tarika (LOCAL DNABERT-S dataset ke saath):
    python s2.py \
        --asv-fasta /kaggle/working/asv_sequences.fasta \
        --asv-counts /kaggle/working/asv_counts.csv \
        --outdir /kaggle/working/step2_output \
        --dnabert-model /kaggle/input/dnabert-s/dnabert-s

Agar DNABERT-S skip karna ho (GPU/internet na ho to):
    python s2.py ... --no-dnabert

Loading mein kya ho raha hai step-by-step dekhna ho (recommended jab tak
sab kaam na kare), --diagnose true lagao:
    python s2.py ... --diagnose true

--diagnose true karne se:
  - path resolution (kaunsa path try ho raha hai, config.json mila ya nahi)
  - transformers/torch/accelerate versions
  - tokenizer load hua ya nahi
  - model load hone ke baad kitne params meta device par hain (0 chahiye)
  - agar meta params bache reh gaye to unke naam bhi list honge
  - GPU par ek real smoke-test forward pass, taaki Triton kernel crash
    (trans_b= issue) 19k+ sequences ke beech mein nahi, turant pata chale
  sab print hoga, chhoti si sample (5 sequences) par pehle try hoga, phir
  poora fallback/crash decide hoga -- taki har baar 19k+ sequences pe wait
  na karna pade sirf loading dekhne ke liye.

Progress bar: ab embedding loop (dnabert_embedding aur kmer_embedding dono)
  tqdm progress bar dikhata hai terminal mein -- live count, % complete,
  elapsed/ETA time, aur seq/sec speed. Ye batata hai script chal rahi hai,
  stuck nahi hai, chahe 19k+ sequences ho ya sirf 5 (diagnose sample).

Is step ke baad ruk kar dekh lena: cluster table mein kitne clusters bane?
Kitne ASVs "-1"/"Unclustered" hue? Reasonable lag raha hai kya? Tabhi
Step 3 (annotation) pe jaana.

-------------------------------------------------------------------------------
NOTE -- DNABERT-S local loading (Kaggle dataset se) -- version history:
-------------------------------------------------------------------------------
v1: HuggingFace Hub se seedha download -- "Tensor on device meta is not on
    the expected device cpu!" aata tha.
v2: local_files_only=True + low_cpu_mem_usage=False + to_empty() fallback --
    to_empty() ne error hide kar diya lekin RANDOM/UNINITIALIZED weights de
    diye (silently wrong embeddings -- worse than the k-mer fallback).
v3: to_empty() hataya. Root cause: DNABERT-S ka custom trust_remote_code=True
    model class (bert_layers.py) transformers==4.28.0 ke against likha gaya
    tha. Naye transformers (4.45+/5.x) ne from_pretrained ka meta-device init
    badal diya, jo is purani custom class ke saath incompatible hai. Fix:
    transformers==4.40.0 + tokenizers>=0.19,<0.20 par pin, aur
    low_cpu_mem_usage=False + torch_dtype=torch.float32 explicitly diya taaki
    accelerate ka low-memory/meta init path bypass ho.
v4 (current): meta-device fix ke baad, GPU (CUDA) par ek NAYA, alag issue
    mila -- DNABERT-S ke andar bundled flash_attn_triton.py Triton 2.x syntax
    (tl.dot(q, k, trans_b=True)) use karta hai. Newer Triton (3.x, jo Kaggle
    ab default deta hai) ne trans_b kwarg hata diya -- "dot() got an
    unexpected keyword argument 'trans_b'". Ye kernel SIRF CUDA par activate
    hota hai (CPU par model plain PyTorch attention use karta hai), isliye
    pehle CPU-only runs mein ye kabhi nahi aaya. Fix do-step: (a) pehle
    config mein attn_impl-jaisa field dhoondo aur "torch" par force karo
    (MosaicBERT-family convention -- DNABERT-S isi family ka fork hai) taaki
    Triton kernel load hi na ho; (b) agar wo field exist nahi karta ya fork
    usko respect nahi karta (hardcoded Triton path), to ek real smoke-test
    forward pass GPU par chala ke crash detect karo aur model ko CPU par
    reload kar do -- guaranteed-safe path, kyunki Triton kernel CPU par kabhi
    trigger hi nahi hota. --dnabert-attn aur --force-cpu-attn flags se ye
    control hota hai.
===============================================================================
"""

import argparse
import logging
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm


# -----------------------------------------------------------------------------
# Arguments -- sab paths user-controlled
# -----------------------------------------------------------------------------

def str2bool(value: str) -> bool:
    """--diagnose true / false / 1 / 0 / yes / no sab accept karta hai."""
    if isinstance(value, bool):
        return value
    v = value.strip().lower()
    if v in ("yes", "true", "t", "1", "y"):
        return True
    if v in ("no", "false", "f", "0", "n"):
        return False
    raise argparse.ArgumentTypeError(f"Boolean value expected, got '{value}'")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="STEP 2: Embed ASVs + cluster with HDBSCAN")
    p.add_argument("--asv-fasta", required=True, type=Path,
                    help="Step 1 ka FASTA output (asv_sequences.fasta)")
    p.add_argument("--asv-counts", required=True, type=Path,
                    help="Step 1 ka counts CSV output (asv_counts.csv)")
    p.add_argument("--outdir", required=True, type=Path,
                    help="Yahan Step 2 ke saare outputs save honge")

    p.add_argument("--no-dnabert", action="store_true",
                    help="DNABERT-S skip karo, k-mer fallback embedding use karo "
                         "(GPU/internet na ho to)")
    p.add_argument("--dnabert-model", default="zhihan1996/DNABERT-S",
                    help="DNABERT-S path -- apna local folder path do (jahan "
                         "config.json, pytorch_model.bin, etc. saari files ek "
                         "hi folder mein hain) -- e.g. /kaggle/input/dnabert-s/dnabert-s. "
                         "Ya HF Hub model-id bhi de sakte ho (zhihan1996/DNABERT-S) "
                         "agar internet download chahiye. Step 3 mein bhi SAME "
                         "path/id use karna.")
    p.add_argument("--kmer-size", type=int, default=6)

    p.add_argument("--min-cluster-size", type=int, default=5,
                    help="HDBSCAN min_cluster_size -- chhota value = zyada, fine clusters")

    p.add_argument("--diagnose", type=str2bool, default=False,
                    metavar="true|false",
                    help="true karo to DNABERT-S loading ka har step verbose "
                         "print hoga (path resolution, versions, meta-device "
                         "param check, GPU smoke test) aur pehle sirf 5 "
                         "sequences par test hoga before poora dataset embed "
                         "karne. Default: false.")
    p.add_argument("--diagnose-sample-size", type=int, default=5,
                    help="--diagnose true hone par kitni sequences pe pehle "
                         "quick test chalana hai (default 5)")

    p.add_argument("--dnabert-attn", choices=["auto", "torch", "flash"], default="auto",
                    help="GPU par DNABERT-S ke bundled Triton flash-attention kernel "
                         "ka 'trans_b' crash (naye Triton versions mein) handle karne "
                         "ka tarika. auto (default) = pehle config.attn_impl='torch' "
                         "override try karo (Triton kernel load hi na ho), phir ek "
                         "real forward pass se verify karo -- agar phir bhi crash "
                         "hua to model CPU par reload karo (guaranteed-safe, kyunki "
                         "Triton kernel CPU par trigger hi nahi hota). torch = sirf "
                         "config override try karo, CPU fallback mat karo (fail-fast "
                         "agar override kaam nahi karta). flash = kuch mat chhedo "
                         "(is repo ke saath GPU par crash karega).")
    p.add_argument("--force-cpu-attn", action="store_true",
                    help="Config-override attempt skip karo, DNABERT-S ko seedha "
                         "CPU par load/run karo. Use karo agar pehle se pata hai "
                         "ki attn_impl override is fork ke liye kaam nahi karta, "
                         "ya sirf guaranteed-safe path chahiye. Dheema hai GPU se "
                         "lekin Triton bug se bilkul bachta hai.")

    return p.parse_args()


# -----------------------------------------------------------------------------
# Diagnostics helpers
# -----------------------------------------------------------------------------

def print_env_versions() -> None:
    print("---- [DIAGNOSE] Environment versions ----")
    try:
        import torch
        print(f"  torch          : {torch.__version__}  (cuda available: {torch.cuda.is_available()})")
    except Exception as exc:
        print(f"  torch          : import FAILED -- {exc}")
    try:
        import transformers
        print(f"  transformers   : {transformers.__version__}")
    except Exception as exc:
        print(f"  transformers   : import FAILED -- {exc}")
    try:
        import tokenizers
        print(f"  tokenizers     : {tokenizers.__version__}")
    except Exception as exc:
        print(f"  tokenizers     : import FAILED -- {exc}")
    try:
        import accelerate
        print(f"  accelerate     : {accelerate.__version__}")
    except Exception as exc:
        print(f"  accelerate     : not importable -- {exc}")
    try:
        import triton
        print(f"  triton         : {triton.__version__}")
    except Exception as exc:
        print(f"  triton         : not importable -- {exc}")
    print(f"  python         : {sys.version.split()[0]}")
    print("------------------------------------------")


def print_path_resolution(model_name: str, model_path: Path, is_local: bool, has_config: bool) -> None:
    print("---- [DIAGNOSE] Path resolution ----")
    print(f"  requested path/id : {model_name}")
    print(f"  Path.exists()      : {model_path.exists()}")
    print(f"  config.json found  : {has_config}")
    print(f"  resolved is_local  : {is_local}")
    if model_path.exists():
        try:
            found = sorted(p.name for p in model_path.iterdir())
            print(f"  files in folder    : {found}")
        except Exception as exc:
            print(f"  could not list folder contents -- {exc}")
    print("-------------------------------------")


def print_meta_param_report(model) -> int:
    """Har parameter ka device check karta hai. Return: kitne meta device par hain.
    0 hona chahiye -- agar 0 se zyada hai, model REAL weights ke saath load
    NAHI hua (bhale hi koi exception na aaye), aur embeddings garbage honge."""
    meta_params = [name for name, p in model.named_parameters() if p.device.type == "meta"]
    total_params = sum(1 for _ in model.named_parameters())
    print("---- [DIAGNOSE] Meta-device parameter check ----")
    print(f"  total parameters       : {total_params}")
    print(f"  parameters on 'meta'   : {len(meta_params)}  (0 = good, real weights loaded)")
    if meta_params:
        preview = meta_params[:10]
        print(f"  first meta param names : {preview}{' ...' if len(meta_params) > 10 else ''}")
        print("  >>> WARNING: model has meta-device params. If you proceed anyway, "
              "embeddings will be computed from UNINITIALIZED/RANDOM weights, "
              "not real DNABERT-S weights. This will silently corrupt your "
              "clustering results. Do NOT use to_empty() to paper over this -- "
              "fix the load instead. <<<")
    print("--------------------------------------------------")
    return len(meta_params)


# -----------------------------------------------------------------------------
# Embedding functions
# -----------------------------------------------------------------------------

def kmer_embedding(sequences: list[str], k: int = 6) -> np.ndarray:
    """Fallback: sequence ko k-mer frequency vector mein convert karta hai.
    Database use nahi hoti -- sirf sequence ka apna internal pattern."""
    bases = "ACGT"
    all_kmers = ["".join(c) for c in product(bases, repeat=k)]
    kmer_index = {kmer: i for i, kmer in enumerate(all_kmers)}

    vectors = np.zeros((len(sequences), len(all_kmers)), dtype=np.float32)
    for row, seq in enumerate(tqdm(sequences, desc="Embedding ASVs (k-mer)", unit="seq")):
        seq = seq.upper()
        total = 0
        for i in range(len(seq) - k + 1):
            idx = kmer_index.get(seq[i:i + k])
            if idx is not None:
                vectors[row, idx] += 1
                total += 1
        if total > 0:
            vectors[row] /= total
    return vectors


def _resolve_local_path(model_name: str, diagnose: bool) -> tuple[Path, bool]:
    model_path = Path(model_name)
    looks_like_path = str(model_name).startswith("/") or str(model_name).startswith(".")
    has_config = (model_path / "config.json").is_file()
    is_local = model_path.exists() or has_config or looks_like_path

    if diagnose:
        print_path_resolution(model_name, model_path, is_local, has_config)

    if looks_like_path and not model_path.exists() and not has_config:
        raise FileNotFoundError(
            f"'{model_name}' local path jaisa dikh raha hai lekin filesystem par "
            f"nahi mila (ya config.json andar nahi hai). Kaggle notebook mein "
            f"'import os; [print(r) for r,d,f in os.walk(\"/kaggle/input\") if "
            f"\"config.json\" in f]' chala ke exact path verify karo, phir "
            f"--dnabert-model se sahi path do."
        )
    return model_path, is_local


def _load_dnabert_model(model_name: str, is_local: bool, diagnose: bool, args=None):
    """Actual tokenizer + model load. transformers==4.40.0 range ke against
    likha gaya hai -- low_cpu_mem_usage=False + torch_dtype=torch.float32
    explicitly diya hai taaki accelerate ka meta-device/low-memory init
    path bypass ho aur real weights load hon (see module docstring v3).

    DNABERT-S ka bundled flash_attn_triton.py Triton 2.x syntax (trans_b=
    kwarg on tl.dot) use karta hai, jo Triton 3.x (Kaggle current) mein
    hata diya gaya hai. Ye kernel sirf CUDA par activate hota hai -- CPU
    par model plain PyTorch attention use karta hai aur bug hi nahi aata
    (see module docstring v4). Do fixes try karte hain, in order:
      1. config.attn_impl = "torch" (ya similar field) -- MosaicBERT-family
         override, agar ye fork respect karta hai to Triton kernel kabhi
         load hi nahi hoga.
      2. force CPU forward pass -- guaranteed to work, dheema hai.
    """
    import torch
    from transformers import AutoModel, AutoTokenizer, AutoConfig

    dnabert_attn = getattr(args, "dnabert_attn", "auto") if args else "auto"
    force_cpu_attn = getattr(args, "force_cpu_attn", False) if args else False

    logging.info("Loading pretrained model: %s (local=%s)", model_name, is_local)

    tokenizer = AutoTokenizer.from_pretrained(
        model_name,
        trust_remote_code=True,
        local_files_only=is_local,
    )
    if diagnose:
        print(f"  [DIAGNOSE] tokenizer loaded OK: {type(tokenizer).__name__}")

    # ---- attempt 1: config-flag override (skip if forcing CPU or user said flash) ----
    config = None
    if not force_cpu_attn and dnabert_attn in ("auto", "torch"):
        try:
            config = AutoConfig.from_pretrained(
                model_name, trust_remote_code=True, local_files_only=is_local,
            )
            applied = []
            for field in ("attn_impl", "use_flash_attn", "flash_attn"):
                if hasattr(config, field):
                    old = getattr(config, field)
                    new = "torch" if field == "attn_impl" else False
                    setattr(config, field, new)
                    applied.append((field, old, new))
            if diagnose:
                if applied:
                    print(f"  [DIAGNOSE] config overrides applied: {applied}")
                else:
                    print("  [DIAGNOSE] no known attn_impl-style field found on config -- "
                          "this fork likely hardcodes the Triton path; CPU fallback "
                          "(if triggered) will be needed.")
        except Exception as exc:
            if diagnose:
                print(f"  [DIAGNOSE] could not load/modify config separately: {exc}")
            config = None

    def _load_model(target_device: str):
        kwargs = dict(
            trust_remote_code=True,
            local_files_only=is_local,
            low_cpu_mem_usage=False,
            torch_dtype=torch.float32,
        )
        if config is not None:
            kwargs["config"] = config
        m = AutoModel.from_pretrained(model_name, **kwargs)
        n_meta = sum(1 for _, p in m.named_parameters() if p.device.type == "meta")
        if n_meta > 0:
            raise RuntimeError(f"{n_meta} parameters on meta device after load.")
        m.eval()
        m.to(target_device)
        return m

    gpu_available = torch.cuda.is_available()
    device = "cuda" if (gpu_available and not force_cpu_attn) else "cpu"

    model = _load_model(device)
    if diagnose:
        n_meta = print_meta_param_report(model)
        if n_meta > 0:
            raise RuntimeError(
                "Meta-device params present after load -- see report above. "
                "Fix the transformers/accelerate version mismatch before retrying."
            )

    # ---- verify with one real forward pass; if Triton still explodes, retry on CPU ----
    if device == "cuda" and dnabert_attn == "auto":
        try:
            with torch.no_grad():
                test_ids = tokenizer("ACGTACGTACGT", return_tensors="pt")["input_ids"].to(device)
                _ = model(test_ids)
            if diagnose:
                print("  [DIAGNOSE] GPU forward-pass smoke test OK -- Triton kernel not "
                      "triggered (or config override successfully avoided it).")
        except Exception as exc:
            if diagnose:
                print(f"  [DIAGNOSE] GPU forward pass failed ({exc}) -- config override "
                      "did not avoid the Triton kernel for this fork. Falling back to CPU.")
            logging.warning(
                "DNABERT-S Triton kernel incompatibility on GPU (%s); retrying on CPU.", exc
            )
            del model
            torch.cuda.empty_cache()
            device = "cpu"
            model = _load_model(device)

    logging.info("Running on device: %s", device)
    if diagnose:
        print(f"  [DIAGNOSE] model final device: {device}")

    return tokenizer, model, device


def dnabert_embedding(sequences: list[str], model_name: str, diagnose: bool = False, args=None) -> np.ndarray:
    """Pretrained DNABERT-S se embedding nikalta hai -- koi training nahi,
    sirf inference (feature extraction).

    model_name ek LOCAL folder path honi chahiye (jaise Kaggle dataset
    /kaggle/input/dnabert-s/dnabert-s) jisme saari DNABERT-S files hon.
    local_files_only ye path check se decide hota hai -- HF Hub ko touch
    nahi karta agar local mila.

    diagnose=True hone par: versions print, path resolution print, model
    load hone ke baad meta-device param count verify hota hai (0 hona
    chahiye), aur GPU par ek real smoke-test forward pass hota hai (Triton
    trans_b crash turant pakadne ke liye) -- agar 0 se zyada meta params hain
    to explicitly RuntimeError raise hota hai instead of silently continuing
    with random weights.

    Loop par tqdm progress bar hai -- terminal mein live % complete,
    elapsed/ETA time, aur seq/sec dikhega, taaki pata chale script chal
    rahi hai, stuck nahi (especially useful for 19k+ sequence runs)."""
    import torch

    if diagnose:
        print_env_versions()

    model_path, is_local = _resolve_local_path(model_name, diagnose)
    tokenizer, model, device = _load_dnabert_model(model_name, is_local, diagnose, args=args)

    vectors = []
    with torch.no_grad():
        progress = tqdm(sequences, desc=f"Embedding ASVs (DNABERT-S, {device})", unit="seq")
        for i, seq in enumerate(progress):
            try:
                input_ids = tokenizer(seq, return_tensors="pt")["input_ids"].to(device)
                if diagnose and i == 0:
                    print(f"  [DIAGNOSE] about to run forward pass on seq[0], "
                          f"len(seq)={len(seq)}, input_ids.shape={tuple(input_ids.shape)}, "
                          f"device={device}")
                hidden = model(input_ids)[0]
                mean_vec = torch.mean(hidden[0], dim=0)
                vectors.append(mean_vec.cpu().numpy())
                if diagnose and i == 0:
                    print(f"  [DIAGNOSE] first sequence embedding shape: {mean_vec.shape} "
                          f"(DNABERT-S real hidden size should be 768 -- if you see "
                          f"4096 anywhere downstream, that's the k-mer fallback, not this)")
            except Exception as exc:
                if diagnose:
                    print(f"  [DIAGNOSE] forward pass FAILED on sequence index {i} "
                          f"(ASV len={len(seq)}): {type(exc).__name__}: {exc!r}")
                raise
    return np.vstack(vectors)


# -----------------------------------------------------------------------------
# Main logic
# -----------------------------------------------------------------------------

def load_asvs(fasta_path: Path, counts_path: Path) -> pd.DataFrame:
    from Bio import SeqIO
    records = list(SeqIO.parse(str(fasta_path), "fasta"))
    seq_df = pd.DataFrame({
        "ASV_ID": [r.id for r in records],
        "Sequence": [str(r.seq) for r in records],
    })
    counts_df = pd.read_csv(counts_path)[["ASV_ID", "Count"]]
    merged = seq_df.merge(counts_df, on="ASV_ID", how="left")
    merged["Count"] = merged["Count"].fillna(0).astype(int)
    return merged


def run_clustering(embeddings: np.ndarray, min_cluster_size: int) -> np.ndarray:
    import hdbscan
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, metric="euclidean")
    return clusterer.fit_predict(embeddings)


def get_embeddings(sequences: list[str], args: argparse.Namespace) -> np.ndarray:
    """DNABERT-S try karta hai (agar --no-dnabert nahi diya), aur --diagnose
    true hone par pehle chhoti sample par verify karta hai before poora
    dataset embed karna -- taaki loading issue 19k+ sequences ke baad pata
    na chale."""
    if args.no_dnabert:
        logging.info("--no-dnabert diya gaya hai, seedha k-mer embedding use kar rahe hain.")
        return kmer_embedding(sequences, k=args.kmer_size)

    if args.diagnose:
        import traceback
        n = min(args.diagnose_sample_size, len(sequences))
        print(f"\n===== [DIAGNOSE] Quick test on {n} sequences before full run =====")
        try:
            sample_emb = dnabert_embedding(sequences[:n], args.dnabert_model, diagnose=True, args=args)
            print(f"===== [DIAGNOSE] Quick test PASSED -- shape {sample_emb.shape} -- "
                  f"proceeding with full {len(sequences)}-sequence embedding =====\n")
        except Exception as exc:
            # exc ka str() kabhi khaali ho sakta hai (e.g. bare `assert` bina
            # message ke, ya kuch C-extension errors) -- isliye sirf {exc}
            # print karna kaafi nahi hai. Poora traceback + type dikhate hain
            # taaki asli jagah pata chale jahan crash hua.
            print(f"===== [DIAGNOSE] Quick test FAILED: {type(exc).__name__}: {exc!r} =====")
            print("---- [DIAGNOSE] Full traceback ----")
            traceback.print_exc()
            print("------------------------------------")
            print("===== [DIAGNOSE] Falling back to k-mer embedding for the FULL run. "
                  "Fix the error above and re-run with --diagnose true to use real "
                  "DNABERT-S embeddings. =====\n")
            return kmer_embedding(sequences, k=args.kmer_size)

    try:
        return dnabert_embedding(sequences, args.dnabert_model, diagnose=False, args=args)
    except Exception as exc:
        logging.warning("DNABERT-S fail hua (%s). k-mer fallback use kar rahe hain.", exc)
        return kmer_embedding(sequences, k=args.kmer_size)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                         datefmt="%H:%M:%S")
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("========================================")
    print("STEP 2: Embedding + Clustering shuru")
    if args.diagnose:
        print("(diagnose mode: ON -- verbose loading diagnostics)")
    if args.force_cpu_attn:
        print("(force-cpu-attn: ON -- DNABERT-S will run on CPU only)")
    print("========================================")

    logging.info("Loading Step 1 output ...")
    asv_df = load_asvs(args.asv_fasta, args.asv_counts)
    logging.info("Loaded %d ASVs", len(asv_df))

    logging.info("Generating embeddings ...")
    sequences = asv_df["Sequence"].tolist()
    embeddings = get_embeddings(sequences, args)

    logging.info("Embeddings shape: %s", embeddings.shape)
    if embeddings.shape[1] == 4096:
        logging.warning("Embedding dimension is 4096 -- this is the k-mer fallback "
                         "(k=6 -> 4^6=4096), NOT DNABERT-S. DNABERT-S real hidden "
                         "size is 768. Re-run with --diagnose true to see why "
                         "DNABERT-S loading failed.")
    np.save(args.outdir / "embeddings.npy", embeddings)

    logging.info("Clustering with HDBSCAN (min_cluster_size=%d) ...", args.min_cluster_size)
    labels = run_clustering(embeddings, args.min_cluster_size)

    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise = int(np.sum(labels == -1))
    print(f"\n== STEP 2 DONE: {n_clusters} clusters bane, {n_noise} ASVs unclustered/noise ==\n")

    asv_df["Cluster_ID"] = [f"Cluster_{lbl}" if lbl != -1 else "Unclustered" for lbl in labels]

    out_csv = args.outdir / "asv_with_clusters.csv"
    asv_df.to_csv(out_csv, index=False)
    logging.info("Saved: %s", out_csv)

    # cluster representatives (highest-abundance ASV per cluster) -- Step 3 ke liye
    representatives = (
        asv_df.sort_values("Count", ascending=False)
              .groupby("Cluster_ID", as_index=False)
              .first()
    )
    rep_fasta = args.outdir / "cluster_representatives.fasta"
    with open(rep_fasta, "w") as fh:
        for _, row in representatives.iterrows():
            fh.write(f">{row['Cluster_ID']}\n{row['Sequence']}\n")
    logging.info("Saved: %s", rep_fasta)

    # representative embeddings bhi save karo taaki Step 3 dobara embed na kare
    rep_id_to_idx = {aid: i for i, aid in enumerate(asv_df["ASV_ID"])}
    rep_embeddings = np.vstack([embeddings[rep_id_to_idx[aid]] for aid in representatives["ASV_ID"]])
    np.save(args.outdir / "representative_embeddings.npy", rep_embeddings)
    representatives[["Cluster_ID", "ASV_ID"]].to_csv(
        args.outdir / "representative_ids.csv", index=False
    )

    print(">>> Ab cluster table check karo, phir Step 3 (annotate_species.py) chalao <<<")


if __name__ == "__main__":
    main()
