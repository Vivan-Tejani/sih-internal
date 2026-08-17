%%writefile step1.py

#!/usr/bin/env python3
"""
===============================================================================
STEP 1 (Python-only alternative) -- FASTQ Preprocessing WITHOUT R/DADA2
===============================================================================
Agar R/DADA2 install nahi ho pa raha, ye script uska kaam approximate karta
hai pure Python mein:
    1. Quality filtering  -- Phred score ke basis pe low-quality reads hataana
    2. Dereplication       -- identical sequences ko group karke count nikalna
                              (DADA2 ka "denoising" utna sophisticated nahi hai,
                              lekin duplicate/near-duplicate collapsing ka kaam
                              yahi karta hai)
    3. Length filtering     -- bahut chhote/lambe reads hataana

NOTE: Ye asli DADA2 ke jitna accurate error-correction nahi karta (DADA2 ek
statistical error model seekhta hai). Lekin hackathon demo/testing ke liye
ye kaafi hai, aur baad mein chaho to R fix karke DADA2 se replace kar sakte ho
-- Step 2/3 scripts ko koi farak nahi padega, wo sirf FASTA + counts CSV
expect karte hain, chahe wo kahin se bhi aaye ho.

INPUT  : Raw FASTQ file
OUTPUT : asv_sequences.fasta  (unique sequences)
         asv_counts.csv       (ASV_ID, Sequence, Count)

Chalane ka tarika:
    python step1_preprocess_python.py \
        --input raw_sample.fastq \
        --outdir step1_output \
        --min-quality 20 \
        --min-length 100
===============================================================================
"""

import argparse
import logging
from collections import Counter
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pure Python FASTQ preprocessing (DADA2 alternative)")
    p.add_argument("--input", required=True, type=Path, help="Raw FASTQ file")
    p.add_argument("--outdir", required=True, type=Path, help="Output directory")
    p.add_argument("--min-quality", type=int, default=20,
                    help="Minimum average Phred quality score to keep a read")
    p.add_argument("--min-length", type=int, default=100,
                    help="Minimum sequence length to keep")
    p.add_argument("--max-length", type=int, default=2000,
                    help="Maximum sequence length to keep")
    p.add_argument("--min-count", type=int, default=2,
                    help="Drop sequences seen fewer than this many times "
                         "(removes likely one-off sequencing errors)")
    return p.parse_args()


def phred_score(qual_char: str) -> int:
    """FASTQ quality character ko Phred score (number) mein convert karta hai."""
    return ord(qual_char) - 33


def average_quality(qual_string: str) -> float:
    scores = [phred_score(c) for c in qual_string]
    return sum(scores) / len(scores) if scores else 0.0


def parse_fastq(path: Path):
    """FASTQ file ko manually parse karta hai (4 lines per record)."""
    with open(path) as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().strip()
            plus = fh.readline()
            qual = fh.readline().strip()
            yield seq, qual


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
                         datefmt="%H:%M:%S")
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    print("========================================")
    print("STEP 1 (Python-only): FASTQ preprocessing")
    print("========================================")

    if not args.input.exists():
        raise FileNotFoundError(f"Input FASTQ nahi mila: {args.input}")

    logging.info("Reading FASTQ: %s", args.input)

    total_reads = 0
    passed_quality = 0
    passed_length = 0
    sequence_counter = Counter()

    for seq, qual in parse_fastq(args.input):
        total_reads += 1
        if not seq or not qual:
            continue

        seq = seq.upper()

        # Step A: quality filter
        avg_q = average_quality(qual)
        if avg_q < args.min_quality:
            continue
        passed_quality += 1

        # Step B: length filter
        if not (args.min_length <= len(seq) <= args.max_length):
            continue
        passed_length += 1

        # Step C: keep only valid ACGT sequences (drop reads with too many Ns)
        if seq.count("N") / len(seq) > 0.05:
            continue

        sequence_counter[seq] += 1

    logging.info("Total reads       : %d", total_reads)
    logging.info("Passed quality     : %d", passed_quality)
    logging.info("Passed length       : %d", passed_length)
    logging.info("Unique sequences (before min-count filter): %d", len(sequence_counter))

    # Step D: dereplication + minimum count filter (removes likely one-off errors)
    filtered = {seq: count for seq, count in sequence_counter.items() if count >= args.min_count}
    logging.info("Unique sequences (after min-count=%d filter): %d", args.min_count, len(filtered))

    if len(filtered) == 0:
        logging.warning("Koi sequence bacha nahi filters ke baad! --min-count ya "
                         "--min-quality kam karke try karo.")

    # sort by abundance, descending
    sorted_seqs = sorted(filtered.items(), key=lambda x: x[1], reverse=True)

    asv_ids = [f"ASV_{i+1:04d}" for i in range(len(sorted_seqs))]
    sequences = [s for s, _ in sorted_seqs]
    counts = [c for _, c in sorted_seqs]

    # write FASTA
    fasta_path = args.outdir / "asv_sequences.fasta"
    with open(fasta_path, "w") as fh:
        for aid, seq in zip(asv_ids, sequences):
            fh.write(f">{aid}\n{seq}\n")

    # write counts CSV
    counts_path = args.outdir / "asv_counts.csv"
    pd.DataFrame({"ASV_ID": asv_ids, "Sequence": sequences, "Count": counts}).to_csv(
        counts_path, index=False
    )

    print(f"\n== STEP 1 DONE: {len(sorted_seqs)} unique sequences ==")
    print(f"Saved: {fasta_path}")
    print(f"Saved: {counts_path}")
    print("\n>>> Ab Step 2 (step2_embed_and_cluster.py) chalao <<<")


if __name__ == "__main__":
    main()