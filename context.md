1. Overview / Home

Summary cards: total ASVs processed, # clusters formed, % confidently labeled, % novel, % flagged
Quick diversity snapshot (Shannon/Simpson at a glance)
Sample metadata (site, depth, date) if multi-sample
"Run pipeline" / upload eDNA data entry point

2. UMAP Explorer (this is your centerpiece)

Interactive UMAP scatter, color-coded exactly like you specced: confidently-labeled / disagreement-flag / novel cluster / singleton / contaminant
Click a cluster → side panel with centroid taxon, confidence, member count, stability score
Toggle layers on/off (e.g. hide contaminants)
This is the "wow" visual for judges — worth the most polish time

3. Taxonomy Table

Sortable/filterable table: cluster ID, predicted taxon, rank, confidence, stability score, reproducibility score
Filter by rank-level threshold, confidence band, novel-only, flagged-only

4. Abundance & Diversity

Taxa × sample abundance matrix (heatmap or stacked bar)
Alpha diversity (Shannon, Simpson, Chao1, Pielou) per sample
Beta diversity (Bray-Curtis/Jaccard) — maybe a small dendrogram or PCoA plot
Rarefaction curves

5. Novel Taxa / Candidates

Ranked list of novel clusters by stability + reproducibility
Each candidate: cluster size, nearest known neighbor (from Stage 7 soft annotation), distance metrics
This is your "discovery" showcase page — good for storytelling

6. Disagreement / QC Review

Stage 5 flagged clusters: centroid label vs member-read spread
Contamination-gated reads (Stage 2) with reasons
Lets you demo the "catches Model 1 being confidently wrong" pitch live

7. Benchmark / Comparison

Side-by-side vs QIIME2/BLAST: runtime, taxa recovered, novel taxa found
This is often what actually wins hackathons — quantified superiority claims need a visual home