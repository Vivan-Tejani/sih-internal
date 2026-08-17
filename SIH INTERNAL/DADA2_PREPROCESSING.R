%%writefile dada2_denoise.R

.libPaths(c("/kaggle/input/datasets/mistridaivya/dada2-rlibs", .libPaths()))
library(dada2)

# ---- Inputs ----
fwd_fastq <- "/kaggle/input/datasets/mistridaivya/edna-trimmed-fastq/edna-trimmed-fastq/ERR3466765_1.fastq"
rev_fastq <- "/kaggle/input/datasets/mistridaivya/edna-trimmed-fastq/edna-trimmed-fastq/ERR3466765_2.fastq"

# ---- Filtered output paths ----
filt_fwd <- "/kaggle/working/ERR3466765_1.filt.fastq.gz"
filt_rev <- "/kaggle/working/ERR3466765_2.filt.fastq.gz"

# ---- Filter: remove N's, additional quality trim ----
cat("Filtering reads (removing N bases, quality trimming)...\n")
filt_out <- filterAndTrim(fwd_fastq, filt_fwd, rev_fastq, filt_rev,
                           maxN = 0, maxEE = c(2, 2), truncQ = 2,
                           rm.phix = TRUE, compress = TRUE, multithread = TRUE)
print(filt_out)

# ---- Learn error rates (on FILTERED files now) ----
cat("Learning forward error model...\n")
errF <- learnErrors(filt_fwd, multithread = TRUE)

cat("Learning reverse error model...\n")
errR <- learnErrors(filt_rev, multithread = TRUE)

# ---- Denoise ----
cat("Denoising forward reads...\n")
dadaF <- dada(filt_fwd, err = errF, multithread = TRUE)

cat("Denoising reverse reads...\n")
dadaR <- dada(filt_rev, err = errR, multithread = TRUE)

# ---- Merge paired reads ----
cat("Merging pairs...\n")
merged <- mergePairs(dadaF, filt_fwd, dadaR, filt_rev, verbose = TRUE)

# ---- Build ASV table ----
seqtab <- makeSequenceTable(merged)
cat("ASV table dimensions:", dim(seqtab), "\n")

# ---- Remove chimeras ----
cat("Removing chimeras...\n")
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)

cat("ASV table after chimera removal:", dim(seqtab_nochim), "\n")
cat("Fraction of reads kept after chimera removal:",
    sum(seqtab_nochim) / sum(seqtab), "\n")

# ---- Save outputs ----
saveRDS(seqtab_nochim, "/kaggle/working/seqtab_nochim_ERR3466765.rds")
write.csv(t(seqtab_nochim), "/kaggle/working/asv_table_ERR3466765.csv")

cat("Done. ASV table saved.\n")