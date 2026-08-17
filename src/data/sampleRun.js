export const SAMPLE_RUN = {
  totalASVs: 19857,
  totalClusters: 597,          // distinct clusters, excludes the "Unclustered" bucket
  totalReads: 157333,
  totalRows: 598,              // 597 clusters + 1 "Unclustered" catch-all row

  statusBreakdown: {
    confidentMatch: { count: 40, pct: 6.7 },
    candidateNovelTaxon: { count: 150, pct: 25.1 },
    lowConfidenceDivergent: { count: 408, pct: 68.2 },
  },

  // Map onto the existing 5 stat cards. NOTE: this dataset only has 3 status
  // categories (Confident match / Candidate novel taxon / Low confidence-divergent).
  // There is no explicit "flagged/disagreement" or "contaminant" category in this
  // CSV — map "Flagged for Review" to lowConfidenceDivergent as the closest proxy,
  // and label it clearly so it's not misrepresented as a QC-flag stage output.
  confidentlyLabeledPct: 6.7,
  novelTaxaPct: 25.1,
  flaggedForReviewPct: 68.2,   // proxy: low-confidence/divergent matches

  diversity: {
    shannon: 4.83,
    simpson: 0.965,
    chao1: 598.0,
    pielou: 0.755,
  },

  similarityRange: { min: 34.23, max: 85.51 },

  topNovelCandidates: [
    { clusterId: "Unclustered", nearestMatch: "Paecilomyces_nostocoides", reads: 15596, similarityPct: 53.17 },
    { clusterId: "Cluster_27", nearestMatch: "Uvigerina_peregrina", reads: 2024, similarityPct: 58.63 },
    { clusterId: "Cluster_289", nearestMatch: "Microplana_sp.", reads: 1860, similarityPct: 44.62 },
    { clusterId: "Cluster_188", nearestMatch: "Aspergillus_niger", reads: 1309, similarityPct: 60.12 },
    { clusterId: "Cluster_392", nearestMatch: "Chytridiales_XX_sp.", reads: 883, similarityPct: 45.27 },
  ],
};
