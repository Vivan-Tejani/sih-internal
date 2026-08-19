import React, { useMemo, useState } from 'react';
import { Info } from 'lucide-react';
import umapData from '../data/umap_explorer_data.json';

const CATEGORY_MAP = {
  "Confident match": { label: "Confidently Labeled", color: "#2dd4bf" },      // teal
  "Candidate novel taxon": { label: "Novel Cluster", color: "#0ea5e9" },      // ocean blue
  "Low confidence / possible divergent": { label: "Disagreement Flag", color: "#fbbf24" }, // sandy amber (kept as-is, warm contrast)
  "Singleton": { label: "Singleton", color: "#5eead4" },                     // seafoam
  "Contaminant": { label: "Contaminant", color: "#fb7185" },                 // coral
};

function computeDiversity(clusters) {
  const abund = clusters.map(c => c.totalReads).filter(n => n > 0);
  const total = abund.reduce((a, b) => a + b, 0);
  const S = abund.length; // observed richness (clusters with reads)

  // Shannon (natural log)
  let shannon = 0;
  abund.forEach(n => {
    const p = n / total;
    shannon += -p * Math.log(p);
  });

  // Simpson's Index (dominance) and inverse (diversity form, 1-D)
  let sumPsq = 0;
  abund.forEach(n => {
    const p = n / total;
    sumPsq += p * p;
  });
  const simpsonDominance = sumPsq;
  const simpsonDiversity = 1 - sumPsq;

  // Pielou's evenness J = H / ln(S)
  const pielou = S > 1 ? shannon / Math.log(S) : 0;

  // Chao1 richness estimator: S_obs + (f1^2) / (2*f2)
  // f1 = singletons (clusters with exactly 1 read), f2 = doubletons (exactly 2 reads)
  const f1 = abund.filter(n => n === 1).length;
  const f2 = abund.filter(n => n === 2).length;
  const chao1 = f2 > 0
    ? S + (f1 * f1) / (2 * f2)
    : S + (f1 * (f1 - 1)) / 2; // bias-corrected form when f2 = 0

  return {
    S, total, shannon, simpsonDominance, simpsonDiversity, pielou, chao1, f1, f2,
  };
}

function MetricCard({ label, value, sub, tooltip }) {
  return (
    <div className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-2 flex flex-col gap-1">
      <div className="flex items-center justify-between">
        <span className="text-xs text-gray-400 uppercase tracking-wide">{label}</span>
        {tooltip && (
          <span title={tooltip}>
            <Info className="w-3.5 h-3.5 text-gray-500" />
          </span>
        )}
      </div>
      <div className="font-mono text-lg font-semibold text-white">{value}</div>
      {sub && <div className="text-xs text-gray-500">{sub}</div>}
    </div>
  );
}

export default function DiversityPage() {
  const clusters = umapData.clusters;
  const metrics = useMemo(() => computeDiversity(clusters), [clusters]);

  // Abundance bar data — top 20 clusters by reads for a readable stacked/bar view
  const topClusters = useMemo(() => {
    return [...clusters]
      .filter(c => c.totalReads > 0)
      .sort((a, b) => b.totalReads - a.totalReads)
      .slice(0, 20);
  }, [clusters]);

  const maxReads = topClusters.length ? topClusters[0].totalReads : 1;

  return (
    <div className="max-w-[95rem] mx-auto px-8 py-6 flex flex-col gap-3">

      {/* Title */}
      <div>
        <h2 className="text-lg font-light tracking-wide text-white">Abundance & Diversity</h2>
        <p className="text-gray-400 mt-1">
          Alpha diversity computed from real per-cluster read abundances &middot; single-sample run
        </p>
      </div>

      {/* Alpha diversity metric cards */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-1">
        <MetricCard
          label="Shannon Index (H)"
          value={metrics.shannon.toFixed(3)}
          sub={`${metrics.S} clusters with reads`}
          tooltip="Natural-log Shannon diversity index — higher means more even, more diverse."
        />
        <MetricCard
          label="Simpson's Diversity (1-D)"
          value={metrics.simpsonDiversity.toFixed(3)}
          sub={`Dominance D = ${metrics.simpsonDominance.toFixed(3)}`}
          tooltip="1 minus the probability two random reads belong to the same cluster."
        />
        <MetricCard
          label="Chao1 Richness"
          value={metrics.chao1.toFixed(1)}
          sub={`f1=${metrics.f1}, f2=${metrics.f2} (singletons/doubletons)`}
          tooltip="Estimated true richness accounting for unseen rare taxa."
        />
        <MetricCard
          label="Pielou's Evenness (J)"
          value={metrics.pielou.toFixed(3)}
          sub="0 = uneven, 1 = perfectly even"
          tooltip="Shannon index normalized by max possible diversity (ln richness)."
        />
      </div>

      {/* Abundance bar chart — top 20 clusters */}
      <div className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-3">
        <h3 className="text-xs font-medium text-white mb-1">Top Clusters by Read Abundance</h3>
        <p className="text-xs text-gray-400 mb-6">Top 20 of {metrics.S} clusters with nonzero reads, out of {metrics.total.toLocaleString()} total reads</p>

        <div className="space-y-2.5">
          {topClusters.map(c => {
            const pct = (c.totalReads / maxReads) * 100;
            const info = CATEGORY_MAP[c.status] || CATEGORY_MAP["Low confidence / possible divergent"];
            return (
              <div key={c.clusterId} className="flex items-center gap-1">
                <div className="w-40 shrink-0 text-xs text-gray-300 truncate font-mono" title={c.clusterId}>
                  {c.clusterId}
                </div>
                <div className="flex-1 bg-[#050a12]/50 rounded-md h-5 relative overflow-hidden">
                  <div
                    className="h-full rounded-md transition-all duration-500"
                    style={{ width: `${pct}%`, backgroundColor: info.color, opacity: 0.75 }}
                  />
                </div>
                <div className="w-24 shrink-0 text-right text-xs font-mono text-gray-300">
                  {c.totalReads.toLocaleString()}
                </div>
                <div className="w-32 shrink-0 text-xs text-gray-500 truncate" title={c.taxon}>
                  {c.taxon.replace(/_/g, ' ')}
                </div>
              </div>
            );
          })}
        </div>
      </div>

    </div>
  );
}