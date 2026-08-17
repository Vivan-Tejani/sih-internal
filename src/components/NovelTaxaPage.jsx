import React, { useMemo, useState } from 'react';
import { Sparkles, Info, ChevronDown } from 'lucide-react';
import umapData from '../data/umap_explorer_data.json';

export default function NovelTaxaPage() {
  const clusters = umapData.clusters;
  const [sortBy, setSortBy] = useState('reads'); // reads | asvs | confidence
  const [expandedId, setExpandedId] = useState(null);

  const novel = useMemo(() => {
    const filtered = clusters.filter(c => c.status === 'Candidate novel taxon');
    const sorted = [...filtered].sort((a, b) => {
      if (sortBy === 'reads') return b.totalReads - a.totalReads;
      if (sortBy === 'asvs') return b.numASVs - a.numASVs;
      return a.similarityPct - b.similarityPct; // lower similarity = more divergent/novel
    });
    return sorted;
  }, [clusters, sortBy]);

  const totalNovelReads = novel.reduce((sum, c) => sum + c.totalReads, 0);

  return (
    <div className="max-w-[95rem] mx-auto px-8 py-6 flex flex-col gap-6">

      {/* Title */}
      <div>
        <h2 className="text-3xl font-light tracking-wide text-white flex items-center gap-3">
          <Sparkles className="w-7 h-7 text-fuchsia-400" />
          Novel Taxa Candidates
        </h2>
        <p className="text-gray-400 mt-1">
          {novel.length} clusters flagged as candidate novel taxa &middot; {totalNovelReads.toLocaleString()} total reads
        </p>
      </div>

      {/* Sort control */}
      <div className="flex items-center gap-3">
        <span className="text-sm text-gray-400">Sort by:</span>
        <div className="flex gap-2">
          {[
            { key: 'reads', label: 'Read abundance' },
            { key: 'asvs', label: 'ASV count' },
            { key: 'confidence', label: 'Most divergent' },
          ].map(opt => (
            <button
              key={opt.key}
              onClick={() => setSortBy(opt.key)}
              className={`text-sm px-3 py-1.5 rounded-lg border transition-colors ${
                sortBy === opt.key
                  ? 'bg-fuchsia-400/15 border-fuchsia-400/40 text-fuchsia-200'
                  : 'bg-white/5 border-white/10 text-gray-400 hover:text-gray-200'
              }`}
            >
              {opt.label}
            </button>
          ))}
        </div>
      </div>

      {/* Candidate cards */}
      <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
        {novel.map((c, idx) => {
          const isExpanded = expandedId === c.clusterId;
          return (
            <div
              key={c.clusterId}
              className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-5 hover:border-fuchsia-400/30 transition-colors cursor-pointer"
              onClick={() => setExpandedId(isExpanded ? null : c.clusterId)}
            >
              <div className="flex items-start justify-between mb-3">
                <div>
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono px-2 py-0.5 rounded bg-fuchsia-400/10 text-fuchsia-300 border border-fuchsia-400/20">
                      #{idx + 1}
                    </span>
                    <span className="font-mono text-sm text-gray-400">{c.clusterId}</span>
                  </div>
                </div>
                <ChevronDown className={`w-4 h-4 text-gray-500 transition-transform ${isExpanded ? 'rotate-180' : ''}`} />
              </div>

              <div className="grid grid-cols-3 gap-2 mb-3">
                <div>
                  <p className="text-xs text-gray-500 uppercase tracking-wide">ASVs</p>
                  <p className="font-mono text-lg text-white">{c.numASVs.toLocaleString()}</p>
                </div>
                <div>
                  <p className="text-xs text-gray-500 uppercase tracking-wide">Reads</p>
                  <p className="font-mono text-lg text-white">{c.totalReads.toLocaleString()}</p>
                </div>
                <div>
                  <p className="text-xs text-gray-500 uppercase tracking-wide">% Sample</p>
                  <p className="font-mono text-lg text-white">{c.pctOfSample}%</p>
                </div>
              </div>

              <div className="border-t border-white/5 pt-3">
                <p className="text-xs text-gray-500 uppercase tracking-wide mb-1">Nearest Known Match</p>
                <p className="text-sm text-gray-200">{c.taxon.replace(/_/g, ' ')}</p>
                <p className="text-xs text-fuchsia-300 font-mono mt-0.5">{c.similarityPct}% similarity</p>
              </div>

              {isExpanded && (
                <div className="border-t border-white/5 mt-3 pt-3">
                  <p className="text-xs text-gray-500 uppercase tracking-wide mb-2">Taxonomic Lineage of Nearest Match</p>
                  <div className="text-xs text-gray-400 flex flex-wrap gap-1 leading-relaxed">
                    {c.lineage.split(';').map((seg, i, arr) => (
                      <React.Fragment key={i}>
                        <span className={i === arr.length - 1 ? 'text-gray-200' : ''}>
                          {seg.replace(/_/g, ' ')}
                        </span>
                        {i < arr.length - 1 && <span className="text-fuchsia-400/40">›</span>}
                      </React.Fragment>
                    ))}
                  </div>
                  <p className="text-xs text-gray-500 mt-2">Representative ASV: <span className="font-mono text-gray-300">{c.representativeASV}</span></p>
                </div>
              )}
            </div>
          );
        })}
        {novel.length === 0 && (
          <p className="text-gray-500 col-span-full text-center py-10">No candidate novel taxa in this run.</p>
        )}
      </div>

    </div>
  );
}