import React, { useState, useMemo } from 'react';
import { Search, ArrowUpDown, ArrowUp, ArrowDown, ChevronDown } from 'lucide-react';
import umapData from '../data/umap_explorer_data.json';

const CATEGORY_MAP = {
  "Confident match": { label: "Confidently Labeled", color: "#22d3ee" },
  "Candidate novel taxon": { label: "Novel Cluster", color: "#e879f9" },
  "Low confidence / possible divergent": { label: "Disagreement Flag", color: "#fbbf24" },
  "Singleton": { label: "Singleton", color: "#2dd4bf" },
  "Contaminant": { label: "Contaminant", color: "#fb7185" },
};

const RANK_LABELS = ['Domain', 'Supergroup/Clade', 'Clade', 'Kingdom', 'Phylum', 'Subphylum', 'Class', 'Order', 'Family', 'Genus', 'Species'];

const getRank = (lineage) => {
  const parts = lineage.split(';');
  return parts[parts.length - 1].replace(/_/g, ' ');
};

export default function TaxonomyTable() {
  const clusters = umapData.clusters;

  const [search, setSearch] = useState('');
  const [statusFilter, setStatusFilter] = useState('all'); // all | novel | flagged
  const [confBand, setConfBand] = useState('all'); // all | high | mid | low
  const [sortKey, setSortKey] = useState('similarityPct');
  const [sortDir, setSortDir] = useState('desc');

  const confBandOf = (pct) => {
    if (pct >= 80) return 'high';
    if (pct >= 65) return 'mid';
    return 'low';
  };

  const filtered = useMemo(() => {
    let rows = clusters;

    if (search.trim()) {
      const q = search.trim().toLowerCase();
      rows = rows.filter(c =>
        c.clusterId.toLowerCase().includes(q) ||
        c.taxon.toLowerCase().includes(q) ||
        c.lineage.toLowerCase().includes(q)
      );
    }

    if (statusFilter === 'novel') {
      rows = rows.filter(c => c.status === 'Candidate novel taxon');
    } else if (statusFilter === 'flagged') {
      rows = rows.filter(c => c.status === 'Low confidence / possible divergent');
    }

    if (confBand !== 'all') {
      rows = rows.filter(c => confBandOf(c.similarityPct) === confBand);
    }

    const sorted = [...rows].sort((a, b) => {
      let av = a[sortKey];
      let bv = b[sortKey];
      if (typeof av === 'string') {
        av = av.toLowerCase();
        bv = bv.toLowerCase();
        return sortDir === 'asc' ? av.localeCompare(bv) : bv.localeCompare(av);
      }
      return sortDir === 'asc' ? av - bv : bv - av;
    });

    return sorted;
  }, [clusters, search, statusFilter, confBand, sortKey, sortDir]);

  const toggleSort = (key) => {
    if (sortKey === key) {
      setSortDir(prev => (prev === 'asc' ? 'desc' : 'asc'));
    } else {
      setSortKey(key);
      setSortDir('desc');
    }
  };

  const SortIcon = ({ colKey }) => {
    if (sortKey !== colKey) return <ArrowUpDown className="w-3.5 h-3.5 opacity-30" />;
    return sortDir === 'asc' ? <ArrowUp className="w-3.5 h-3.5 text-bio-cyan" /> : <ArrowDown className="w-3.5 h-3.5 text-bio-cyan" />;
  };

  const columns = [
    { key: 'clusterId', label: 'Cluster ID' },
    { key: 'taxon', label: 'Predicted Taxon' },
    { key: 'rank', label: 'Rank', sortable: false },
    { key: 'similarityPct', label: 'Confidence' },
    { key: 'status', label: 'Status' },
    { key: 'stability', label: 'Stability', sortable: false },
    { key: 'reproducibility', label: 'Reproducibility', sortable: false },
    { key: 'numASVs', label: 'ASVs' },
    { key: 'totalReads', label: 'Reads' },
  ];

  return (
    <div className="max-w-[95rem] mx-auto px-8 py-6 h-[calc(100vh-80px)] flex flex-col gap-1">

      {/* Title */}
      <div>
        <h2 className="text-lg font-light tracking-wide text-white">Taxonomy Table</h2>
        <p className="text-gray-400 mt-1">{filtered.length} of {clusters.length} clusters shown</p>
      </div>

      {/* Filters bar */}
      <div className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-2 flex flex-col md:flex-row gap-1 md:items-center">
        <div className="relative flex-1 min-w-[220px]">
          <Search className="w-4 h-4 text-gray-500 absolute left-3 top-1/2 -translate-y-1/2" />
          <input
            type="text"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            placeholder="Search cluster ID, taxon, or lineage..."
            className="w-full h-8 bg-[#050a12]/60 border border-white/10 rounded-lg pl-9 pr-3 text-xs text-gray-200 placeholder-gray-500 focus:outline-none focus:border-bio-cyan/50"
          />
        </div>

        <div className="relative flex-1">
          <select
            value={statusFilter}
            onChange={(e) => setStatusFilter(e.target.value)}
            className="w-full h-8 appearance-none bg-[#050a12]/60 border border-white/10 rounded-lg pl-3 pr-9 text-xs text-gray-200 focus:outline-none focus:border-bio-cyan/50"
          >
            <option value="all">All statuses</option>
            <option value="novel">Novel only</option>
            <option value="flagged">Flagged only</option>
          </select>
          <ChevronDown className="w-4 h-4 text-gray-500 absolute right-3 top-1/2 -translate-y-1/2 pointer-events-none" />
        </div>

        <div className="relative flex-1">
          <select
            value={confBand}
            onChange={(e) => setConfBand(e.target.value)}
            className="w-full h-8 appearance-none bg-[#050a12]/60 border border-white/10 rounded-lg pl-3 pr-9 text-xs text-gray-200 focus:outline-none focus:border-bio-cyan/50"
          >
            <option value="all">All confidence bands</option>
            <option value="high">High (&ge;80%)</option>
            <option value="mid">Mid (65-79%)</option>
            <option value="low">Low (&lt;65%)</option>
          </select>
          <ChevronDown className="w-4 h-4 text-gray-500 absolute right-3 top-1/2 -translate-y-1/2 pointer-events-none" />
        </div>
      </div>

      {/* Table */}
      <div className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl overflow-hidden shadow-lg flex-1 min-h-0 flex flex-col">
        <div className="overflow-x-auto overflow-y-auto flex-1 min-h-0">
          <table className="w-full text-xs">
            <thead className="sticky top-0 bg-[#0a1420]/95 backdrop-blur-md z-10">
              <tr className="border-b border-white/10">
                {columns.map(col => (
                  <th
                    key={col.key}
                    onClick={() => col.sortable !== false && toggleSort(col.key)}
                    className={`text-left px-4 py-2 text-xs uppercase tracking-wide text-gray-400 font-medium whitespace-nowrap ${
                      col.sortable !== false ? 'cursor-pointer hover:text-gray-200 select-none' : ''
                    }`}
                  >
                    <div className="flex items-center gap-1.5">
                      {col.label}
                      {col.sortable !== false && <SortIcon colKey={col.key} />}
                    </div>
                  </th>
                ))}
              </tr>
            </thead>
            <tbody>
              {filtered.map((c) => {
                const info = CATEGORY_MAP[c.status] || CATEGORY_MAP["Low confidence / possible divergent"];
                return (
                  <tr key={c.clusterId} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                    <td className="px-4 py-2 font-mono text-gray-300 whitespace-nowrap">{c.clusterId}</td>
                    <td className="px-4 py-2 text-white whitespace-nowrap">{c.taxon.replace(/_/g, ' ')}</td>
                    <td className="px-4 py-2 text-gray-400 whitespace-nowrap">{getRank(c.lineage)}</td>
                    <td className="px-4 py-2 font-mono whitespace-nowrap" style={{ color: info.color }}>
                      {c.similarityPct}%
                    </td>
                    <td className="px-4 py-2 whitespace-nowrap">
                      <span
                        className="text-xs px-2 py-1 rounded-full border"
                        style={{ color: info.color, borderColor: `${info.color}40`, backgroundColor: `${info.color}15` }}
                      >
                        {info.label}
                      </span>
                    </td>
                    <td className="px-4 py-2 text-gray-500 whitespace-nowrap">
                      <span title="Not yet computed by pipeline">N/A</span>
                    </td>
                    <td className="px-4 py-2 text-gray-500 whitespace-nowrap">
                      <span title="Not yet computed by pipeline">N/A</span>
                    </td>
                    <td className="px-4 py-2 font-mono text-gray-300 whitespace-nowrap">{c.numASVs.toLocaleString()}</td>
                    <td className="px-4 py-2 font-mono text-gray-300 whitespace-nowrap">{c.totalReads.toLocaleString()}</td>
                  </tr>
                );
              })}
              {filtered.length === 0 && (
                <tr>
                  <td colSpan={columns.length} className="px-4 py-10 text-center text-gray-500">
                    No clusters match the current filters.
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>
      </div>
    </div>
  );
}