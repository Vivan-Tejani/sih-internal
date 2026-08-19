import React, { useState, useRef, useMemo, useEffect } from 'react';
import { X, Target, Info, MousePointerClick, RefreshCcw } from 'lucide-react';
import umapData from '../data/umap_explorer_data.json';

const CATEGORY_MAP = {
  "Confident match": { label: "Confidently Labeled", color: "#22d3ee" }, // cyan
  "Candidate novel taxon": { label: "Novel Cluster", color: "#e879f9" }, // violet
  "Low confidence / possible divergent": { label: "Disagreement Flag", color: "#fbbf24" }, // amber
  "Singleton": { label: "Singleton", color: "#2dd4bf" }, // teal
  "Contaminant": { label: "Contaminant", color: "#fb7185" }, // rose
};

export default function UmapExplorer() {
  const containerRef = useRef(null);
  
  // Legend toggle state
  const [visibleCategories, setVisibleCategories] = useState({
    "Confident match": true,
    "Candidate novel taxon": true,
    "Low confidence / possible divergent": true,
    "Singleton": true,
    "Contaminant": true,
  });

  // Pan / Zoom state
  const [transform, setTransform] = useState({ x: 0, y: 0, scale: 1 });
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState({ x: 0, y: 0 });

  // Hover & Selection state
  const [hoveredCluster, setHoveredCluster] = useState(null);
  const [selectedCluster, setSelectedCluster] = useState(null);

  // SVG Data dimensions
  const meta = umapData.meta;
  const padding = Math.max(meta.xMax - meta.xMin, meta.yMax - meta.yMin) * 0.1;
  const viewBox = {
    minX: meta.xMin - padding,
    maxX: meta.xMax + padding,
    minY: meta.yMin - padding,
    maxY: meta.yMax + padding,
  };
  const width = viewBox.maxX - viewBox.minX;
  const height = viewBox.maxY - viewBox.minY;

  // Render variables
  const clusters = umapData.clusters;
  
  // Calculate size scale (log scale)
  const maxAsv = Math.max(...clusters.map(c => c.numASVs));
  const minAsv = Math.min(...clusters.map(c => c.numASVs));
  
  const getRadius = (numASVs) => {
    // Log scale clamping between 0.3 and 2.2 in view coordinates
    // We scale this depending on view box width so it maps to roughly 3-22px on screen
    const minR = width * 0.003; 
    const maxR = width * 0.022;
    
    if (numASVs <= 1) return minR;
    const logVal = Math.log(numASVs) / Math.log(maxAsv);
    return minR + logVal * (maxR - minR);
  };

  const handleWheel = (e) => {
    e.preventDefault();
    if (!containerRef.current) return;
    
    const scaleFactor = e.deltaY > 0 ? 0.9 : 1.1;
    setTransform(prev => ({
      ...prev,
      scale: Math.min(Math.max(0.5, prev.scale * scaleFactor), 10)
    }));
  };

  const handleMouseDown = (e) => {
    // Only left click
    if (e.button !== 0) return;
    setIsDragging(true);
    setDragStart({ x: e.clientX - transform.x, y: e.clientY - transform.y });
  };

  const handleMouseMove = (e) => {
    if (!isDragging) return;
    setTransform(prev => ({
      ...prev,
      x: e.clientX - dragStart.x,
      y: e.clientY - dragStart.y,
    }));
  };

  const handleMouseUp = () => {
    setIsDragging(false);
  };

  useEffect(() => {
    if (isDragging) {
      window.addEventListener('mouseup', handleMouseUp);
      window.addEventListener('mousemove', handleMouseMove);
    } else {
      window.removeEventListener('mouseup', handleMouseUp);
      window.removeEventListener('mousemove', handleMouseMove);
    }
    return () => {
      window.removeEventListener('mouseup', handleMouseUp);
      window.removeEventListener('mousemove', handleMouseMove);
    };
  }, [isDragging, dragStart]);

  const resetView = () => {
    setTransform({ x: 0, y: 0, scale: 1 });
  };

  const toggleCategory = (category) => {
    setVisibleCategories(prev => ({
      ...prev,
      [category]: !prev[category]
    }));
  };

  const formatLineage = (lineage) => {
    return lineage.replace(/_/g, ' ').split(';');
  };

  return (
    <div className="max-w-[90rem] mx-auto px-8 py-6 h-[calc(100vh-80px)] flex flex-col gap-1">
      
      {/* Title */}
      <div className="shrink-0">
        <h2 className="text-lg font-light tracking-wide text-white">UMAP Explorer</h2>
        <p className="text-gray-400 mt-1">{meta.totalClusters} clusters &middot; real embedding-space projection</p>
      </div>

      <div className="flex flex-col lg:flex-row gap-2 flex-1 min-h-0">
        
        {/* Main Plot Area */}
        <div className="relative flex-1 bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl overflow-hidden shadow-lg flex items-center justify-center">
          
          {/* Controls */}
          <div className="absolute top-4 right-4 z-10 flex gap-1">
            <button 
              onClick={resetView}
              className="p-2 bg-[#0a1420]/80 hover:bg-[#0a1420] text-gray-300 hover:text-white border border-white/10 rounded-lg backdrop-blur-md transition-colors"
              title="Reset View"
            >
              <RefreshCcw className="w-4 h-4" />
            </button>
          </div>

          {/* SVG Map Container */}
          <div 
            ref={containerRef}
            className="w-full h-full cursor-grab active:cursor-grabbing outline-none"
            onWheel={handleWheel}
            onMouseDown={handleMouseDown}
          >
            <div 
              className="w-full h-full origin-center"
              style={{
                transform: `translate(${transform.x}px, ${transform.y}px) scale(${transform.scale})`,
                transition: isDragging ? 'none' : 'transform 0.1s ease-out'
              }}
            >
              <svg 
                viewBox={`${viewBox.minX} ${viewBox.minY} ${width} ${height}`} 
                className="w-full h-full overflow-visible"
              >
                {/* Y-axis is inverted here: y={-c.y} because standard scatter math goes up, svg goes down */}
                {clusters.filter(c => visibleCategories[c.status]).map(c => {
                  const r = getRadius(c.numASVs);
                  const isHovered = hoveredCluster?.clusterId === c.clusterId;
                  const isSelected = selectedCluster?.clusterId === c.clusterId;
                  const colorInfo = CATEGORY_MAP[c.status] || CATEGORY_MAP["Low confidence / possible divergent"];

                  return (
                    <circle
                      key={c.clusterId}
                      cx={c.x}
                      cy={-c.y}
                      r={isSelected || isHovered ? r * 1.5 : r}
                      fill={colorInfo.color}
                      fillOpacity={isHovered || isSelected ? 1 : 0.7}
                      stroke={colorInfo.color}
                      strokeWidth={isSelected ? width * 0.005 : width * 0.001}
                      className="cursor-pointer transition-all duration-200"
                      style={{
                        filter: isSelected ? `drop-shadow(0 0 4px ${colorInfo.color})` : 'none',
                        transformOrigin: `${c.x}px ${-c.y}px`,
                      }}
                      onMouseEnter={(e) => {
                        const rect = containerRef.current.getBoundingClientRect();
                        setHoveredCluster({
                          ...c,
                          mouseX: e.clientX - rect.left,
                          mouseY: e.clientY - rect.top
                        });
                      }}
                      onMouseLeave={() => setHoveredCluster(null)}
                      onClick={(e) => {
                        e.stopPropagation();
                        setSelectedCluster(c);
                      }}
                    />
                  );
                })}
              </svg>
            </div>
          </div>

          {/* Tooltip */}
          {hoveredCluster && !isDragging && (
            <div 
              className="absolute pointer-events-none bg-[#050a12]/90 backdrop-blur-md border border-bio-cyan/30 p-2 rounded-lg shadow-xl text-xs z-20"
              style={{
                left: hoveredCluster.mouseX + 15,
                top: hoveredCluster.mouseY + 15,
              }}
            >
              <p className="font-mono text-bio-cyan font-semibold">{hoveredCluster.clusterId}</p>
              <p className="text-white mt-1">{hoveredCluster.taxon.replace(/_/g, ' ')}</p>
              <p className="text-gray-400 text-xs mt-1">Sim: {hoveredCluster.similarityPct}%</p>
            </div>
          )}

          {/* Legend Panel (overlay inside map, top left) */}
          <div className="absolute top-4 left-4 z-10 bg-[#0a1420]/80 backdrop-blur-md border border-white/10 rounded-xl p-2 shadow-lg min-w-[260px]">
            <h3 className="text-xs font-medium text-white mb-3 tracking-wide">Category Layers</h3>
            <div className="space-y-3">
              {Object.keys(CATEGORY_MAP).map(status => {
                const count = meta.statusCounts[status] || 0;
                const info = CATEGORY_MAP[status];
                const isVisible = visibleCategories[status];
                
                return (
                  <label key={status} className="flex items-center gap-1 cursor-pointer group">
                    <div className="relative flex items-center justify-center">
                      <input 
                        type="checkbox" 
                        className="sr-only"
                        checked={isVisible}
                        onChange={() => toggleCategory(status)}
                      />
                      <div className={`w-4 h-4 rounded border flex items-center justify-center transition-colors ${
                        isVisible ? 'bg-transparent border-transparent' : 'border-gray-500 bg-transparent'
                      }`}>
                        {isVisible && (
                          <div className="w-full h-full rounded-sm" style={{ backgroundColor: info.color }} />
                        )}
                      </div>
                    </div>
                    <div className={`text-xs transition-colors flex-1 flex justify-between items-center ${isVisible ? 'text-gray-200' : 'text-gray-500 line-through'}`}>
                      <span>{info.label}</span>
                      <span className="font-mono text-xs opacity-70">{count}</span>
                    </div>
                  </label>
                );
              })}
            </div>
          </div>
        </div>

        {/* Side Panel (Right) */}
        <div className="w-full lg:w-[350px] shrink-0 bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-3 shadow-lg flex flex-col relative overflow-hidden transition-all duration-300">
          {!selectedCluster ? (
            <div className="flex-1 flex flex-col items-center justify-center text-center opacity-50">
              <MousePointerClick className="w-12 h-12 text-bio-cyan mb-4 opacity-50" />
              <p className="text-xs text-white font-medium">No Cluster Selected</p>
              <p className="text-xs text-gray-400 mt-2">Click a cluster point on the map to view detailed metadata and lineage.</p>
            </div>
          ) : (
            <div className="animate-in fade-in slide-in-from-right-4 duration-300 h-full flex flex-col">
              <button 
                onClick={() => setSelectedCluster(null)}
                className="absolute top-4 right-4 p-1.5 text-gray-400 hover:text-white hover:bg-white/10 rounded-lg transition-colors"
              >
                <X className="w-5 h-5" />
              </button>

              <div className="mb-6 pr-8">
                <p className="font-mono text-xs text-gray-400 mb-1">{selectedCluster.clusterId}</p>
                <h3 className="text-base font-semibold text-white tracking-tight break-words">
                  {selectedCluster.taxon.replace(/_/g, ' ')}
                </h3>
              </div>

              <div className="bg-[#050a12]/50 border border-white/5 rounded-xl p-2 mb-6">
                <div className="flex justify-between items-end mb-1">
                  <span className="text-xs text-gray-500 uppercase tracking-wide">Confidence</span>
                  <span className="text-xs font-mono text-white" style={{ color: CATEGORY_MAP[selectedCluster.status]?.color }}>
                    {selectedCluster.similarityPct}%
                  </span>
                </div>
                <div className="text-xs text-gray-300 mb-4">{CATEGORY_MAP[selectedCluster.status]?.label}</div>
                
                <div className="space-y-3">
                  <div className="flex justify-between items-center border-t border-white/5 pt-3">
                    <span className="text-xs text-gray-400">ASVs in cluster</span>
                    <span className="font-mono text-white">{selectedCluster.numASVs.toLocaleString()}</span>
                  </div>
                  <div className="flex justify-between items-center border-t border-white/5 pt-3">
                    <span className="text-xs text-gray-400">Total reads</span>
                    <span className="font-mono text-white">{selectedCluster.totalReads.toLocaleString()}</span>
                  </div>
                  <div className="flex justify-between items-center border-t border-white/5 pt-3">
                    <span className="text-xs text-gray-400">% of sample</span>
                    <span className="font-mono text-white">{selectedCluster.pctOfSample}%</span>
                  </div>
                </div>
              </div>

              <div className="mb-6">
                <p className="text-xs text-gray-500 uppercase tracking-wide mb-2">Representative ASV</p>
                <p className="text-xs font-mono text-gray-300 bg-white/5 py-1 px-2 rounded inline-block">
                  {selectedCluster.representativeASV}
                </p>
              </div>

              <div className="mb-6 flex-1">
                <p className="text-xs text-gray-500 uppercase tracking-wide mb-2">Taxonomic Lineage</p>
                <div className="text-xs text-gray-300 flex flex-wrap gap-1 leading-relaxed">
                  {formatLineage(selectedCluster.lineage).map((segment, i, arr) => (
                    <React.Fragment key={i}>
                      <span className={i === arr.length - 1 ? "text-white font-medium" : "opacity-80"}>
                        {segment}
                      </span>
                      {i < arr.length - 1 && <span className="text-bio-cyan/50 mx-1">›</span>}
                    </React.Fragment>
                  ))}
                </div>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}