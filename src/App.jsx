import React, { useState, useRef, useEffect } from 'react';
import { 
  UploadCloud, 
  Dna, 
  Layers, 
  CheckCircle2, 
  Sparkles, 
  AlertTriangle, 
  X,
  Play,
  Database
} from 'lucide-react';
import { SAMPLE_RUN } from './data/sampleRun';
import UmapExplorer from './components/UmapExplorer';
import TaxonomyTable from './components/TaxonomyTable';
import DiversityPage from './components/DiversityPage';
import NovelTaxaPage from './components/NovelTaxaPage';

export default function App() {
  const [file, setFile] = useState(null);
  const [isDragging, setIsDragging] = useState(false);
  const [isProcessing, setIsProcessing] = useState(false);
  const [hasRun, setHasRun] = useState(false);
  const [statsAnimated, setStatsAnimated] = useState(false);
  const [isSampleData, setIsSampleData] = useState(false);
  const [processingMessage, setProcessingMessage] = useState('');
  const [currentPage, setCurrentPage] = useState('overview');

  const fileInputRef = useRef(null);

  const mockStats = {
    totalAsvs: 12847,
    clusters: 342,
    confidentlyLabeled: 10072,
    novelTaxa: 1824,
    flagged: 951
  };

  const currentStats = hasRun ? (
    isSampleData ? {
      totalAsvs: SAMPLE_RUN.totalASVs,
      clusters: SAMPLE_RUN.totalClusters,
      confidentlyLabeled: SAMPLE_RUN.statusBreakdown.confidentMatch.count,
      novelTaxa: SAMPLE_RUN.statusBreakdown.candidateNovelTaxon.count,
      flagged: SAMPLE_RUN.statusBreakdown.lowConfidenceDivergent.count
    } : mockStats
  ) : {
    totalAsvs: 0,
    clusters: 0,
    confidentlyLabeled: 0,
    novelTaxa: 0,
    flagged: 0
  };

  const handleDragOver = (e) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = (e) => {
    e.preventDefault();
    setIsDragging(false);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setIsDragging(false);
    if (e.dataTransfer.files && e.dataTransfer.files[0]) {
      validateAndSetFile(e.dataTransfer.files[0]);
    }
  };

  const handleFileChange = (e) => {
    if (e.target.files && e.target.files[0]) {
      validateAndSetFile(e.target.files[0]);
    }
  };

  const validateAndSetFile = (uploadedFile) => {
    // Simple validation
    const validExtensions = ['.fa', '.fasta', '.csv'];
    const fileName = uploadedFile.name.toLowerCase();
    if (validExtensions.some(ext => fileName.endsWith(ext))) {
      setFile(uploadedFile);
      setHasRun(false); // reset state on new file
      setStatsAnimated(false);
      setIsSampleData(false);
    } else {
      alert("Unsupported file format. Please upload .fasta, .fa, or .csv");
    }
  };

  const removeFile = () => {
    setFile(null);
    setHasRun(false);
    setStatsAnimated(false);
    setIsSampleData(false);
    if (fileInputRef.current) {
      fileInputRef.current.value = "";
    }
  };

  const runPipeline = () => {
    if (!file) return;
    setIsProcessing(true);
    setProcessingMessage('Clustering ASVs & Assigning Taxonomy...');
    setIsSampleData(false);
    
    // Simulate processing delay
    setTimeout(() => {
      setIsProcessing(false);
      setHasRun(true);
      
      // Trigger stat animation after a short delay
      setTimeout(() => {
        setStatsAnimated(true);
      }, 100);
    }, 2500);
  };

  const loadSampleRun = () => {
    setIsProcessing(true);
    setProcessingMessage('Loading sample run...');
    setIsSampleData(true);
    setHasRun(false);
    setStatsAnimated(false);
    
    setTimeout(() => {
      setIsProcessing(false);
      setHasRun(true);
      
      setTimeout(() => {
        setStatsAnimated(true);
      }, 100);
    }, 2500);
  };

  // Helper for displaying file size
  const formatBytes = (bytes, decimals = 2) => {
    if (!+bytes) return '0 Bytes';
    const k = 1024;
    const dm = decimals < 0 ? 0 : decimals;
    const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return `${parseFloat((bytes / Math.pow(k, i)).toFixed(dm))} ${sizes[i]}`;
  };

  return (
    <div className="relative min-h-screen text-gray-200 font-sans selection:bg-bio-cyan/30">
      {/* Background Video */}
      <video
        autoPlay
        muted
        loop
        playsInline
        className="fixed inset-0 w-full h-full object-cover -z-20"
      >
        <source src="/deep-sea-bg.mp4" type="video/mp4" />
      </video>

      {/* Dark Gradient Scrim */}
      <div className="fixed inset-0 bg-gradient-to-t from-[#050a12] via-[#050a12]/80 to-[#050a12]/40 -z-10 pointer-events-none" />

      {/* Header */}
      <header className="sticky top-0 z-50 bg-[#050a12]/50 backdrop-blur-md border-b border-bio-cyan/20 px-8 py-4 flex justify-between items-center">
        <div className="flex items-center gap-4">
          <Dna className="w-8 h-8 text-bio-cyan drop-shadow-[0_0_10px_rgba(34,211,238,0.5)]" />
          <div>
            <h1 className="text-xl font-bold tracking-wider text-white">eDNA ASV PIPELINE</h1>
            <p className="text-xs text-bio-teal/80 uppercase tracking-widest">Deep-Sea Taxonomic Discovery</p>
          </div>
        </div>
        <nav className="hidden md:flex gap-6 text-lg font-medium tracking-wide">
          <button 
            onClick={() => setCurrentPage('overview')}
            className={`transition-colors ${currentPage === 'overview' ? 'text-bio-cyan drop-shadow-[0_0_8px_rgba(34,211,238,0.5)]' : 'text-gray-400 hover:text-gray-200'}`}
          >
            Overview
          </button>
          <button 
            onClick={() => setCurrentPage('umap')}
            className={`transition-colors ${currentPage === 'umap' ? 'text-bio-cyan drop-shadow-[0_0_8px_rgba(34,211,238,0.5)]' : 'text-gray-400 hover:text-gray-200'}`}
          >
            UMAP Explorer
          </button>
          <button 
            onClick={() => setCurrentPage('taxonomy')}
            className={`transition-colors ${currentPage === 'taxonomy' ? 'text-bio-cyan drop-shadow-[0_0_8px_rgba(34,211,238,0.5)]' : 'text-gray-400 hover:text-gray-200'}`}
          >
            Taxonomy
          </button>
          <button 
            onClick={() => setCurrentPage('diversity')}
            className={`transition-colors ${currentPage === 'diversity' ? 'text-bio-cyan drop-shadow-[0_0_8px_rgba(34,211,238,0.5)]' : 'text-gray-400 hover:text-gray-200'}`}
          >
            Diversity
          </button>
          <button 
            onClick={() => setCurrentPage('novel')}
            className={`transition-colors ${currentPage === 'novel' ? 'text-bio-cyan drop-shadow-[0_0_8px_rgba(34,211,238,0.5)]' : 'text-gray-400 hover:text-gray-200'}`}
          >
            Novel Taxa
          </button>
        </nav>
      </header>

      {/* Main Content */}
      {currentPage === 'overview' ? (
      <main className="max-w-7xl mx-auto px-8 py-12 space-y-12">
        
        {/* Hero / Upload Section */}
        <section className="bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-8 shadow-[0_8px_32px_rgba(0,0,0,0.5)] max-w-3xl mx-auto">
          <h2 className="text-3xl font-light tracking-wide text-white mb-6">Upload Sequence Data</h2>
          
          <div 
            className={`border-2 border-dashed rounded-xl p-10 text-center transition-all duration-300 ${
              isDragging 
                ? 'border-bio-cyan bg-bio-cyan/10 shadow-[0_0_20px_rgba(34,211,238,0.2)]' 
                : 'border-bio-cyan/30 hover:border-bio-cyan/60 hover:bg-[#050a12]/40'
            }`}
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
          >
            {!file ? (
              <div className="flex flex-col items-center gap-4 cursor-pointer" onClick={() => fileInputRef.current?.click()}>
                <div className="p-4 bg-bio-cyan/10 rounded-full">
                  <UploadCloud className="w-10 h-10 text-bio-cyan" />
                </div>
                <div>
                  <p className="text-lg font-medium text-gray-200">Drag & drop or click to browse</p>
                  <p className="text-sm text-gray-400 mt-2">Supports FASTA (.fa, .fasta) and CSV ASV tables</p>
                </div>
              </div>
            ) : (
              <div className="flex items-center justify-between bg-[#050a12]/80 p-4 rounded-lg border border-bio-teal/30">
                <div className="flex items-center gap-4">
                  <div className="p-2 bg-bio-teal/20 rounded">
                    <UploadCloud className="w-6 h-6 text-bio-teal" />
                  </div>
                  <div className="text-left">
                    <p className="font-medium text-white truncate max-w-[200px] sm:max-w-xs">{file.name}</p>
                    <p className="text-xs text-gray-400">{formatBytes(file.size)}</p>
                  </div>
                </div>
                <button 
                  onClick={removeFile}
                  className="p-2 text-gray-400 hover:text-red-400 transition-colors"
                  disabled={isProcessing}
                >
                  <X className="w-5 h-5" />
                </button>
              </div>
            )}
            <input 
              type="file" 
              ref={fileInputRef} 
              className="hidden" 
              accept=".fasta,.fa,.csv"
              onChange={handleFileChange}
            />
          </div>

          <div className="mt-8 flex flex-col md:flex-row md:items-center justify-between gap-4">
            <div className="flex items-center gap-2">
              {hasRun && isSampleData ? (
                <>
                  <div className="w-2.5 h-2.5 rounded-full bg-bio-cyan"></div>
                  <p className="text-sm text-bio-cyan">Sample run loaded — {SAMPLE_RUN.totalRows} clusters analyzed.</p>
                </>
              ) : (
                <>
                  <div className={`w-2.5 h-2.5 rounded-full ${!file && !hasRun ? 'bg-gray-600' : 'bg-bio-cyan animate-pulse'}`}></div>
                  <p className="text-sm text-gray-400">
                    {hasRun ? 'Pipeline complete.' : file ? 'Pipeline ready. Awaiting execution.' : 'Awaiting data input.'}
                  </p>
                </>
              )}
            </div>
            
            <div className="flex flex-col items-end gap-2">
              <div className="flex gap-4 items-center">
                <button 
                  onClick={loadSampleRun}
                  disabled={isProcessing}
                  className="flex items-center gap-2 px-4 py-2 text-sm text-gray-400 hover:text-bio-cyan transition-colors"
                >
                  <Database className="w-4 h-4" /> Load Sample Run
                </button>

                <button 
                  onClick={runPipeline}
                  disabled={!file || isProcessing || (hasRun && !isSampleData)}
                  className={`flex items-center gap-2 px-8 py-3 rounded-lg font-medium tracking-wide transition-all duration-300 ${
                    !file 
                      ? 'bg-gray-800 text-gray-500 cursor-not-allowed' 
                      : isProcessing && !isSampleData
                        ? 'bg-bio-cyan/50 text-white cursor-wait animate-pulse'
                        : hasRun && !isSampleData
                          ? 'bg-green-500/20 text-green-400 border border-green-500/50'
                          : 'bg-gradient-to-r from-bio-cyan to-bio-teal text-[#050a12] shadow-[0_0_20px_rgba(34,211,238,0.4)] hover:shadow-[0_0_30px_rgba(34,211,238,0.6)] hover:scale-105'
                  }`}
                >
                  {isProcessing && !isSampleData ? (
                    <>Processing...</>
                  ) : hasRun && !isSampleData ? (
                    <>
                      <CheckCircle2 className="w-5 h-5" /> Pipeline Complete
                    </>
                  ) : (
                    <>
                      <Play className="w-5 h-5" /> Run Pipeline
                    </>
                  )}
                </button>
              </div>
              {!file && <p className="text-xs text-gray-500 pr-4">No data? Explore a real pipeline output ({SAMPLE_RUN.totalASVs.toLocaleString()} ASVs, {SAMPLE_RUN.totalClusters} clusters)</p>}
            </div>
          </div>
          
          {/* Fake Progress Bar */}
          {isProcessing && (
            <div className="mt-6">
              <div className="h-1.5 w-full bg-gray-800 rounded-full overflow-hidden">
                <div className="h-full bg-bio-cyan w-full animate-[progress_2.5s_ease-in-out]" style={{ transformOrigin: 'left' }} />
              </div>
              <p className="text-xs text-bio-cyan mt-2 text-right animate-pulse">{processingMessage}</p>
            </div>
          )}
        </section>

        {/* Stat Cards */}
        <section className={`grid grid-cols-2 md:grid-cols-5 gap-4 transition-all duration-700 ${statsAnimated ? 'opacity-100 translate-y-0' : 'opacity-90 translate-y-2'}`}>
          <StatCard 
            icon={<Layers className="w-5 h-5 text-bio-cyan" />} 
            value={currentStats.totalAsvs} 
            label="Total ASVs" 
            active={statsAnimated}
          />
          <StatCard 
            icon={<Dna className="w-5 h-5 text-bio-teal" />} 
            value={currentStats.clusters} 
            label="Clusters Formed" 
            active={statsAnimated}
          />
          <StatCard 
            icon={<CheckCircle2 className="w-5 h-5 text-green-400" />} 
            value={currentStats.confidentlyLabeled} 
            label="Confidently Labeled" 
            active={statsAnimated}
          />
          <StatCard 
            icon={<Sparkles className="w-5 h-5 text-bio-violet" />} 
            value={currentStats.novelTaxa} 
            label="Novel Taxa Detected" 
            active={statsAnimated}
            glowColor="rgba(167,139,250,0.3)"
          />
          <StatCard 
            icon={<AlertTriangle className="w-5 h-5 text-yellow-400" />} 
            value={currentStats.flagged} 
            label="Flagged for Review" 
            active={statsAnimated}
            glowColor="rgba(250,204,21,0.2)"
          />
        </section>

        <p className="text-xs text-gray-500/70 text-center max-w-4xl mx-auto">
          Status categories reflect this pipeline's 3-tier confidence system (Confident / Divergent / Novel candidate). Disagreement-flag and contaminant-gating stages are computed separately in the QC Review module.
        </p>

        {/* Bottom Row: Diversity & Metadata */}
        <div className="grid md:grid-cols-3 gap-8">
          
          {/* Diversity Snapshot */}
          <section className="md:col-span-1 bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-6 shadow-lg">
            <h3 className="text-lg font-medium text-white mb-6">Diversity Snapshot</h3>
            
            <div className="space-y-6">
              <div>
                <div className="flex justify-between items-end mb-2">
                  <span className="text-sm text-gray-400">Shannon Index</span>
                  <span className="font-mono text-xl text-bio-cyan">{hasRun ? (isSampleData ? SAMPLE_RUN.diversity.shannon : '3.42') : '-.--'}</span>
                </div>
                <div className="h-2 bg-gray-800 rounded-full overflow-hidden">
                  <div 
                    className="h-full bg-gradient-to-r from-bio-teal to-bio-cyan transition-all duration-1000 ease-out" 
                    style={{ width: hasRun ? (isSampleData ? `${(SAMPLE_RUN.diversity.shannon / 5) * 100}%` : '75%') : '0%' }}
                  />
                </div>
              </div>
              
              <div>
                <div className="flex justify-between items-end mb-2">
                  <span className="text-sm text-gray-400">Simpson Index</span>
                  <span className="font-mono text-xl text-bio-teal">{hasRun ? (isSampleData ? SAMPLE_RUN.diversity.simpson : '0.87') : '-.--'}</span>
                </div>
                <div className="h-2 bg-gray-800 rounded-full overflow-hidden">
                  <div 
                    className="h-full bg-gradient-to-r from-bio-cyan to-bio-violet transition-all duration-1000 ease-out delay-150" 
                    style={{ width: hasRun ? (isSampleData ? `${SAMPLE_RUN.diversity.simpson * 100}%` : '87%') : '0%' }}
                  />
                </div>
              </div>
            </div>
          </section>

          {/* Sample Metadata Table */}
          <section className="md:col-span-2 bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-6 shadow-lg overflow-x-auto">
            <h3 className="text-lg font-medium text-white mb-4">Sample Metadata</h3>
            
            <table className="w-full text-left text-sm whitespace-nowrap">
              <thead>
                <tr className="border-b border-bio-cyan/20 text-gray-400">
                  <th className="pb-3 font-medium">Sample ID</th>
                  <th className="pb-3 font-medium">Site</th>
                  <th className="pb-3 font-medium">Depth (m)</th>
                  <th className="pb-3 font-medium">Collection Date</th>
                </tr>
              </thead>
              <tbody className={`transition-opacity duration-500 ${hasRun ? 'opacity-100' : 'opacity-30'}`}>
                <tr className="border-b border-gray-800/50 hover:bg-[#050a12]/50 transition-colors">
                  <td className="py-3 font-mono text-bio-cyan">DS-001</td>
                  <td className="py-3 text-gray-300">Mariana Trench North</td>
                  <td className="py-3 text-gray-300">4200</td>
                  <td className="py-3 text-gray-400">2026-03-12</td>
                </tr>
                <tr className="border-b border-gray-800/50 hover:bg-[#050a12]/50 transition-colors">
                  <td className="py-3 font-mono text-bio-cyan">DS-002</td>
                  <td className="py-3 text-gray-300">Mariana Trench South</td>
                  <td className="py-3 text-gray-300">4350</td>
                  <td className="py-3 text-gray-400">2026-03-14</td>
                </tr>
                <tr className="hover:bg-[#050a12]/50 transition-colors">
                  <td className="py-3 font-mono text-bio-cyan">DS-003</td>
                  <td className="py-3 text-gray-300">Abyssal Plain Alpha</td>
                  <td className="py-3 text-gray-300">3800</td>
                  <td className="py-3 text-gray-400">2026-03-18</td>
                </tr>
              </tbody>
            </table>
            {!hasRun && (
               <div className="absolute inset-0 bg-[#0a1420]/50 backdrop-blur-[2px] flex items-center justify-center rounded-2xl">
                 <p className="text-sm text-bio-cyan">Awaiting data...</p>
               </div>
            )}
          </section>

        </div>
      </main>
      ) : currentPage === 'umap' ? (
        <UmapExplorer />
      ) : currentPage === 'taxonomy' ? (
        <TaxonomyTable />
      ) : currentPage === 'diversity' ? (
        <DiversityPage />
      ) : (
        <NovelTaxaPage />
      )}
      
      <style dangerouslySetInnerHTML={{__html: `
        @keyframes progress {
          0% { transform: scaleX(0); }
          50% { transform: scaleX(0.7); }
          100% { transform: scaleX(1); }
        }
      `}} />
    </div>
  );
}

function StatCard({ icon, value, label, suffix = '', active, glowColor = 'rgba(34,211,238,0.3)' }) {
  // Simple count up effect
  const [displayValue, setDisplayValue] = useState(0);

  useEffect(() => {
    if (active && value > 0) {
      let startTimestamp = null;
      const duration = 1500; // ms
      
      const step = (timestamp) => {
        if (!startTimestamp) startTimestamp = timestamp;
        const progress = Math.min((timestamp - startTimestamp) / duration, 1);
        // easeOutQuart
        const easeProgress = 1 - Math.pow(1 - progress, 4);
        
        // Handle decimals if suffix is %
        if (suffix === '%') {
           setDisplayValue((easeProgress * value).toFixed(1));
        } else {
           setDisplayValue(Math.floor(easeProgress * value));
        }

        if (progress < 1) {
          window.requestAnimationFrame(step);
        } else {
          setDisplayValue(value);
        }
      };
      
      window.requestAnimationFrame(step);
    } else if (!active) {
      setDisplayValue(0);
    }
  }, [active, value, suffix]);

  return (
    <div className={`bg-white/5 backdrop-blur-xl border border-white/10 rounded-2xl p-5 flex flex-col justify-between transition-all duration-500 ${active ? 'shadow-lg' : ''}`}
         style={{ boxShadow: active ? `0 0 20px ${glowColor}` : 'none' }}>
      <div className="flex items-center justify-between mb-4">
        {icon}
      </div>
      <div>
        <div className="font-mono text-2xl font-semibold text-white tracking-tight mb-1">
          {displayValue.toLocaleString()}{suffix}
        </div>
        <div className="text-xs text-gray-400 font-medium tracking-wide uppercase">
          {label}
        </div>
      </div>
    </div>
  );
}