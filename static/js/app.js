/* app.js — CRISPR Guide RNA Finder  Phase 4 */

// ── State ──────────────────────────────────────────────────────────────────
let allGuides   = [];
let seqStats    = {};
let currentSort = "score";
let currentSeq  = "";
let mlAvailable = false;

// ── Examples ───────────────────────────────────────────────────────────────
const EXAMPLES = {
  short:  "ATCGATCGTGCAGGCTACGGTAGCTATCGATCGATCGGGGCTACGATCGATCAGG",
  medium: "ATCGATCGTGCAGGCTACGGTAGCTATCGATCGATCGGGGCTACGATCGATCAGGCTAGGCTAGGATCGATCGATCGTTGCAGGATCGGCTATCGATCGATCGATCGGGGCTACGATCGATCAGGCTAGGCTAGGATCGATCGATCGTTGCAGGATCGGCT",
  rich:   "GCTAGGCTACGGTAGGCTAGGCTATGCAGGCTACGGTAGGCTAGGCTATCGATCGGGGCTAGGCTATCGATCGGGGCTAGCTAGGCTATCGATCGGGG",
  fasta:  ">BRCA1_sample\nATCGATCGTGCAGGCTACGGTAGCTATCGATCGATCGGGGCTACGATCGATCAGGCTAGGCTAGGATCGATCGATCGTTGCAGGATCGGCTATCGATCGATCGATCGGGG",
};
function loadExample(k) { document.getElementById("dna-input").value = EXAMPLES[k]; clearError(); }

// ── Tabs ───────────────────────────────────────────────────────────────────
document.querySelectorAll(".itab").forEach(btn => {
  btn.addEventListener("click", () => {
    const tab = btn.dataset.tab;
    document.querySelectorAll(".itab").forEach(b => b.classList.remove("active"));
    document.querySelectorAll(".itab-panel").forEach(p => p.classList.add("hidden"));
    btn.classList.add("active");
    document.getElementById("tab-" + tab).classList.remove("hidden");
  });
});

// ── FASTA upload ───────────────────────────────────────────────────────────
const uploadZone = document.getElementById("upload-zone");
const fastaFile  = document.getElementById("fasta-file");
uploadZone.addEventListener("click", () => fastaFile.click());
fastaFile.addEventListener("change", e => { if (e.target.files[0]) showFileName(e.target.files[0].name); });
uploadZone.addEventListener("dragover",  e => { e.preventDefault(); uploadZone.classList.add("drag-over"); });
uploadZone.addEventListener("dragleave", () => uploadZone.classList.remove("drag-over"));
uploadZone.addEventListener("drop", e => {
  e.preventDefault(); uploadZone.classList.remove("drag-over");
  const file = e.dataTransfer.files[0];
  if (file) { fastaFile.files = e.dataTransfer.files; showFileName(file.name); }
});
function showFileName(name) {
  const el = document.getElementById("upload-name");
  el.textContent = "📄 " + name; el.classList.remove("hidden");
}

// ── Analysis ───────────────────────────────────────────────────────────────
async function runAnalysis() {
  const activeTab = document.querySelector(".itab.active").dataset.tab;
  if (activeTab === "upload" && fastaFile.files[0]) await runFastaUpload();
  else await runSequenceAnalysis();
}

async function runSequenceAnalysis() {
  const raw = document.getElementById("dna-input").value.trim();
  if (!raw) { showError("Please enter a DNA sequence."); return; }
  setLoading(true); clearError();

  try {
    const res  = await fetch("/analyze", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        sequence:              raw,
        guide_len:             parseInt(document.getElementById("guide-len").value),
        strand:                document.getElementById("strand").value,
        min_gc:                parseFloat(document.getElementById("min-gc").value),
        max_gc:                parseFloat(document.getElementById("max-gc").value),
        min_score:             parseInt(document.getElementById("min-score").value),
        exclude_high_offtarget: document.getElementById("exclude-offtarget").checked,
        use_ml:                document.getElementById("use-ml")?.checked ?? true,
      }),
    });
    const data = await res.json();
    if (!data.success) { showError(data.error); return; }
    allGuides   = data.guides;
    seqStats    = data.stats;
    mlAvailable = data.ml_available;
    currentSeq  = raw;
    renderResults();
  } catch { showError("Could not reach the server. Is app.py running?"); }
  finally  { setLoading(false); }
}

async function runFastaUpload() {
  const file = fastaFile.files[0];
  if (!file) { showError("Please select a FASTA file."); return; }
  setLoading(true); clearError();

  const fd = new FormData();
  fd.append("file",      file);
  fd.append("guide_len", document.getElementById("guide-len").value);
  fd.append("strand",    document.getElementById("strand").value);
  fd.append("min_score", document.getElementById("min-score").value);
  fd.append("use_ml",    document.getElementById("use-ml")?.checked ? "true" : "false");

  try {
    const res  = await fetch("/upload-fasta", { method: "POST", body: fd });
    const data = await res.json();
    if (!data.success) { showError(data.error); return; }
    allGuides   = data.records.flatMap(r => r.guides);
    seqStats    = data.records[0]?.stats || {};
    mlAvailable = data.ml_available;
    currentSeq  = "";
    renderResults();
  } catch { showError("Upload failed. Is app.py running?"); }
  finally  { setLoading(false); }
}

// ── Render ─────────────────────────────────────────────────────────────────
function renderResults() {
  document.getElementById("empty-state").classList.add("hidden");
  document.getElementById("results-content").classList.remove("hidden");
  renderStats(); renderGuideList(); renderVisualizer();
}

function renderStats() {
  const fwd   = allGuides.filter(g => g.strand === "+").length;
  const top   = allGuides[0]?.score ?? 0;
  const avgGC = allGuides.length ? Math.round(allGuides.reduce((s,g) => s+g.gc_percent,0)/allGuides.length) : 0;
  const avgML = mlAvailable && allGuides.length
    ? (allGuides.reduce((s,g) => s+g.ml_efficiency,0)/allGuides.length).toFixed(2)
    : "N/A";

  document.getElementById("stats-bar").innerHTML = `
    <div class="stat-card"><div class="stat-label">Guides</div><div class="stat-value">${allGuides.length}</div><div class="stat-sub">${fwd}+ / ${allGuides.length-fwd}−</div></div>
    <div class="stat-card"><div class="stat-label">Top score</div><div class="stat-value">${top}</div><div class="stat-sub">/ 80 max</div></div>
    <div class="stat-card"><div class="stat-label">Avg GC%</div><div class="stat-value">${avgGC}%</div><div class="stat-sub">opt. 40–70%</div></div>
    <div class="stat-card"><div class="stat-label">Avg ML eff.</div><div class="stat-value">${avgML}</div><div class="stat-sub">0–1 scale</div></div>
  `;
}

function renderGuideList() {
  const sorted = [...allGuides].sort((a,b) => {
    if (currentSort === "score")      return b.score - a.score;
    if (currentSort === "gc")         return b.gc_percent - a.gc_percent;
    if (currentSort === "pos")        return a.position - b.position;
    if (currentSort === "offtarget")  return a.off_target_hits - b.off_target_hits;
    if (currentSort === "ml")         return b.ml_efficiency - a.ml_efficiency;
    return 0;
  });

  const otClass = { Low:"ot-low", Medium:"ot-medium", High:"ot-high" };
  const strandClass = { "+":"strand-fwd", "-":"strand-rev" };

  document.getElementById("guide-list").innerHTML = sorted.map((g, i) => `
    <div class="guide-card">
      <div class="guide-top">
        <div class="guide-rank">${i+1}</div>
        <div class="guide-seq">${g.guide}</div>
        <div class="quality-badge quality-${g.quality}">${g.quality} · ${g.score}</div>
        ${mlAvailable ? `<div class="ml-badge ml-${g.ml_label}">ML: ${g.ml_label}</div>` : ""}
      </div>
      <div class="guide-meta">
        <span>PAM: <span class="val pam-chip">${g.pam}</span></span>
        <span>Pos: <span class="val">${g.position}</span></span>
        <span>GC: <span class="val">${g.gc_percent}%</span></span>
        <span>Strand: <span class="${strandClass[g.strand]||''}">${g.strand==='+'?'(+)':'(−)'}</span></span>
        <span>Off-target: <span class="${otClass[g.off_target_risk]||''}">${g.off_target_hits} · ${g.off_target_risk}</span></span>
      </div>
      <div class="score-bar-wrap">
        <span class="score-label">Rule score</span>
        <div class="score-bar-bg"><div class="score-bar-fill" style="width:${(g.score/80)*100}%"></div></div>
        <span class="score-label">${g.score}/80</span>
      </div>
      ${mlAvailable ? `
      <div class="ml-bar-wrap">
        <span class="ml-label">ML eff.</span>
        <div class="ml-bar-bg"><div class="ml-bar-fill" style="width:${g.ml_efficiency*100}%"></div></div>
        <span class="score-label">${(g.ml_efficiency*100).toFixed(0)}%</span>
      </div>` : ""}
    </div>
  `).join("");
}

function renderVisualizer() {
  if (!currentSeq) { document.getElementById("vis-seq").textContent = "Upload mode — paste a sequence to see visualization."; return; }
  const seq = currentSeq.split("\n").filter(l=>!l.startsWith(">")).join("").replace(/\s/g,"").toUpperCase().slice(0,300);
  const n   = seq.length;
  const marker = new Array(n).fill("");
  allGuides.filter(g=>g.strand==="+").slice(0,5).forEach(g => {
    const s = g.position-1, e = s+g.guide.length;
    for (let i=s; i<Math.min(e,n); i++) marker[i]="guide";
    for (let i=e; i<Math.min(e+3,n); i++) marker[i]="pam";
  });
  let html="";
  for (let i=0;i<n;i++){
    const nt=seq[i];
    if      (marker[i]==="guide") html+=`<span class="nt-guide">${nt}</span>`;
    else if (marker[i]==="pam")   html+=`<span class="nt-pam">${nt}</span>`;
    else                          html+=nt;
    if ((i+1)%60===0) html+="\n";
  }
  document.getElementById("vis-seq").innerHTML=html;
}

// ── Sort ───────────────────────────────────────────────────────────────────
function sortGuides(btn) {
  document.querySelectorAll(".chip").forEach(c=>c.classList.remove("active"));
  btn.classList.add("active");
  currentSort = btn.dataset.sort;
  if (allGuides.length) renderGuideList();
}

// ── Export ─────────────────────────────────────────────────────────────────
function exportCSV() {
  if (!allGuides.length) return;
  const header = "rank,guide,pam,position,strand,gc_percent,score,quality,off_target_hits,off_target_risk,ml_efficiency,ml_label";
  const rows   = allGuides.map((g,i)=>`${i+1},${g.guide},${g.pam},${g.position},${g.strand},${g.gc_percent},${g.score},${g.quality},${g.off_target_hits},${g.off_target_risk},${g.ml_efficiency},${g.ml_label}`);
  downloadFile([header,...rows].join("\n"), "crispr_guides.csv", "text/csv");
}
function exportJSON() {
  if (!allGuides.length) return;
  downloadFile(JSON.stringify({generated:new Date().toISOString(),stats:seqStats,guides:allGuides},null,2), "crispr_guides.json", "application/json");
}
function downloadFile(content, filename, type) {
  const a = Object.assign(document.createElement("a"), {href: URL.createObjectURL(new Blob([content],{type})), download: filename});
  a.click(); URL.revokeObjectURL(a.href);
}

// ── UI helpers ─────────────────────────────────────────────────────────────
function setLoading(on) {
  document.querySelector(".btn-text").classList.toggle("hidden",on);
  document.querySelector(".btn-loading").classList.toggle("hidden",!on);
  document.getElementById("analyze-btn").disabled=on;
}
function showError(msg) { const el=document.getElementById("input-error"); el.textContent="⚠ "+msg; el.classList.remove("hidden"); }
function clearError()   { document.getElementById("input-error").classList.add("hidden"); }
