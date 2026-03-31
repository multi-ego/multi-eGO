import { useState, useCallback, useRef } from "react";
import { parseMoleculeFile } from "../utils/parseMolecule.js";

// ─── pdb2gmx option definitions ───────────────────────────────────────────────

const WATER_MODELS = [
  { value: "none", label: "None (no water)" },
];

const CHAINSEP_OPTIONS = [
  { value: "id_or_ter", label: "id_or_ter — new chain on ID change or TER record (recommended)" },
  { value: "id", label: "id — new chain on chain ID change only" },
  { value: "ter", label: "ter — new chain on TER record only" },
  { value: "id_and_ter", label: "id_and_ter — new chain only when both change" },
  { value: "interactive", label: "interactive — ask at runtime" },
];

const DEFAULT_FORM = {
  ff: "multi-ego-basic",
  water: "none",
  outputGro: "conf.gro",
  outputTop: "topol.top",
  outputPosre: "posre.itp",
  chainsep: "id_or_ter",
  ignh: true,
  missing: false,
  mergeAll: false,
  renum: false,
};

function buildCommand(filename, form, parsed) {
  const parts = ["gmx pdb2gmx"];
  parts.push(`-f ${filename}`);
  parts.push(`-o ${form.outputGro}`);
  parts.push(`-p ${form.outputTop}`);
  parts.push(`-i ${form.outputPosre}`);
  parts.push(`-ff ${form.ff}`);
  parts.push(`-water ${form.water}`);
  parts.push(`-chainsep ${form.chainsep}`);
  if (form.ignh) parts.push("-ignh");
  if (form.missing) parts.push("-missing");
  if (form.mergeAll) parts.push("-merge all");
  if (form.renum) parts.push("-renum");
  return parts.join(" \\\n    ");
}

// ─── sub-components ───────────────────────────────────────────────────────────

function SectionTitle({ children }) {
  return (
    <h2 className="text-lg font-semibold text-white border-b border-gray-800 pb-2 mb-4">{children}</h2>
  );
}

function Toggle({ label, description, checked, onChange }) {
  return (
    <label className="flex items-start gap-3 cursor-pointer group">
      <div className="mt-0.5 shrink-0">
        <div
          className={`w-10 h-6 rounded-full transition cursor-pointer ${checked ? "bg-brand-600" : "bg-gray-700"} relative`}
          onClick={() => onChange(!checked)}
        >
          <div
            className={`absolute top-1 h-4 w-4 rounded-full bg-white shadow transition-transform ${
              checked ? "translate-x-5" : "translate-x-1"
            }`}
          />
          <input
            type="checkbox"
            className="sr-only"
            checked={checked}
            onChange={(e) => onChange(e.target.checked)}
          />
        </div>
      </div>
      <div>
        <p className="text-sm font-medium text-gray-200 group-hover:text-white transition">{label}</p>
        {description && <p className="text-xs text-gray-500 mt-0.5">{description}</p>}
      </div>
    </label>
  );
}

function MoleculeTable({ parsed }) {
  const totalResidues = parsed.chains.reduce((n, c) => n + c.residues.length, 0);
  return (
    <div className="space-y-4">
      <div className="grid grid-cols-3 gap-3">
        <div className="card text-center">
          <p className="text-2xl font-bold text-brand-400">{parsed.chains.length}</p>
          <p className="text-xs text-gray-500 mt-1">Chain{parsed.chains.length !== 1 ? "s" : ""}</p>
        </div>
        <div className="card text-center">
          <p className="text-2xl font-bold text-brand-400">{totalResidues}</p>
          <p className="text-xs text-gray-500 mt-1">Residues</p>
        </div>
        <div className="card text-center">
          <p className="text-2xl font-bold text-brand-400">{parsed.atomCount.toLocaleString()}</p>
          <p className="text-xs text-gray-500 mt-1">Atoms</p>
        </div>
      </div>

      {parsed.chains.map((chain, ci) => (
        <div key={ci} className="card">
          <p className="text-sm font-semibold text-gray-300 mb-2">
            Chain <span className="font-mono text-brand-400">{chain.id}</span>
            <span className="ml-2 text-xs text-gray-500">({chain.residues.length} residues)</span>
          </p>
          <div className="flex flex-wrap gap-1 max-h-36 overflow-y-auto">
            {chain.residues.map((r) => (
              <span
                key={`${r.number}:${r.name}`}
                className="inline-block rounded bg-gray-800 px-1.5 py-0.5 font-mono text-xs text-gray-300"
              >
                {r.name}
                <span className="text-gray-600 ml-0.5">{r.number}</span>
              </span>
            ))}
          </div>
        </div>
      ))}

      {parsed.warnings.length > 0 && (
        <div className="space-y-2">
          {parsed.warnings.map((w, i) => (
            <div key={i} className="flex gap-2 rounded-lg border border-yellow-800 bg-yellow-950/40 p-3 text-sm text-yellow-300">
              <span className="shrink-0">⚠️</span>
              <span>{w}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

function DropZone({ onFile }) {
  const [dragging, setDragging] = useState(false);
  const inputRef = useRef();

  const handleFile = (file) => {
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (e) => onFile(file.name, e.target.result);
    reader.readAsText(file);
  };

  const onDrop = (e) => {
    e.preventDefault();
    setDragging(false);
    handleFile(e.dataTransfer.files[0]);
  };

  return (
    <div
      onDragOver={(e) => { e.preventDefault(); setDragging(true); }}
      onDragLeave={() => setDragging(false)}
      onDrop={onDrop}
      onClick={() => inputRef.current.click()}
      className={`flex flex-col items-center justify-center gap-3 rounded-xl border-2 border-dashed p-12 cursor-pointer transition
        ${dragging ? "border-brand-500 bg-brand-950/30" : "border-gray-700 bg-gray-900 hover:border-gray-500"}`}
    >
      <span className="text-4xl">📂</span>
      <p className="text-sm font-medium text-gray-300">Drop a PDB or GRO file here</p>
      <p className="text-xs text-gray-600">or click to browse</p>
      <input
        ref={inputRef}
        type="file"
        accept=".pdb,.ent,.gro"
        className="sr-only"
        onChange={(e) => handleFile(e.target.files[0])}
      />
    </div>
  );
}

// ─── main page ────────────────────────────────────────────────────────────────

export default function SetupAssistant() {
  const [file, setFile] = useState(null); // { name, content }
  const [parsed, setParsed] = useState(null);
  const [parseError, setParseError] = useState(null);
  const [form, setForm] = useState(DEFAULT_FORM);
  const [copied, setCopied] = useState(false);

  const set = useCallback((key, val) => setForm((f) => ({ ...f, [key]: val })), []);

  const handleFile = (name, content) => {
    setParseError(null);
    try {
      const result = parseMoleculeFile(name, content);
      setFile({ name, content });
      setParsed(result);
      // auto-suggest merge when single chain detected from multi-chain PDB
      if (result.format === "pdb" && result.chains.length > 1) {
        setForm((f) => ({ ...f, mergeAll: false }));
      }
    } catch (err) {
      setParseError(err.message);
      setFile(null);
      setParsed(null);
    }
  };

  const reset = () => {
    setFile(null);
    setParsed(null);
    setParseError(null);
  };

  const command = file ? buildCommand(file.name, form, parsed) : "";

  const copy = () => {
    navigator.clipboard.writeText(command);
    setCopied(true);
    setTimeout(() => setCopied(false), 2000);
  };

  return (
    <div className="mx-auto max-w-6xl px-6 py-12">
      <div className="mb-10">
        <h1 className="text-4xl font-bold text-white mb-2">Setup Assistant</h1>
        <p className="text-gray-400">
          Upload a PDB or GRO file to analyse the system and generate the{" "}
          <code className="text-gray-200">pdb2gmx</code> command for your multi-eGO topology.
        </p>
      </div>

      <div className="grid gap-8 lg:grid-cols-2">
        {/* ── LEFT column ── */}
        <div className="space-y-6">

          {/* File upload / summary */}
          {!file ? (
            <div>
              <SectionTitle>1. Upload structure file</SectionTitle>
              <DropZone onFile={handleFile} />
              {parseError && (
                <div className="mt-3 rounded-lg border border-red-800 bg-red-950/40 p-3 text-sm text-red-300">
                  ❌ {parseError}
                </div>
              )}
            </div>
          ) : (
            <div>
              <div className="flex items-center justify-between mb-3">
                <SectionTitle>1. Structure file</SectionTitle>
                <button onClick={reset} className="text-xs text-gray-500 hover:text-gray-300 transition">
                  ✕ Remove
                </button>
              </div>
              <div className="card mb-4 flex items-center gap-3">
                <span className="text-2xl">{parsed.format === "pdb" ? "🧬" : "📄"}</span>
                <div>
                  <p className="font-medium text-white text-sm">{file.name}</p>
                  {parsed.title && <p className="text-xs text-gray-500 mt-0.5">{parsed.title}</p>}
                </div>
              </div>
              <MoleculeTable parsed={parsed} />
            </div>
          )}

          {/* pdb2gmx options */}
          <div className={file ? "" : "opacity-40 pointer-events-none"}>
            <SectionTitle>2. pdb2gmx options</SectionTitle>

            <div className="card space-y-5">
              {/* Force field */}
              <div>
                <label className="label">Force field (-ff)</label>
                <p className="text-xs text-gray-500 mb-1">
                  Use <code className="text-gray-400">multi-ego-basic</code> for multi-eGO workflows. The{" "}
                  <code className="text-gray-400">multi-ego-basic.ff</code> directory must be present in your
                  working directory.
                </p>
                <input
                  type="text"
                  value={form.ff}
                  onChange={(e) => set("ff", e.target.value)}
                  className="input-field"
                />
              </div>

              {/* Chain separation */}
              <div>
                <label className="label">Chain separation (-chainsep)</label>
                <select
                  value={form.chainsep}
                  onChange={(e) => set("chainsep", e.target.value)}
                  className="input-field"
                >
                  {CHAINSEP_OPTIONS.map((o) => (
                    <option key={o.value} value={o.value}>{o.label}</option>
                  ))}
                </select>
              </div>

              {/* Output filenames */}
              <div>
                <label className="label">Output files</label>
                <div className="grid grid-cols-3 gap-2">
                  {[
                    { key: "outputGro", placeholder: "conf.gro" },
                    { key: "outputTop", placeholder: "topol.top" },
                    { key: "outputPosre", placeholder: "posre.itp" },
                  ].map(({ key, placeholder }) => (
                    <input
                      key={key}
                      type="text"
                      value={form[key]}
                      placeholder={placeholder}
                      onChange={(e) => set(key, e.target.value)}
                      className="input-field font-mono text-xs"
                    />
                  ))}
                </div>
                <p className="text-xs text-gray-600 mt-1">.gro &nbsp;·&nbsp; .top &nbsp;·&nbsp; posre.itp</p>
              </div>

              {/* Flags */}
              <div className="space-y-3 pt-1">
                <Toggle
                  label="-ignh — ignore hydrogens"
                  description="Remove existing hydrogen atoms and let pdb2gmx re-add them. Recommended for most PDB files."
                  checked={form.ignh}
                  onChange={(v) => set("ignh", v)}
                />
                <Toggle
                  label="-merge all — merge chains"
                  description="Treat all chains as a single molecule. Use for homo-oligomers or when all chains share the same topology."
                  checked={form.mergeAll}
                  onChange={(v) => set("mergeAll", v)}
                />
                <Toggle
                  label="-missing — allow missing atoms"
                  description="Continue even when atoms are missing from residues. Use with caution."
                  checked={form.missing}
                  onChange={(v) => set("missing", v)}
                />
                <Toggle
                  label="-renum — renumber residues"
                  description="Renumber residues sequentially starting from 1."
                  checked={form.renum}
                  onChange={(v) => set("renum", v)}
                />
              </div>
            </div>
          </div>
        </div>

        {/* ── RIGHT column: command output ── */}
        <div className="lg:sticky lg:top-20 lg:self-start space-y-6">
          <div>
            <div className="flex items-center justify-between mb-3">
              <h2 className="text-sm font-semibold text-gray-400 uppercase tracking-wider">
                Generated command
              </h2>
              {file && (
                <button onClick={copy} className="btn-secondary text-xs py-1.5 px-3">
                  {copied ? "✓ Copied!" : "Copy"}
                </button>
              )}
            </div>
            <pre className={`code-block text-xs leading-relaxed min-h-32 ${!file ? "opacity-30" : ""}`}>
              {file
                ? command
                : "# Upload a structure file to generate the command"}
            </pre>
          </div>

          {/* Next steps */}
          <div className="card space-y-4">
            <h3 className="font-semibold text-white">Next steps</h3>
            <ol className="space-y-3 text-sm text-gray-400 list-none">
              {[
                {
                  n: "1",
                  text: "Run the command above. pdb2gmx will interactively ask you to choose terminal groups for each chain.",
                },
                {
                  n: "2",
                  text: (
                    <>
                      Verify the output <code className="text-gray-300">topol.top</code> and{" "}
                      <code className="text-gray-300">conf.gro</code> look correct. Edit <code 
                      className="text-gray-300">topol.top</code> [moleculetype] name 
                      and correspondig [molecules] name to reflect your system.
                    </>
                  ),
                },
                {
                  n: "3",
                  text: "Generate training data by MD simulation, AI, single structure or reweighin and extract contact matrices using cmdata.",
                },
                {
                  n: "4",
                  text: (
                    <>
                      Use the{" "}
                      <a href="../config" className="text-brand-400 hover:underline">
                        Config Builder
                      </a>{" "}
                      to set up the prior and then the production runs.
                    </>
                  ),
                },
              ].map(({ n, text }) => (
                <li key={n} className="flex gap-3">
                  <span className="shrink-0 flex h-5 w-5 items-center justify-center rounded-full bg-brand-900 text-brand-400 text-xs font-bold">
                    {n}
                  </span>
                  <span>{text}</span>
                </li>
              ))}
            </ol>
          </div>

          {/* pdb2gmx docs link */}
          <p className="text-xs text-gray-600">
            Reference:{" "}
            <a
              href="https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html"
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-500 hover:underline"
            >
              gmx pdb2gmx documentation ↗
            </a>
          </p>
        </div>
      </div>
    </div>
  );
}
