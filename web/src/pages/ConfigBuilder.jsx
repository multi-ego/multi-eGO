import { useState, useCallback } from "react";
import yaml from "js-yaml";

// ─── default symmetries present in most example configs ──────────────────────
const DEFAULT_SYMMETRIES = [
  // Side-chain chemical equivalences
  "ARG NH1 NH2",
  "ASP OD1 OD2",
  "GLU OE1 OE2",
  "PHE CD1 CD2",
  "PHE CE1 CE2",
  "TYR CD1 CD2",
  "TYR CE1 CE2",
  // Backbone O1/O2 terminal equivalences for all standard amino acids
  "ALA O1 O2",
  "ARG O1 O2",
  "ASN O1 O2",
  "ASP O1 O2",
  "CYS O1 O2",
  "GLN O1 O2",
  "GLU O1 O2",
  "GLY O1 O2",
  "HIS O1 O2",
  "ILE O1 O2",
  "LEU O1 O2",
  "LYS O1 O2",
  "MET O1 O2",
  "PHE O1 O2",
  "PRO O1 O2",
  "SER O1 O2",
  "THR O1 O2",
  "TRP O1 O2",
  "TYR O1 O2",
  "VAL O1 O2",
];

const MATRIX_OPTIONS = ["intramat_1_1", "intermat_1_1", "intramat_1_2", "intermat_1_2"];

const emptyRef = () => ({
  id: crypto.randomUUID(),
  reference: "",
  train: "",
  matrix: "intramat_1_1",
  epsilon: 0.25,
});

function buildYaml(form) {
  const doc = [];

  doc.push({ system: form.system });
  doc.push({ egos: form.egos });

  if (form.p_to_learn !== 0.9995) doc.push({ p_to_learn: form.p_to_learn });
  if (form.epsilon_min !== 0.07) doc.push({ epsilon_min: form.epsilon_min });
  if (form.learn_tolerance !== 0.01) doc.push({ learn_tolerance: form.learn_tolerance });
  if (form.no_header) doc.push("no_header");
  if (form.force_split) doc.push("force_split");
  if (form.single_molecule) doc.push("single_molecule");
  if (form.explicit_name) doc.push({ explicit_name: form.explicit_name });

  const syms = form.symmetry
    .split("\n")
    .map((s) => s.trim())
    .filter(Boolean);
  if (syms.length) doc.push({ symmetry: syms });

  if (form.egos === "production" && form.inputRefs.length) {
    const refs = form.inputRefs.map(({ reference, train, matrix, epsilon }) => ({
      reference,
      train,
      matrix,
      epsilon,
    }));
    doc.push({ input_refs: refs });
  }

  return "---\n" + yaml.dump(doc, { lineWidth: -1, quotingType: '"' });
}

// ─── sub-components ───────────────────────────────────────────────────────────

function SectionTitle({ children }) {
  return <h2 className="text-lg font-semibold text-white border-b border-gray-800 pb-2 mb-4">{children}</h2>;
}

function Collapsible({ title, children, defaultOpen = false }) {
  const [open, setOpen] = useState(defaultOpen);
  return (
    <div className="card">
      <button
        type="button"
        onClick={() => setOpen((o) => !o)}
        className="flex w-full items-center justify-between text-left"
      >
        <h2 className="text-lg font-semibold text-white">{title}</h2>
        <svg
          className={`h-5 w-5 text-gray-400 transition-transform ${open ? "rotate-180" : ""}`}
          fill="none"
          viewBox="0 0 24 24"
          stroke="currentColor"
          strokeWidth={2}
        >
          <path strokeLinecap="round" strokeLinejoin="round" d="M19 9l-7 7-7-7" />
        </svg>
      </button>
      {open && <div className="mt-4 border-t border-gray-800 pt-4">{children}</div>}
    </div>
  );
}

function Toggle({ label, description, checked, onChange }) {
  return (
    <label className="flex items-start gap-3 cursor-pointer group">
      <div className="mt-0.5">
        <input type="checkbox" className="sr-only" checked={checked} onChange={(e) => onChange(e.target.checked)} />
        <div
          className={`w-10 h-6 rounded-full transition ${checked ? "bg-brand-600" : "bg-gray-700"} relative`}
          onClick={() => onChange(!checked)}
        >
          <div
            className={`absolute top-1 h-4 w-4 rounded-full bg-white shadow transition-transform ${
              checked ? "translate-x-5" : "translate-x-1"
            }`}
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

function NumberInput({ label, description, value, onChange, min, max, step }) {
  return (
    <div>
      <label className="label">{label}</label>
      {description && <p className="text-xs text-gray-500 mb-1">{description}</p>}
      <input
        type="number"
        value={value}
        min={min}
        max={max}
        step={step ?? 0.001}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="input-field w-40"
      />
    </div>
  );
}

function InputRefRow({ ref_: r, index, onChange, onRemove }) {
  const update = (key, val) => onChange(index, { ...r, [key]: val });
  return (
    <div className="card relative">
      <button
        onClick={() => onRemove(index)}
        className="absolute right-3 top-3 text-gray-600 hover:text-red-400 transition text-lg leading-none"
        aria-label="Remove"
      >
        ×
      </button>
      <p className="text-xs font-mono text-brand-400 mb-3">input_refs[{index}]</p>
      <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-4">
        <div>
          <label className="label">reference</label>
          <input
            type="text"
            value={r.reference}
            placeholder="reference"
            onChange={(e) => update("reference", e.target.value)}
            className="input-field"
          />
        </div>
        <div>
          <label className="label">train</label>
          <input
            type="text"
            value={r.train}
            placeholder="native_MD"
            onChange={(e) => update("train", e.target.value)}
            className="input-field"
          />
        </div>
        <div>
          <label className="label">matrix</label>
          <select value={r.matrix} onChange={(e) => update("matrix", e.target.value)} className="input-field">
            {MATRIX_OPTIONS.map((m) => (
              <option key={m} value={m}>
                {m}
              </option>
            ))}
          </select>
        </div>
        <div>
          <label className="label">epsilon (kJ/mol)</label>
          <input
            type="number"
            value={r.epsilon}
            min={0.001}
            max={1}
            step={0.001}
            onChange={(e) => update("epsilon", parseFloat(e.target.value))}
            className="input-field"
          />
        </div>
      </div>
    </div>
  );
}

// ─── main page ────────────────────────────────────────────────────────────────

export default function ConfigBuilder() {
  const [form, setForm] = useState({
    system: "",
    egos: "production",
    p_to_learn: 0.9995,
    epsilon_min: 0.07,
    learn_tolerance: 0.01,
    no_header: false,
    force_split: false,
    single_molecule: false,
    explicit_name: "",
    symmetry: DEFAULT_SYMMETRIES.join("\n"),
    inputRefs: [emptyRef()],
  });

  const set = useCallback((key, val) => setForm((f) => ({ ...f, [key]: val })), []);

  const addRef = () => setForm((f) => ({ ...f, inputRefs: [...f.inputRefs, emptyRef()] }));
  const removeRef = (i) =>
    setForm((f) => ({ ...f, inputRefs: f.inputRefs.filter((_, idx) => idx !== i) }));
  const updateRef = (i, updated) =>
    setForm((f) => ({ ...f, inputRefs: f.inputRefs.map((r, idx) => (idx === i ? updated : r)) }));

  const yamlOutput = buildYaml(form);

  const download = () => {
    const blob = new Blob([yamlOutput], { type: "text/yaml" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "config.yml";
    a.click();
    URL.revokeObjectURL(url);
  };

  const copyToClipboard = () => navigator.clipboard.writeText(yamlOutput);

  return (
    <div className="mx-auto max-w-6xl px-6 py-12">
      <div className="mb-10">
        <h1 className="text-4xl font-bold text-white mb-2">Config Builder</h1>
        <p className="text-gray-400">
          Fill in the form to generate a <code className="text-gray-200">config.yml</code> for multi-eGO.
          The YAML preview updates in real time.
        </p>
      </div>

      <div className="grid gap-8 lg:grid-cols-2">
        {/* ── LEFT: form ── */}
        <div className="space-y-8">
          {/* System */}
          <div className="card space-y-4">
            <SectionTitle>System</SectionTitle>
            <div>
              <label className="label">System name</label>
              <p className="text-xs text-gray-500 mb-1">
                Must match the folder name inside <code className="text-gray-400">inputs/</code>
              </p>
              <input
                type="text"
                value={form.system}
                placeholder="e.g. GB1, ttrref, AB42"
                onChange={(e) => set("system", e.target.value)}
                className="input-field"
              />
            </div>
            <div>
              <label className="label">Mode (egos)</label>
              <div className="flex gap-3 mt-1">
                {["mg", "production"].map((opt) => (
                  <button
                    key={opt}
                    onClick={() => set("egos", opt)}
                    className={`rounded-lg border px-4 py-2 text-sm font-medium transition ${
                      form.egos === opt
                        ? "border-brand-500 bg-brand-900 text-brand-300"
                        : "border-gray-700 bg-gray-800 text-gray-400 hover:border-gray-600"
                    }`}
                  >
                    {opt}
                  </button>
                ))}
              </div>
              <p className="text-xs text-gray-500 mt-2">
                {form.egos === "mg"
                  ? "Generates a molten-globule prior force field."
                  : "Generates a production force field from reference + training simulations."}
              </p>
            </div>
          </div>

          {/* input_refs */}
          {form.egos === "production" && (
            <div className="space-y-4">
              <div className="flex items-center justify-between">
                <SectionTitle>input_refs</SectionTitle>
                <button onClick={addRef} className="btn-secondary text-xs py-1.5 px-3">
                  + Add reference
                </button>
              </div>
              {form.inputRefs.map((r, i) => (
                <InputRefRow key={r.id} ref_={r} index={i} onChange={updateRef} onRemove={removeRef} />
              ))}
            </div>
          )}

          {/* Symmetry */}
          <Collapsible title="Symmetry">
            <p className="text-xs text-gray-500 mb-3">
              One entry per line, format: <code className="text-gray-400">RESNAME ATOM1 ATOM2</code>
            </p>
            <textarea
              value={form.symmetry}
              onChange={(e) => set("symmetry", e.target.value)}
              rows={10}
              className="input-field font-mono text-xs resize-y"
              placeholder={"ARG NH1 NH2\nASP OD1 OD2"}
            />
          </Collapsible>

          {/* Flags */}
          <Collapsible title="Flags">
            <div className="space-y-4">
              <Toggle
                label="no_header"
                description="Remove headers from output files"
                checked={form.no_header}
                onChange={(v) => set("no_header", v)}
              />
              <Toggle
                label="force_split"
                description="Split inter and intra-molecular interactions in output files"
                checked={form.force_split}
                onChange={(v) => set("force_split", v)}
              />
              <Toggle
                label="single_molecule"
                description="Enable optimisations for single-molecule simulations"
                checked={form.single_molecule}
                onChange={(v) => set("single_molecule", v)}
              />
            </div>
          </Collapsible>

          {/* Advanced parameters */}
          <Collapsible title="Advanced parameters">
            <div className="space-y-5">
              <NumberInput
                label="p_to_learn"
                description="Fraction of training contacts to learn (default: 0.9995)"
                value={form.p_to_learn}
                onChange={(v) => set("p_to_learn", v)}
                min={0.9}
                max={1}
                step={0.0001}
              />
              <NumberInput
                label="epsilon_min (kJ/mol)"
                description="Minimum meaningful interaction energy (default: 0.07)"
                value={form.epsilon_min}
                onChange={(v) => set("epsilon_min", v)}
                min={0.001}
                max={1}
              />
              <NumberInput
                label="learn_tolerance"
                description="Relative deviation from default to set new c6/c12 (default: 0.01)"
                value={form.learn_tolerance}
                onChange={(v) => set("learn_tolerance", v)}
                min={0}
                max={1}
              />
              <div>
                <label className="label">explicit_name</label>
                <p className="text-xs text-gray-500 mb-1">Optional name for the output directory</p>
                <input
                  type="text"
                  value={form.explicit_name}
                  placeholder="leave blank for auto"
                  onChange={(e) => set("explicit_name", e.target.value)}
                  className="input-field"
                />
              </div>
            </div>
          </Collapsible>
        </div>

        {/* ── RIGHT: YAML preview ── */}
        <div className="lg:sticky lg:top-20 lg:self-start space-y-3">
          <div className="flex items-center justify-between">
            <h2 className="text-sm font-semibold text-gray-400 uppercase tracking-wider">config.yml preview</h2>
            <div className="flex gap-2">
              <button onClick={copyToClipboard} className="btn-secondary text-xs py-1.5 px-3">
                Copy
              </button>
              <button onClick={download} className="btn-primary text-xs py-1.5 px-3">
                Download
              </button>
            </div>
          </div>
          <pre className="code-block min-h-64 text-xs leading-relaxed whitespace-pre overflow-auto max-h-[80vh]">
            {yamlOutput}
          </pre>
          <p className="text-xs text-gray-600">
            Save as <code>inputs/{form.system || "&lt;system&gt;"}/config.yml</code> and run:{" "}
            <code className="text-gray-400">python multiego.py --config inputs/{form.system || "&lt;system&gt;"}/config.yml</code>
          </p>
        </div>
      </div>
    </div>
  );
}
