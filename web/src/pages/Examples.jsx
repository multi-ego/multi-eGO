import { useState } from "react";
import { Link } from "react-router-dom";

// ─── example data ─────────────────────────────────────────────────────────────

const EXAMPLES = [
  {
    id: "gpref",
    title: "G-protein B1 domain",
    subtitle: "Our benchmark for protein folding and dynamics",
    tags: ["1 training set", "1 reference", "singe molecule", "folded", "intramolecular"],
    description:
      "The simplest and most common setup: one training simulation, one reference simulation, " +
      "one molecule in the box. GB1 is a fast-folding protein widely used as a benchmark. " +
      "The training data was generated with CHARMM22* and contact matrices cover intra-molecular contacts only.",
    highlights: [
      { label: "Training FF", value: "CHARMM22*" },
      { label: "Contacts",    value: "intramat_1_1 (intra-molecular)" },
      { label: "ε",           value: "0.350 kJ/mol" },
    ],
    tree: `inputs/gpref/
├── topol.top
├── config.yml
├── md_ensemble/           ← training simulation
│   ├── topol.top
│   ├── intramat_1_1.ndx.h5
│   └── charmm22st.ff/
└── reference/             ← mg reference simulation
    ├── topol.top
    ├── intramat_1_1.ndx.h5
    ├── ffnonbonded.itp
    └── topol_mego.top`,
    config: `---
- system: gpref
- egos: production
- symmetry:
  - ARG NH1 NH2
  - ASP OD1 OD2
  - GLU OE1 OE2
  - PHE CD1 CD2
  - PHE CE1 CE2
  - TYR CD1 CD2
  - TYR CE1 CE2
  - CTER O1 O2
- input_refs:
  - reference: reference
    train: md_ensemble
    matrix: intramat_1_1
    epsilon: 0.350`,
    command: "mego --config inputs/gpref/config.yml",
  },
  {
    id: "abetaref",
    title: "Aβ42",
    subtitle: "Our benchmark for intrinsically disordered peptides dynamics",
    tags: ["single molecule", "1 training set", "1 reference", "IDP", "intramolecular"],
    description:
      "Amyloid-β is an intrinsically disordered peptide (IDP) prone to aggregation. " +
      "This example demonstrates that multi-eGO is not limited to folded proteins: " +
      "the mg prior and learning procedure work equally well for disordered systems. " +
      "The training data was produced with CHARMM22* and the folder uses a non-standard " +
      "reference directory name (ref instead of reference), illustrating that the name " +
      "only needs to match the config file.",
    highlights: [
      { label: "Training FF", value: "CHARMM22*" },
      { label: "Contacts",    value: "intramat_1_1 (intra-molecular)" },
      { label: "ε",           value: "0.325 kJ/mol" },
    ],
    tree: `inputs/abetaref/
├── topol.top
├── config.yml
├── native_MD/             ← training simulation
│   ├── topol.top
│   ├── intramat_1_1.ndx.h5
│   └── charmm27.ff/
└── ref/                   ← mg reference simulation (custom name)
    ├── topol.top
    ├── intramat_1_1.ndx.h5
    ├── ffnonbonded.itp
    └── topol_mego.top`,
    config: `---
- system: abetaref
- egos: production
- symmetry:
  - ARG NH1 NH2
  - ASP OD1 OD2
  - GLU OE1 OE2
  - PHE CD1 CD2
  - PHE CE1 CE2
  - TYR CD1 CD2
  - TYR CE1 CE2
  - CTER O1 O2
- input_refs:
  - reference: ref
    train: native_MD
    matrix: intramat_1_1
    epsilon: 0.325`,
    command: "mego --config inputs/abetaref/config.yml",
  },
  {
    id: "ttrref",
    title: "TTR peptide",
    subtitle: "Our benchmark for peptides self-assembly",
    tags: ["2 training sets", "2 references", "multi-state", "intra+inter molecular"],
    description:
      "A more advanced setup: a disordered peptide trained from two  simulations — mononomer and fibril state — " +
      "using two different force fields. Both intra-molecular and inter-molecular contacts " +
      "are learned, with intra-molecular resulting from the combination of the two different states. " +
      "Two reference state, the molten globule for intramolecular interactions and a zero prior for intermolecular ones.",
    highlights: [
      { label: "Training FFs", value: "AMBER99SB-disp (native_MD), CHARMM22* (fibril)" },
      { label: "Contacts",     value: "intramat_1_1 + intermat_1_1 (intra & inter-molecular)" },
      { label: "ε",            value: "0.25 kJ/mol (all refs)" },
    ],
    tree: `inputs/ttrref/
├── topol.top
├── config.yml
├── native_MD/             ← training set 1 (native simulation)
│   ├── topol.top
│   ├── intramat_1_1.ndx.h5
│   └── amber99SBdisp.ff/
├── fibril/                ← training set 2 (fibril simulation)
│   ├── topol.top
│   ├── intramat_1_1.ndx.h5
│   ├── intermat_1_1.ndx.h5
│   └── charmm22st.ff/
└── reference/             ← reference simulation
    ├── topol.top
    ├── intramat_1_1.ndx.h5 ← mg reference simulation
    ├── intermat_1_1.ndx.h5 ← zero reference
    ├── ffnonbonded.itp
    └── topol_mego.top`,
    config: `---
- system: ttrref
- egos: production
- symmetry:
  - ARG NH1 NH2
  - ASP OD1 OD2
  - GLU OE1 OE2
  - PHE CD1 CD2
  - PHE CE1 CE2
  - TYR CD1 CD2
  - TYR CE1 CE2
  - CTER O1 O2
- input_refs:
  - reference: reference
    train: native_MD
    matrix: intramat_1_1
    epsilon: 0.25
  - reference: reference
    train: fibril
    matrix: intramat_1_1
    epsilon: 0.25
  - reference: reference
    train: fibril
    matrix: intermat_1_1
    epsilon: 0.25`,
    command: "mego --config inputs/ttrref/config.yml",
  },
  {
    id: "lyso-bnz",
    title: "Lysozyme + Benzene",
    subtitle: "Our benchmark for multi-domains protein and small molecule binding",
    tags: ["2 molecules", "1 training set", "3 references", "multi-domains", "protein-ligand", "intra+inter molecular"],
    description:
      "A protein–ligand complex: hen-egg-white lysozyme (LYZ) with benzene (BNZ) as a probe ligand. " +
      "This example introduces two concepts not present in the single-molecule examples. " +
      "First, Lysozyme is treated as a two domains protein, so the same single training is learned in two steps: " +
      "Each domain is learned on top of the standard molten globule prior, while iterdomain contacts are learned " +
      "on a prior resulting from the meGO production simulation obtained after the first step (i.e. a meGO simulation " +
      "where the domains are folded but inter-domain contacts are defined as molten globule. " +
      "Third, in this case we also introduce inter-molecular contacts between the protein and the ligand. " +
      "These are learned from a training simulation with a BNZ in the pocket and 5 BNZ in solution. " +
      "The intermolecular prior is obtained by a meGO simulation with the same box of the training and 5 BNZ molecules " +
      "described with the correct internal geometry and molten globule intermolecular interactions. " +
      "The benzene ring symmetry (all six ring atoms interchangeable) is declared explicitly in the config.",
    highlights: [
      { label: "Training FF",  value: "DES-Amber" },
      { label: "Contacts",     value: "intramat_1_1 + intermat_1_2 (intra & inter-molecular)" },
      { label: "ε (LYZ intra)", value: "0.28 kJ/mol" },
      { label: "ε (LYZ-BNZ inter)", value: "0.53 kJ/mol" },
    ],
    tree: `inputs/lyso-bnz_ref/
├── topol.top
├── topol_BNZ.itp
├── config.yml
├── mg/                      ← MG prior for LYZ intra + same training concentration prior for LYZ-BNZ inter
│   ├── intramat_1_1.ndx.h5
│   ├── intermat_1_2.ndx.h5
│   ├── ffnonbonded.itp
│   ├── top_BNZ.itp
│   └── topol_mego.top
├── mg_id/                   ← multi-domain prior for LYZ
│   ├── intramat_1_1.ndx.h5
│   ├── ffnonbonded.itp
│   ├── top_BNZ.itp
│   └── topol_mego.top
└── training/                ← training simulation (DES-Amber FF)
    ├── topol.top
    ├── topol_BNZ.itp
    ├── topol_Lyso.itp
    ├── intramat_1_1.ndx.h5
    ├── intermat_1_2.ndx.h5
    └── des-amber.ff/`,
    config: `---
- system: lyso-bnz_ref
- egos: production
- no_header
- symmetry:
  - ASP OD1 OD2
  - GLU OE1 OE2
  - PHE CD1 CD2
  - PHE CE1 CE2
  - TYR CD1 CD2
  - TYR CE1 CE2
  - ARG NH1 NH2
  - CTER O1 O2
  - BNZ CD1 CD2 CE1 CE2 CZ CG
- input_refs:
  # LYZ intra-domain — standard MG prior
  - reference: mg
    train: training
    matrix: intramat_1_1
    epsilon: 0.28
  # LYZ inter-domain — multi-domain prior
  - reference: mg_id
    train: training
    matrix: intramat_1_1
    epsilon: 1.0
  # LYZ-BNZ intermolecular — cross-matrix prior
  - reference: mg
    train: training
    matrix: intermat_1_2
    epsilon: 0.53`,
    command: "mego --config inputs/lyso-bnz_ref/config.yml",
  },
];

// ─── sub-components ───────────────────────────────────────────────────────────

function Tag({ children }) {
  return (
    <span className="inline-flex items-center rounded-full border border-brand-800 bg-brand-950/50 px-2.5 py-0.5 text-xs font-medium text-brand-300">
      {children}
    </span>
  );
}

function Highlight({ label, value }) {
  return (
    <div className="rounded-lg border border-gray-800 bg-gray-900/50 px-4 py-3">
      <p className="mb-1 text-xs font-medium uppercase tracking-wide text-gray-500">{label}</p>
      <p className="text-sm text-gray-200">{value}</p>
    </div>
  );
}

function CodeBlock({ children, label }) {
  const [copied, setCopied] = useState(false);

  function copy() {
    navigator.clipboard.writeText(children).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    });
  }

  return (
    <div>
      {label && <p className="label mb-1">{label}</p>}
      <div className="relative">
        <pre className="code-block overflow-x-auto pr-16 text-sm leading-relaxed">{children}</pre>
        <button
          onClick={copy}
          className="absolute right-3 top-3 rounded border border-gray-700 bg-gray-800 px-2 py-1 text-xs text-gray-400 transition hover:border-gray-600 hover:text-white"
        >
          {copied ? "Copied!" : "Copy"}
        </button>
      </div>
    </div>
  );
}

// ─── main page ────────────────────────────────────────────────────────────────

export default function Examples() {
  const [active, setActive] = useState(EXAMPLES[0].id);
  const example = EXAMPLES.find((e) => e.id === active);

  return (
    <div className="mx-auto max-w-6xl space-y-10 px-6 py-12">
      {/* Header */}
      <div>
        <h1 className="section-heading mb-3">Examples</h1>
        <p className="max-w-2xl text-gray-400">
          Four ready-to-run test cases from the repository, covering the most common multi-<em>e</em>GO
          workflows — from a simple single-molecule protein to a multi-chain system learned from
          multiple training simulations simultaneously. Each example lives under{" "}
          <code className="rounded bg-gray-800 px-1.5 py-0.5 text-xs text-brand-300">
            tests/test_inputs/
          </code>{" "}
          and can be run directly with{" "}
          <code className="rounded bg-gray-800 px-1.5 py-0.5 text-xs text-brand-300">
            mego --config &lt;path&gt;/config.yml
          </code>
          .
        </p>
      </div>

      {/* Tab selector */}
      <div className="flex gap-2 border-b border-gray-800">
        {EXAMPLES.map((e) => (
          <button
            key={e.id}
            onClick={() => setActive(e.id)}
            className={`-mb-px border-b-2 px-4 py-2.5 text-sm font-medium transition ${
              active === e.id
                ? "border-brand-400 text-brand-400"
                : "border-transparent text-gray-400 hover:text-white"
            }`}
          >
            {e.title}
          </button>
        ))}
      </div>

      {/* Content */}
      <div key={example.id} className="space-y-8">
        {/* Title + tags */}
        <div>
          <div className="mb-2 flex flex-wrap items-center gap-2">
            <h2 className="text-2xl font-bold text-white">{example.title}</h2>
            <span className="text-gray-500">·</span>
            <span className="text-gray-400">{example.subtitle}</span>
          </div>
          <div className="flex flex-wrap gap-2">
            {example.tags.map((t) => (
              <Tag key={t}>{t}</Tag>
            ))}
          </div>
        </div>

        {/* Description */}
        <p className="leading-relaxed text-gray-300">{example.description}</p>

        {/* Highlights grid */}
        <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
          {example.highlights.map((h) => (
            <Highlight key={h.label} label={h.label} value={h.value} />
          ))}
        </div>

        {/* Directory tree + config side by side on wide screens */}
        <div className="grid gap-6 lg:grid-cols-2">
          <CodeBlock label="Input directory structure">{example.tree}</CodeBlock>
          <CodeBlock label="config.yml">{example.config}</CodeBlock>
        </div>

        {/* Run command */}
        <div>
          <p className="label mb-1">Run this example</p>
          <div className="flex flex-wrap items-center gap-4">
            <div className="flex-1">
              <CodeBlock>{example.command}</CodeBlock>
            </div>
            <Link
              to="/config"
              className="btn-secondary whitespace-nowrap"
            >
              Open in Config Builder →
            </Link>
          </div>
        </div>
      </div>

      {/* Progression note */}
      <div className="rounded-xl border border-gray-800 bg-gray-900/40 p-6">
        <h3 className="mb-2 font-semibold text-white">Progression of complexity</h3>
        <p className="text-sm text-gray-400">
          The four examples are ordered by complexity. <strong className="text-gray-200">GB1</strong> is
          the canonical starting point: one molecule, one training simulation, intra-molecular contacts
          only. <strong className="text-gray-200">Aβ</strong> shows the same setup applied to a
          disordered peptide, and illustrates that the reference folder can have any name as long as
          it matches the config. <strong className="text-gray-200">TTR</strong> adds inter-molecular
          contacts and two independent training simulations run with different force fields — a
          more general multi-<em>e</em>GO workflow. <strong className="text-gray-200">Lysozyme + Benzene</strong> introduces
          multi-domain proteins and protein–ligand binding.
        </p>
      </div>
    </div>
  );
}
