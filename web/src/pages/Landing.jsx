import { Link } from "react-router-dom";

const steps = [
  {
    number: "01",
    title: "Prepare the topology",
    description:
      "Build a GROMACS topology using the multi-ego-basic.ff force field and run a training simulation with your all-atom force field.",
  },
  {
    number: "02",
    title: "Generate the mg prior",
    description:
      "Create a molten-globule (mg) baseline force field from physico-chemical and statistical properties. Run a short reference simulation.",
  },
  {
    number: "03",
    title: "Learn and produce",
    description:
      "Extract contact probabilities from the reference simulation, then run multi-eGO to learn pairwise Lennard-Jones interactions and produce the production force field.",
  },
];

const features = [
  {
    icon: "⚛️",
    title: "Atomic resolution",
    description:
      "Full all-atom Lennard-Jones parametrisation — no coarse-graining, full structural detail retained.",
  },
  {
    icon: "📊",
    title: "Data-driven",
    description:
      "Learns directly from contact probability distributions extracted from MD trajectories, structures, or AI ensembles.",
  },
  {
    icon: "🔁",
    title: "GROMACS compatible",
    description:
      "Output is standard ffnonbonded.itp and topol_mego.top — drop-in for any GROMACS workflow.",
  },
  {
    icon: "⚡",
    title: "Fast simulations",
    description:
      "Orders of magnitude faster than all-atom force fields while reproducing folding, aggregation and binding thermodynamics.",
  },
  {
    icon: "🧩",
    title: "Multi-reference",
    description:
      "Combine multiple training simulations and reference states in a single config file for complex systems.",
  },
  {
    icon: "🔬",
    title: "Validated",
    description:
      "Benchmarked on protein folding, amyloid aggregation, ligand binding and intrinsically disordered proteins.",
  },
];

const installSteps = [
  { label: "Conda", code: "conda env create -f conda/environment.yml\nconda activate meGO" },
  { label: "pip", code: "pip install -r requirements.txt" },
];

export default function Landing() {
  return (
    <div className="space-y-24 pb-24">
      {/* Hero */}
      <section className="relative overflow-hidden border-b border-gray-800 bg-gradient-to-b from-gray-900 to-gray-950">
        <div className="mx-auto max-w-6xl px-6 py-24 text-center">
          <div className="mb-6 inline-flex items-center rounded-full border border-brand-800 bg-brand-950/50 px-4 py-1 text-xs font-medium text-brand-300">
            v1.0 — GPL v3 Open Source
          </div>
          <h1 className="mb-6 text-5xl font-bold tracking-tight text-white sm:text-6xl">
            Multi-<em>e</em>GO
          </h1>
          <p className="mx-auto mb-10 max-w-2xl text-lg text-gray-400">
            Data-driven, atomic-resolution force fields for molecular dynamics simulations. Learn pairwise
            Lennard-Jones interactions directly from contact probability distributions.
          </p>
          <div className="flex flex-wrap justify-center gap-4">
            <a
              href="https://github.com/multi-ego/multi-eGO"
              target="_blank"
              rel="noopener noreferrer"
              className="btn-primary"
            >
              <svg className="h-4 w-4" fill="currentColor" viewBox="0 0 24 24">
                <path
                  fillRule="evenodd"
                  clipRule="evenodd"
                  d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z"
                />
              </svg>
              View on GitHub
            </a>
            <Link to="/config" className="btn-secondary">
              Build a config file →
            </Link>
          </div>
          {/* Badges */}
          <div className="mt-10 flex flex-wrap justify-center gap-3">
            <img src="https://img.shields.io/badge/Version-1.0-blue" alt="Version" />
            <img src="https://img.shields.io/badge/License-GPL%20v3-blue.svg" alt="License" />
            <img src="https://img.shields.io/badge/code%20style-black-000000.svg" alt="Code style: black" />
          </div>
        </div>
      </section>

      {/* How it works */}
      <section className="mx-auto max-w-6xl px-6">
        <h2 className="section-heading mb-4 text-center">How it works</h2>
        <p className="mx-auto mb-12 max-w-2xl text-center text-gray-400">
          Multi-<em>e</em>GO constructs an atomistic Lennard-Jones force field following a Bayesian
          learning framework in three stages.
        </p>
        <div className="grid gap-6 sm:grid-cols-3">
          {steps.map((step) => (
            <div key={step.number} className="card relative">
              <span className="mb-4 block font-mono text-4xl font-bold text-brand-800">{step.number}</span>
              <h3 className="mb-2 text-lg font-semibold text-white">{step.title}</h3>
              <p className="text-sm text-gray-400">{step.description}</p>
            </div>
          ))}
        </div>
      </section>

      {/* Features */}
      <section className="mx-auto max-w-6xl px-6">
        <h2 className="section-heading mb-12 text-center">Features</h2>
        <div className="grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
          {features.map((f) => (
            <div key={f.title} className="card">
              <span className="mb-3 block text-3xl">{f.icon}</span>
              <h3 className="mb-2 font-semibold text-white">{f.title}</h3>
              <p className="text-sm text-gray-400">{f.description}</p>
            </div>
          ))}
        </div>
      </section>

      {/* Installation */}
      <section className="mx-auto max-w-6xl px-6">
        <h2 className="section-heading mb-4">Installation</h2>
        <p className="mb-8 text-gray-400">
          Requires Python ≥ 3.10. Choose your preferred environment manager:
        </p>
        <div className="grid gap-6 sm:grid-cols-2">
          {installSteps.map((s) => (
            <div key={s.label}>
              <p className="label mb-2">{s.label}</p>
              <pre className="code-block whitespace-pre">{s.code}</pre>
            </div>
          ))}
        </div>
        <p className="mt-6 text-sm text-gray-500">
          The <code className="text-gray-300">cmdata</code> trajectory analysis tool requires a separate
          installation — see the{" "}
          <a
            href="https://github.com/multi-ego/multi-eGO/blob/main/tools/cmdata/README.md"
            className="text-brand-400 hover:underline"
            target="_blank"
            rel="noopener noreferrer"
          >
            cmdata README
          </a>{" "}
          for instructions.
        </p>
      </section>

      {/* Quick start */}
      <section className="mx-auto max-w-6xl px-6">
        <h2 className="section-heading mb-4">Quick start</h2>
        <p className="mb-6 text-gray-400">
          Generate a molten-globule prior, then a production force field:
        </p>
        <div className="space-y-4">
          <div>
            <p className="label">1. Generate the mg prior</p>
            <pre className="code-block">python multiego.py --system GB1 --egos mg</pre>
          </div>
          <div>
            <p className="label">2. Generate the production force field</p>
            <pre className="code-block">
              python multiego.py --system GB1 --egos production --train md_monomer --epsilon 0.3
            </pre>
          </div>
          <div>
            <p className="label">3. Or use a config file for multiple training sets</p>
            <pre className="code-block">python multiego.py --config inputs/my_system/config.yml</pre>
          </div>
        </div>
        <div className="mt-6">
          <Link to="/config" className="btn-primary">
            Build your config file →
          </Link>
        </div>
      </section>
    </div>
  );
}
