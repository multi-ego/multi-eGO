import { Link } from "react-router-dom";

const features = [
  {
    icon: "⚛️",
    title: "Atomic resolution",
    description:
      "No coarse-graining, full structural detail retained.",
  },
  {
    icon: "📊",
    title: "Data-driven & Bayesian",
    description:
      "Reweight prior information by learning contact probability distributions from any combination of MD trajectories, structures, or AI ensembles.",
  },
  {
    icon: "🧩",
    title: "Multi-reference",
    description:
      "Combine multiple training simulations and reference states to describe complex systems and processes.",
  },
  {
    icon: "🔁",
    title: "GROMACS compatible",
    description:
      "Output is ready for any GROMACS workflow.",
  },
  {
    icon: "⚡",
    title: "Fast simulations and hypothesis testing",
    description:
      "Orders of magnitude faster than all-atom force fields allow running in-silico experiments.",
  },
  {
    icon: "🔬",
    title: "Validated",
    description:
      "Benchmarked on protein folding, amyloid aggregation, ligand binding and intrinsically disordered proteins.",
  },
];

export default function Landing() {
  return (
    <div className="space-y-24 pb-24">
      {/* Hero */}
      <section className="relative overflow-hidden border-b border-gray-800 bg-gradient-to-b from-gray-900 to-gray-950">
        <div className="mx-auto max-w-6xl px-6 py-24 text-center">
          <div className="mb-6 inline-flex items-center rounded-full border border-brand-800 bg-brand-950/50 px-4 py-1 text-xs font-medium text-brand-300">
            v beta.6 — GPL v3 Open Source
          </div>
          <h1 className="mb-6 text-5xl font-bold tracking-tight text-white sm:text-6xl">
            Multi-<em>e</em>GO
          </h1>
          <p className="mx-auto mb-10 max-w-2xl text-lg text-gray-400">
            Data-driven, atomic-resolution force fields for molecular dynamics simulations.
          </p>
        </div>
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

      {/* Quick start */}
      <section className="mx-auto max-w-6xl px-6">
        <h2 className="section-heading mb-4">Quick start</h2>
        <div className="mt-6 flex flex-wrap gap-3">
          <Link to="/examples" className="btn-secondary">
            Examples →
          </Link>
          <Link to="/simulation" className="btn-secondary">
            Run simulation →
          </Link>
        </div>
      </section>
    </div>
  );
}
