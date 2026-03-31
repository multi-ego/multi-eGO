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
            <Link to="/install" className="btn-secondary">
              Installation →
            </Link>
            <Link to="/setup" className="btn-secondary">
              Setup assistant →
            </Link>
            <Link to="/config" className="btn-secondary">
              Config builder →
            </Link>
            <Link to="/examples" className="btn-secondary">
              Examples →
            </Link>
            <Link to="/simulation" className="btn-secondary">
              Run simulation →
            </Link>
          </div>
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
