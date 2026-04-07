import { Link } from "react-router-dom";

const COLAB_URL =
  "https://colab.research.google.com/github/multi-ego/multi-eGO/blob/main/tools/colab/run_mg_gromacs.ipynb";

const OPTIONS = [
  {
    label: "Option A",
    badge: "~2 min",
    title: "Install via apt",
    description:
      "Installs the Ubuntu-packaged GROMACS binary with a single apt command. " +
      "Covers all steps of the mg workflow and is the fastest way to get started.",
    pros: ["Ready in ~2 minutes", "No compilation required", "Full mg workflow supported"],
    cons: ["Older GROMACS version", "No multi-eGO-specific patches"],
  },
  {
    label: "Option B",
    badge: "~20 min",
    title: "Compile from source",
    description:
      "Builds the multi-eGO GROMACS fork from the release-2023 branch. " +
      "Compilation takes ~20 minutes on Colab CPUs but produces the exact binary " +
      "used during multi-eGO development, with CUDA GPU support when a GPU runtime is selected.",
    pros: ["multi-eGO GROMACS fork", "CUDA GPU support", "Protocol fidelity"],
    cons: ["~20 min compile time"],
  },
];

const STEPS = [
  {
    number: "01",
    title: "Install GROMACS",
    description:
      "Choose between a quick apt install (~2 min) or compiling the multi-eGO GROMACS fork " +
      "from source (~20 min). Only one option needs to be run.",
  },
  {
    number: "02",
    title: "Install multi-eGO",
    description:
      "The multi-eGO package is cloned from GitHub and installed with pip. " +
      "GMXLIB is configured automatically so GROMACS finds the multi-eGO force-field files.",
  },
  {
    number: "03",
    title: "Upload PDB file",
    description:
      "Upload the PDB file of your protein. The notebook derives the system name from the filename.",
  },
  {
    number: "04",
    title: "Generate topology (pdb2gmx)",
    description:
      "gmx pdb2gmx builds a GROMACS topology from the input structure using the multi-eGO-basic force field.",
  },
  {
    number: "05",
    title: "Generate mg force field (mego)",
    description:
      "mego --egos mg reads the AMBER topology and writes topol_mego.top and ffnonbonded.itp " +
      "with the molten-globule C6/C12 non-bonded parameters.",
  },
  {
    number: "06",
    title: "Energy minimisation",
    description:
      "Steepest-descent minimisation with the mg force field removes bad contacts " +
      "in the starting structure before dynamics begin.",
  },
  {
    number: "07",
    title: "NVT production run",
    description:
      "Langevin dynamics at the chosen temperature and timestep. GPU is used automatically " +
      "when a Colab GPU runtime is active. Produces run.xtc and run.tpr.",
  },
  {
    number: "08",
    title: "Energy analysis",
    description:
      "Potential energy and temperature are extracted from the GROMACS energy file " +
      "and plotted so you can verify the simulation is well-behaved.",
  },
  {
    number: "09",
    title: "Download outputs",
    description:
      "run.xtc and run.tpr are downloaded to your computer — everything needed " +
      "to extract contact histograms with cmdata.",
  },
];

const PREREQS = [
  {
    file: "protein.pdb",
    note: "PDB coordinates of the protein (no water, no ions)",
  },
];

export default function Simulation() {
  return (
    <div className="mx-auto max-w-6xl space-y-16 px-6 py-12">
      {/* Header */}
      <div className="space-y-4">
        <h1 className="section-heading">mg Reference Simulation</h1>
        <p className="max-w-2xl text-gray-400">
          Run the complete molten-globule reference simulation pipeline directly in your
          browser — no local GROMACS installation required. The notebook runs on Google Colab
          and covers every setup step from a raw PDB file to a production trajectory.
        </p>
      </div>

      {/* Colab CTA */}
      <div className="rounded-xl border border-brand-800 bg-brand-950/30 p-8 text-center">
        <p className="mb-2 text-sm font-medium uppercase tracking-widest text-brand-400">
          Run in your browser
        </p>
        <h2 className="mb-6 text-2xl font-bold text-white">
          Multi-<em>e</em>GO mg simulation · GROMACS notebook
        </h2>
        <a
          href={COLAB_URL}
          target="_blank"
          rel="noopener noreferrer"
          className="inline-flex items-center gap-3 rounded-lg bg-[#F9AB00] px-6 py-3 text-sm font-semibold text-gray-900 transition hover:bg-[#F9AB00]/90"
        >
          <img
            src="https://colab.research.google.com/assets/colab-badge.svg"
            alt="Open in Colab"
            className="h-5"
          />
          Open in Google Colab
        </a>
        <p className="mt-4 text-xs text-gray-500">
          Free GPU available · requires a Google account · no local installation needed
        </p>
      </div>

      {/* Two installation options */}
      <div>
        <h2 className="section-heading mb-2">GROMACS installation options</h2>
        <p className="mb-6 text-gray-400">
          The notebook offers two ways to install GROMACS — run exactly one of them.
        </p>
        <div className="grid gap-6 lg:grid-cols-2">
          {OPTIONS.map((opt) => (
            <div key={opt.label} className="card space-y-3">
              <div className="flex items-center gap-3">
                <span className="rounded bg-brand-900 px-2 py-0.5 text-xs font-bold text-brand-300">
                  {opt.label}
                </span>
                <span className="rounded bg-gray-800 px-2 py-0.5 text-xs text-gray-400">
                  {opt.badge}
                </span>
                <h3 className="font-semibold text-white">{opt.title}</h3>
              </div>
              <p className="text-sm text-gray-400">{opt.description}</p>
              <ul className="space-y-1 text-sm">
                {opt.pros.map((p) => (
                  <li key={p} className="flex gap-2 text-green-400">
                    <span>✓</span>
                    <span>{p}</span>
                  </li>
                ))}
                {opt.cons.map((c) => (
                  <li key={c} className="flex gap-2 text-gray-500">
                    <span>–</span>
                    <span>{c}</span>
                  </li>
                ))}
              </ul>
            </div>
          ))}
        </div>
      </div>

      {/* Prerequisites */}
      <div>
        <h2 className="section-heading mb-6">Prerequisites</h2>
        <p className="mb-4 text-gray-400">
          All you need is a PDB file of the protein you want to simulate. The notebook handles
          topology generation, force-field parameterisation, and simulation setup automatically.
        </p>
        <div className="overflow-x-auto rounded-xl border border-gray-800">
          <table className="w-full text-sm">
            <thead>
              <tr className="border-b border-gray-800 bg-gray-900/60">
                <th className="px-4 py-3 text-left font-medium text-gray-400">File</th>
                <th className="px-4 py-3 text-left font-medium text-gray-400">Description</th>
              </tr>
            </thead>
            <tbody>
              {PREREQS.map((p) => (
                <tr key={p.file} className="border-b border-gray-800/50 last:border-0">
                  <td className="px-4 py-3">
                    <code className="rounded bg-gray-800 px-1.5 py-0.5 text-xs text-brand-300">
                      {p.file}
                    </code>
                  </td>
                  <td className="px-4 py-3 text-gray-400">{p.note}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <div className="mt-4 flex gap-3">
          <Link to="/setup" className="btn-secondary text-sm">
            Setup assistant →
          </Link>
          <Link to="/config" className="btn-secondary text-sm">
            Config builder →
          </Link>
        </div>
      </div>

      {/* What the notebook does */}
      <div>
        <h2 className="section-heading mb-6">What the notebook does</h2>
        <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {STEPS.map((step) => (
            <div key={step.number} className="card">
              <span className="mb-3 block font-mono text-3xl font-bold text-brand-800">
                {step.number}
              </span>
              <h3 className="mb-2 font-semibold text-white">{step.title}</h3>
              <p className="text-sm text-gray-400">{step.description}</p>
            </div>
          ))}
        </div>
      </div>

      {/* Output & next steps */}
      <div className="grid gap-6 lg:grid-cols-2">
        <div className="card">
          <h3 className="mb-3 font-semibold text-white">Output files</h3>
          <ul className="space-y-2 text-sm text-gray-400">
            {[
              ["run.xtc", "Production trajectory in GROMACS XTC format"],
              ["run.tpr", "GROMACS run-input file (required by cmdata)"],
              ["run.edr", "Binary energy file"],
              ["energy_plot.png", "Potential energy and temperature plots"],
              ["topol_mego.top", "GROMACS topology with mg force field"],
              ["ffnonbonded.itp", "Non-bonded C6/C12 parameters"],
            ].map(([f, desc]) => (
              <li key={f} className="flex gap-3">
                <code className="shrink-0 rounded bg-gray-800 px-1.5 py-0.5 text-xs text-brand-300">
                  {f}
                </code>
                <span>{desc}</span>
              </li>
            ))}
          </ul>
        </div>

        <div className="card">
          <h3 className="mb-3 font-semibold text-white">Next steps</h3>
          <ol className="space-y-3 text-sm text-gray-400">
            <li className="flex gap-3">
              <span className="font-mono text-brand-400">1.</span>
              <span>
                Extract contact histograms from the trajectory with{" "}
                <code className="rounded bg-gray-800 px-1 text-xs text-brand-300">cmdata</code>
                {" "}(requires both{" "}
                <code className="rounded bg-gray-800 px-1 text-xs text-brand-300">run.xtc</code>
                {" "}and{" "}
                <code className="rounded bg-gray-800 px-1 text-xs text-brand-300">run.tpr</code>).
              </span>
            </li>
            <li className="flex gap-3">
              <span className="font-mono text-brand-400">2.</span>
              <span>
                Convert histograms to contact matrices with{" "}
                <code className="rounded bg-gray-800 px-1 text-xs text-brand-300">make_mat.py</code>.
              </span>
            </li>
            <li className="flex gap-3">
              <span className="font-mono text-brand-400">3.</span>
              <span>
                Generate the production force field with{" "}
                <code className="rounded bg-gray-800 px-1 text-xs text-brand-300">
                  mego --config config.yml
                </code>
                .
              </span>
            </li>
          </ol>
        </div>
      </div>
    </div>
  );
}
