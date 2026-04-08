const COLAB_MG_URL =
  "https://colab.research.google.com/github/multi-ego/multi-eGO/blob/main/tools/colab/run_mg_gromacs.ipynb";

const COLAB_EXAMPLES_URL =
  "https://colab.research.google.com/github/multi-ego/multi-eGO/blob/main/tools/colab/run_examples_gromacs.ipynb";

const NOTEBOOKS = [
  {
    url: COLAB_MG_URL,
    tag: "Your own protein",
    title: "mg reference simulation",
    description:
      "Start from any PDB or GRO file and run the complete molten-globule reference " +
      "simulation pipeline — topology, mg force field, energy minimisation, and NVT " +
      "production — ending with a trajectory ready for cmdata.",
    steps: "9 automated steps · upload / test input / FASTA sequence",
  },
  {
    url: COLAB_EXAMPLES_URL,
    tag: "Bundled examples",
    title: "Production simulation of example systems",
    description:
      "Select one of the multi-eGO test systems (GB1, Aβ42, TTR peptide, " +
      "Lysozyme+benzene). All contact matrices are pre-computed — the notebook " +
      "generates the production force field with mego and runs the simulation immediately.",
    steps: "5 automated steps · no file upload needed",
  },
];

const OUTPUT_FILES = [
  ["run.xtc", "Production trajectory in GROMACS XTC format"],
  ["run.tpr", "GROMACS run-input file (required by cmdata)"],
  ["run.edr", "Binary energy file"],
  ["energy_plot.png", "Potential energy and temperature plots"],
  ["topol_mego.top", "GROMACS topology with mg force field"],
  ["ffnonbonded.itp", "Non-bonded C6/C12 parameters"],
];

const OPTIONS = [
  {
    label: "Option A",
    badge: "~2 min",
    title: "Install via apt",
    description:
      "Installs the Ubuntu-packaged GROMACS binary with a single apt command. " +
      "Covers all steps of the mg workflow and is the fastest way to get started.",
    pros: ["Ready in ~2 minutes", "No compilation required", "Full mg workflow supported"],
    cons: ["Not optimized, slow simulations", "Not compatible with CMDATA from multi-eGO"],
  },
  {
    label: "Option B",
    badge: "~20 min",
    title: "Compile from source",
    description:
      "Builds GROMACS release-2023 branch with system FFTW3 (installed via apt). " +
      "Compilation takes ~20 minutes on Colab CPUs and produces an optimized binary with CUDA GPU support when a GPU runtime is selected.",
    pros: ["CUDA GPU support, fast", "Compatible with CMDATA"],
    cons: ["~20 min compile time"],
  },
];

const STEPS = [
  {
    number: "01",
    title: "Install GROMACS",
    description:
      "Choose between a quick apt install (~2 min) or compiling GROMACS " +
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
    title: "Upload/select/generate a PDB/GRO file",
    description:
      "Upload the PDB/GRO file of your protein. One can also use a provided example or generate a linear polymer.",
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
      "mego --egos mg reads the GROMACS topology and writes topol_mego.top and ffnonbonded.itp " +
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

export default function Simulation() {
  return (
    <div className="mx-auto max-w-6xl space-y-16 px-6 py-12">
      {/* Header */}
      <div className="space-y-4">
        <h1 className="section-heading">Run a Multi-<em>e</em>GO Simulation</h1>
        <p className="max-w-2xl text-gray-400">
          Run multi-eGO simulations directly in your browser — no local GROMACS installation
          required. Two Google Colab notebooks are available: one for your own protein starting
          from a raw PDB or GRO file, and one for the bundled example systems with pre-computed
          contact matrices ready to simulate immediately.
        </p>
      </div>

      {/* Colab notebooks */}
      <div>
        <p className="mb-2 text-center text-sm font-medium uppercase tracking-widest text-brand-400">
          Run in your browser
        </p>
        <p className="mb-6 text-center text-xs text-gray-500">
          Free GPU available · requires a Google account · no local installation needed
        </p>
        <div className="grid gap-6 lg:grid-cols-2">
          {NOTEBOOKS.map((nb) => (
            <div
              key={nb.url}
              className="flex flex-col justify-between rounded-xl border border-brand-800 bg-brand-950/30 p-8"
            >
              <div>
                <p className="mb-1 text-xs font-semibold uppercase tracking-widest text-brand-400">
                  {nb.tag}
                </p>
                <h2 className="mb-3 text-xl font-bold text-white">
                  Multi-<em>e</em>GO · {nb.title}
                </h2>
                <p className="mb-4 text-sm text-gray-400">{nb.description}</p>
                <p className="mb-6 text-xs text-gray-500">{nb.steps}</p>
              </div>
              <a
                href={nb.url}
                target="_blank"
                rel="noopener noreferrer"
                className="inline-flex items-center justify-center gap-3 rounded-lg bg-[#F9AB00] px-6 py-3 text-sm font-semibold text-gray-900 transition hover:bg-[#F9AB00]/90"
              >
                <img
                  src="https://colab.research.google.com/assets/colab-badge.svg"
                  alt="Open in Colab"
                  className="h-5"
                />
                Open in Google Colab
              </a>
            </div>
          ))}
        </div>
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

      {/* What the mg notebook does */}
      <div>
        <h2 className="section-heading mb-1">What the mg simulation notebook does</h2>
        <p className="mb-6 text-sm text-gray-500">
          Applies to the <em>mg reference simulation</em> notebook (your own protein).
        </p>
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

      {/* What the examples notebook does */}
      <div>
        <h2 className="section-heading mb-1">What the examples notebook does</h2>
        <p className="mb-6 text-sm text-gray-500">
          Applies to the <em>production simulation of example systems</em> notebook.
        </p>
        <div className="grid gap-4 sm:grid-cols-2 lg:grid-cols-3">
          {[
            {
              number: "01",
              title: "Install GROMACS + multi-eGO",
              description:
                "Same two-option GROMACS installation as the mg notebook, followed by " +
                "cloning and pip-installing the multi-eGO package.",
            },
            {
              number: "02",
              title: "Select an example system",
              description:
                "Choose from GB1, Aβ42, TTR peptide, or " +
                "Lysozyme+benzene. The config, reference topology, and starting structure " +
                "are set automatically.",
            },
            {
              number: "03",
              title: "Generate production force field (mego)",
              description:
                "mego --config reads the pre-computed contact matrices and reference " +
                "topology bundled in the repository and writes topol_mego.top and " +
                "ffnonbonded.itp.",
            },
            {
              number: "04",
              title: "Energy minimisation",
              description:
                "Steepest-descent minimisation with the production multi-eGO force field " +
                "removes bad contacts in the starting structure.",
            },
            {
              number: "05",
              title: "NVT production run",
              description:
                "Langevin dynamics at the chosen temperature. GPU is used automatically " +
                "when a Colab GPU runtime is active. Produces run.xtc and run.tpr.",
            },
          ].map((step) => (
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
            {OUTPUT_FILES.map(([f, desc]) => (
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
