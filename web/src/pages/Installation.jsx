export default function Installation() {
  return (
    <div className="mx-auto max-w-4xl px-6 py-12 space-y-12">
      <div>
        <h1 className="text-4xl font-bold text-white mb-2">Installation</h1>
        <p className="text-gray-400">
          Multi-<em>e</em>GO requires Python ≥ 3.10 and GROMACS &gt; 2022.
        </p>
      </div>

      {/* Python environment */}
      <section className="space-y-6">
        <h2 className="text-2xl font-semibold text-white border-b border-gray-800 pb-2">
          Python environment
        </h2>
        <p className="text-gray-400">Choose your preferred environment manager:</p>

        <div>
          <p className="label mb-2">Conda (recommended)</p>
          <pre className="code-block whitespace-pre">{`conda env create -f conda/environment.yml
conda activate meGO`}</pre>
        </div>

        <div>
          <p className="label mb-2">pip</p>
          <pre className="code-block">pip install -r requirements.txt</pre>
        </div>
      </section>

      {/* GROMACS */}
      <section className="space-y-4">
        <h2 className="text-2xl font-semibold text-white border-b border-gray-800 pb-2">GROMACS</h2>
        <p className="text-gray-400">
          GROMACS &gt; 2022 is required to run simulations. It must be compiled from source with the
          following CMake flags to enable the legacy and nblib APIs needed by multi-eGO:
        </p>
        <pre className="code-block whitespace-pre">{`cmake .. \\
    -DGMX_INSTALL_LEGACY_API=ON \\
    -DBUILD_SHARED_LIBS=ON \\
    -DGMX_INSTALL_NBLIB_API=ON`}</pre>
        <p className="text-sm text-gray-500">
          Refer to the official guide for the full build procedure (toolchain, GPU support, MPI, etc.):
        </p>
        <a
          href="https://manual.gromacs.org/current/install-guide/index.html"
          target="_blank"
          rel="noopener noreferrer"
          className="inline-flex items-center gap-2 text-brand-400 hover:underline text-sm"
        >
          GROMACS installation guide ↗
        </a>
      </section>

      {/* cmdata */}
      <section className="space-y-4">
        <h2 className="text-2xl font-semibold text-white border-b border-gray-800 pb-2">
          cmdata — contact matrix extraction
        </h2>
        <p className="text-gray-400">
          <code className="text-gray-200">cmdata</code> is the trajectory analysis tool used to extract
          contact histograms from MD simulations. It requires recompiling GROMACS from source with the
          cmdata patch applied.
        </p>
        <div className="card space-y-3 text-sm text-gray-400">
          <p>
            1. Clone the repository and navigate to the cmdata tool directory:
          </p>
          <pre className="code-block">{`git clone https://github.com/multi-ego/multi-eGO.git
cd multi-eGO/tools/cmdata`}</pre>
          <p>2. Follow the instructions in the README to patch and recompile GROMACS:</p>
          <a
            href="https://github.com/multi-ego/multi-eGO/blob/main/tools/cmdata/README.md"
            target="_blank"
            rel="noopener noreferrer"
            className="text-brand-400 hover:underline"
          >
            cmdata README ↗
          </a>
        </div>
      </section>

      {/* Force field */}
      <section className="space-y-4">
        <h2 className="text-2xl font-semibold text-white border-b border-gray-800 pb-2">
          multi-ego-basic force field
        </h2>
        <p className="text-gray-400">
          The <code className="text-gray-200">multi-ego-basic.ff</code> directory is included in the
          repository root. GROMACS needs to find it when preparing your topology. The recommended approach
          is to set the <code className="text-gray-200">GMXLIB</code> environment variable to the repository
          root so GROMACS can always locate it, regardless of your working directory:
        </p>
        <pre className="code-block">{`# Add to your shell profile (~/.bashrc, ~/.zshrc, etc.)
export GMXLIB=/path/to/multi-eGO`}</pre>
        <p className="text-gray-400">
          Alternatively, run <code className="text-gray-200">pdb2gmx</code> from the repository root, or
          copy <code className="text-gray-200">multi-ego-basic.ff</code> to your working directory or your
          GROMACS force field path.
        </p>
        <pre className="code-block">{`# With GMXLIB set, this works from any directory
gmx pdb2gmx -f my_protein.pdb -ff multi-ego-basic -water none`}</pre>
      </section>

      {/* Quick verify */}
      <section className="space-y-4">
        <h2 className="text-2xl font-semibold text-white border-b border-gray-800 pb-2">
          Verify the installation
        </h2>
        <pre className="code-block">{`# Should print usage information
python multiego.py --help

# Run the test suite
pytest tests/`}</pre>
      </section>
    </div>
  );
}
