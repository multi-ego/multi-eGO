export default function Footer() {
  return (
    <footer className="border-t border-gray-800 bg-gray-950">
      <div className="mx-auto max-w-6xl px-6 py-8 flex flex-col sm:flex-row items-center justify-between gap-4 text-sm text-gray-500">
        <p>
          Multi-<em>e</em>GO — data-driven force fields for molecular dynamics
        </p>
        <div className="flex gap-4">
          <a
            href="https://github.com/multi-ego/multi-eGO"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-gray-300 transition"
          >
            GitHub
          </a>
          <a
            href="https://github.com/multi-ego/multi-eGO/blob/main/LICENSE"
            target="_blank"
            rel="noopener noreferrer"
            className="hover:text-gray-300 transition"
          >
            GPL v3
          </a>
        </div>
      </div>
    </footer>
  );
}
