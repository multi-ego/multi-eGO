import { useState } from "react";
import { Link, NavLink, useLocation } from "react-router-dom";

const links = [
  { to: "/", label: "Home", end: true },
  { to: "/install", label: "Installation" },
  { to: "/setup", label: "Setup" },
  { to: "/config", label: "Config Builder" },
  { to: "/examples", label: "Examples" },
  { to: "/simulation", label: "Simulation" },
];

const GitHubIcon = () => (
  <svg className="h-5 w-5" fill="currentColor" viewBox="0 0 24 24" aria-hidden="true">
    <path
      fillRule="evenodd"
      clipRule="evenodd"
      d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z"
    />
  </svg>
);

export default function Navbar() {
  const [open, setOpen] = useState(false);
  const { pathname } = useLocation();

  // Close the mobile menu whenever the route changes
  // (useLocation re-renders the component on navigation)
  const close = () => setOpen(false);

  return (
    <header className="sticky top-0 z-50 border-b border-gray-800 bg-gray-950/90 backdrop-blur">
      <nav className="mx-auto max-w-6xl px-6 py-3">
        {/* ── Top bar ────────────────────────────────────── */}
        <div className="flex items-center justify-between">
          <Link to="/" onClick={close} className="font-bold text-white text-lg">
            <span className="text-brand-400">Multi-<em>e</em>GO</span>
          </Link>

          {/* Desktop links — hidden below md */}
          <div className="hidden md:flex items-center gap-6 text-sm font-medium">
            {links.map(({ to, label, end }) => (
              <NavLink
                key={to}
                to={to}
                end={end}
                className={({ isActive }) =>
                  isActive ? "text-brand-400" : "text-gray-400 hover:text-white transition"
                }
              >
                {label}
              </NavLink>
            ))}
            <a
              href="https://github.com/multi-ego/multi-eGO"
              target="_blank"
              rel="noopener noreferrer"
              className="text-gray-400 hover:text-white transition"
              aria-label="GitHub repository"
            >
              <GitHubIcon />
            </a>
          </div>

          {/* Hamburger button — visible below md */}
          <button
            className="md:hidden flex items-center justify-center rounded p-2 text-gray-400 hover:text-white hover:bg-gray-800 transition"
            aria-label={open ? "Close menu" : "Open menu"}
            aria-expanded={open}
            onClick={() => setOpen((o) => !o)}
          >
            {open ? (
              /* X icon */
              <svg className="h-5 w-5" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" d="M6 18L18 6M6 6l12 12" />
              </svg>
            ) : (
              /* Hamburger icon */
              <svg className="h-5 w-5" fill="none" stroke="currentColor" strokeWidth={2} viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" d="M4 6h16M4 12h16M4 18h16" />
              </svg>
            )}
          </button>
        </div>

        {/* ── Mobile dropdown — visible below md when open ─ */}
        {open && (
          <div className="md:hidden mt-3 flex flex-col gap-1 border-t border-gray-800 pt-3 text-sm font-medium">
            {links.map(({ to, label, end }) => (
              <NavLink
                key={to}
                to={to}
                end={end}
                onClick={close}
                className={({ isActive }) =>
                  "rounded px-3 py-2 transition " +
                  (isActive
                    ? "bg-gray-800 text-brand-400"
                    : "text-gray-400 hover:bg-gray-800 hover:text-white")
                }
              >
                {label}
              </NavLink>
            ))}
            <a
              href="https://github.com/multi-ego/multi-eGO"
              target="_blank"
              rel="noopener noreferrer"
              onClick={close}
              className="flex items-center gap-2 rounded px-3 py-2 text-gray-400 hover:bg-gray-800 hover:text-white transition"
            >
              <GitHubIcon />
              <span>GitHub</span>
            </a>
          </div>
        )}
      </nav>
    </header>
  );
}
