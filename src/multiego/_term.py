"""
_term.py — lightweight terminal-output helpers for multi-eGO.

All colour/style codes are suppressed automatically when stdout is not a TTY
(e.g. when output is piped or redirected) or when the ``NO_COLOR`` environment
variable is set, following the https://no-color.org convention.

Public API
----------
header(msg)      — top-level pipeline step  (bold)
sub(msg)         — second-level step        (normal, indented 2 spaces)
path(p)          — file/directory path      (dim, indented 4 spaces)
timing(s)        — "done in X.XX s" line    (dim green, indented 4 spaces)
section_timing(s)— "done in X.XX s" for top-level sections (dim green, 2 spaces)
note(msg)        — informational note       (dim, indented 4 spaces)
warn(msg)        — warning                  (yellow)
error(msg)       — error                    (bold red)
success(msg)     — success / final line     (bold green)
rule(msg)        — full-width rule with centred label
spinner(msg)     — context manager: shows an animated spinner while a block runs
"""

import itertools
import os
import sys
import threading
import time

# ---------------------------------------------------------------------------
# Colour/style enable check
# ---------------------------------------------------------------------------


def _colours_enabled() -> bool:
    if os.environ.get("NO_COLOR"):
        return False
    return hasattr(sys.stdout, "isatty") and sys.stdout.isatty()


_USE_COLOR = _colours_enabled()

# ANSI codes
_RESET = "\033[0m"
_BOLD = "\033[1m"
_DIM = "\033[2m"
_RED = "\033[31m"
_GREEN = "\033[32m"
_YELLOW = "\033[33m"
_CYAN = "\033[36m"


def _c(*codes: str) -> str:
    """Return the concatenation of ANSI codes, or '' when colours are off."""
    return "".join(codes) if _USE_COLOR else ""


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def header(msg: str) -> None:
    """Top-level pipeline step, e.g. '- Processing contact matrices'."""
    print(f"{_c(_BOLD)}- {msg}{_c(_RESET)}")


def sub(msg: str) -> None:
    """Second-level step, indented two spaces."""
    print(f"  - {msg}")


def path(p: str) -> None:
    """File or directory path, indented four spaces, dim."""
    print(f"{_c(_DIM)}    {p}{_c(_RESET)}")


def timing(seconds: float, indent: int = 4) -> None:
    """'done in X.XX s' line at the given indent level."""
    pad = " " * indent
    print(f"{_c(_DIM, _GREEN)}{pad}done in {seconds:.2f} s{_c(_RESET)}")


def section_timing(seconds: float) -> None:
    """'done in X.XX s' at two-space indent (for top-level sections)."""
    timing(seconds, indent=2)


def note(msg: str) -> None:
    """Informational note, indented four spaces, dim."""
    print(f"{_c(_DIM)}    {msg}{_c(_RESET)}")


def warn(msg: str) -> None:
    """Warning line in yellow."""
    print(f"{_c(_YELLOW)}  WARNING: {msg}{_c(_RESET)}")


def error(msg: str) -> None:
    """Error line in bold red."""
    print(f"{_c(_BOLD, _RED)}ERROR: {msg}{_c(_RESET)}")


def success(msg: str) -> None:
    """Success / final summary line in bold green."""
    print(f"{_c(_BOLD, _GREEN)}{msg}{_c(_RESET)}")


def rule(msg: str = "", width: int = 60) -> None:
    """Print a full-width horizontal rule with an optional centred label."""
    if msg:
        pad = max(0, width - len(msg) - 2)
        left = pad // 2
        right = pad - left
        line = "─" * left + f" {msg} " + "─" * right
    else:
        line = "─" * width
    print(f"{_c(_DIM)}{line}{_c(_RESET)}")


# ---------------------------------------------------------------------------
# Spinner context manager
# ---------------------------------------------------------------------------

_SPINNER_FRAMES = ("⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏")


class spinner:
    """Context manager that shows an animated spinner while work is in progress.

    Usage::

        with _term.spinner("Loading contact matrix"):
            heavy_computation()

    The spinner is suppressed silently when stdout is not a TTY so that log
    files are not cluttered with ``\\r`` characters.
    """

    def __init__(self, msg: str, interval: float = 0.08) -> None:
        self._msg = msg
        self._interval = interval
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None

    def _spin(self) -> None:
        for frame in itertools.cycle(_SPINNER_FRAMES):
            if self._stop_event.is_set():
                break
            sys.stdout.write(f"\r{_c(_CYAN)}{frame}{_c(_RESET)}  {self._msg} ")
            sys.stdout.flush()
            time.sleep(self._interval)

    def __enter__(self) -> "spinner":
        if _USE_COLOR:
            self._thread = threading.Thread(target=self._spin, daemon=True)
            self._thread.start()
        else:
            # Non-TTY: just print the label once so log files stay readable
            print(f"  {self._msg} ...")
        return self

    def __exit__(self, *_) -> None:
        if _USE_COLOR:
            self._stop_event.set()
            if self._thread:
                self._thread.join()
            # Erase the spinner line
            sys.stdout.write("\r" + " " * (len(self._msg) + 6) + "\r")
            sys.stdout.flush()
