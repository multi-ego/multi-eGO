#!/usr/bin/env python3
"""
Update the reference test outputs in tests/test_outputs/.

Reads test_cases.txt, runs multiego for every active test case, then copies
the generated output directories into test_outputs/ to serve as the new
reference for future test runs.

Usage (from any directory):
    python tests/update_testing_data.py
"""

import os
import sys
import shutil
import subprocess
import yaml

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
MEGO_ROOT = os.path.abspath(os.path.join(TEST_DIR, os.pardir))
TEST_CASES = os.path.join(TEST_DIR, "test_cases.txt")


def _system_from_config(config_path):
    with open(config_path) as f:
        data = yaml.safe_load(f)
    for item in data:
        if isinstance(item, dict) and "system" in item:
            return item["system"]
    raise ValueError(f"No 'system' key found in {config_path}")


def read_test_cases():
    """Return (commands, systems) for every active line in test_cases.txt."""
    commands, systems = [], []
    with open(TEST_CASES) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.replace("TEST_ROOT", TEST_DIR) for p in line.split()]
            if "--system" in parts:
                system = parts[parts.index("--system") + 1]
            elif "--config" in parts:
                config_path = parts[parts.index("--config") + 1]
                system = _system_from_config(config_path)
            else:
                continue
            commands.append(parts)
            systems.append(system)
    return commands, systems


def main():
    commands, systems = read_test_cases()
    unique_systems = sorted(set(systems))

    # 1. Clear existing output and reference directories for active systems only
    print("Clearing output directories...")
    for system in unique_systems:
        for path in (
            os.path.join(MEGO_ROOT, "outputs", system),
            os.path.join(TEST_DIR, "test_outputs", system),
        ):
            if os.path.exists(path):
                print(f"  rm -rf {path}")
                shutil.rmtree(path)

    # 2. Run each test case
    print("\nRunning test cases...")
    for cmd, system in zip(commands, systems):
        # --system runs need --inputs_dir; --config runs derive it automatically
        if "--system" in cmd and "--inputs_dir" not in cmd:
            cmd = cmd + ["--inputs_dir", os.path.join(TEST_DIR, "test_inputs")]

        print(f"\n  {' '.join(cmd)}")
        result = subprocess.run([sys.executable, os.path.join(MEGO_ROOT, "multiego.py"), *cmd])
        if result.returncode != 0:
            print(f"ERROR: command failed (exit {result.returncode})", file=sys.stderr)
            sys.exit(1)

    # 3. Copy outputs → test_outputs (no stdout parsing needed: we know the structure)
    print("\nCopying outputs to test_outputs/...")
    for system in unique_systems:
        src_system = os.path.join(MEGO_ROOT, "outputs", system)
        dst_system = os.path.join(TEST_DIR, "test_outputs", system)
        if not os.path.exists(src_system):
            print(f"WARNING: no output found for {system}", file=sys.stderr)
            continue
        os.makedirs(dst_system, exist_ok=True)
        for case_dir in sorted(os.listdir(src_system)):
            src = os.path.join(src_system, case_dir)
            dst = os.path.join(dst_system, case_dir)
            if os.path.isdir(src):
                print(f"  {src} -> {dst}")
                shutil.copytree(src, dst)

    print("\nFinished updating test outputs!")


if __name__ == "__main__":
    main()
