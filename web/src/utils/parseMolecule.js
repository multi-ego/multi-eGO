/**
 * Parse a PDB or GRO file and return a summary of chains and residues.
 *
 * Returns:
 *   {
 *     format: "pdb" | "gro",
 *     title: string,
 *     chains: [{ id, residues: [{ name, number }] }],
 *     atomCount: number,
 *     warnings: string[],
 *   }
 */

export function parseMoleculeFile(filename, content) {
  const ext = filename.split(".").pop().toLowerCase();
  if (ext === "pdb" || ext === "ent") return parsePDB(content);
  if (ext === "gro") return parseGRO(content);
  throw new Error(`Unsupported file type: .${ext}. Please upload a .pdb or .gro file.`);
}

function parsePDB(content) {
  const lines = content.split("\n");
  const warnings = [];
  let title = "";
  const chainMap = new Map(); // chainId → Map(resKey → {name, number})
  let atomCount = 0;
  const hetatmResidues = new Set();

  for (const line of lines) {
    const rec = line.substring(0, 6).trim();

    if (rec === "TITLE") {
      title = line.substring(10).trim();
    }

    if (rec === "ATOM" || rec === "HETATM") {
      atomCount++;
      const chainId = line[21] ?? " ";
      const resName = line.substring(17, 20).trim();
      const resNum = parseInt(line.substring(22, 26).trim(), 10);

      if (rec === "HETATM") hetatmResidues.add(`${chainId}:${resNum}:${resName}`);

      if (!chainMap.has(chainId)) chainMap.set(chainId, new Map());
      const resKey = `${resNum}:${resName}`;
      if (!chainMap.get(chainId).has(resKey)) {
        chainMap.get(chainId).set(resKey, { name: resName, number: resNum });
      }
    }
  }

  if (hetatmResidues.size > 0) {
    warnings.push(
      `Found ${hetatmResidues.size} HETATM residue(s). These are typically ligands, waters, or modified residues. pdb2gmx may not recognise them — you may need to remove or handle them separately.`
    );
  }

  const chains = Array.from(chainMap.entries()).map(([id, resMap]) => ({
    id: id.trim() || "—",
    residues: Array.from(resMap.values()).sort((a, b) => a.number - b.number),
  }));

  if (chains.length === 0) warnings.push("No ATOM or HETATM records found.");

  return { format: "pdb", title, chains, atomCount, warnings };
}

function parseGRO(content) {
  const lines = content.split("\n").filter((l) => l.length > 0);
  const warnings = [];

  const title = lines[0]?.trim() ?? "";
  const atomCount = parseInt(lines[1]?.trim(), 10) || 0;

  // GRO has no chain concept — group by residue number changes
  const residues = [];
  const seen = new Set();

  for (let i = 2; i < 2 + atomCount && i < lines.length; i++) {
    const line = lines[i];
    if (line.length < 20) continue;
    const resNum = parseInt(line.substring(0, 5).trim(), 10);
    const resName = line.substring(5, 10).trim();
    const key = `${resNum}:${resName}`;
    if (!seen.has(key)) {
      seen.add(key);
      residues.push({ name: resName, number: resNum });
    }
  }

  if (residues.length === 0) warnings.push("No residues could be parsed from this GRO file.");

  warnings.push("GRO files have no chain information. If your system has multiple chains, use a PDB file for better analysis.");

  const chains = [{ id: "—", residues }];
  return { format: "gro", title, chains, atomCount, warnings };
}
