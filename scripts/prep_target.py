#!/usr/bin/env python3
"""Stage 1 - target preparation, end to end.

    python3 -I scripts/prep_target.py          # then:  pymol data/target/view.pml

Steps
  1. cache 1CX8 (apo), 3S9L (TfR1 + transferrin), and 11QC (competition reference)
  2. read domain boundaries from SCOP via PDBe       -> apical = the PA domain
  3. extract the apical domain, single chain, ATOM records only
  4. align mouse TfR1 (Q62351) to the human span
  5. build the avoid list: apical residues that touch transferrin and/or the other TfR1 chain
  6. rank conserved surface patches in biological assemblies; write hotspots and views

Everything lands in data/ (git-ignored). Standard library only; needs network on first run.
"""
import json
import math
import csv
import hashlib
import itertools
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW, OUT = ROOT / "data/raw", ROOT / "data/target"

CHAIN = "C"          # copy of 1CX8 to use: best ordered of the 8 (lowest mean B-factor)
CUTOFF = 5.0         # Angstrom contact definition
TF_SITE = (619, 670)  # transferrin docking site in the helical domain, from the docs
SCOP_NAMES = {"53210": "protease_like", "52026": "apical", "47673": "helical"}
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
         "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


# ---------- helpers ----------
def fetch(url, dest=None):
    if dest and dest.exists():
        return dest.read_bytes()
    data = urllib.request.urlopen(url, timeout=60).read()
    if dest:
        dest.write_bytes(data)
    return data


def uniprot(acc):
    return "".join(fetch(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", RAW / f"{acc}.fasta").decode().splitlines()[1:])


def read_atoms(path):
    """{(chain, resnum): {"resn": str, "xyz": [(x, y, z), ...], "ca": (x, y, z)|None}}; ATOM only, primary altloc."""
    res = {}
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or line[16] not in " A":
            continue
        r = res.setdefault((line[21], int(line[22:26])), {"resn": line[17:20], "xyz": [], "names": [], "elements": [], "ca": None})
        xyz = tuple(float(line[i:i + 8]) for i in (30, 38, 46))
        r["xyz"].append(xyz)
        r["names"].append(line[12:16].strip())
        r["elements"].append(line[76:78].strip() or line[12:16].strip()[0])
        if line[12:16].strip() == "CA":
            r["ca"] = xyz
    return res


def touching(res, chain, span, other_chains, cutoff=CUTOFF):
    """Residues of `chain` within `span` having any atom within cutoff of any atom of `other_chains`."""
    others = [x for (c, _), r in res.items() if c in other_chains for x in r["xyz"]]
    return {n for (c, n), r in res.items() if c == chain and span[0] <= n <= span[1]
            and any(math.dist(a, b) < cutoff for a in r["xyz"] for b in others)}


def align(a, b, match=2, mismatch=-1, gap=-4, with_options=False):
    """Global (Needleman-Wunsch) alignment, linear gap penalty."""
    n, m = len(a), len(b)
    s = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        s[i][0] = i * gap
    for j in range(1, m + 1):
        s[0][j] = j * gap
    sub = lambda i, j: match if a[i - 1] == b[j - 1] else mismatch
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i][j] = max(s[i - 1][j - 1] + sub(i, j), s[i - 1][j] + gap, s[i][j - 1] + gap)
    ra, rb, i, j = [], [], n, m
    while i or j:
        if i and j and s[i][j] == s[i - 1][j - 1] + sub(i, j):
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif i and s[i][j] == s[i - 1][j] + gap:
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1
    result = "".join(reversed(ra)), "".join(reversed(rb))
    if not with_options:
        return result
    # Traverse the optimal traceback DAG, without enumerating exponentially many
    # alignments. A chosen traceback must not hide equally scoring gap placements.
    options, pending, seen = defaultdict(set), [(n, m)], set()
    while pending:
        i, j = pending.pop()
        if (i, j) in seen:
            continue
        seen.add((i, j))
        if i and j and s[i][j] == s[i - 1][j - 1] + sub(i, j):
            options[i].add((j, b[j - 1])); pending.append((i - 1, j - 1))
        if i and s[i][j] == s[i - 1][j] + gap:
            options[i].add((None, "-")); pending.append((i - 1, j))
        if j and s[i][j] == s[i][j - 1] + gap:
            pending.append((i, j - 1))
    return *result, options


# ---------- steps ----------
def step_download():
    RAW.mkdir(parents=True, exist_ok=True)
    for pdb in ("1CX8", "3S9L", "11QC"):
        dest = RAW / f"{pdb}.pdb"
        if not dest.exists():
            fetch(f"https://files.rcsb.org/download/{pdb}.pdb", dest)
            print(f"  downloaded {pdb}")


def step_domains():
    d = json.loads(fetch("https://www.ebi.ac.uk/pdbe/api/mappings/1cx8", RAW / "1cx8_mappings.json"))["1cx8"]["SCOP"]
    dom = defaultdict(list)
    for sunid, v in d.items():
        for m in v["mappings"]:
            if m["chain_id"] == CHAIN:   # author numbering = numbering in the PDB file
                dom[SCOP_NAMES[sunid]].append((m["start"]["author_residue_number"], m["end"]["author_residue_number"]))
    dom = {k: sorted(v) for k, v in dom.items()}
    (a, b), = dom["apical"]
    assert not (a <= TF_SITE[1] and TF_SITE[0] <= b), "apical range overlaps the transferrin site"
    for k, v in dom.items():
        print(f"  {k:14s} {v}")
    return dom


def step_extract(dom, apo):
    (start, end), = dom["apical"]
    keep = [n for (c, n) in sorted(apo) if c == CHAIN and start <= n <= end]
    assert keep == list(range(start, end + 1)), "gap or missing residues in apical range"
    lines = [l for l in (RAW / "1CX8.pdb").read_text().splitlines()
             if l.startswith("ATOM") and l[21] == CHAIN and start <= int(l[22:26]) <= end and l[16] in " A"]
    (OUT / "human_apical.pdb").write_text("\n".join(lines) + "\nTER\nEND\n")
    seq = "".join(THREE[apo[(CHAIN, n)]["resn"]] for n in keep)
    (OUT / "human_apical.fasta").write_text(f">TfR1_human_{start}-{end}_1CX8_{CHAIN}\n{seq}\n")
    ca = [apo[(CHAIN, n)]["ca"] for n in keep]
    breaks = [n for n, p, q in zip(keep, ca, ca[1:]) if math.dist(p, q) > 4.2]
    print(f"  {len(keep)} residues {start}-{end}, {len(lines)} atoms, chain breaks: {breaks or 'none'}, "
          f"end-to-end {math.dist(ca[0], ca[-1]):.1f} A")
    return seq


def step_mouse(start, end, human):
    human_full, mouse_full = uniprot("P02786"), uniprot("Q62351")
    assert human_full[start - 1:end] == human, "1CX8 numbering != UniProt P02786 numbering"
    ah, am, options = align(human_full, mouse_full, with_options=True)
    hpos, mpos, rows, mapping = 0, 0, [], []
    for h, m in zip(ah, am):
        hpos += h != "-"
        mpos += m != "-"
        if start <= hpos <= end and not (h == "-" and hpos == end):
            rows.append((hpos if h != "-" else None, h, m))
            mapping.append({"human_pos": hpos if h != "-" else "", "human_aa": h,
                            "mouse_pos": mpos if m != "-" else "", "mouse_aa": m,
                            "unique_optimal_mapping": h != "-" and len(options[hpos]) == 1})
    mouse = "".join(m for _, _, m in rows if m != "-")
    (OUT / "mouse_apical.fasta").write_text(f">TfR1_mouse_Q62351_aligned_to_human_{start}-{end}\n{mouse}\n")
    same = sum(h == m for _, h, m in rows)
    diff = [f"{h}{p}{m}" for p, h, m in rows if "-" not in (h, m) and h != m]
    txt = [f"human {start}-{end} vs mouse Q62351: {same}/{len(rows)} identical ({100 * same / len(rows):.0f}%)",
           "differences (human residue, 1CX8 position, mouse residue): " + " ".join(diff), "",
           "human_pos  hum  mouse"]
    txt += [f"{p or '-':>9}  {h}    {m}   {'' if h == m else '*'}" for p, h, m in rows]
    (OUT / "human_mouse_apical_alignment.txt").write_text("\n".join(txt) + "\n")
    with (OUT / "residue_mapping.tsv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(mapping[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(mapping)
    ambiguity = {str(p): [{"mouse_pos": n, "mouse_aa": aa} for n, aa in sorted(options[p], key=lambda x: x[0] or 0)]
                 for p in range(start, end + 1) if len(options[p]) > 1}
    (OUT / "alignment_ambiguities.json").write_text(json.dumps(ambiguity, indent=2) + "\n")
    print(f"  mouse {len(mouse)} aa, {same}/{len(rows)} identical, {len(diff)} substitutions; "
          f"gaps at {[(p, h, m) for p, h, m in rows if '-' in (h, m)]}")
    return {p: (h, m) for p, h, m in rows if p}


def step_avoid(dom, apo, pairs):
    """Contacts in BOTH biological assemblies, never between unrelated crystal copies."""
    cx = read_atoms(RAW / "3S9L.pdb")
    (start, end), = dom["apical"]
    # numbering check: 3S9L uses the same author numbering as 1CX8
    same_seq = all(cx[("A", n)]["resn"] == apo[(CHAIN, n)]["resn"] for n in range(start, end + 1))
    assert same_seq, "3S9L and 1CX8 apical sequences/numbering differ"
    seam, tf = set(), set()
    for assembly_id, receptor, ligand in ((1, "A", "C"), (2, "B", "D")):
        assembly = biological_assembly(RAW / "3S9L.pdb", assembly_id)
        receptors = [c for c in {c for c, _ in assembly} if c.startswith(receptor + ":")]
        ligands = [c for c in {c for c, _ in assembly} if c.startswith(ligand + ":")]
        for me in receptors:
            seam |= touching(assembly, me, (start, end), set(receptors) - {me})
            tf |= touching(assembly, me, (start, end), ligands)
    both = seam & tf
    cat = {n: ("both" if n in both else "seam" if n in seam else "transferrin")
           for n in sorted(seam | tf)}
    out = {"source": "PDB 3S9L, REMARK 350 biological assemblies 1 and 2, 5 A heavy-atom contacts, UniProt numbering",
           "seam_with_other_tfr1_chain": sorted(seam), "touch_transferrin": sorted(tf),
           "both": sorted(both), "category": {str(n): c for n, c in cat.items()},
           "mouse_substitution": {str(n): pairs[n][0] + str(n) + pairs[n][1] for n in cat if pairs[n][0] != pairs[n][1]}}
    (OUT / "avoid_residues.json").write_text(json.dumps(out, indent=1))
    print(f"  seam {len(seam)}, transferrin {len(tf)}, both {len(both)}, avoid total {len(cat)} of {end - start + 1}")
    return cat


def biological_assembly(path, assembly_id):
    """Apply the PDB's REMARK 350 BIOMT operations; retain distinct copy IDs."""
    active, chains, operations = False, [], []
    for line in path.read_text().splitlines():
        if not line.startswith("REMARK 350"):
            continue
        if "BIOMOLECULE:" in line:
            active = int(line.split(":")[-1]) == assembly_id
        elif active and "APPLY THE FOLLOWING TO CHAINS:" in line:
            chains = [c.strip() for c in line.split(":")[-1].split(",")]
        elif active and "BIOMT" in line:
            fields = line.split()
            if fields[2] == "BIOMT1":
                operations.append((chains[:], fields[3], []))
            operations[-1][2].append(tuple(map(float, fields[4:8])))
    assert operations, f"No assembly {assembly_id} in {path}"
    residues, out = read_atoms(path), {}
    for chains, op, matrix in operations:
        assert len(matrix) == 3
        def transform(x):
            return tuple(sum(row[i] * x[i] for i in range(3)) + row[3] for row in matrix)
        for (c, n), r in residues.items():
            if c in chains:
                out[(f"{c}:{op}", n)] = {**r, "xyz": [transform(x) for x in r["xyz"]],
                                        "ca": transform(r["ca"]) if r["ca"] else None}
    return out


def sidechain_sasa(residues, chain, span, samples=240):
    """Shrake-Rupley SASA (A^2), 1.4 A water probe, Fibonacci sphere.

    All protein atoms in the biological assembly occlude the queried side chains.
    Approximate geometric exposure, NOT binding energy or a pKa calculation.
    """
    radii = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}
    atoms, grid = [], defaultdict(list)
    cell = lambda x: tuple(math.floor(v / 7.0) for v in x)
    for key, r in residues.items():
        for xyz, name, element in zip(r["xyz"], r["names"], r["elements"]):
            if element == "H":
                continue
            atom = (xyz, radii[element] + 1.4, key, name)
            grid[cell(xyz)].append(len(atoms)); atoms.append(atom)
    sphere = []
    for k in range(samples):
        z = 1 - 2 * (k + 0.5) / samples
        phi = k * math.pi * (3 - math.sqrt(5))
        rho = math.sqrt(1 - z * z)
        sphere.append((rho * math.cos(phi), rho * math.sin(phi), z))
    area = defaultdict(float)
    for i, (xyz, radius, (c, n), name) in enumerate(atoms):
        if c != chain or not span[0] <= n <= span[1] or name in {"N", "CA", "C", "O", "OXT"}:
            continue
        bucket = cell(xyz)
        neighbors = []
        for dx, dy, dz in itertools.product((-1, 0, 1), repeat=3):
            for j in grid.get((bucket[0] + dx, bucket[1] + dy, bucket[2] + dz), []):
                other, r, _, _ = atoms[j]
                if i != j and math.dist(xyz, other) < radius + r:
                    neighbors.append((other, r * r))
        exposed = 0
        for unit in sphere:
            x, y, z = (xyz[k] + radius * unit[k] for k in range(3))
            if not any((x - q[0])**2 + (y - q[1])**2 + (z - q[2])**2 < r2 for q, r2 in neighbors):
                exposed += 1
        area[n] += 4 * math.pi * radius**2 * exposed / samples
    return area


def minimum_distance(residue, others):
    return min((math.dist(a, b) for a in residue["xyz"] for b in others), default=999.0)


def step_hotspots(dom, apo, pairs, cat, samples=240):
    """Find conserved, exposed candidate patches; no claim of experimental hotspots."""
    (start, end), = dom["apical"]
    contexts = [("1CX8_C", biological_assembly(RAW / "1CX8.pdb", 2), "C:1"),
                ("3S9L_A", biological_assembly(RAW / "3S9L.pdb", 1), "A:1"),
                ("3S9L_B", biological_assembly(RAW / "3S9L.pdb", 2), "B:1"),
                ("11QC_A", biological_assembly(RAW / "11QC.pdb", 1), "A:1"),
                ("11QC_B", biological_assembly(RAW / "11QC.pdb", 1), "B:1")]
    glycans = [tuple(float(line[i:i + 8]) for i in (30, 38, 46))
               for line in (RAW / "11QC.pdb").read_text().splitlines()
               if line.startswith("HETATM") and line[17:20] == "NAG" and line[16] in " A"]
    assert glycans, "Expected resolved NAG glycans in the competition reference"
    exposures, clearances = {}, {}
    for label, residues, chain in contexts:
        for n in range(start, end + 1):
            assert (chain, n) in residues, f"Missing residue {label}:{n}"
            assert THREE[residues[(chain, n)]["resn"]] == pairs[n][0]
        print(f"  surface exposure in {label}", flush=True)
        exposures[label] = sidechain_sasa(residues, chain, (start, end), samples)
        # Anything deleted by truncation: other chains AND the rest of this chain.
        omitted = [x for (c, n), r in residues.items()
                   if c != chain or not start <= n <= end for x in r["xyz"]]
        clearances[label] = {n: minimum_distance(residues[(chain, n)], omitted)
                             for n in range(start, end + 1)}
    mouse_map = {int(r["human_pos"]): r for r in
                 csv.DictReader((OUT / "residue_mapping.tsv").open(), delimiter="\t") if r["human_pos"]}
    rows = []
    for n in range(start, end + 1):
        ca = apo[(CHAIN, n)]["ca"]
        neighbors = [p for p in pairs if math.dist(ca, apo[(CHAIN, p)]["ca"]) <= 10]
        conservation = sum(pairs[p][0] == pairs[p][1] for p in neighbors) / len(neighbors)
        exposure = min(s.get(n, 0) for s in exposures.values())
        clearance = min(s[n] for s in clearances.values())
        # N251 and N317 have attached glycans in 11QC (LINK records).
        # A C-alpha buffer is a conservative flag, not a complete glycan envelope.
        glycan_distance = min(math.dist(ca, apo[(CHAIN, p)]["ca"]) for p in (251, 317))
        resolved_glycan_distance = min(minimum_distance(residues[(chain, n)], glycans)
                                      for label, residues, chain in contexts if label.startswith("11QC"))
        end_distance = min(math.dist(ca, apo[(CHAIN, p)]["ca"]) for p in (start, end))
        reasons = []
        if pairs[n][0] != pairs[n][1]: reasons.append("not_conserved")
        if mouse_map[n]["unique_optimal_mapping"] != "True": reasons.append("ambiguous_alignment")
        if exposure < 10: reasons.append("sidechain_sasa_below_10")
        if clearance < 5: reasons.append("within_5A_of_omitted_protein")
        if glycan_distance < 12: reasons.append("within_12A_CA_of_glycosylation_site")
        if resolved_glycan_distance < 5: reasons.append("within_5A_of_resolved_glycan")
        if end_distance < 8: reasons.append("within_8A_CA_of_cut_terminus")
        if n in cat: reasons.append("biological_contact_avoid")
        row = {"human_pos": n, "human_aa": pairs[n][0], "mouse_pos": mouse_map[n]["mouse_pos"],
               "mouse_aa": pairs[n][1], "conserved": pairs[n][0] == pairs[n][1],
               "min_sidechain_sasa_A2": round(exposure, 2), "min_omitted_distance_A": round(clearance, 2),
               "glycan_site_CA_distance_A": round(glycan_distance, 2),
               "min_resolved_glycan_distance_A": round(resolved_glycan_distance, 2),
               "local_conservation_10A": round(conservation, 3),
               "eligible": not reasons, "exclusion_reasons": ";".join(reasons)}
        row.update({f"sasa_{label}_A2": round(s.get(n, 0), 2) for label, s in exposures.items()})
        rows.append(row)
    with (OUT / "hotspot_candidates.tsv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    candidates = [r for r in rows if r["eligible"]]
    patches = []
    for triple in itertools.combinations(candidates, 3):
        positions = [r["human_pos"] for r in triple]
        aas = {r["human_aa"] for r in triple}
        # A hydrophobic anchor and a cationic sidechain are design hypotheses.
        if not aas & set("FWYILVM") or not aas & set("KR"):
            continue
        diameter = max(math.dist(residues[(chain, a)]["ca"], residues[(chain, b)]["ca"])
                       for _, residues, chain in contexts for a, b in itertools.combinations(positions, 2))
        if diameter > 12:
            continue
        score = sum(min(r["min_sidechain_sasa_A2"], 100) * r["local_conservation_10A"] for r in triple)
        patches.append({"human_positions": positions, "residues": [r["human_aa"] + str(r["human_pos"]) for r in triple],
                        "mouse_positions": [int(r["mouse_pos"]) for r in triple],
                        "heuristic_score": round(score, 2), "max_CA_diameter_A": round(diameter, 2),
                        "rfdiffusion_hotspot_res": [f"{CHAIN}{p}" for p in positions],
                        "crop_1based_positions": [p - start + 1 for p in positions]})
    patches.sort(key=lambda x: (-x["heuristic_score"], x["human_positions"]))
    selected = []
    for patch in patches:
        if all(len(set(patch["human_positions"]) & set(p["human_positions"])) <= 1 for p in selected):
            selected.append(patch)
        if len(selected) == 3: break
    assert selected, "No patches pass the documented criteria; review instead of inventing hotspots"
    output = {"status": "provisional_structure_derived_candidates_not_experimentally_validated",
              "not_the_unnamed_three": "The competition does not enumerate a complete nine-position set. These are independently derived patches.",
              "numbering": "human UniProt/PDB author numbering; chain C in human_apical.pdb; mouse UniProt Q62351",
              "method": {"sasa": "Shrake-Rupley sidechain only, 1.4 A probe, protein occlusion", "sphere_samples": samples,
                         "min_sasa_A2": 10, "min_omitted_protein_distance_A": 5,
                         "min_glycan_site_CA_distance_A": 12, "min_cut_terminus_CA_distance_A": 8,
                         "min_resolved_NAG_distance_A": 5,
                         "max_patch_CA_diameter_A": 12, "neighborhood_CA_radius_A": 10,
                         "score": "sum(min(minimum_sidechain_SASA,100) * local_conservation); heuristic, not affinity",
                         "chemistry": "at least one F/W/Y/I/L/V/M and one K/R per triplet"},
              "eligible_residues": len(candidates), "patches": selected,
              "input_sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(RAW.iterdir()) if p.is_file()},
              "limitations": ["Mouse conservation is sequence-only; local gap placement is alignment-dependent.",
                              "Glycan buffers do not model full glycan conformations.",
                              "Hotspots do not constrain the whole binder footprint; check full complexes.",
                              "Histidines must be placed on the binder after a complex exists; no pH prediction is made."]}
    (OUT / "hotspots.json").write_text(json.dumps(output, indent=2) + "\n")
    for i, patch in enumerate(selected, 1):
        print(f"  candidate patch {i}: {', '.join(patch['residues'])}")
    pml = [f"load {OUT}/human_apical.pdb, hotspot_target", "hide everything", "show cartoon, hotspot_target", "color grey80, hotspot_target"]
    for i, patch in enumerate(selected, 1):
        name = f"candidate_patch_{i}"
        pml += [f"select {name}, hotspot_target and resi {'+'.join(map(str, patch['human_positions']))}",
                f"show sticks, {name}", f"color {['teal', 'orange', 'violet'][i-1]}, {name}",
                f"label {name} and name CA, resn + resi"]
    pml += ["orient hotspot_target", "deselect"]
    (OUT / "hotspots.pml").write_text("\n".join(pml) + "\n")
    return output


def glycosylation_sites(path, chain):
    """Asn residue numbers carrying a NAG, read from the PDB LINK records."""
    return sorted({int(l[22:26]) for l in path.read_text().splitlines()
                   if l.startswith("LINK") and l[17:20] == "ASN" and l[21] == chain and l[47:50] == "NAG"})


def step_view(dom, cat):
    sel = lambda rs: " or ".join(f"resi {a}-{b}" for a, b in rs)
    pick = lambda c: "+".join(str(n) for n, k in cat.items() if k == c) or "0"
    (start, end), = dom["apical"]
    spans = [s for v in dom.values() for s in v]
    (OUT / "domains.json").write_text(json.dumps(
        {"source": "SCOP via PDBe for 1CX8; transferrin site from the project docs; glycosylation from LINK records",
         "structure": "1CX8", "chain": CHAIN,
         "domains": {k: [list(s) for s in v] for k, v in dom.items()},
         "transferrin_site": list(TF_SITE),
         "glycosylation_sites": glycosylation_sites(RAW / "1CX8.pdb", CHAIN),
         "termini": [min(a for a, _ in spans), max(b for _, b in spans)]}, indent=1) + "\n")
    pml = f"""# generated by scripts/prep_target.py - regenerate, do not edit
# PyMOL has no trailing comments, so comments sit on their own lines.
reinitialize
bg_color white
set cartoon_fancy_helices, 1
set sphere_scale, 0.5
set label_size, 20
set label_color, black

# object 1: the whole chain, coloured by domain
load {RAW}/1CX8.pdb, tfr1
remove tfr1 and not chain {CHAIN}
remove hetatm
hide everything, tfr1
show cartoon, tfr1
select protease_like, tfr1 and ({sel(dom['protease_like'])})
select apical, tfr1 and ({sel(dom['apical'])})
select helical, tfr1 and ({sel(dom['helical'])})
color grey70, protease_like
color teal, apical
color slate, helical
select transferrin_site, tfr1 and resi {TF_SITE[0]}-{TF_SITE[1]}
show sticks, transferrin_site
color orange, transferrin_site
label tfr1 and name CA and resi 122, "N"
label tfr1 and name CA and resi 760, "C"

# object 2: the isolated apical domain with the avoid list (hidden at first)
load {OUT}/human_apical.pdb, apical_only
hide everything, apical_only
show cartoon, apical_only
color teal, apical_only
select seam_only, apical_only and resi {pick('seam')}
select tf_only, apical_only and resi {pick('transferrin')}
select seam_and_tf, apical_only and resi {pick('both')}
show sticks, seam_only or tf_only or seam_and_tf
color magenta, seam_only
color orange, tf_only
color red, seam_and_tf
disable apical_only
deselect
orient tfr1
zoom tfr1

# to switch to the isolated domain, type:
#   disable tfr1
#   enable apical_only
#   orient apical_only
# legend: magenta = touches the other TfR1 chain, orange = touches transferrin, red = both
"""
    (OUT / "view.pml").write_text(pml)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print("[1/6] inputs");                  step_download()
    print("[2/6] domains (SCOP via PDBe)"); dom = step_domains()
    apo = read_atoms(RAW / "1CX8.pdb")
    (start, end), = dom["apical"]
    print("[3/6] extract apical domain");   human = step_extract(dom, apo)
    print("[4/6] mouse alignment");         pairs = step_mouse(start, end, human)
    print("[5/6] avoid list");              cat = step_avoid(dom, apo, pairs)
    print("[6/6] candidate hotspots and views")
    step_hotspots(dom, apo, pairs, cat)
    step_view(dom, cat)
    print(f"\ndone -> {OUT}\nview it:  pymol {OUT / 'view.pml'}")


if __name__ == "__main__":
    main()
