# Stage 1 — target preparation

What the target-prep code does, what it writes, and how to run it.

Stage 1 turns three public PDB entries and two UniProt sequences into the design target handed to stage 2: a single-chain apical-domain PDB, a span-matched mouse sequence, a list of residues a binder must avoid, and a ranked set of provisional hotspot candidates. Everything is derived; nothing here is experimentally validated. For *why* the recommended residues were chosen, with citations, see [the hotspot evidence report](task6-hotspot-analysis.md).

Status as of 28 October 2026: steps 1–6 computationally complete and reproducible. Binding, cross-species activity and pH-dependent release are untested.

## Quick start

From the repository root:

```sh
python3 -I scripts/prep_target.py                     # build the target (a few minutes)
python3 -I -m unittest discover -s tests -v           # 7 regression checks
pymol data/target/hotspots.pml                        # look at the candidate patches
python3 -I scripts/visualize_target.py --open         # render three views, open the session
```

Python standard library only — no pip install. The first run needs network access to download the public inputs; afterwards everything runs offline from the cache in `data/raw/`. The whole `data/` tree is git-ignored, so a fresh clone starts empty and rebuilds it.

Re-running is safe and idempotent: a re-run from the populated cache reproduces every output byte-for-byte.

## The three files

| File | Role |
|---|---|
| [scripts/prep_target.py](../scripts/prep_target.py) | The pipeline. Downloads inputs, carves out the target, aligns mouse, computes the avoid list and the hotspot candidates, writes PyMOL views. |
| [scripts/visualize_target.py](../scripts/visualize_target.py) | Rendering only. Drives an installed PyMOL to produce three images and an interactive session from what `prep_target.py` already wrote. |
| [tests/test_target_prep.py](../tests/test_target_prep.py) | Regression checks on the parts that were wrong before or are easy to get wrong silently. |

`visualize_target.py` and the tests both read the pipeline's outputs, so run `prep_target.py` first.

---

## What `prep_target.py` does

Six steps, printed as `[n/6]` as it goes.

### 1 · Fetch inputs

Caches three PDB entries in `data/raw/`:

| Entry | What it is | Used for |
|---|---|---|
| `1CX8` | Human TfR1 ectodomain, apo | The target itself. Chain C, biological assembly 2 (chains C/D). |
| `3S9L` | TfR1 bound to transferrin | The avoid list, and a second conformation for exposure screening. |
| `11QC` | Newer human TfR1 ectodomain dimer | A third conformation, and the only entry with resolved apical glycans. |

Also caches UniProt `P02786` (human) and `Q62351` (mouse), and the PDBe SCOP domain mapping for 1CX8.

Inputs are never silently refreshed. The SHA-256 of every cached input is recorded in `hotspots.json`, so results can be tied to the exact bytes they came from. To pick up an upstream revision, delete the cached file deliberately.

### 2 · Domain boundaries

Reads SCOP domain assignments for 1CX8 chain C from the PDBe mappings API rather than hard-coding them:

```
protease-like   122–189 and 383–608
apical / PA     190–382
helical         609–760
```

Asserts that the apical range does not overlap the transferrin docking site at 619–670. These boundaries are one database's definition; other sources draw them slightly differently.

### 3 · Extract the apical domain

Writes chain C residues 190–382 as ATOM records only, primary altloc only, keeping the original author numbering. 193 residues, 1469 atoms, no chain breaks, termini 5.3 Å apart — the apical domain is an insertion whose two ends emerge side by side, so the cut is clean.

Dropping heteroatoms is an input-format choice for the generator. It does **not** mean glycans are absent from the assay receptor; see step 6.

### 4 · Mouse alignment

Global Needleman–Wunsch of full-length human P02786 against mouse Q62351 (match +2, mismatch −1, linear gap −4), then the apical span is taken from the alignment. The human side is asserted to match the 1CX8 crop exactly, which also confirms that PDB numbering equals UniProt numbering.

Result: **144 identical columns out of 194** (74%), 48 substitutions, two gap columns.

Mouse positions come from the alignment, never from adding a constant offset. Where several equally scoring alignments exist, the code walks the whole optimal traceback graph instead of picking one path, and records every alternative. Human positions **196 and 210–212** are ambiguous under this scoring; ambiguous positions are disqualified as hotspot anchors. This is a sequence-level uncertainty check, not a mouse structural validation.

### 5 · The avoid list

Residues a binder must not bury, measured as 5 Å heavy-atom contacts in **reconstructed biological assemblies** of 3S9L — assembly 1 from chains A/C and assembly 2 from B/D, each with its own `REMARK 350` transform applied.

This is a correction to an earlier version that compared asymmetric-unit chains A and B. Those are not each other's physiological partners, so the old 21-residue list was measuring crystal packing. It is superseded.

| | Old (asymmetric unit) | Current (biological assemblies) |
|---|---:|---:|
| Contacts with the other TfR1 subunit | 14 | 11 |
| Contacts with transferrin | 10 | 0 |
| Union | 21 | **11** |

The 11 positions are **312, 313, 314, 316, 317, 318, 319, 320, 322, 323, 324**.

Zero transferrin contacts is a statement about a 5 Å cutoff, not about clearance. The apical domain's closest approach to transferrin in these assemblies is **6.4–6.8 Å** — just outside the cutoff. A binder with any reach can obstruct transferrin from an apical anchor. Check complete complexes.

### 6 · Hotspot candidates and views

Scores all 193 apical positions across **five receptor copies**: 1CX8 chain C, both 3S9L receptors in their own reconstructed assemblies, and both chains of the 11QC dimer. A position is eligible only if it passes every filter in every copy:

| Filter | Threshold |
|---|---|
| Identical in mouse, with an unambiguous alignment mapping | required |
| Side-chain SASA (Shrake–Rupley, 1.4 Å probe, 240 sphere points) | ≥ 10 Å² |
| Distance to protein atoms removed by the truncation | ≥ 5 Å |
| Cα distance to the glycosylation sites N251/N317 | ≥ 12 Å |
| Distance to resolved NAG atoms in 11QC | ≥ 5 Å |
| Cα distance to either artificial cut terminus | ≥ 8 Å |
| Membership in the step-5 avoid list | excluded |

30 of 193 positions survive. Triplets are then formed from eligible residues, keeping only those within a 12 Å maximum Cα diameter in every copy and containing at least one hydrophobic/aromatic residue (F/W/Y/I/L/V/M) and one positive residue (K/R). The chemistry rule is a design hypothesis aimed at later histidine work, not a general rule for good binders.

Ranking is explicitly heuristic — the sum of each residue's minimum side-chain SASA (capped at 100 Å²) weighted by local sequence identity within 10 Å. **It is not an affinity estimate.** The top three triplets sharing at most one residue are kept:

| Rank | Human | Mouse | Score | Max Cα diameter |
|---|---|---|---:|---:|
| 1 | **F237 / K240 / K241** | F239 / K242 / K243 | 89.87 | 9.45 Å |
| 2 | V291 / A293 / R339 | V293 / A295 / R341 | 84.49 | 10.50 Å |
| 3 | K240 / D242 / I262 | K242 / D244 / I264 | 69.93 | 11.84 Å |

Rank 3 overlaps rank 1, so these are not three independent epitopes.

---

## Outputs

Everything lands in `data/target/` (git-ignored).

| File | Contents |
|---|---|
| `human_apical.pdb` | **The design target.** 1CX8 chain C, residues 190–382, original numbering. |
| `human_apical.fasta` | Matching human sequence, 193 aa. |
| `mouse_apical.fasta` | Span-matched mouse sequence. |
| `human_mouse_apical_alignment.txt` | Readable alignment, one optimal traceback. |
| `residue_mapping.tsv` | Human ↔ mouse numbering per position, with an `unique_optimal_mapping` flag. |
| `alignment_ambiguities.json` | Every alternative mapping among equally scoring alignments. |
| `avoid_residues.json` | The biological-assembly contact list, split into seam / transferrin / both. |
| `hotspot_candidates.tsv` | All 193 residues: per-copy SASA, clearances, local conservation, eligibility and the reason for each exclusion. |
| `hotspots.json` | The selected patches, all three numbering systems, every threshold, input hashes and stated limitations. |
| `domains.json` | Domain boundaries, the transferrin site, glycosylation sites and chain termini — the metadata the viewer draws from. |
| `view.pml` | Domain colouring on the full chain, plus the avoid list on the isolated domain. |
| `hotspots.pml` | The three candidate triplets, coloured and labelled. |

`hotspot_candidates.tsv` is the file to open when asking "why was position *n* rejected?" — the `exclusion_reasons` column names every failed filter.

### Handing the target to the generator

`hotspots.json` carries the same patch in three numbering systems because downstream tools disagree about which one they want:

```text
'ppi.hotspot_res=[C237,C240,C241]'
```

That is the RFdiffusion argument fragment for `human_apical.pdb` as written — chain C, original human numbering. It is a fragment, not a launch command. If a hosted wrapper renumbers the crop to 1–193 the positions become **48, 51, 52**; if it renames the chain, remap that too. The mouse numbers are **239, 242, 243**. Check the chain and numbering the generator actually received before launching a large run, and run each candidate triplet separately.

---

## Visualizing it

`scripts/visualize_target.py` renders what the pipeline produced. It needs an installed PyMOL executable but no Python packages; PyMOL does the rendering in its own interpreter.

```sh
python3 -I scripts/visualize_target.py --open                          # render five views and open
python3 -I scripts/visualize_target.py                                 # render only
python3 -I scripts/visualize_target.py --pymol /opt/homebrew/bin/pymol # explicit PyMOL
pymol data/visualization/tfr1.pse                                      # reopen a prior result
```

No network needed. Everything it draws comes from `domains.json`, `avoid_residues.json` and `hotspots.json`, so the views follow the pipeline rather than hard-coded residue numbers. Before rendering it checks that the patches are three residues each on chain C, that the apical object has exactly 193 Cα atoms, that the receptor contains only chains C and D, that resolved sugars survived the cleanup, and that each selected residue appears exactly once per object. Failure is reported against `data/visualization/render.log` rather than a half-written image.

### The five scenes

Switch with PyMOL's scene buttons, or `scene 01_overview, recall`. Each scene displays a one-line caption on recall explaining what you are looking at.

| Scene | Shows |
|---|---|
| `01_overview` | The whole ectodomain dimer, chain C coloured by domain. |
| `02_domains` | Chain C alone: three domains, termini, transferrin site. |
| `03_constraints` | The apical domain in context, with everything a binder must avoid. |
| `04_target` | `human_apical.pdb` exactly as the generator receives it. |
| `05_hotspots` | All three candidate patches, close up. |

**Colour key, consistent across every scene:**

| Colour | Meaning |
|---|---|
| Grey | Protease-like domain (122–189, 383–608). |
| Teal | Apical domain (190–382) — the design target. |
| Slate | Helical domain (609–760). |
| Orange cartoon | Transferrin docking site, 619–670. Blocking it stops iron delivery body-wide. |
| Translucent blue | The partner subunit, chain D. |
| Magenta sticks | Apical residues buried against the partner subunit — the avoid list. |
| Yellow | Glycosylation sites N251/N317/N727 and their resolved NAG sugars. |
| Gold / orange / pink | The three residues of a candidate patch, in order. |

`01_overview` is 1CX8 biological assembly 2 — chains C and D in their deposited positions, no extra transform. It is the extracellular part only: no membrane, no intracellular tail. The samarium ions are dropped as a crystallization additive, but the N-linked sugars are kept, because they constrain where a binder can actually reach.

`03_constraints` is the scene worth studying before launching a generation run. It puts the three things that disqualify a surface — the dimer seam, the glycans and the rest of the receptor — next to the patch being proposed.

`05_hotspots` draws all three patches from `hotspots.json`; scenes 01–04 carry only the primary one, so the labels stay readable. Sphere size is set for legibility and is not a molecular surface. Colours identify residues and encode nothing measured. These are provisional anchors, and they are not histidines — histidine substitutions go on the designed binder later.

Outputs in `data/visualization/`: one 1400 × 1000 PNG per scene, `tfr1.pse` (the session with all five scenes and the named selections), `render_target.py` (the generated PyMOL driver, including the exact input paths), and `render.log` / `viewer.log`.

The session also carries named selections — `apical_domain`, `helical`, `protease_like`, `transferrin_site`, `dimer_seam`, `glycans`, `glycosites`, `partner_chain`, `patch_1`–`patch_3` — so you can toggle any of them from the right-hand panel or use them in your own commands.

Re-rendering overwrites only these files; it never touches the target or the hotspot rankings. If a re-run of `prep_target.py` changes the patches, re-render — labels and captions follow the metadata automatically — and update the residue names quoted in this document.

Because `data/` is git-ignored, these images exist only on the machine that generated them. They are not committed documentation assets.

### Troubleshooting

- **Missing input** — run `python3 -I scripts/prep_target.py` first.
- **PyMOL not found** — install it, or pass `--pymol /path/to/pymol`.
- **Rendering failed** — the real error is in `data/visualization/render.log`.
- **Images appear but no window does** — check `viewer.log`, then open `tfr1.pse` directly. The GUI needs a desktop display.

---

## The tests

```sh
python3 -I -m unittest discover -s tests -v
```

Seven checks, all passing, run offline from the cache. They cover the failure modes that would otherwise be silent:

- Both 3S9L biological-assembly transforms, against the literal coordinates in `REMARK 350`.
- The crystal-packing bug itself: residue 211 contacts asymmetric-unit chain B but *not* the real assembly partner. This test fails if anyone reintroduces the old comparison.
- SASA against an analytic single-atom sphere area, and that a neighbouring chain reduces it.
- That the alignment code exposes *all* optimal mappings, not one traceback.
- That F237/K240/K241 really are identical and uniquely mapped in both UniProt sequences.

---

## What this does not establish

The pipeline is reproducible and its geometry is checked. None of that makes the candidates correct.

- **The hotspots are a hypothesis.** They are structure-derived candidates, not measured energetic hotspots, and not a recovery of the competition's unnamed three conserved positions — no complete nine-position set was found in the sources checked.
- **The thresholds are permissive.** A stricter trial (20 Å² exposure, 8 Å clearance) produced no qualifying triplet at all, so the recommendation depends on the current settings. Raising SASA sampling from 240 to 960 points changed nothing, which checks numerical sampling, not biology.
- **The primary patch is tight.** F237 is only modestly exposed (14.09 Å² minimum), and all three anchors sit 6–7 Å from protein atoms the truncation removed. A binder approaching from the wrong direction clashes while still contacting valid anchors.
- **Conservation is sequence-only.** No mouse structure was predicted. Cross-species binding depends on the whole interface, not three shared positions.
- **Glycans are only buffered, not modelled.** 11QC resolves NAG at N251 and N317 in the apical domain; the filters keep a margin around them but no glycan conformation is simulated.
- **The apical domain is not unoccupied.** Ferritin is a natural ligand, and transferrin passes within 6.4–6.8 Å.
- **Nothing about pH has been tested.** Conserved K240/K241 are *possible* partners for binder-side histidines; that is not a pH switch. Histidines are placed on the binder after a complex exists, and the S142/G142 allele difference between these structures and the competition reference sits outside the apical crop but matters for full-ectodomain sequence models.

Every design that comes out of stage 2 must be checked against the full receptor dimer, the rest of its own chain, transferrin and the glycan-bearing regions before the pipeline scales up.

> `data/target/apical_avoid.png` predates the biological-assembly correction and depicts the superseded 21-residue list. It is not regenerated by the current code and can be deleted; use `view.pml` instead.
