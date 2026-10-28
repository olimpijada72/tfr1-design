# Task 6: provisional TfR1 hotspot selection

28 October 2026. Stage 1, task 6 — not pipeline stage 6 (final design selection).

**Decision:** use **human F237, K240 and K241** as the first *structure-derived candidate patch* for a pilot generation run. They correspond to **mouse F239, K242 and K243**. This finishes the computational candidate-selection deliverable; it does not establish binding affinity, mouse cross-reactivity, or pH-dependent release.

These are **not** a recovered list of the competition's “three conserved positions.” No complete, unambiguous nine-position list was located in the sources checked. The task has therefore been completed using the structural-candidate route described in [the target-prep guide](target-prep.md).

## What the research changes

The [competition page](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1) does say that a particular apical patch differs at six of nine positions. Its detailed description names six substitutions plus a deletion, but does not enumerate the complete nine-position set or the conserved three. That statement cannot justify treating the entire apical domain as having only three conserved binding opportunities. Its recommended epitope is a preference, not a requirement.

The project's actual alignment contains **144 identical columns out of 194 (74.2%)**, covering 193 residues of each species, with 48 substitutions and two gap columns. The mouse positions must be obtained by alignment, not by blindly adding an offset. The script now records all possible mappings on equally scoring optimal alignment paths. Human positions 196 and 210–212 are ambiguous under its scoring scheme; none of the recommended anchors is ambiguous. The page's Y211→P214/deletion-at-L212 convention and the old script's Y211→D213/deletion-at-V210 convention are different optimal gap placements, not evidence of different mouse sequences.

Two different choices had been conflated:

1. **Target hotspots:** TfR1 residues used to guide where the generated binder attaches. TfR1 remains unchanged.
2. **Histidine placements:** residues in the subsequently designed **binder** that may be replaced or redesigned to create a pH-dependent interaction.

Published computational work found that simply increasing histidine content often gave little pH sensitivity. Productive mechanisms included disrupting acceptor interactions and creating repulsion near positive residues, as well as destabilizing buried networks within the binder. Their effects depend on the environment and on the free and bound states. Much of that evidence compares pH 7.4 with 5.4; it does not establish release at this project's pH 6.0. See [Computational design of pH-sensitive binders](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621740/).

TfR-specific experimental work also supports histidine engineering as a possibility, but concerns an existing mouse-binding nanobody and measurements at pH 5.5. Its successful substitutions cannot be transferred by residue number into an unrelated de novo binder. See [the mouse TfR nanobody study](https://link.springer.com/article/10.1186/s12987-023-00462-z).

## The structural bug that had to be fixed first

The old script compared chain A against chain B in the asymmetric unit of **3S9L**. Those are not each other's physiological dimer partners. The PDB's `REMARK 350` specifies:

- Assembly 1: chains A (TfR1) and C (transferrin), plus transformed copies of both. The second operation is `(x,y,z) → (231.676−x,y,−z)`.
- Assembly 2: chains B and D, plus transformed copies. The second operation is `(x,y,z) → (y,x,−126.71625−z)`.

Reconstructing both assemblies changes the 5 Å heavy-atom contact result:

| Analysis | Apical dimer contacts | Apical transferrin contacts | Union |
|---|---:|---:|---:|
| Old asymmetric-unit comparison | 14 | 10 | 21 |
| Correct biological assemblies | 11 | 0 | 11 |

The corrected dimer-contact positions are **312, 313, 314, 316, 317, 318, 319, 320, 322, 323, 324**. The old list is superseded. Zero direct apical–transferrin contacts at this cutoff does **not** establish that a binder attached there cannot obstruct transferrin: the apical domain's closest approach to transferrin in these assemblies is 6.4–6.8 Å, just outside the cutoff.

The assembly definitions are in the downloaded [3S9L record](https://www.rcsb.org/structure/3S9L); the transformed-coordinate calculations are reproducible in `prep_target.py`. The [original structural paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3156180/) also documents conformational changes around the apical arm, supporting the use of both apo and transferrin-bound structures rather than one conformation.

The newer competition reference, [11QC](https://www.rcsb.org/structure/11QC), supplies an additional apo dimer. Its coordinates and LINK records show NAG attached at **N251 and N317** in the apical domain, and N727 elsewhere. Removing HETATM records from the generator input does not make these glycans disappear from the assay receptor. Full glycan conformations remain unmodelled. Also, the apical domain is not biologically unoccupied: ferritin is a natural ligand, as discussed in [this TfR binder study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12626770/).

## How the candidates were derived

The script evaluates all 193 apical positions, using five human receptor copies:

- 1CX8 chain C, in biological assembly 2 (C/D).
- 3S9L receptor A and receptor B, each in its own reconstructed biological assembly.
- Both receptor chains of 11QC's biological dimer.

Each position must have an identical, uniquely mapped mouse residue. It must have at least **10 Å² of exposed side-chain surface in every copy**, computed with a standard Shrake–Rupley geometric calculation using a 1.4 Å solvent probe and 240 sphere points per atom. All protein atoms in the assembly occlude the surface. This is an exposure calculation, not a binding-energy calculation.

Additional filters require:

- At least 5 Å from any protein atoms omitted by apical truncation, across all five contexts. This includes other receptor domains, the dimer partner and transferrin where present.
- At least 12 Å between the residue's Cα and the N251/N317 Cα positions in 1CX8 C, and at least 5 Å from resolved NAG atoms in either 11QC receptor copy. These are screening buffers, not a complete glycan envelope.
- At least 8 Å Cα distance from both artificial cut termini in 1CX8 C.
- No membership in the corrected contact avoid list.

Triplets must contain a hydrophobic/aromatic candidate (F/W/Y/I/L/V/M) and a positive candidate (K/R), and have maximum pairwise Cα distance at most 12 Å in every context. This chemistry preference deliberately explores interfaces potentially useful for subsequent histidine design; it is not a universal rule for good binders.

Ranking is explicitly heuristic: sum each residue's minimum side-chain SASA (capped at 100 Å²), weighted by sequence identity among apical neighbours within 10 Å Cα distance in 1CX8 C. It does not estimate KD or energetic hotspot strength. The output selects up to three triplets sharing at most one residue; they need not be three spatially independent epitopes.

## First candidate patch

| Human residue | Mouse residue | Minimum side-chain SASA, Å² | Minimum distance to omitted protein, Å | Local identity within 10 Å |
|---|---|---:|---:|---:|
| F237 | F239 | 14.09 | 6.96 | 95.8% |
| K240 | K242 | 39.25 | 6.47 | 91.3% |
| K241 | K243 | 45.29 | 7.17 | 89.5% |

Maximum Cα diameter is **9.45 Å** across the five copies. All three sites are conserved without alignment ambiguity. F237 offers a partially exposed aromatic side chain, while K240/K241 provide potential positive partners for future binder-side histidines. Neighbouring acidic residues, including D242 and E244, could instead stabilize a protonated histidine; the eventual geometry matters.

**Material limitations:** F237 is only modestly exposed, and all three anchors are within about 6–7 Å of omitted protein atoms in at least one context. A binder approaching from the wrong direction can therefore clash despite contacting valid anchors. Glycans are another footprint-level risk. The patch is a new design hypothesis, not an experimentally demonstrated functional epitope.

A stricter trial requiring 20 Å² side-chain exposure and 8 Å clearance produced no triplet satisfying the chemistry and diameter conditions. Thus the recommendation depends on reasonable but permissive screening thresholds. Increasing SASA sampling from 240 to 960 points retained the same first patch and the same two alternatives; this checks numerical sampling, not biological validity.

The two secondary candidate sets are **V291/A293/R339** (mouse V293/A295/R341) and **K240/D242/I262** (mouse K242/D244/I264). The first has a more divergent neighbourhood and only 5.65 Å clearance at V291. The second overlaps the primary region and includes an acidic anchor. Both are lower-priority hypotheses, not demonstrated pH-switch sites.

## Generator handoff and remaining validation

For the current `data/target/human_apical.pdb`, original human residue numbers and **chain C** are retained. The upstream RFdiffusion hotspot override for the primary candidate is:

```text
'ppi.hotspot_res=[C237,C240,C241]'
```

This follows [RFdiffusion's input-PDB chain/residue convention](https://github.com/RosettaCommons/RFdiffusion#binder-design). It is an argument fragment, not a full launch command. A hosted wrapper may require a different schema. If it renumbers the crop to 1–193, the positions become **48, 51, 52**; if it renames the target chain, remap that too. The mouse positions are **239, 242, 243**, not the human numbers.

Run each candidate set separately. Start with a pilot before allocating thousands of generations. For each resulting complex, inspect the full binder footprint against the full receptor dimer, the remainder of the same receptor chain, transferrin and glycan-bearing regions. Compare its complete interface against mouse, and reject reliance on species-specific neighbouring contacts. The old structures' S142/G142 allele difference is outside the cropped apical target, but matters when preparing full-ectodomain sequence models.

Only after a plausible cross-species complex exists should binder histidines be placed and evaluated in the neutral and protonated environments. Retain neutral-pH binding controls. Structural confidence scores do not establish affinity or a pH 6.0 release response; those outcomes remain experimental questions.

## Reproduction

[The target-prep guide](target-prep.md) documents the scripts, every output file and how to run them. In short:

```sh
python3 -I scripts/prep_target.py
python3 -I -m unittest discover -s tests -v
pymol data/target/hotspots.pml
```

Python standard library only, offline from the cache after the first run, with input SHA-256 hashes recorded in `hotspots.json`. Seven regression tests cover both assembly transformations, the crystal-contact failure mode, analytical SASA/occlusion, alignment ties and the recommended human/mouse identities.

No binder generation, energetic hotspot measurement, mouse structure prediction or pH validation has been performed in task 6.
