# TfR1 Binder Design

De novo design of a pH-switchable protein binder for human transferrin receptor 1, for the [Anthropic × Adaptyv 2026 competition, challenge 04](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1).

**Deadline:** Sunday 1 November 2026, 23:59 AoE.

---

## The challenge

Design a small protein, from scratch and entirely on a computer, that binds **TfR1** — a receptor on the surface of brain blood vessels whose day job is pulling iron into the cell. A laboratory then synthesises the top designs and measures whether they work.

Why anyone wants this: the brain is walled off from the bloodstream, and a typical therapeutic antibody reaches roughly a third of one per cent of its blood concentration in brain tissue. TfR1 already carries cargo across that barrier, so a binder that grabs it can smuggle a drug through. Several drugs on the market already work this way.

The binder has to satisfy three requirements at once.

| # | Requirement | Why it is there |
|---|---|---|
| ① | Bind human TfR1 at pH 7.4, K<sub>D</sub> under 1 µM | pH 7.4 is blood. This is the basic job. |
| ② | Also bind mouse TfR1 at pH 7.4 | Drugs are tested in mice before humans. A human-only binder cannot be tested at all. |
| ③ | **Release** human TfR1 at pH 6.0 | pH 6.0 means we are inside the cell. Letting go there frees the cargo to continue into the brain instead of being carried back out. |

Requirement ③ is the hard one. Only **histidine** changes its charge between pH 7.4 and pH 6.0, so it is the only switch chemistry offers — and the generative models we design with have no representation of pH at all.

## The goal

Submit **twenty designs** that bind both species tightly, carrying a deliberate and documented spread of histidine placements.

That is the success criterion even if none of them switch, because it is the version of this experiment that produces an interpretable result either way. We are engineering requirements ① and ②, and making an informed bet on ③.

---

## Architecture

Six stages. Generation and filtering run on hosted endpoints; the pH switch and the final selection are done locally, by hand.

```mermaid
flowchart TD
    S1["1 · Target prep<br/>truncate TfR1 to the apical domain"]
    S2["2 · Generation<br/>RFdiffusion → ProteinMPNN → Boltz-2<br/>thousands of candidates"]
    S3["3 · Affinity filter<br/>ipSAE rescore<br/>hundreds remain"]
    S4["4 · Cross-species gate<br/>re-fold against mouse TfR1<br/>50–100 survive"]
    S5["5 · Install the pH switch<br/>graft histidines into the interface"]
    S6["6 · Choose twenty<br/>a portfolio, not a ranking"]
    OUT(["20 designs submitted"])

    S1 --> S2 --> S3 --> S4 --> S5
    S5 --> RA["Route A · by hand<br/>baseline, always available"]
    S5 --> RB["Route B · Baker his_ph_interface<br/>PRIORITY"]
    S5 --> RC["Route C · Proton-PottsMPNN<br/>wild card"]
    RA --> S6
    RB --> S6
    RC --> S6
    S6 --> OUT

    classDef pri fill:#dbeeeb,stroke:#0b766c,stroke-width:3px,color:#14181d
    classDef wild fill:#fae7d9,stroke:#b4460d,stroke-dasharray:5 3,color:#14181d
    classDef switch fill:#fae7d9,stroke:#b4460d,stroke-width:2px,color:#14181d
    class RB pri
    class RC wild
    class S5 switch
```

Every stage exists because of a specific biological fact:

| Stage | The biology that demands it |
|---|---|
| **1 · Target prep** | TfR1 has one patch we may bind (apical) and one we must not (the transferrin site, residues 619–670 — block it and you cause anaemia). Cut away everything else. |
| **2 · Generation** | Proteins bind by surfaces that fit. Make thousands of shapes and see which ones do. Binder length is capped at 65–90 residues so the interface stays breakable. |
| **3 · Affinity filter** | Requirement ① sets a floor of 1 µM. ipSAE scores the interface we designed rather than the target we already understood. |
| **4 · Cross-species gate** | Human and mouse share only 3 of 9 key positions on the apical patch. Keep what grips those three. Expect to lose half to three quarters of the pool. |
| **5 · pH switch** | Histidine is the only residue that changes charge between pH 7.4 and 6.0, and burying a charge where no water reaches it breaks an interface. |
| **6 · Choose twenty** | Requirement ③ can only be scored approximately, so spread the bet rather than ranking on it. Includes 2 unmodified controls. |

### Stage 5 routes

The pH switch is the least certain part of the project, so it has three independent routes with a clear order of preference.

| Route | What it is | Standing |
|---|---|---|
| **A · By hand** | Rank positions by burial and adjacency, mutate, re-fold, repair breakages with ArcRefine | Baseline — zero install risk, happens regardless |
| **B · [`baker-laboratory/ph_binder_paper`](https://github.com/baker-laboratory/ph_binder_paper)** | Published, wet-lab-validated pipeline that places histidines **and scores pH sensitivity directly** | **Priority** |
| **C · [`christian-creator/ProtonPottsMPNN`](https://github.com/christian-creator/ProtonPottsMPNN)** | Inverse-folding model with protonation states in its alphabet | Wild card — never on the critical path |

Routes B and C are alternatives to *stage 5 only*. Both take the same input stage 4 already produces — a binder-on-target complex, binder as chain A, target as chain B — which is what makes them safe to attempt in parallel.

---

## Documents

| Document | What it covers |
|---|---|
| [The TfR1 Problem](docs/tfr1-biological-doc.html) &nbsp;·&nbsp; [md](docs/tfr1-biological-doc.md) | The biology from zero. Assumes no background. Start here if you are new. |
| [The TfR1 Pipeline](docs/tfr1-architecture.html) &nbsp;·&nbsp; [md](docs/tfr1-architecture.md) | The six stages in detail, with the biology behind each one. The working plan. |

The HTML versions are the originals and are better to read — they carry diagrams, an interactive pH simulator and clickable glossary terms. The Markdown versions exist for diffing, grepping and review.

## Contributing

The useful background is not biology but pipeline engineering, GPU orchestration and careful evaluation design. If you can run batch jobs on rented GPUs, reason about overfitting to a scoring function, or build the bookkeeping that tracks thousands of candidates down to twenty, you can contribute from day one.

Work is delegated only when it is off the critical path, has a defined input/output contract, and has a do-nothing fallback that is the plan we already have. Route B above is the clearest example.
