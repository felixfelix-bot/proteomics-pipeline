# Lydia's STRING Network Approach vs. Our Pipeline: Deep-Dive Comparison

**Author:** Hermes subagent (delegated analysis)
**Date:** 2026-07-06
**Scope:** Compare `R/05_string_network.R`, `R/16_string_network_targeted.R`,
`R/17_crac_string_network.R`, `R/19_gsea.R`, `R/20_pathway_network.R`.
**Repo:** `~/proteomics-pipeline` (TRIP4/ASCC proteomics study, Dr. Aruna)

---

## TL;DR — The Core Confusion, Resolved

There are **two fundamentally different STRING workflows** living in this repo,
and they are easy to conflate:

| Approach | Gene input to STRING | Where | Matches Lydia? |
|----------|----------------------|-------|----------------|
| **A. Direct PPI of significant proteins** | The significant proteins themselves (from limma/DE results) | `R/05`, `R/16`, `R/17`, `R/19` (Part C) | **YES — this is Lydia** |
| **B. PPI of GO-pathway genes** | Genes belonging to enriched GO terms (after `enrichGO`) | `R/20` (and the GO half of `R/16`) | **NO — this is *our* invention** |

**Dr. Aruna's clarification is correct and decisive:**
> Lydia takes significant proteins → checks physical STRING interactions between
> them directly. It's NOT GO pathways → STRING.

`R/20_pathway_network.R` does GO enrichment **first**, then takes the pathway
gene lists and maps those to STRING. That is the opposite order and answers a
different question. `R/05`, `R/16`, `R/17`, and the network half of `R/19` all
do the Lydia-style direct mapping. So **the Lydia approach already exists in
the repo** — it just lives in scripts other than `R/20`.

---

## 1. What Lydia's Approach Actually Does (`R/05_string_network.R`)

This file is the faithful adaptation of Lydia's original script
**`20260521_stats_cutoffs_TurboID_mapToStringDBnetwork.R`** (referenced in the
header comment at `R/05` line 6, but the original `.R` file is **not present in
the repo** — only the adapted copy). The filename itself tells you the whole
philosophy: *"stats cutoffs → TurboID → map to STRING DB network"*.

### Step-by-step (matching the numbered STEP comments in the file)

| Step | What happens | Code anchor (`R/05`) |
|------|--------------|----------------------|
| **Init** | Connect to STRINGdb v12.0, human (taxon 9606), score ≥ 400 (medium confidence) | lines 67–72 |
| **1** | **Map the raw DE result data frame** to STRING IDs via `string_db$map(df, "gene", ...)`. Input = the experiment's *full* results table (gene, log2FC, padj). No GO step. | lines 97–112 |
| **2** | Define **"seed"** proteins using Lydia's two-pronged high-confidence rule: (a) `log2FC > 2 AND -log10(padj) > 6` **OR** (b) `log2FC > 7 AND -log10(padj) > 2`. Need ≥2 seeds. | lines 125–144 |
| **3** | **Network expansion (1-hop):** for every seed, call `string_db$get_neighbors(seed)`. Union seeds + neighbors → `expanded_ids`. | lines 153–161 |
| **4** | Pull all STRING interactions among the expanded set via `get_interactions(expanded_ids)`. | lines 163–167 |
| **5** | **Filter to "core" edges** — keep only interactions where at least one endpoint is a seed. This focuses the graph on the seeds' local neighborhood. | lines 169–181 |
| **6** | Build an **undirected `igraph`** from the core edge list. | lines 183–194 |
| **7** | Annotate vertices with `gene_name`, `log2FC`, `padj`, and an `is_seed` flag from the mapped table. | lines 196–211 |
| **8** | For each node, count how many of its neighbors are seeds → `core_links` (a connectivity-to-seed score). | lines 213–223 |
| **9** | Flag `inNetwork` membership back onto the original data frame (used for volcano overlay). | lines 225–235 |
| **10** | Build a **candidate interactor table**: non-seed nodes in the network, ranked by `core_links` desc, then `|log2FC|` desc. These are "proteins we may have missed." | lines 237–254 |
| **11** | **Plot the network** (Fruchterman–Reingold, seeds red & big, candidates teal & small, label seeds + high-connectivity nodes). | lines 256–295 |
| **12** | **Volcano plot with STRING overlay**: color significant vs. `inNetwork` vs. `high` seed. | lines 297–344 |

### What this is, conceptually

```
DE results (all proteins) ──► STRINGdb$map() ──► pick SEEDS (strict cutoffs)
                                                       │
                                                       ▼
                                              get_neighbors(seeds)  ◄── 1-hop expansion
                                                       │
                                                       ▼
                                         get_interactions(seeds ∪ neighbors)
                                                       │
                                                       ▼
                                       keep edges touching ≥1 seed  (core subgraph)
                                                       │
                                                       ▼
                                              igraph + plot + volcano overlay
```

**The defining properties of Lydia's approach:**
1. **Input = significant/DE proteins directly** — never GO-term gene sets.
2. **Seed definition is statistical** (log2FC + p-adj cutoffs), not ontology-based.
3. **Network is grown from the seeds outward by 1 hop** — the goal is to
   *discover* physical interactors of the confident hits, including proteins
   that were **not** detected in the mass-spec experiment.
4. **Output is a physical PPI network map** + a ranked candidate list +
   volcano overlay. The biological question is *"do my top hits sit in a
   connected physical interaction web, and who else is in that web?"*

---

## 2. Script-by-Script Comparison

### 2.1 `R/16_string_network_targeted.R` — Lydia-style, **closest match**

This is an almost line-for-line refactor of `R/05` ("Adapted from Lydia's
`05_string_network.R`", lines 16 & 44). It keeps Lydia's exact pipeline:
map full DE → seeds (same two-pronged cutoff) → 1-hop neighbor expansion →
core-edge filter → igraph → candidate ranking.

**What it adds beyond Lydia:**
- A **network-membership classification** of significant genes into three
  buckets: `seed_high_confidence`, `in_network`, `not_in_network` (lines 129–147).
  This split is then fed into downstream GO enrichment **as a *secondary*
  analysis** (lines 236–326), i.e. GO is run *on the network output*, not
  before it.
- Git-hashed, sanitized filenames; PDF + PNG dual output.
- A pinned experiment list (`TRIP4_TurboID_vs_WT`).

**Verdict:** `R/16` **is Lydia's approach**, with extra reporting glued on the
back. The STRING-network portion is faithful. The GO half is a follow-on that
does **not** feed back into the network — so it does not violate Lydia's logic.

### 2.2 `R/17_crac_string_network.R` — Lydia-style for the CRAC dataset

Same Lydia logic, applied to the FLAG-TRIP4 CRAC RNA-interactome CSV
(different column names: `external_gene_name`, `logFC`, `FDR`).

**Differences from `R/05`/`R/16`:**
- No 1-hop neighbor expansion. It maps the **significant CRAC genes themselves**
  and pulls `get_interactions(mapped_ids)` directly (lines 142–150). The network
  contains **only detected proteins** — no STRING-discovered "missing" neighbors.
- Seeds defined with the same Lydia cutoffs but expressed as
  `padj < 1e-6 & log2FC > 2` OR `log2FC > 7` (lines 86–94) — logically identical
  to `R/05` (since `-log10(1e-6) ≈ 6`), just written differently.
- Output: network PNG+PDF + degree-ranked candidate table merged with log2FC/padj.

**Verdict:** Pure Lydia-style direct PPI, slightly *simpler* than `R/05` because
it skips neighbor expansion. Still "significant proteins → STRING."

### 2.3 `R/19_gsea.R` — Two-part: GSEA (Lydia-adjacent) + STRING maps (Lydia-style)

This file has **three parts**:
- **Part A** (`run_gsea`): GSEA via `gseGO` on ranked gene lists
  (`sign(log2FC) * -log10(padj)`), all three ontologies, all experiments.
  This is *ranked-list enrichment* — Lydia-flavored analytics, not a network.
- **Part B**: GSEA on the CRAC dataset (same helper).
- **Part C** (`build_string_network`): **this is the Lydia network** —
  `get_significant_genes(df)` → STRING `$map()` → `$get_interactions()` → igraph
  → plot with seed coloring (Lydia cutoffs at lines 198–199). Capped at 200 sig
  genes for readability.

**Verdict:** Part C **is** Lydia's direct-mapping approach (no GO-first step).
The GSEA parts are separate enrichment products that ship in the same script.
Note Part C **does not do the 1-hop neighbor expansion** that `R/05` does — it
only shows interactions among the detected significant proteins themselves,
like `R/17`.

### 2.4 `R/20_pathway_network.R` — **THE OFFENDING ONE (GO-first)**

This is the script that **does not match Lydia**. Flow:

| Step | Code (lines) |
|------|--------------|
| For each experiment, get significant genes via `get_significant_genes` | 219 |
| Run `enrichGO` for BP/MF/CC with `universe` background | 232–246 |
| Take **top 3 enriched GO terms** by p.adjust | 257 |
| For each top term, `strsplit(geneID, "/")` → that GO term's gene list | 262 |
| **Now** map *those GO-term genes* to STRING, get interactions, plot | `build_go_network`, 56–213 |

**This is the inverse of Lydia.** The genes entering STRING are *members of an
enriched Gene Ontology term*, not the significant proteins themselves. The
biological question becomes *"is there a PPI network among the genes annotated
to GO:000XXXX?"* — which is a pathway-coherence check, not a
"do-my-hits-connect-physically" check.

Why it drifts:
1. STRING membership is determined by GO annotation, so you can only ever see
   proteins that GO *already* grouped together — the network becomes a
   tautological restatement of the GO result, not new information.
2. STRING-discovered "missing" interactors (the whole point of Lydia's neighbor
   expansion) never appear, because the input is a closed GO gene set.
3. Seed definitions (Lydia's log2FC/p-adj cutoffs) are absent — coloring is by
   logFC direction of the GO-term genes instead.

### 2.5 Side-by-side matrix

| Property | `R/05` (Lydia) | `R/16` | `R/17` | `R/19` Part C | `R/20` |
|----------|:--:|:--:|:--:|:--:|:--:|
| Input to STRING = significant proteins directly | ✅ | ✅ | ✅ | ✅ | ❌ (GO genes) |
| Lydia seed cutoffs (logFC>2&-log10p>6 ‖ logFC>7) | ✅ | ✅ | ✅ | ✅ | ❌ |
| 1-hop neighbor expansion (`get_neighbors`) | ✅ | ✅ | ❌ | ❌ | ❌ |
| Core-edge filter (keep edges touching a seed) | ✅ | ✅ | ❌ | ❌ | ❌ |
| GO enrichment runs **before** STRING | ❌ | ❌ | ❌ | ❌ | ✅ |
| GO enrichment runs **after** STRING (secondary) | ❌ | ✅ | ❌ | (separate part) | ❌ |
| Vertex color = seed vs candidate | ✅ | ✅ | ✅ | ✅ | ❌ (up/down logFC) |
| Candidate-interactor ranking by `core_links` | ✅ | ✅ | ❌ (degree) | ❌ (degree) | ❌ |
| Volcano overlay with network membership | ✅ | ❌ | ❌ | ❌ | ❌ |

**Bottom line:** `R/05` and `R/16` are the most faithful Lydia implementations
(full pipeline incl. neighbor expansion). `R/17` and `R/19` Part C are Lydia in
spirit but simplified (no neighbor expansion — only show interactions among
detected proteins). `R/20` is **not Lydia at all**.

---

## 3. Where They Diverge — and Why It Matters

### 3.1 The order of operations is reversed in `R/20`

```
Lydia (R/05/16/17/19C):   Significant proteins ─► STRING PPI map
Our R/20:                  Significant proteins ─► enrichGO ─► GO genes ─► STRING PPI map
```

In `R/20`, the GO step is an upstream **filter and re-labeler**. By the time
proteins reach STRING, they are no longer "the significant hits" — they are
"the members of enriched GO term #N". This means:

- A significant protein that does **not** fall in a top-3 GO term never enters
  the network, even if it has strong STRING edges to other hits.
- A non-significant protein that happens to be annotated to a top GO term
  **does** enter the network (STRING `$map` only needs a recognized symbol).
  → The network can contain proteins that were *not* significant in the
  experiment.

Lydia's approach has neither problem: the STRING node set is exactly the
significant hits (plus discovered neighbors), gated by the seed cutoffs.

### 3.2 `R/20` loses the "candidate interactor discovery" feature

The single most useful thing Lydia's pipeline does is **find proteins you
didn't detect** that are physically connected to the proteins you did detect
(via `get_neighbors`). That is impossible in `R/20` because the input is a
closed GO gene set with no expansion step.

### 3.3 `R/20` answers a different scientific question

- Lydia (`R/05`): *"Do my high-confidence TRIP4-proximal proteins form a
  physical interaction network, and what other proteins sit in that network?"*
- `R/20`: *"For each enriched GO term, do the annotated genes have STRING
  edges among them?"* — i.e. a pathway-internal-coherence visualization.

Both are legitimate analyses, but only the first is "Lydia's approach," and
only the first answers the question Dr. Aruna actually asked for.

### 3.4 `R/17` and `R/19` Part C are simplified Lydia (still correct, less rich)

They skip the 1-hop neighbor expansion and the core-edge filter, so the network
shows **only edges among detected significant proteins**. If two significant
proteins are bridged only through an undetected third protein, that bridge will
not appear. For sparse datasets (CRAC), this can make a network look emptier
than Lydia's version would. Worth noting but **not wrong** — they still map
significant proteins directly, which is the essence of Lydia's method.

---

## 4. What Needs to Change to Match Lydia *Exactly*

If the goal is "produce Lydia-style STRING maps for every experiment,"
`R/20_pathway_network.R` should either be **deprecated/relabelled** or
**rewritten** to follow the `R/05`/`R/16` recipe. Concretely:

### Option A (recommended): promote `R/16` to the canonical network script
`R/16_string_network_targeted.R` already implements Lydia faithfully. To make
it the single source of truth:
1. Expand `STRING_EXPERIMENTS` (line 213) from one entry to the full target
   list (the `TARGET_NAMES` vector in `R/19` lines 286–310 is a good template).
2. Keep the secondary GO-by-membership block (lines 236–326) but clearly
   document it as a **follow-on**, not part of the network construction.
3. **Retire or rename `R/20`** to something like `R/20_go_term_coherence.R`
   and update the Makefile target (`make pathway-network`) and README so users
   don't mistake it for the Lydia network.

### Option B: rewrite `R/20` to Lydia logic
Replace `analyze_experiment_pathways()` with the `R/05` pipeline:
1. Drop the `enrichGO` call (lines 232–246).
2. Drop the "top 3 GO terms" loop (lines 256–278).
3. Replace `build_go_network()` with `R/05`'s `analyze_string_network()`:
   map the full DE df → define seeds (Lydia cutoffs) → `get_neighbors` →
   core-edge filter → igraph → candidate ranking.
4. Re-add the volcano overlay from `R/05` Step 12.
5. Optionally keep a *separate* GO section that runs **after** and uses the
   network's `in_network` gene set as input (mirroring `R/16` lines 257–326).

### Concrete checklist of behavioral changes for `R/20`
- [ ] Remove `enrichGO` from the network path; STRING input becomes
      `get_significant_genes(df)` (or the full `df` for `$map`).
- [ ] Re-introduce Lydia seed cutoffs:
      `(log2FC > 2 & -log10(padj) > 6) | (log2FC > 7 & -log10(padj) > 2)`.
- [ ] Add `string_db$get_neighbors(seed_string_ids)` and union with seeds.
- [ ] Add the core-edge filter
      `int$from %in% seed_ids | int$to %in% seed_ids`.
- [ ] Switch node coloring from up/down-logFC to seed-vs-candidate.
- [ ] Add the candidate-interactor table ranked by `core_links`.
- [ ] Add the volcano overlay with `inNetwork` / `high` categories.

### What to *keep* from the current `R/20`
- The nice up/down logFC node coloring and size-by-significance styling is
  genuinely useful **as an optional annotation layer** on the Lydia network —
  it can be applied to the seed-expanded graph without changing the logic.
- The PDF+PNG dual output and git-hashed filenames.

---

## 5. Repo-Wide Lydia Footprint

Searches for `Lydia`/`lydia` hit **39 lines across 13 files**. Summary of
where Lydia's influence is claimed:

| File | Lydia connection |
|------|------------------|
| `R/05_string_network.R` | Direct adaptation of `20260521_stats_cutoffs_TurboID_mapToStringDBnetwork.R` |
| `R/16_string_network_targeted.R` | Adapted from `R/05` (Lydia pipeline) |
| `R/17_crac_string_network.R` | Uses "Lydia's criteria" for seeds |
| `R/19_gsea.R` | Labelled "(all experiments — Lydia style)"; network part uses Lydia seed cutoffs |
| `R/02_volcano_plots.R` | Lydia-style multi-category volcano (`plot_lydia_volcano`) |
| `R/06_gene_families.R` | Based on Lydia's `20260424_Overlap_TurboID_wOthers.R` gene-family highlighting |
| `R/07_overlap_analysis.R` | Based on Lydia's `20260424_Overlap_TurboID_wOthers.R` |
| `R/12_go_network_volcano.R` | "Lydia-style" GO-network volcano (color, no text labels) |
| `R/01_config.R` | `LOG2FC_CUTOFF = 0.5` "matching Lydia's scripts"; legacy color scheme |
| `README.md`, `Makefile`, `run_all.R`, `docs/ARCHITECTURE.md` | Documentation references |

**Two original Lydia script filenames are cited but NOT in the repo:**
1. `20260521_stats_cutoffs_TurboID_mapToStringDBnetwork.R` — the STRING-network
   original (referenced only in `R/05` header).
2. `20260424_Overlap_TurboID_wOthers.R` — the overlap/gene-family original
   (referenced in `R/06`, `R/07`).

If the originals are needed for ground-truth verification, they would need to
be obtained from Lydia directly; they are not recoverable from this repo.

---

## 6. Recommendations (prioritized)

1. **Immediate (clarity):** Update `R/20`'s header and the `make pathway-network`
   help text to state plainly that it is *"STRING networks of GO-term genes,
   not Lydia-style direct-hit networks."* This alone prevents the
   misunderstanding recurring. (~5 min.)
2. **Short term (correctness for Dr. Aruna's request):** Point Dr. Aruna at
   `R/16` (and `R/19` Part C) for the Lydia-style maps she wants. If she needs
   maps for *all* experiments, expand `R/16`'s `STRING_EXPERIMENTS` list
   (Option A above). (~30 min.)
3. **Medium term (consolidation):** Either retire `R/20` or rewrite it to
   Lydia logic (Option B). Keep one canonical Lydia network script, not three
   half-overlapping ones (`R/05`, `R/16`, `R/19` Part C).
4. **Nice-to-have:** Add the 1-hop neighbor expansion back to `R/17` and
   `R/19` Part C so they recover Lydia's candidate-discovery feature.

---

## 7. One-Line Summary for the Parent Agent

`R/05`, `R/16`, `R/17`, and Part C of `R/19` **already implement Lydia's
approach** (significant proteins → STRING directly, with Lydia's seed cutoffs);
`R/20_pathway_network.R` is the **only** script that does GO-enrichment-first
and is therefore the one that does *not* match Lydia. Dr. Aruna's STRING
network request should be served by `R/16` (most faithful, full pipeline) or
`R/19` Part C (simpler, no neighbor expansion) — not by `R/20`.
