# nsea / mnsea Visualization Implementation Checklist

## Goal

Turn the design spec in `nsea-mnsea-visualization-methods-spec.md` into an execution-ready checklist for `enrichplot`.

This checklist has two top-level goals:

1. **Complete the inherited `gseaResult`-style plotting workflow for `nseaResult` and `mnseaResult`.**
2. **Implement a first batch of new mechanism-oriented methods that are specific to network enrichment results.**

## Working Rule

- **First complete the old workflow.**
- **Then add new methods.**
- **Prefer helper-first implementation over plot-first patching.**
- **Each batch should end with tests, NEWS, and a clean commit boundary.**

## Batch 0: Baseline Confirmation

### Goal

Confirm what is already done and what remains missing, so later batches do not re-open closed work.

### Checklist

- [ ] Reconfirm current `mnseaResult` support in:
  - `dotplot()`
  - `heatplot()`
  - `cnetplot()`
  - `emapplot()`
  - `gseaplot()`
  - `ridgeplot()`
  - `upsetplot()`
  - `ssplot()`
- [ ] Reconfirm existing helper contracts:
  - `fortify.mnseaResult()`
  - `fortify_mnsea_contribution()`
  - `fortify_mnsea_subnetwork()`
- [ ] Reconfirm current gaps:
  - `gseaplot2()`
  - `gsearank()`
  - `hplot()`
  - `treeplot()`
  - `barplot()` compatibility confirmation

### Verification

- [ ] `tests/testthat/test-mnsea-helpers.R` remains green
- [ ] No regression in existing `mnsea` helper outputs

## Batch 1: Complete the Running-Score Family

### Goal

Finish the plots that should naturally follow from `gseaplot.mnseaResult()`.

### Target Functions

- [ ] `gseaplot2()` for `nseaResult`
- [ ] `gseaplot2()` for `mnseaResult`
- [ ] `gsearank()` for `nseaResult`
- [ ] `gsearank()` for `mnseaResult`
- [ ] `hplot()` for `nseaResult`
- [ ] `hplot()` for `mnseaResult`

### Implementation Notes

- [ ] Reuse `gsInfo.mnseaResult()` where possible
- [ ] Do not invent a second ranked-score pipeline
- [ ] For `mnseaResult`, keep the same score-space rule:
  - `layer = NULL` -> `collapsed_scores`
  - `layer = "<single-layer>"` -> `layer_scores[[layer]]`
- [ ] Reject ambiguous multi-layer running-score requests unless the plot truly supports them

### Likely Files

- [ ] `R/gseaplot.R`
- [ ] `man/gseaplot.Rd`
- [ ] `tests/testthat/test-mnsea-helpers.R`

### Tests

- [ ] `gseaplot2()` smoke test for `nseaResult`
- [ ] `gseaplot2()` smoke test for `mnseaResult`
- [ ] `gseaplot2()` stable `geneSetID` resolution test
- [ ] `gseaplot2()` single-layer score-space test
- [ ] `gsearank()` smoke tests for `nsea` and `mnsea`
- [ ] `hplot()` smoke tests for `nsea` and `mnsea`
- [ ] boundary test for invalid `layer`

### Exit Criteria

- [ ] Running-score family methods are usable end-to-end for both `nseaResult` and `mnseaResult`
- [ ] No duplicated score-selection logic remains

## Batch 2: Complete the Similarity / Tree Family

### Goal

Finish the plots that depend on pathway similarity and term clustering.

### Target Functions

- [ ] `treeplot()` for `nseaResult`
- [ ] `treeplot()` for `mnseaResult`
- [ ] `barplot()` compatibility confirmation for `nseaResult`
- [ ] `barplot()` compatibility confirmation for `mnseaResult`

### Implementation Notes

- [ ] Reuse current `emapplot` / `ssplot` similarity logic where possible
- [ ] Do not create a second, incompatible similarity definition for `mnsea`
- [ ] If needed, stabilize an internal helper for `mnsea` pairwise similarity before wiring `treeplot()`

### Likely Files

- [ ] `R/treeplot.R`
- [ ] `R/barplot.R` or tests only
- [ ] `man/treeplot.Rd`
- [ ] `man/barplot.Rd` if signature or docs change
- [ ] `tests/testthat/test-mnsea-helpers.R`

### Tests

- [ ] `treeplot()` smoke test for `nsea`
- [ ] `treeplot()` smoke test for `mnsea`
- [ ] `treeplot()` no-precomputed-termsim path for `mnsea`
- [ ] single-pathway / two-pathway boundary tests
- [ ] `barplot()` compatibility checks for `nsea` and `mnsea`

### Exit Criteria

- [ ] Existing `gseaResult`-style workflow is effectively complete for `nseaResult` / `mnseaResult`
- [ ] Remaining gaps are no longer “missing old plots”, only “new method work”

## Batch 3: Stabilize Mechanism-Oriented Helper Layer

### Goal

Create the reusable summary helpers needed by new methods.

### Target Helpers

- [ ] `compute_rewiring_score()`
- [ ] `classify_mechanism_state()`
- [ ] `summarize_nsea_mechanism()`
- [ ] `extract_rewiring_features()`

### Minimal Contracts

#### `compute_rewiring_score()`

- [ ] Accept `nseaResult` or `mnseaResult`
- [ ] Return one numeric rewiring score per pathway
- [ ] Use a transparent rule, not an opaque heuristic blob

#### `classify_mechanism_state()`

- [ ] Use enrichment shift + rewiring score
- [ ] Return one of:
  - `conserved`
  - `rewired`
  - `context_specific`
  - `contradictory`

#### `summarize_nsea_mechanism()`

- [ ] Return a term-level summary table with at least:
  - `ID`
  - `Description`
  - `NES`
  - `p.adjust`
  - `leading_edge_size`
  - `leading_edge_overlap`
  - `rewiring_score`
  - `centrality_shift`
  - `mechanism_class`

#### `extract_rewiring_features()`

- [ ] Return a feature-level comparison table with at least:
  - `Feature`
  - `score`
  - `abs_score`
  - `sign`
  - `status`

### Likely Files

- [ ] `R/nsea-mechanism-helpers.R` or equivalent new helper file
- [ ] roxygen docs for helper contracts
- [ ] `tests/testthat/test-nsea-mechanism-helpers.R`

### Tests

- [ ] contract tests for each helper
- [ ] deterministic classification tests
- [ ] empty / single-term / single-layer boundary tests
- [ ] invalid-layer and invalid-pathway tests

### Exit Criteria

- [ ] All new mechanism plots can consume the same helper layer
- [ ] Helper outputs are explicit enough to test without plotting

## Batch 4: Implement `phaseplot()`

### Goal

Deliver the first new method with clear `nsea` / `mnsea` identity.

### Method Definition

- [ ] X-axis = enrichment shift
- [ ] Y-axis = rewiring score
- [ ] Size = leading-edge size or overlap
- [ ] Color = significance or mechanism class

### Work Items

- [ ] Add generic to `R/AllGenerics.R`
- [ ] Create `R/phaseplot.R`
- [ ] Define methods for:
  - `nseaResult`
  - `mnseaResult`
- [ ] Reuse `summarize_nsea_mechanism()`
- [ ] Make default axis labels explicit and readable

### Likely Files

- [ ] `R/AllGenerics.R`
- [ ] `R/phaseplot.R`
- [ ] `man/phaseplot.Rd`
- [ ] `tests/testthat/test-phaseplot.R`

### Tests

- [ ] smoke test for `nseaResult`
- [ ] smoke test for `mnseaResult`
- [ ] mechanism class mapping test
- [ ] size/color semantic tests
- [ ] empty result and one-term boundary tests

### Exit Criteria

- [ ] `phaseplot()` provides information not already available from `dotplot()` / `emapplot()`
- [ ] The default plot already separates conserved vs rewired patterns in a readable way

## Batch 5: Implement `rewireplot()`

### Goal

Deliver the first pathway-specific mechanism explanation plot.

### Method Definition

- [ ] Nodes = leading-edge or pathway-driving features
- [ ] Edges = pathway-specific subnetwork
- [ ] Feature status displayed as:
  - `shared`
  - `gained`
  - `lost`
  - `shifted`
- [ ] Optional display of coupling-mediated edges for `mnsea`

### Work Items

- [ ] Add generic to `R/AllGenerics.R`
- [ ] Create `R/rewireplot.R`
- [ ] Reuse:
  - `fortify_mnsea_subnetwork()`
  - `extract_rewiring_features()`
- [ ] Define stable pathway selection rules
- [ ] Decide whether to facet by layer or color by layer for the first version

### Likely Files

- [ ] `R/AllGenerics.R`
- [ ] `R/rewireplot.R`
- [ ] `man/rewireplot.Rd`
- [ ] `tests/testthat/test-rewireplot.R`

### Tests

- [ ] smoke test for `mnseaResult`
- [ ] stable `pathway_id` resolution test
- [ ] feature-status mapping test
- [ ] no-edge / no-coupling boundary tests

### Exit Criteria

- [ ] `rewireplot()` answers “same pathway name, same mechanism or not?”
- [ ] pathway-specific rewiring evidence is readable without extra manual preprocessing

## Batch 6: Implement `consensusmap()`

### Goal

Add a multi-network / multi-layer overview plot for mechanism agreement and disagreement.

### Work Items

- [ ] Add generic and method file
- [ ] Build a term × context summary matrix
- [ ] Show enrichment strength and topology consistency together
- [ ] Attach one mechanism class per pathway

### Dependencies

- [ ] `summarize_nsea_mechanism()`
- [ ] `classify_mechanism_state()`

### Exit Criteria

- [ ] Users can quickly identify conserved, rewired, and context-specific pathways

## Batch 7: Implement `mechanismflow()`

### Goal

Add an evolution-style view for pathway state transitions across conditions or layers.

### Work Items

- [ ] Add generic and method file
- [ ] Define state transitions between contexts
- [ ] Choose a first rendering strategy:
  - sankey
  - river
- [ ] Keep the first version focused on `mnseaResult` if needed

### Exit Criteria

- [ ] The plot makes pathway state transitions easier to read than side-by-side `NES` comparisons

## Shared Verification Rules

### For Every Batch

- [ ] Add or update focused tests
- [ ] Run targeted verification first
- [ ] Update `NEWS.md`
- [ ] Regenerate docs if roxygen changes
- [ ] Keep commit boundaries narrow and readable

### Test Priority

- [ ] contract tests before visual tests
- [ ] helper tests before plot tests
- [ ] edge-case tests before “pretty plot” checks

## Suggested Commit Boundaries

### Commit Group A: Complete Old Plots

- [ ] `feat: add nsea gseaplot2 family support`
- [ ] `feat: add mnsea gseaplot2 family support`
- [ ] `feat: add nsea mnsea treeplot support`
- [ ] `test: confirm nsea mnsea barplot compatibility`

### Commit Group B: Add New Helper Layer

- [ ] `feat: add nsea mechanism helper summaries`

### Commit Group C: Add New Methods

- [ ] `feat: add phaseplot for nsea and mnsea`
- [ ] `feat: add rewireplot for mnsea`
- [ ] `feat: add consensusmap for network mechanism overview`
- [ ] `feat: add mechanismflow for pathway state transitions`

## Recommended Immediate Next Step

If implementation resumes now, the next coding batch should be:

1. **`gseaplot2()`**
2. **`gsearank()` / `hplot()`**
3. **then `treeplot()`**

Only after that should the work move to:

4. **`phaseplot()`**
5. **`rewireplot()`**

This keeps the development path clean:

- first complete compatibility
- then establish new method identity
