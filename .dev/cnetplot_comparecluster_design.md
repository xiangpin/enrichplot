# compareCluster cnetplot Design Draft

## Goal
Update `enrichplot::cnetplot.compareClusterResult()` so that `compareCluster` pies behave consistently with the newer `ggtangle::cnetplot()` sizing model, while keeping the current visual style and minimizing API breakage.

The immediate goals are:
- Make `categorySizeBy` actually affect category pie size for `compareClusterResult`.
- Clarify how pie colors should be customized.
- Align `enrichplot` documentation with the current `ggtangle` semantics.

## Current State

### Data Flow
- `cnetplot.enrichResult()` and `cnetplot.gseaResult()` convert the enrichment object to gene sets, attach numeric term-level columns as attributes, and pass `categorySizeBy` through to `ggtangle::cnetplot()`.
- `cnetplot.compareClusterResult()` takes a different route:
  - `tidy_compareCluster()` produces a long data frame.
  - Gene sets are reconstructed from `Description`.
  - `ggtangle::cnetplot()` is called with `size_category = 0`, `size_item = 0`, and `node_label = "none"` to generate only the network layout.
  - `add_node_pie()` overlays category and gene nodes with `scatterpie::geom_scatterpie()`.

### Consequences
- The visible category nodes in `compareClusterResult` are not the base `ggtangle` points. They are pies added later in `add_node_pie()`.
- `categorySizeBy` is currently declared in the `compareClusterResult` method signature, but it is not passed to `ggtangle::cnetplot()` and it is not used inside `add_node_pie()`.
- Category pie radius is currently derived from summed `Count` values per term, scaled by `size_category`.
- Pie fill colors are currently determined by the `Cluster` columns passed to `scatterpie`, not by `color_category`.

## Problem Statement
The current implementation has three mismatches.

### 1. API Mismatch
`categorySizeBy` is exposed for `compareClusterResult`, but the argument has no effect on the category pie radius.

### 2. Semantic Mismatch
`ggtangle::cnetplot()` now supports expression/formula-based category sizing such as `~itemNum` and `~ -log10(p.adjust)`, but the `compareClusterResult` path still uses a hard-coded size computation based on total `Count`.

### 3. Documentation Mismatch
The `enrichplot` docs still describe `categorySizeBy` using the older wording ("itemNum", "pvalue", "p.adjust", "qvalue" or numeric vector), which no longer reflects the current `ggtangle` interface.

## Design Principles
- Preserve the current `compareCluster` visual style based on pies.
- Keep the change local to `enrichplot`; do not require new behavior from `ggtangle`.
- Reuse `ggtangle` sizing semantics as much as practical.
- Prefer minimal API expansion.
- Avoid silently changing gene-node pie behavior unless explicitly intended.

## Proposed Changes

### 1. Separate Layout Generation from Pie Radius Computation
Keep the current two-step approach:
- Step 1: use `ggtangle::cnetplot()` only to compute layout and edges.
- Step 2: use `add_node_pie()` to draw visible nodes.

No change is needed to the basic rendering architecture.

### 2. Introduce a Term-Level Size Evaluation Step
Add an internal helper that computes one numeric size value per category term before pies are drawn.

Suggested helper:
- `compute_comparecluster_category_size(d, categorySizeBy)`

Responsibilities:
- Start from the `tidy_compareCluster()` output.
- Build one row per `Description`.
- Preserve numeric term-level columns needed for evaluation, such as:
  - `Count` summary
  - `pvalue`
  - `p.adjust`
  - `qvalue`
  - any other numeric columns that are constant within a term
- Add an explicit `itemNum` column representing the number of unique genes in each term.
- Evaluate `categorySizeBy` against this term-level data.

### 3. Adopt Formula/Expression Semantics for `compareClusterResult`
Normalize `categorySizeBy` in `cnetplot.compareClusterResult()` to match `ggtangle::cnetplot.list()` behavior:
- Default to `~itemNum`.
- Accept bare expressions such as `p.adjust`.
- Accept formulas such as `~ -log10(p.adjust)`.
- Accept scalar numeric values, which should recycle to all categories.

This avoids having `compareClusterResult` behave differently from `enrichResult`.

### 4. Use Evaluated Size Values as Pie Radius Inputs
Modify `add_node_pie()` to accept precomputed category size values.

Suggested signature:
- `add_node_pie(p, d, pie = "equal", category_scale = 1, item_scale = 1, category_size = NULL)`

Behavior:
- If `category_size` is `NULL`, preserve existing `Count`-based behavior for backward compatibility.
- If `category_size` is provided, map it to category `Description` and use it as the base value for pie radius.

Important detail:
- The current code multiplies by `category_scale` twice: once when computing `dd$pathway_size`, and once again in `aes(r = pathway_size * category_scale)`.
- The redesign should make scaling single-source and explicit.

Recommended rule:
- Compute a normalized base radius from `category_size`.
- Apply `category_scale` exactly once when converting the normalized value to plotted radius.

### 5. Keep Gene Pie Size Behavior Unchanged
For this change, keep gene-node pie size fixed as it is today.

Reason:
- The issue is specifically about category node sizing.
- Changing item-node sizing in the same patch would expand scope and complicate validation.

### 6. Clarify Color Customization Strategy
Do not overload `color_category` to mean pie slice fill colors.

Recommended behavior:
- Keep `color_category` semantics unchanged for non-pie category nodes.
- Document that `compareCluster` pie slice colors are controlled by the fill scale.
- In examples and docs, show customization via `scale_fill_manual()`.

Optional future enhancement:
- Add a new parameter such as `pie_colors = NULL` to apply a manual fill scale internally.

This should not be part of the minimal fix unless a stronger convenience API is desired.

## Implementation Sketch

### `cnetplot.compareClusterResult()`
- Change default `categorySizeBy` from `NULL` to `~itemNum`.
- Compute `category_size <- compute_comparecluster_category_size(d, categorySizeBy)`.
- Keep the existing `ggtangle::cnetplot()` call for layout generation.
- Pass `category_size` into `add_node_pie()`.

### `compute_comparecluster_category_size()`
- Build a term-level summary table keyed by `Description`.
- Add `itemNum` as unique gene count per term.
- Summarize `Count` as total term count across clusters.
- For columns like `p.adjust`, `pvalue`, `qvalue`, use one representative value per term if they are invariant within term.
- If a requested variable is not available, raise a clear error that names the missing field.
- Validate that the evaluated result is numeric, length 1 or `n_terms`, and contains no `NA`.

### `add_node_pie()`
- Use `category_size` instead of hard-coded `pathway_size` when supplied.
- Normalize category radii in a stable way so legends remain interpretable.
- Update the legend labeller accordingly.

## Open Design Choice

### How should pie legend labels behave when `categorySizeBy` is not `Count`?
There are two reasonable options.

Option A: keep a numeric radius legend only
- The legend shows the numeric values produced by `categorySizeBy`.
- Best for expressions like `-log10(p.adjust)`.
- Most faithful to the actual radius mapping.

Option B: keep the current gene-count-style legend when possible
- Use count labels only when the size source is effectively count-based.
- Switch to numeric-value labels for all other cases.

Recommendation:
- Use Option A as the general rule.
- If `categorySizeBy` is exactly `~itemNum` or equivalent count-based behavior, labels can remain count-like.

## Backward Compatibility
- Existing plots without `categorySizeBy` should continue to work.
- Existing visual style should remain close to current behavior when using the default.
- Existing code that adds `scale_fill_manual()` after `cnetplot(compareClusterResult)` should continue to work.

Potential visible changes:
- Default pie radii may shift slightly if the new default is implemented through normalized `itemNum` rather than the current summed `Count` rule.
- Pie legend labels may change if they are tied to the new size source.

## Documentation Changes

### `R/cnetplot.R` roxygen
- Update the `categorySizeBy` description for all relevant methods to match current expression/formula semantics.
- For `compareClusterResult`, explicitly state that the argument controls category pie size.
- Add examples using:
  - `categorySizeBy = ~itemNum`
  - `categorySizeBy = ~ -log10(p.adjust)`
  - `scale_fill_manual()` for cluster pie colors

### `man/cnetplot.Rd`
- Regenerate documentation after roxygen updates.

### NEWS
- Add an entry describing:
  - `categorySizeBy` now works for `compareClusterResult` category pies
  - pie color customization should use fill scales
  - improved consistency with `ggtangle::cnetplot()`

## Verification Plan

### Manual Verification
Create a script under `.dev`, for example:
- `.dev/test_cnetplot_comparecluster_size.R`

Check at least the following cases:
1. Default `compareCluster` plot with no explicit `categorySizeBy`
2. `categorySizeBy = ~itemNum`
3. `categorySizeBy = ~ -log10(p.adjust)`
4. `categorySizeBy = 2`
5. `scale_fill_manual()` on cluster pies
6. GSEA-style `compareCluster` result if available through `core_enrichment`

### Behavior Checks
- Category pie radii change when `categorySizeBy` changes.
- Gene pie radii remain unchanged.
- Edge layout is unaffected.
- Labels still align with category nodes.
- The pie legend remains readable and correctly reflects the size source.

### Error Handling Checks
- Unknown variable in `categorySizeBy` gives a clear error.
- Non-numeric evaluation gives a clear error.
- `NA` output from the size expression gives a clear error.

## Non-Goals
- Redesigning the overall `compareCluster` cnetplot appearance.
- Changing gene pie size semantics.
- Adding a new high-level palette API unless needed after the minimal fix.
- Refactoring `enrichResult` and `gseaResult` methods beyond documentation alignment.

## Recommended Implementation Order
1. Add the internal term-level size computation helper.
2. Thread computed `category_size` into `add_node_pie()`.
3. Fix radius scaling so `category_scale` is applied exactly once.
4. Update roxygen and examples.
5. Add a small manual test script under `.dev`.
