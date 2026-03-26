# Manhattan Plot Implementation Plan

## Goal
Implement a `manhattanplot()` function in the `enrichplot` package to visualize enrichment results as a Manhattan plot, drawing conceptual inspiration from `gprofiler2::gostplot` (Figure 1 of f1000research paper 9-709). The plot will display enriched terms across the X-axis and their significance (-log10 P-value) on the Y-axis.

## Proposed Changes

### R/AllGenerics.R
- **[MODIFY] R/AllGenerics.R**
  Add the generic function `manhattanplot` to support S4 method dispatch.
  ```R
  #' Manhattan plot for enrichment result
  #'
  #' @title manhattanplot
  #' @rdname manhattanplot
  #' @param x enrichment result.
  #' @param ... additional parameters.
  #' @return ggplot object
  #' @export
  setGeneric("manhattanplot", function(x, ...) {
      standardGeneric("manhattanplot")
  })
  ```

### R/manhattanplot.R
- **[NEW] R/manhattanplot.R**
  Create a new S4 method implementation file to house `manhattanplot`.
  - Define methods for `enrichResult`, `gseaResult`, `compareClusterResult`, `enrichResultList`, and `gseaResultList` (similar to `dotplot.R`).
  - **Visualization Logic**: 
    - Use `fortify()` to convert the `x` enrichment object into a `df` data frame.
    - **X-axis**: Enriched terms (`Description`), ordered or grouped by category (e.g. `ONTOLOGY`, `Cluster`, or semantic similarity if `pairwise_termsim` was used). By default, they will follow the base `fortify()` order or `orderBy` parameter. Note: standard Manhattan plots space the points evenly or group by semantic similarities.
    - **Y-axis**: Computed significance, representing `-log10(p.adjust)` (or `-log10(pvalue)` if specified).
    - **Size**: Point size maps to the gene `Count`.
    - **Color**: Point color maps to the grouping variable (e.g., `ONTOLOGY`, dataset `Cluster`, or default to `p.adjust` gradient if un-grouped).
    - **Styling**: Integrate standard `enrichplot` configurations, including `theme_dose(font.size)` and `enrichplot_point_shape`. Incorporate an optional horizontal dashed line referencing significant thresholds (e.g. `yintercept = -log10(0.05)`).
    - Return the resulting `ggplot` object.

## Verification Plan

### Manual Verification
Verification will be performed via manual inspection with standard DOSE/clusterProfiler data.

A script (`.dev/test_manhattan.R`) will be composed to iteratively test:
1. **Single Enrichment**: `DOSE::enrichDO()` to ensure un-grouped Manhattan plotting handles axes and coloring correctly.
2. **Compared Clusters**: `clusterProfiler::compareCluster()` to verify grouped coloring, facet behaviors (if implemented), and spacing.
3. Check visual properties:
   - Make sure X-axis labels are readable or omitted effectively (Manhattan plots often omit strict x-axis text in favor of broad grouping or interactive labels).
   - Confirm Y-axis correctly shows the log scale of significance.
   - Point sizes scale predictably by `Count`.
