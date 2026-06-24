test_that("duplicate term descriptions get stable display labels", {
    x <- mock_enrich_result()

    mapping <- enrichplot:::get_term_mapping(x)

    expect_equal(mapping$ID, c("T1", "T2"))
    expect_equal(mapping$label, c("dup [T1]", "dup [T2]"))
})

test_that("extract_geneSets keeps ID selection and carries display labels", {
    x <- mock_enrich_result()

    gene_sets <- enrichplot:::extract_geneSets(x, c("T2", "dup [T1]"))

    expect_equal(names(gene_sets), c("T2", "T1"))
    expect_equal(
        attr(gene_sets, "term_labels"),
        c(T2 = "dup [T2]", T1 = "dup [T1]")
    )
})

test_that("heatplot(showTop) fails early without foldChange", {
    x <- mock_enrich_result()

    expect_error(
        heatplot(x, showCategory = c("T1", "T2"), showTop = 1),
        "`showTop` requires `foldChange`."
    )
})

test_that("heatplot uses disambiguated labels for duplicate descriptions", {
    x <- mock_enrich_result()

    p <- heatplot(
        x,
        showCategory = c("T1", "T2"),
        showTop = 1,
        foldChange = mock_foldchange()
    )

    expect_s3_class(p, "ggplot")
    expect_setequal(unique(p$data$categoryID), c("dup [T1]", "dup [T2]"))
    expect_equal(
        anyDuplicated(ggplot2::ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()),
        0L
    )
})

test_that("pairwise_termsim uses stable labels when descriptions repeat", {
    x <- mock_enrich_result()

    y <- pairwise_termsim(x, method = "JC", showCategory = c("T1", "T2"))

    expect_equal(rownames(y@termsim), c("dup [T1]", "dup [T2]"))
    expect_equal(colnames(y@termsim), c("dup [T1]", "dup [T2]"))
})

test_that("get_enrichplot_color expands two custom colors to three safely", {
    old <- getOption("enrichplot.colours")
    on.exit(options(enrichplot.colours = old), add = TRUE)

    options(enrichplot.colours = c("#111111", "#222222"))

    expect_equal(
        enrichplot:::get_enrichplot_color(3),
        c("#111111", "white", "#222222")
    )
})

test_that("cnetplot smoke test works for compareClusterResult", {
    x <- mock_comparecluster_result()

    p <- cnetplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
})

test_that("emapplot smoke test works for compareClusterResult", {
    x <- mock_comparecluster_result()
    x <- pairwise_termsim(x, method = "JC", showCategory = 2)

    p <- emapplot(x, showCategory = 2)

    expect_s3_class(p, "ggplot")
})
