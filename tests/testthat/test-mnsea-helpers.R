test_that("fortify_mnsea_contribution standardizes pathway tables", {
    x <- mock_mnsea_result()

    df <- fortify_mnsea_contribution(x, level = "pathway")

    expect_true(all(
        c("ID", "Description", "layer", "contribution", "share", "n_feature") %in%
            colnames(df)
    ))
    expect_true(nrow(df) > 0)
})

test_that("fortify_mnsea_subnetwork standardizes nodes and edges", {
    x <- mock_mnsea_result()

    subnet <- fortify_mnsea_subnetwork(x, pathway_id = "T1")

    expect_true(all(
        c("node", "label", "node_type", "layer", "abs_score", "sign") %in%
            colnames(subnet$nodes)
    ))
    expect_true(all(
        c("from", "to", "edge_type", "abs_weight") %in%
            colnames(subnet$edges)
    ))
    expect_true(any(subnet$edges$edge_type %in% c("intra", "coupling")))
})

test_that("dotplot works for mnseaResult contribution view", {
    x <- mock_mnsea_result()

    p <- dotplot(
        x,
        x = "share",
        color = "contribution",
        showCategory = 2
    )

    expect_s3_class(p, "ggplot")
    expect_true(all(c("layer", "share", "contribution") %in% colnames(p$data)))
})

test_that("heatplot works for mnseaResult pathway contribution view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, showCategory = 2, value = "contribution")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Description", "layer", "fill_value") %in% colnames(p$data)))
    expect_setequal(unique(as.character(p$data$Description)), c("Pathway 1", "Pathway 2"))
})

test_that("heatplot works for mnseaResult pathway-specific feature view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, pathway_id = "T1", showTop = 2, value = "score")

    expect_s3_class(p, "ggplot")
    expect_true(all(c("Feature", "layer", "fill_value") %in% colnames(p$data)))
    expect_lte(length(unique(as.character(p$data$Feature))), 2)
    expect_true(all(p$data$ID == "T1"))
})

test_that("heatplot rejects incompatible mnseaResult value semantics", {
    x <- mock_mnsea_result()

    expect_error(
        heatplot(x, value = "score"),
        "When `pathway_id` is NULL, `value` must be `share` or `contribution`."
    )
    expect_error(
        heatplot(x, pathway_id = "T1", value = "share"),
        "When `pathway_id` is provided, `value` must be `score` or `abs_score`."
    )
})
