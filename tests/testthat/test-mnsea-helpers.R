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
    expect_true("pathway" %in% subnet$nodes$node_type)
    expect_true(any(subnet$edges$edge_type %in% c("membership", "intra", "coupling")))
})

test_that("mnsea helpers use a stable default pathway when pathway_id is NULL", {
    x <- mock_mnsea_result()

    feature_df <- fortify_mnsea_contribution(x, level = "feature")
    subnet <- fortify_mnsea_subnetwork(x)

    expect_true(all(feature_df$ID == "T1"))
    expect_true(all(subnet$nodes$pathway_id == "T1"))
    expect_true(all(subnet$pathway$ID == "T1"))
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

test_that("cnetplot uses the default mnsea pathway when pathway_id is NULL", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, node_label = "category")

    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$pathway_id == "T1"))
})

test_that("heatplot uses the default mnsea pathway for feature view", {
    x <- mock_mnsea_result()

    p <- heatplot(x, showTop = 2, value = "score")

    expect_s3_class(p, "ggplot")
    expect_true(all(p$data$ID == "T1"))
    expect_lte(length(unique(as.character(p$data$Feature))), 2)
})

test_that("heatplot rejects incompatible mnseaResult value semantics", {
    x <- mock_mnsea_result()

    expect_error(
        heatplot(x, pathway_id = "T1", value = "share"),
        "When `pathway_id` is provided, `value` must be `score` or `abs_score`."
    )
})

test_that("cnetplot works for mnseaResult subnetworks", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, node_label = "category")

    expect_s3_class(p, "ggplot")
    expect_true("pathway" %in% p$data$node_type)
})

test_that("cnetplot supports share and exclusive mnsea labels", {
    x <- mock_mnsea_result()

    p_share <- cnetplot(x, pathway_id = "T1", node_label = "share")
    p_exclusive <- cnetplot(x, pathway_id = "T1", node_label = "exclusive")

    share_labels <- unique(unlist(lapply(p_share$layers, function(layer) layer$data$label)))
    exclusive_labels <- unique(unlist(lapply(p_exclusive$layers, function(layer) layer$data$label)))

    expect_true("g1" %in% share_labels || "g2" %in% share_labels)
    expect_true(length(exclusive_labels) == 0 || "g1" %in% exclusive_labels || "g2" %in% exclusive_labels)
})

test_that("cnetplot default mnsea labels keep pathway labels and deduplicate features", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "all")
    label_data <- p$layers[[4]]$data

    expect_true("Pathway 1" %in% label_data$label)
    feature_labels <- label_data$label[label_data$node_type == "feature"]
    expect_equal(anyDuplicated(feature_labels), 0L)
})

test_that("cnetplot item labels use one representative label per feature", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "item")
    item_labels <- p$layers[[4]]$data$label

    expect_equal(anyDuplicated(item_labels), 0L)
})

test_that("cnetplot encodes mnsea edge semantics and honors size_edge", {
    x <- mock_mnsea_result()

    p_small <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, size_edge = 0.5)
    p_large <- cnetplot(x, pathway_id = "T1", include_couplings = TRUE, size_edge = 2)
    p_no_coupling <- cnetplot(x, pathway_id = "T1", include_couplings = FALSE, size_edge = 2)

    edge_small <- ggplot2::ggplot_build(p_small)$data[[1]]
    edge_large <- ggplot2::ggplot_build(p_large)$data[[1]]
    edge_no_coupling <- ggplot2::ggplot_build(p_no_coupling)$data[[1]]

    expect_gt(mean(edge_large$linewidth, na.rm = TRUE), mean(edge_small$linewidth, na.rm = TRUE))
    expect_gte(length(unique(edge_large$linetype)), 2)
    expect_equal(length(unique(edge_no_coupling$linetype)), 1)
})

test_that("cnetplot encodes mnsea feature sign with node colors", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")
    feature_nodes <- ggplot2::ggplot_build(p)$data[[2]]

    expect_gte(length(unique(feature_nodes$colour)), 2)
})

test_that("cnetplot distinguishes mnsea pathway and feature nodes with shapes", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")
    built <- ggplot2::ggplot_build(p)$data

    expect_true(all(built[[2]]$shape == 21))
    expect_true(all(built[[3]]$shape == 23))
})

test_that("cnetplot uses explicit mnsea legend titles", {
    x <- mock_mnsea_result()

    p <- cnetplot(x, pathway_id = "T1", node_label = "none")

    expect_equal(p$scales$get_scales("linetype")$name, "Edge type")
    expect_equal(p$scales$get_scales("shape")$name, "Node type")
    expect_equal(p$scales$get_scales("fill")$name, "Layer")
    expect_equal(p$scales$get_scales("colour")$name, "Feature sign")
    expect_equal(p$scales$get_scales("size")$name, "Feature magnitude")
})
