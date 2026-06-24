#' Plot-facing helpers for `mnseaResult`
#'
#' These helpers convert explanation-ready outputs from `enrichit` into stable
#' data structures that are easier to consume in `enrichplot`.
#'
#' @param x A `mnseaResult` object.
#' @param level One of `"pathway"` or `"feature"`.
#' @param pathway_id Optional pathway ID.
#' @param include_couplings Logical, whether to keep inter-layer coupling edges.
#' @param include_isolated Logical, whether to keep nodes without retained edges.
#' @param ... Additional parameters reserved for future extensions.
#'
#' @return `fortify_mnsea_contribution()` returns a `data.frame`.
#' `fortify_mnsea_subnetwork()` returns a list containing standardized `nodes`
#' and `edges` tables together with the selected pathway metadata.
#' @importFrom enrichit extract_mnsea_subnetwork get_mnsea_contribution
#' @export
fortify_mnsea_contribution <- function(
    x,
    level = c("pathway", "feature"),
    pathway_id = NULL,
    ...
) {
    level <- match.arg(level)
    df <- get_mnsea_contribution(x, pathway_id = pathway_id, level = level)
    df <- as.data.frame(df, stringsAsFactors = FALSE)

    if (nrow(df) == 0) {
        return(df)
    }

    if (level == "pathway") {
        keep <- c("ID", "Description", "layer", "contribution", "share", "n_feature")
        df <- df[, intersect(keep, colnames(df)), drop = FALSE]
        df$Description <- as.character(df$Description)
        df$layer <- as.character(df$layer)
        df <- df[order(df$Description, -df$contribution, df$layer), , drop = FALSE]
        rownames(df) <- NULL
        return(df)
    }

    df$Description <- as.character(df$Description)
    df$layer <- as.character(df$layer)
    df$Feature <- as.character(df$Feature)
    df$sign <- ifelse(
        df$score < 0,
        "suppressed",
        ifelse(df$score > 0, "activated", "neutral")
    )
    df$pathway_id <- as.character(df$ID)
    df$pathway_description <- as.character(df$Description)
    df <- df[order(-df$abs_score, df$layer, df$Feature), , drop = FALSE]
    rownames(df) <- NULL
    df
}

#' @rdname fortify_mnsea_contribution
#' @export
fortify_mnsea_subnetwork <- function(
    x,
    pathway_id = NULL,
    include_couplings = TRUE,
    include_isolated = FALSE,
    ...
) {
    subnet <- extract_mnsea_subnetwork(
        x,
        pathway_id = pathway_id,
        include_couplings = include_couplings,
        include_isolated = include_isolated
    )

    if (!is.null(subnet$pathway) && nrow(subnet$pathway) > 0) {
        pathway_id <- as.character(subnet$pathway$ID[1])
        pathway_description <- as.character(subnet$pathway$Description[1])
    } else {
        pathway_id <- NA_character_
        pathway_description <- NA_character_
    }

    nodes <- as.data.frame(subnet$nodes, stringsAsFactors = FALSE)
    edges <- as.data.frame(subnet$edges, stringsAsFactors = FALSE)

    if (nrow(nodes) > 0) {
        nodes$layer <- as.character(nodes$layer)
        nodes$Feature <- as.character(nodes$Feature)
        nodes$node <- nodes$node_key
        nodes$label <- nodes$Feature
        nodes$node_type <- "feature"
        if (!"abs_score" %in% colnames(nodes) && "score" %in% colnames(nodes)) {
            nodes$abs_score <- abs(nodes$score)
        }
        nodes$sign <- ifelse(
            nodes$score < 0,
            "suppressed",
            ifelse(nodes$score > 0, "activated", "neutral")
        )
        nodes$pathway_id <- pathway_id
        nodes$pathway_description <- pathway_description
    }

    if (nrow(edges) > 0) {
        edges$edge_type <- as.character(edges$edge_type)
        edges$abs_weight <- abs(edges$weight)
    }

    subnet$nodes <- nodes
    subnet$edges <- edges
    subnet
}

.result_data <- function(x) {
    if (isS4(x) && "result" %in% methods::slotNames(x)) {
        return(as.data.frame(x@result, stringsAsFactors = FALSE))
    }
    as.data.frame(x, stringsAsFactors = FALSE)
}

.rank_mnsea_terms <- function(df, by = c("p.adjust", "NES", "contribution", "share")) {
    by <- match.arg(by)

    if (nrow(df) == 0) {
        return(df[, c("ID", "Description"), drop = FALSE])
    }

    if (by %in% c("contribution", "share")) {
        rank_df <- stats::aggregate(
            df[[by]],
            by = list(ID = df$ID, Description = df$Description),
            FUN = max,
            na.rm = TRUE
        )
        colnames(rank_df)[3] <- ".rank"
        rank_df <- rank_df[order(rank_df$.rank, decreasing = TRUE), , drop = FALSE]
        return(rank_df)
    }

    keep <- c("ID", "Description", by)
    rank_df <- unique(df[, keep, drop = FALSE])
    if (by == "NES") {
        rank_df$.rank <- abs(rank_df$NES)
        rank_df <- rank_df[order(rank_df$.rank, decreasing = TRUE), , drop = FALSE]
    } else {
        rank_df$.rank <- rank_df[[by]]
        rank_df <- rank_df[order(rank_df$.rank, decreasing = FALSE), , drop = FALSE]
    }
    rank_df
}
