#' @rdname ssplot
#' @exportMethod ssplot
setMethod(
    "ssplot",
    signature(x = "enrichResult"),
    function(x, showCategory = 30, ...) {
        ssplot.enrichResult(x, showCategory = showCategory, ...)
    }
)

#' @rdname ssplot
#' @exportMethod ssplot
setMethod(
    "ssplot",
    signature(x = "gseaResult"),
    function(x, showCategory = 30, ...) {
        ssplot.enrichResult(x, showCategory = showCategory, ...)
    }
)

#' @rdname ssplot
#' @exportMethod ssplot
setMethod(
    "ssplot",
    signature(x = "compareClusterResult"),
    function(x, showCategory = 30, ...) {
        ssplot.compareClusterResult(x, showCategory = showCategory, ...)
    }
)


#' @rdname ssplot
#' @param drfun The function used for dimension reduction,
#' e.g. `stats::cmdscale` (the default), `vegan::metaMDS`, or `ape::pcoa`.
#' @param dr.params list, the parameters of `tidydr::dr`.
#' @param ... additional parameters
#'
#' additional parameters can refer the following parameters.
#'     \itemize{
#'       \item \code{layout} igraph layout function for node positioning
#'       \item \code{color} Variable that used to color enriched terms, e.g. 'pvalue','p.adjust' or 'qvalue'.
#'       \item \code{size_category} relative size of the categories
#'       \item \code{min_edge} The minimum similarity threshold for whether
#'         two nodes are connected, should between 0 and 1, default value is 0.2.
#'       \item \code{color_edge} color of the network edge
#'       \item \code{size_edge} relative size of edge width.
#'       \item \code{node_label} Select which labels to be displayed,
#'         one of 'category', 'group', 'all' and 'none'.
#'       \item \code{node_label_size} size of node label, default is 5.
#'       \item \code{pie} one of 'equal' or 'Count' to set the slice ratio of the pies (for `compareClusterResult` only).
#'       \item \code{label_format} a numeric value sets wrap length, alternatively a custom function to format axis labels.
#'       \item \code{clusterFunction} function of Clustering method, such as stats::kmeans(the default),
#'         cluster::clara, cluster::fanny or cluster::pam.
#'       \item \code{nWords} Numeric, the number of words in the cluster tags, the default value is 4.
#'       \item \code{nCluster} Numeric, the number of clusters,
#'         the default value is square root of the number of nodes.
#'     }
#'
#' additional parameters can refer the emapplot function: \link{emapplot}.
#' @importFrom tidydr theme_dr
ssplot.enrichResult <- function(
    x,
    showCategory = 30,
    drfun = NULL,
    dr.params = list(),
    #group = TRUE,
    node_label = "group",
    ...
) {
    if (is.null(drfun)) {
        drfun = stats::cmdscale
        dr.params = list(eig = TRUE)
    }
    if (is.character(drfun)) {
        drfun <- eval(parse(text = drfun))
    }

    drResult <- get_drResult(
        x = x,
        showCategory = showCategory,
        drfun = drfun,
        dr.params = dr.params
    )
    coords <- drResult$drdata[, c(1, 2)]
    colnames(coords) <- c("x", "y")
    rownames(coords) <- attr(drResult$data, "Labels")
    p <- emapplot(
        x = x,
        showCategory = showCategory,
        #group = group,
        node_label = node_label,
        ...
    )

    ## Set axis label according to drfun
    p <- adj_axis(p = p, drResult = drResult)

    p + theme_dr()
}


#' @rdname ssplot
#' @param pie one of 'equal' or 'Count' to set the slice ratio of the pies
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 coord_equal
# @param cex_pie2axis It is used to adjust the relative size of the pie chart on the coordinate axis,
# the default value is 0.0125.
#' @importFrom stats setNames
ssplot.compareClusterResult <- function(
    x,
    showCategory = 30,
    #split = NULL,
    pie = "equal",
    drfun = NULL,
    #cex_pie2axis = 0.0125,
    dr.params = list(),
    node_label = "group",
    ...
) {
    if (is.null(drfun)) {
        drfun = stats::cmdscale
        dr.params = list(eig = TRUE)
    }

    if (is.character(drfun)) {
        drfun <- eval(parse(text = drfun))
    }
    split = NULL
    drResult <- get_drResult(
        x = x,
        showCategory = showCategory,
        split = split,
        pie = pie,
        drfun = drfun,
        dr.params = dr.params
    )
    coords <- drResult$drdata[, c(1, 2)]
    colnames(coords) <- c("x", "y")
    rownames(coords) <- attr(drResult$data, "Labels")
    p <- emapplot(
        x,
        showCategory = showCategory,
        coords = coords,
        split = split,
        pie = pie,
        #with_edge = with_edge,
        #cex_pie2axis = cex_pie2axis,
        #group = group,
        node_label = node_label,
        ...
    )
    ## Set axis label according to the method parameter
    p <- adj_axis(p = p, drResult = drResult)

    p + theme_dr()
}


#' Get a distance matrix
#'
#' @param x enrichment result.
#' @param showCategory number of enriched terms to display.
#' @param split separate result by 'category' variable.
#' @param pie proportion of clusters in the pie chart.
#' @noRd
build_dist <- function(x, showCategory, split = NULL, pie = NULL) {
    sim = get_pairwise_sim(
        x = x,
        showCategory = showCategory,
        split = split,
        pie = pie
    )

    # ensure symmetry
    if (!isSymmetric(sim)) {
        sim <- (sim + t(sim)) / 2
    }

    # clamp to [0,1]
    sim[is.na(sim)] <- 0
    sim <- pmin(pmax(sim, 0), 1)

    # avoid exact 1 for off-diagonal entries (some DR methods may fail)
    eps <- .Machine$double.eps
    diag(sim) <- 1
    offdiag_idx <- row(sim) != col(sim)
    sim[offdiag_idx & sim >= 1] <- 1 - eps

    stats::as.dist(1 - sim)
}


#' Get a similarity matrix
#'
#' @param x enrichment result.
#' @param showCategory number of enriched terms to display.
#' @param split separate result by 'category' variable.
#' @param pie proportion of clusters in the pie chart.
#' @noRd
get_pairwise_sim <- function(x, showCategory, split = NULL, pie = NULL) {
    if (inherits(x, "compareClusterResult")) {
        ## Optimized fortify call for large datasets
        y <- fortify(
            model = x,
            showCategory = showCategory,
            includeAll = TRUE,
            split = split
        )
        y$Cluster <- sub("\n.*", "", y$Cluster)
        
        ## Optimized pie category preparation
        pie_data <- prepare_pie_category(y, pie = pie)
        keep <- rownames(pie_data)
    } else {
        n <- update_n(x, showCategory)
        if (is.numeric(n)) {
            keep <- seq_len(min(n, nrow(x@result)))
        } else {
            keep <- match(n, rownames(x@termsim))
            keep <- keep[!is.na(keep)]
        }
    }
    
    if (length(keep) == 0) {
        yulab_abort("no enriched term found (no rows selected by showCategory).",
                        class = "no_terms_error")
    }
    
    ## Optimized termsim filling
    fill_termsim(x, keep)
}


#' Adjust axis label according to the dimension reduction method
#'
#' @param p ggplot2 object
#' @param drs dimension reduction result
#' @noRd
adj_axis <- function(p, drResult) {
    title = NULL
    eigenvalue <- drResult$eigenvalue
    if (!is.null(eigenvalue) && length(eigenvalue) >= 2) {
        total <- sum(eigenvalue)
        if (total == 0) {
            total <- 1
        }

        xlab = paste0(
            "Dimension1 (",
            format(100 * eigenvalue[1] / total, digits = 4, "%)")
        )
        ylab = paste0(
            "Dimension2 (",
            format(100 * eigenvalue[2] / sum(eigenvalue), digits = 4, "%)")
        )
    } else {
        xlab = "Dimension1"
        ylab = "Dimension2"
        if (!is.null(drResult$stress)) {
            title <- paste0("stress = ", drResult$stress)
        }
    }
    p <- p + labs(x = xlab, y = ylab, title = title)
    return(p)
}

#' Get the result of dimension reduction
#'
#' @param x enrichment result.
#' @param showCategory number of enriched terms to display.
#' @param split separate result by 'category' variable.
#' @param pie proportion of clusters in the pie chart.
#' @param drfun The function used for dimension reduction.
#' @param dr.params list, the parameters of tidydr::dr.
#' @importFrom rlang check_installed
#' @noRd
get_drResult <- function(
    x,
    showCategory,
    split = NULL,
    pie = NULL,
    drfun,
    dr.params
) {
    ## Input validation
    check_input(x, arg_name = "x")
    check_input(showCategory, arg_name = "showCategory")
    
    ## Optimized distance matrix building
    distance_mat <- build_dist(
        x = x,
        showCategory = showCategory,
        split = split,
        pie = pie
    )
    
    check_installed('tidydr', 'for `get_drResult()`')
    
    ## Optimized error handling
    drResult <- tryCatch({
        do.call(tidydr::dr, c(list(data = distance_mat, fun = drfun), dr.params))
    }, error = function(e) {
        yulab_warn("dimensionality reduction failed with provided drfun; falling back to stats::cmdscale",
                       class = "dr_fallback_warning")
        
        tryCatch({
            tidydr::dr(distance_mat, stats::cmdscale, eig = TRUE)
        }, error = function(e2) {
            yulab_abort("dimensionality reduction failed (both provided method and fallback)",
                           class = "dr_failure_error")
        })
    })

    if (is.null(drResult$drdata)) {
        yulab_warn("Wrong drfun parameter or unsupported dimensionality reduction method; using stats::cmdscale",
                       class = "dr_parameter_warning")
        drResult <- tidydr::dr(distance_mat, stats::cmdscale, eig = TRUE)
    }
    
    return(drResult)
}
