#' @rdname ridgeplot
#' @exportMethod ridgeplot
setMethod(
    "ridgeplot",
    signature(x = "gseaResult"),
    function(
        x,
        showCategory = 30,
        fill = "p.adjust",
        core_enrichment = TRUE,
        label_format = 30,
        ...
    ) {
        ridgeplot.gseaResult(
            x,
            showCategory = showCategory,
            fill = fill,
            core_enrichment = core_enrichment,
            label_format = label_format,
            ...
        )
    }
)


#' @rdname ridgeplot
#' @param orderBy The order of the Y-axis
#' @param decreasing logical. Should the orderBy order be increasing or decreasing?
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_x_reverse
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom rlang check_installed
#' @importFrom yulab.utils yulab_abort
#' @importFrom yulab.utils yulab_warn
#' @author Guangchuang Yu
ridgeplot.gseaResult <- function(
    x,
    showCategory = 30,
    fill = "p.adjust",
    core_enrichment = TRUE,
    label_format = 30,
    orderBy = "NES",
    decreasing = FALSE
) {
    ## Input validation with better error messages
    check_input(x, type = "gseaResult", arg_name = "x")
    
    if (!fill %in% colnames(x@result)) {
        yulab_abort(paste0("'", fill, "' variable not available in result"), 
                        class = "missing_column_error")
    }

    ## geom_density_ridges <- get_fun_from_pkg('ggridges', 'geom_density_ridges')
    if (orderBy != 'NES' && !orderBy %in% colnames(x@result)) {
        yulab_warn('wrong orderBy parameter; set to default `orderBy = "NES"`',
                     class = "parameter_warning")
        orderBy <- "NES"
    }
    
    ## Optimized category selection
    if (inherits(showCategory, 'numeric')) {
        selected <- seq_len(min(showCategory, nrow(x@result)))
    } else if (inherits(showCategory, "character")) {
        ii <- match(showCategory, x@result$Description)
        if (all(is.na(ii))) {
            ii <- match(showCategory, x@result$ID)
        }
        ii <- ii[!is.na(ii)]
        if (length(ii) == 0) {
            yulab_warn("No matching categories found, using first 10",
                          class = "category_warning")
            ii <- seq_len(min(10, nrow(x@result)))
        }
        selected <- x@result[ii, "ID"]
    } else {
        yulab_warn("showCategory should be a number of pathways or a vector of selected pathways",
                       class = "parameter_warning")
        selected <- seq_len(min(10, nrow(x@result)))
    }

    ## Optimized gene set extraction
    if (core_enrichment) {
        gs2id <- geneInCategory(x)[selected]
    } else {
        gs2id <- x@geneSets[names(x@geneSets) %in% selected]
    }

    ## Optimized gene name mapping
    if (x@readable && length(x@gene2Symbol) > 0) {
        gene_names <- names(x@geneList)
        symbol_match <- match(gene_names, names(x@gene2Symbol))
        valid_matches <- !is.na(symbol_match)
        names(x@geneList)[valid_matches] <- x@gene2Symbol[symbol_match[valid_matches]]
    }

    ## Vectorized data preparation
    gs2val <- lapply(gs2id, function(id) {
        res <- x@geneList[id]
        res[!is.na(res)]
    })

    nn <- names(gs2val)
    i <- match(nn, x$ID)
    nn <- x$Description[i]

    ## Optimized ordering
    order_values <- x@result[[orderBy]][i]
    j <- order(order_values, decreasing = decreasing)
    
    ## Efficient data frame construction
    len <- lengths(gs2val)
    total_len <- sum(len)
    
    gs2val.df <- data.frame(
        category = rep(nn, times = len),
        color = rep(x[i, fill], times = len),
        value = unlist(gs2val, use.names = FALSE)
    )

    colnames(gs2val.df)[2] <- fill
    gs2val.df$category <- factor(gs2val.df$category, levels = nn[j])

    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    check_installed('ggridges', 'for `ridgeplot()`.')

    ggplot(
        gs2val.df,
        aes(x = .data[["value"]], y = .data[["category"]], fill = .data[[fill]])
    ) +
        ggridges::geom_density_ridges() +
        set_enrichplot_color(type = "fill", name = fill, transform = 'log10') +
        scale_y_discrete(labels = label_func) +
        xlab(NULL) +
        ylab(NULL) +
        theme_dose()
}
