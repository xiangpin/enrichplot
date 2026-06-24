#' @rdname heatplot
#' @exportMethod heatplot
setMethod(
    "heatplot",
    signature(x = "enrichResult"),
    function(x, showCategory = 30, ...) {
        heatplot.enrichResult(x, showCategory, ...)
    }
)

#' @rdname heatplot
#' @exportMethod heatplot
setMethod(
    "heatplot",
    signature(x = "gseaResult"),
    function(x, showCategory = 30, ...) {
        heatplot.enrichResult(x, showCategory, ...)
    }
)


#' @rdname heatplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom rlang check_installed
#' @param showTop number of top genes ranked by `abs(foldChange) * frequency`
#' to be shown in the heatmap, default NULL means all genes are shown
#' @param label_format a numeric value sets wrap length, alternatively a
#' custom function to format axis labels.
#' by default wraps names longer than 30 characters
#' @param symbol symbol of the nodes, one of "rect" (the default) or "dot"
#' @param pvalue pvalue of genes
#' @author Guangchuang Yu
heatplot.enrichResult <- function(
    x,
    showCategory = 30,
    showTop = NULL,
    symbol = "rect",
    foldChange = NULL,
    pvalue = NULL,
    label_format = 30
) {
    symbol <- match.arg(symbol, c("rect", "dot"))
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    n <- update_n(x, showCategory)
    geneSets <- extract_geneSets(x, n)
    if (!is.null(showTop) && showTop > 0) {
        if (is.null(foldChange)) {
            yulab.utils::yulab_abort(
                "`showTop` requires `foldChange`.",
                class = "missing_foldchange_error"
            )
        }
        nfreq <- table(unlist(geneSets))
        nfc <- nfreq * abs(foldChange[names(nfreq)])
        topgenes <- head(names(sort(nfc, decreasing = TRUE)), showTop)
        geneSets <- lapply(geneSets, function(s) intersect(s, topgenes))
        geneSets <- set_geneSet_labels(geneSets, get_geneSet_labels(geneSets))
    }
    foldChange <- fc_readable(x, foldChange)
    pvalue <- fc_readable(x, pvalue)
    d <- list2df(geneSets)
    d$categoryID <- get_geneSet_labels(geneSets)[as.character(d$categoryID)]
    if (!is.null(foldChange)) {
        d$foldChange <- foldChange[as.character(d$Gene)]
    }

    if (!is.null(pvalue)) {
        d$pvalue <- pvalue[as.character(d$Gene)]
    }

    p <- ggplot(d, aes(x = .data$Gene, y = .data$categoryID))

    if (symbol == "rect") {
        p <- p + geom_tile(color = 'white')
    }

    get_dotp <- function(p, foldChange, pvalue) {
        if (is.null(foldChange) && is.null(pvalue)) {
            p <- p +
                geom_point(
                    color = 'black',
                    shape = 21,
                    fill = "black",
                    size = 5
                )
            return(p)
        }
        if (!is.null(foldChange) && !is.null(pvalue)) {
            p <- p + geom_point(color = 'black', shape = 21)
            return(p)
        }

        if (is.null(foldChange)) {
            p <- p + geom_point(color = 'black', shape = 21, fill = "black")
        } else {
            p <- p + geom_point(color = 'black', shape = 21, size = 5)
        }

        return(p)
    }

    # copy from https://stackoverflow.com/questions/11053899/how-to-get-a-reversed-log10-scale-in-ggplot2
    reverselog_trans <- function(base = exp(1)) {
        trans <- function(x) -log(x, base)

        check_installed('scales', 'for `heatplot()`.')

        inv <- function(x) base^(-x)
        scales::trans_new(
            paste0("reverselog-", format(base)),
            trans,
            inv,
            scales::log_breaks(base = base),
            domain = c(1e-100, Inf)
        )
    }

    if (symbol == "dot") {
        p <- get_dotp(p, foldChange, pvalue)
        ## only dot need size(pvalue) parameter
        if (!is.null(pvalue)) {
            p <- p +
                aes(size = .data$pvalue) +
                scale_size_continuous(
                    range = c(3, 8),
                    trans = reverselog_trans(10)
                )
        }
    }

    if (!is.null(foldChange)) {
        p <- p +
            aes(fill = !!sym('foldChange')) +
            set_enrichplot_color(
                colors = get_enrichplot_color(3),
                type = "fill",
                reverse = FALSE,
                transform = 'identity'
            )
    }

    p +
        xlab(NULL) +
        ylab(NULL) +
        theme_minimal() +
        scale_y_discrete(labels = label_func) +
        theme(
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 60, hjust = 1)
        )
}
