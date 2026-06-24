#' Color utility functions for enrichplot package
#'
#' This file contains all color-related helper functions

#' Get default enrichplot colors
#'
#' @param n number of colors (2 or 3)
#' @return color vector
#' @export
get_enrichplot_color <- function(n = 2) {
    colors <- getOption("enrichplot.colours")
    if (!is.null(colors)) {
        colors <- as.character(colors)
        if (length(colors) < 2) {
            rlang::abort(
                "`options(enrichplot.colours = ...)` must provide at least 2 colors."
            )
        }
        if (n == 2) {
            return(colors[seq_len(2)])
        }
        if (n == 3) {
            if (length(colors) == 2) {
                return(c(colors[1], "white", colors[2]))
            }
            return(colors[seq_len(3)])
        }
        return(colors)
    }

    if (n != 2 && n != 3) {
        rlang::abort("'n' should be 2 or 3", .arg = "n")
    }

    colors = c("#e06663", "#327eba")
    if (n == 2) {
        return(colors)
    }

    if (n == 3) return(c(colors[1], "white", colors[2]))
}

#' Helper function to set color scale for enrichplot
#'
#' @param colors user provided color vector
#' @param type one of 'color', 'colour' or 'fill'
#' @param name name of the color legend
#' @param .fun force to use user provided color scale function
#' @param reverse whether reverse the color scheme
#' @param transform transform the color scale
#' @param ... additional parameters
#' @return a color scale
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 scale_color_continuous
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 guide_colorbar
#' @importFrom yulab.utils yulab_abort
#' @importFrom yulab.utils yulab_warn
#' @importFrom yulab.utils check_input
#' @author Guangchuang Yu
#' @export
set_enrichplot_color <- function(
    colors = get_enrichplot_color(2),
    type = "color",
    name = NULL,
    .fun = NULL,
    reverse = TRUE,
    transform = 'identity',
    ...
) {
    type <- match.arg(type, c("color", "colour", "fill"))
    if (reverse) {
        colors = rev(colors)
    }
    n <- length(colors)
    ## Input validation
    check_input(colors, type = "character", min_length = 2, arg_name = "colors")
    check_input(type, type = "character", arg_name = "type")
    
    if (n < 2) {
        yulab_abort("'colors' should be of length >= 2", class = "color_length_error")
    } else if (n == 2) {
        params <- list(low = colors[1], high = colors[2])
        fn_suffix <- "continuous"
    } else if (n == 3) {
        params <- list(low = colors[1], mid = colors[2], high = colors[3])
        fn_suffix <- "gradient2"
    } else {
        params <- list(colors = colors)
        fn_suffix <- "gradientn"
    }

    if (!is.null(.fun)) {
        if (n == 3) {
            fn_type <- which_scale_fun(.fun)
            if (fn_type == "gradientn") {
                params <- list(colors = colors)
            } else {
                params <- list(
                    low = colors[1],
                    mid = colors[2], 
                    high = colors[3]
                )
            }
        }
    } else {
        fn <- sprintf("scale_%s_%s", type, fn_suffix)
        .fun <- getFromNamespace(fn, "ggplot2")
    }

    params$guide <- guide_colorbar(reverse = reverse, order = 1)
    params$name <- name
    params$transform <- transform

    params <- modifyList(params, list(...))

    do.call(.fun, params)
}

#' Determine which scale function to use
#'
#' @param .fun function to check
#' @return scale function type
#' @noRd
which_scale_fun <- function(.fun) {
    params <- args(.fun) |> as.list() |> names()
    if ("colours" %in% params) {
        return("gradientn")
    }
    if ("mid" %in% params) {
        return("gradient2")
    }
    return("continuous")
}

#' Create color palette for continuous data
#'
#' @param colors colors of length >=2
#' @return color vector
#' @importFrom rlang check_installed
#' @importFrom yulab.utils check_input
#' @export
#' @examples
#' color_palette(c("red", "yellow", "green"))
#' @author guangchuang yu
color_palette <- function(colors) {
    ## Check input validity
    yulab.utils::check_input(colors, type = "character", min_length = 2, arg_name = "colors")
    
    rlang::check_installed('grDevices', 'for `color_palette()`.')
    grDevices::colorRampPalette(colors)(n = 299)
}

#' Predefined color palettes
enrichplot_point_shape <- ggfun:::enrichplot_point_shape
sig_palette <- color_palette(c("red", "yellow", "blue"))
heatmap_palette <- color_palette(c("red", "yellow", "green"))
