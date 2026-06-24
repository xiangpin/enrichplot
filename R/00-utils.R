#' Internal package helpers loaded early
#'
#' @noRd
require_suggested <- function(package, feature) {
    for (pkg in as.character(package)) {
        rlang::check_installed(pkg, feature)
    }
    invisible(TRUE)
}
