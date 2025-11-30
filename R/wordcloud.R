##' Use wordcloud algorithm to get group tags
##'
##' @param cluster a cluster name
##' @param node_data the data section of the ggplot object,
##' which contains clustering information.
##' @param nWords the number of words in the cluster tags
#' @importFrom utils head
##' @noRd
get_wordcloud <- function(cluster, node_data, nWords = 4) {
    cluster_terms <- node_data$name[node_data$color2 == cluster]    

    if (length(cluster_terms) == 0) {
        return(cluster)
    }

    words <- cluster_terms |>
        tolower() |>
        gsub(" in ", " ", x = _) |>
        gsub(" [0-9]+ ", " ", x = _) |>
        gsub("^[0-9]+ ", "", x = _) |>
        gsub(" [0-9]+$", "", x = _) |>
        gsub(" [a-z] ", " ", x = _) |>
        gsub("^[a-z] ", "", x = _) |>
        gsub(" [a-z]$", "", x = _) |>
        gsub(" / ", " ", x = _) |>
        gsub(" and ", " ", x = _) |>
        gsub(" of ", " ", x = _) |>
        gsub(",", " ", x = _) |>
        gsub(" - ", " ", x = _) |>
        gsub("\\s+", " ", x = _) |> # multiple spaces to single space
        trimws() # remove leading/trailing whitespace

    # Split into words and calculate frequencies
    all_words <- unlist(strsplit(words, "\\s+"))

    if (length(all_words) == 0) {
        return(cluster)
    }

    word_freq <- table(all_words)
    word_freq <- word_freq[order(word_freq, decreasing = TRUE)]

    # Remove common stop words
    stop_words <- c(
        "the",
        "and",
        "for",
        "with",
        "via",
        "by",
        "to",
        "a",
        "an",
        "in",
        "of",
        "on",
        "at"
    )
    meaningful_words <- names(word_freq)[
        !tolower(names(word_freq)) %in% stop_words
    ]

    # Get top nWords meaningful words
    if (length(meaningful_words) > 0) {
        top_words <- head(meaningful_words, nWords)

        # Consider word position for ordering (optional enhancement)
        word_positions <- calculate_word_positions(cluster_terms, top_words)
        if (!is.null(word_positions)) {
            top_words <- word_positions
        }

        return(paste(top_words, collapse = " "))
    } else {
        # Fallback: use most frequent words regardless
        top_words <- head(names(word_freq), nWords)
        return(paste(top_words, collapse = " "))
    }
}


#' Calculate word positions to improve label ordering
#'
#' @param terms vector of terms
#' @param candidate_words candidate words for the label
#' @return ordered words based on position
#' @noRd
calculate_word_positions <- function(terms, candidate_words) {
    if (length(terms) == 0 || length(candidate_words) == 0) {
        return(NULL)
    }

    # Split all terms into words
    all_term_words <- strsplit(tolower(terms), "\\s+")

    # Calculate average position for each candidate word
    word_ranks <- list()

    for (word in candidate_words) {
        positions <- c()
        for (term_words in all_term_words) {
            word_idx <- which(term_words == tolower(word))
            if (length(word_idx) > 0) {
                positions <- c(positions, word_idx)
            }
        }
        if (length(positions) > 0) {
            word_ranks[[word]] <- mean(positions)
        } else {
            word_ranks[[word]] <- Inf # Word not found in any term
        }
    }

    # Order words by their average position
    if (length(word_ranks) > 0) {
        sorted_words <- names(sort(unlist(word_ranks)))
        return(sorted_words)
    }

    return(NULL)
}
