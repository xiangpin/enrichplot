loadNamespace("DOSE")

mock_enrich_result <- function() {
    result <- data.frame(
        ID = c("T1", "T2"),
        Description = c("dup", "dup"),
        GeneRatio = c("1/2", "1/2"),
        BgRatio = c("1/10", "1/10"),
        pvalue = c(0.01, 0.02),
        p.adjust = c(0.02, 0.03),
        qvalue = c(0.02, 0.03),
        geneID = c("g1/g2", "g2/g3"),
        Count = c(2L, 2L),
        stringsAsFactors = FALSE
    )
    rownames(result) <- result$ID

    methods::new(
        "enrichResult",
        result = result,
        pvalueCutoff = 1,
        pAdjustMethod = "BH",
        qvalueCutoff = 1,
        organism = "mock",
        ontology = "mock",
        gene = c("g1", "g2", "g3"),
        keytype = "UNKNOWN",
        universe = character(),
        gene2Symbol = character(),
        geneSets = list(
            T1 = c("g1", "g2"),
            T2 = c("g2", "g3")
        ),
        readable = FALSE,
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list()
    )
}

mock_foldchange <- function() {
    c(g1 = 1, g2 = 2, g3 = 3)
}
