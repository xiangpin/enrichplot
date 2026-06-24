loadNamespace("DOSE")
loadNamespace("clusterProfiler")

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

mock_comparecluster_result <- function() {
    gc <- list(
        A = c("1", "2", "3"),
        B = c("2", "3", "4")
    )
    df <- data.frame(
        Cluster = c("A", "A", "B", "B"),
        ID = c("T1", "T2", "T1", "T2"),
        Description = c("dup", "other", "dup", "other"),
        geneID = c("1/2", "2/3", "2/3", "3/4"),
        Count = c(2L, 2L, 2L, 2L),
        p.adjust = c(0.01, 0.02, 0.03, 0.04),
        stringsAsFactors = FALSE
    )

    methods::new(
        "compareClusterResult",
        compareClusterResult = df,
        geneClusters = gc,
        fun = "mock",
        gene2Symbol = character(),
        keytype = "UNKNOWN",
        readable = FALSE,
        .call = quote(mock()),
        termsim = matrix(0, 0, 0),
        method = "",
        dr = list(),
        organism = "mock"
    )
}
