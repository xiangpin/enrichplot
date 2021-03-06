##' @rdname emapplot_cluster
##' @exportMethod emapplot_cluster
setMethod("emapplot_cluster", signature(x = "enrichResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        emapplot_cluster.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname emapplot_cluster
##' @exportMethod emapplot_cluster
setMethod("emapplot_cluster", signature(x = "gseaResult"),
    function(x, showCategory = 30, color = "p.adjust", label_format = 30, ...) {
        emapplot_cluster.enrichResult(x, showCategory = showCategory,
            color = color, label_format = label_format, ...)
    })

##' @rdname emapplot_cluster
##' @exportMethod emapplot_cluster
setMethod("emapplot_cluster", signature(x = "compareClusterResult"),
    function(x, showCategory = 30, color = "p.adjust", ...) {
        emapplot_cluster.compareClusterResult(x, showCategory = showCategory,
            color=color, label_format = label_format, ...)
    })


##' @rdname emapplot_cluster
##' @param with_edge Logical, if TRUE (the default), draw the edges of the network diagram.
##' @param cex_line Scale of line width
##' @param nWords Numeric, the number of words in the cluster tags, the default value is 4.
##' @param nCluster Numeric, the number of clusters, 
##' the default value is square root of the number of nodes.
##' @param split Separate result by 'category' variable
##' @param min_edge The minimum similarity threshold for whether 
##' two nodes are connected, should between 0 and 1, default value is 0.2.
##' @param cex_label_group Numeric, scale of group labels size, the default value is 1.
##' @param label_style style of group label, one of "shadowtext" and "ggforce".
##' @param group_legend Logical, if TRUE, the grouping legend will be displayed.
##' The default is FALSE
##' @param cex_category Numeric, indicating the size by which plotting category
##' nodes should be scaled relative to the default.
##' @param label_format a numeric value sets wrap length, alternatively a
##' custom function to format axis labels.
##' @param repel whether to correct the position of the label. Defaults to FALSE.
##' @param shadowtext a logical value, whether to use shadow font. Defaults to TRUE.
##' @param ... Additional parameters used to set the position of the group label.
##' When the parameter repel is set to TRUE, additional parameters will take effect.
##' 
##' additional parameters can refer the following parameters.
##'     \itemize{
##'        \item \code{force} Force of repulsion between overlapping text labels. Defaults to 1. 
##'        \item \code{nudge_x, nudge_y} Horizontal and vertical adjustments to nudge 
##'         the starting position of each text label. 
##'        \item \code{direction} "both", "x", or "y" – direction in which to adjust position of labels.
##'     }
##'
##' @importFrom igraph layout_with_fr
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 scale_color_discrete
##' @importFrom ggplot2 scale_size_continuous
##' @importFrom ggplot2 scale_fill_discrete
##' @importFrom ggplot2 geom_text
##' @importFrom stats kmeans
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_node_point
##' @importFrom ggraph geom_edge_link
##' @importFrom DOSE geneInCategory
##' @importFrom GOSemSim godata
##' @importFrom shadowtext geom_shadowtext
emapplot_cluster.enrichResult <- function(x, showCategory = 30,
                                          color = "p.adjust", cex_line = 1,
                                          with_edge = TRUE,
                                          nWords = 4, nCluster = NULL,
                                          split = NULL, min_edge = 0.2,
                                          cex_label_group = 1, 
                                          label_style = "shadowtext", 
                                          group_legend = FALSE, cex_category = 1, 
                                          label_format = 30, repel = FALSE, 
                                          shadowtext = TRUE, ...){
                                          

    has_pairsim(x)
    label_group <- 3
    n <- update_n(x, showCategory)
    y <- as.data.frame(x)

    g <- get_igraph(x=x, y=y, n=n, color=color, cex_line=cex_line,
                    min_edge=min_edge)
    if(n == 1) {
        return(ggraph(g, "tree") + geom_node_point(color="red", size=5) +
                 geom_node_text(aes_(label=~name)))
    }
    edgee <- igraph::get.edgelist(g)
 
    edge_w <- E(g)$weight
    set.seed(123)
    lw <- layout_with_fr(g, weights=edge_w)

    p <- ggraph::ggraph(g, layout=lw)
    # cluster_label1 <- lapply(clusters, function(i){i[order(y[i, "pvalue"])[1]]})

    ## Using k-means clustering to group
    pdata2 <- p$data
    dat <- data.frame(x = pdata2$x, y = pdata2$y)
    colnames(pdata2)[5] <- "color2"

    if(is.null(nCluster)){
        pdata2$color <- kmeans(dat, ceiling(sqrt(nrow(dat))))$cluster
    } else {
        if(nCluster > nrow(dat)) nCluster <- nrow(dat)
        pdata2$color <- kmeans(dat, nCluster)$cluster
    }

    goid <- y$ID
    cluster_color <- unique(pdata2$color)
    clusters <- lapply(cluster_color, function(i){goid[which(pdata2$color == i)]})
    cluster_label <- sapply(cluster_color, get_wordcloud, pdata2 = pdata2,
                            nWords=nWords)
    names(cluster_label) <- cluster_color
    pdata2$color <- cluster_label[as.character(pdata2$color)]
    p$data <- pdata2
    ## Take the location of each group's center nodes as the location of the label
    # label_func <- default_labeller(label_format)
    # if (is.function(label_format)) {
    #     label_func <- label_format
    # }

    # label_x <- stats::aggregate(x ~ color, pdata2, mean)
    # label_y <- stats::aggregate(y ~ color, pdata2, mean)
    # label_location <- data.frame(x = label_x$x, y = label_y$y,
    #                              # label = label_x$color)
    #                              label = label_func(label_x$color))
    label_location <- get_label_location(pdata2, label_format)
    ## Adjust the label position up and down to avoid overlap
    # rownames(label_location) <- label_location$label
    # label_location <- adjust_location(label_location, x_adjust, y_adjust)
    ## use spread.labs
    # label_location$y <- TeachingDemos::spread.labs(x = label_location$y, mindiff = cex_label_group*y_adjust)
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color")

    if (with_edge) {
        p <-  p +  ggraph::geom_edge_link(alpha = .8,
                       aes_(width =~ I(width*cex_line)), colour='darkgrey')
    }

    if (label_style == "shadowtext") {
        p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
                     fill =~ color), show.legend = show_legend)
    } else {
        p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
                     fill =~ color, label =~ color), show.legend = show_legend)
    }

    if (group_legend) p <- p + scale_fill_discrete(name = "groups")
    p <- p + ggnewscale::new_scale_fill() +
        geom_point(shape = 21, aes_(x =~ x, y =~ y, fill =~ color2,
                                    size =~ size)) +
        scale_size_continuous(name = "number of genes",
                              range = c(3, 8) * cex_category) +
        scale_fill_continuous(low = "red", high = "blue", name = color,
                              guide = guide_colorbar(reverse = TRUE)) + 
        theme(legend.title = element_text(size = 10),
                   legend.text  = element_text(size = 10)) +
        theme(panel.background = element_blank())
    if (label_style == "ggforce") return(p)  
    add_group_label(repel = repel, shadowtext = shadowtext, p = p, 
        label_location = label_location, label_group = label_group, 
        cex_label_group = cex_label_group, ...)                 
}



##' @rdname emapplot_cluster
##' @importFrom igraph E
##' @importFrom ggplot2 aes_
##' @importFrom ggplot2 coord_equal
##' @importFrom ggraph ggraph
##' @importFrom ggraph geom_edge_link
##' @importFrom scatterpie geom_scatterpie
##' @importFrom scatterpie geom_scatterpie_legend
##' @importClassesFrom DOSE compareClusterResult
##' @importFrom stats setNames
##' @param pie Proportion of clusters in the pie chart, one of
##' 'equal' (default) and 'Count'.
##' @param legend_n Number of circle in legend, the default value is 5.
emapplot_cluster.compareClusterResult <- function(x, showCategory = 30,
    color = "p.adjust", cex_line = 1, with_edge = TRUE,
    nWords = 4, nCluster = NULL, split = NULL, min_edge = 0.2,
    cex_label_group = 1, pie = "equal", legend_n = 5,
    cex_category = 1, label_style = "shadowtext", group_legend = FALSE, 
    label_format = 30, repel = FALSE, shadowtext = TRUE, ...){

    has_pairsim(x)
    label_group <- 3
    y <- get_selected_category(showCategory, x, split)  
    ## Data structure transformation, combining the same ID (Description) genes
    y_union <- merge_compareClusterResult(y)
    
    geneSets <- setNames(strsplit(as.character(y_union$geneID), "/",
                                  fixed = TRUE), y_union$ID) 

                                  
    g <- build_emap_graph(y=y_union,geneSets=geneSets,color=color,
        cex_line=cex_line, min_edge=min_edge, pair_sim = x@termsim, 
        method = x@method)
    p <- build_ggraph(y = y, g = g, y_union = y_union, cex_category = cex_category,
        pie = pie, layout = "nicely")
    if (is.null(dim(y)) | nrow(y) == 1 | is.null(dim(y_union)) | nrow(y_union) == 1)
        return(p)

    ## then add the pie plot
    ## Get the matrix data for the pie plot
    ID_Cluster_mat <- prepare_pie_category(y, pie=pie)

    # Start the cluster diagram
    edge_w <- E(g)$weight
    set.seed(123)
    lw <- layout_with_fr(g, weights=edge_w)
    p <- ggraph(g, layout=lw)

    ## Using k-means clustering to group
    pdata2 <- p$data
    dat <- data.frame(x = pdata2$x, y = pdata2$y)
    colnames(pdata2)[5] <- "color2"
    
    if (is.null(nCluster)){
        pdata2$color <- kmeans(dat, floor(sqrt(nrow(dat))))$cluster
    } else {
        if (nCluster > nrow(dat)) nCluster <- nrow(dat)
        pdata2$color <- kmeans(dat, nCluster)$cluster
    }

    goid <- y_union$ID
    cluster_color <- unique(pdata2$color)
    clusters <- lapply(cluster_color, function(i){goid[which(pdata2$color == i)]})
    cluster_label <- sapply(cluster_color,  get_wordcloud, pdata2 = pdata2,
                            nWords=nWords)
    names(cluster_label) <- cluster_color
    pdata2$color <- cluster_label[as.character(pdata2$color)]
    p$data <- pdata2

    #plot the edge
    #get the X-coordinate and y-coordinate of pies
    pdata2 <- p$data

    desc <- y_union$Description[match(rownames(ID_Cluster_mat),
                                      y_union$Description)]
    i <- match(desc, pdata2$name)

    ID_Cluster_mat$x <- pdata2$x[i]
    ID_Cluster_mat$y <- pdata2$y[i]

    #Change the radius value to fit the pie plot
    radius <- NULL
    ID_Cluster_mat$radius <- sqrt(pdata2$size[i] / sum(pdata2$size) * cex_category) 

    x_loc1 <- min(ID_Cluster_mat$x)
    y_loc1 <- min(ID_Cluster_mat$y)

    ## Take the location of each group's center nodes as the location of the label
    label_func <- default_labeller(label_format)
    if (is.function(label_format)) {
        label_func <- label_format
    }

    # label_x <- stats::aggregate(x ~ color, pdata2, mean)
    # label_y <- stats::aggregate(y ~ color, pdata2, mean)
    # label_location <- data.frame(x = label_x$x, y = label_y$y,
    #                              # label = label_x$color)
    #                              label = label_func(label_x$color))
    label_location <- get_label_location(pdata2, label_format)
    ## Adjust the label position up and down to avoid overlap
    # rownames(label_location) <- label_location$label
    # label_location <- adjust_location(label_location, x_adjust, y_adjust)
    # label_location$y <- TeachingDemos::spread.labs(x = label_location$y, mindiff = cex_label_group*y_adjust)
    show_legend <- c(group_legend, FALSE)
    names(show_legend) <- c("fill", "color")

    if(with_edge) {
        p <-  p +  geom_edge_link(alpha = .8, aes_(width =~ I(width*cex_line)),
                                  colour='darkgrey')
    }

    if(label_style == "shadowtext") {
        p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
            fill =~ color), show.legend = show_legend)
    } else {
        p <- p + ggforce::geom_mark_ellipse(aes_(x =~ x, y =~ y, color =~ color,
            fill =~ color, label =~ color), show.legend = show_legend)
    }

    if(group_legend) p <- p + scale_fill_discrete(name = "groups")

    p <- p + ggnewscale::new_scale_fill() +
        geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=ID_Cluster_mat,
            cols=colnames(ID_Cluster_mat)[1:(ncol(ID_Cluster_mat)-3)],color=NA) +
        coord_equal()+
        geom_scatterpie_legend(ID_Cluster_mat$radius, x=x_loc1, y=y_loc1,
            n = legend_n,
            labeller = function(x) round(sum(pdata2$size) * x^2 / cex_category)) + 
        theme(legend.title = element_text(size = 10),
              legend.text  = element_text(size = 10)) +
        theme(panel.background = element_blank())
    if (label_style == "ggforce") return(p)
    add_group_label(repel = repel, shadowtext = shadowtext, p = p, 
        label_location = label_location, label_group = label_group, 
        cex_label_group = cex_label_group, ...) 
}



























