#' @import dplyr tibble ggplot2 ggraph tidygraph pals
NULL

#' Removes legend from plot.
#'
#' @export
no_legend <- function(){
    theme(
        legend.position = 'none'
    )
}

#' Removes margins from plot.
#'
#' @export
no_margin <- function(){
    theme(
        plot.margin = margin(0,0,0,0,unit='lines')
    )
}

#' Rangeframe scales.
scale_axis_rangeframe <- function(){
    guides(x='axis_truncated', y='axis_truncated')
}

#' Theme rangeframe.
theme_rangeframe <- function(){
    theme(
        axis.line = element_line(colour='black', lineend='round', size=0.3),
        axis.ticks = element_line(size=0.3),
        panel.border = element_blank()
    )
}

#' Removes x axis text.
#'
#' @export
no_x_text <- function(){
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )
}

#' Removes y axis text.
#'
#' @export
no_y_text <- function(){
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    )
}

#' Compute UMAP embedding
get_umap <- function(
    x,
    ...
){
    if (ncol(x)>100){
        pca_mat <- irlba::prcomp_irlba(x, n=n_pcs)$x
        rownames(pca_mat) <- rownames(x)
        x <- as.matrix(pca_mat)
    }
    umap_tbl <- uwot::umap(x, ...) %>%
        {colnames(.) <- c('UMAP_1', 'UMAP_2'); .} %>%
        as_tibble(rownames='gene')
    return(umap_tbl)
}


#' Plot goodness-of-fit metrics.
#'
#' @import ggpointdensity patchwork
#'
#' @param object An object.
#' @param network Name of the network to use.
#' @param point_size Float indicating the point size.
#'
#' @return A ggplot2 object.
#'
#' @rdname plot_gof
#' @export
#' @method plot_gof SeuratPlus
plot_gof.SeuratPlus <- function(
    object,
    network = 'glm_network',
    point_size = 0.5
){

    module_params <- NetworkModules(object, network=network)@params

    gof <- gof(object, network=network) %>%
        filter(rsq<=1, rsq>=0) %>%
        mutate(nice=rsq>module_params$rsq_thresh&nvariables>module_params$nvar_thresh)

    p1 <- ggplot(gof, aes(rsq, nvariables, alpha=nice)) +
        geom_pointdensity(size=point_size, shape=16) +
        geom_hline(yintercept=module_params$nvar_thresh, size=0.5, color='darkgrey', linetype='dashed') +
        geom_vline(xintercept=module_params$rsq_thresh, size=0.5, color='darkgrey', linetype='dashed') +
        scale_color_gradientn(colors=rev(pals::brewer.reds(100))) +
        scale_alpha_discrete(range=c(0.5,1)) +
        scale_y_continuous(
            trans=scales::pseudo_log_trans(base = 10),
            breaks=c(0, 1, 10, 100, 1000, 10000, 100000)
        ) +
        scale_x_continuous(breaks=seq(0,1,0.2)) +
        theme_bw() +
        no_legend() +
        labs(x=expression('Explained variance'~(R**2)), y='# variables in model') +
        theme(
            plot.margin = unit(c(0,0,0,0), 'line'),
            strip.text = element_blank()
        )

    p2 <- ggplot(gof, aes(rsq)) +
        geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
        theme_void() +
        no_legend()

    p3 <- ggplot(gof, aes(nvariables)) +
        geom_histogram(fill='darkgray', bins=20, color='black', size=0.2) +
        scale_x_continuous(
            trans=scales::pseudo_log_trans(base = 10),
            breaks=c(0, 1, 10, 100, 1000, 10000, 100000)
        ) +
        theme_void() +
        coord_flip() +
        no_legend()

    layout <- '
    AAAA#
    BBBBC
    BBBBC
    '
    p_out <- p2 + p1 + p3 + plot_layout(design = layout) & no_margin()
    return(p_out)
}


#' Plot module metrics number of genes, number of peaks and number of TFs per gene.
#'
#' @import patchwork ggh4x
#'
#' @param object An object.
#' @param network Name of the network to use.
#'
#' @return A ggplot2 object.
#'
#' @rdname plot_module_metrics
#' @export
#' @method plot_module_metrics SeuratPlus
plot_module_metrics.SeuratPlus <- function(
    object,
    network = 'glm_network'
){
    modules <- NetworkModules(object, network=network)@meta
    plot_df <- modules %>%
        distinct(target, n_regions)

    p1 <- ggplot(plot_df, aes(1, n_regions)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        theme_rangeframe() + scale_axis_rangeframe() +
        no_x_text() +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y='# peaks')


    plot_df <- modules %>%
        distinct(target, n_tfs)

    p2 <- ggplot(plot_df, aes(1, n_tfs)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        theme_rangeframe() + scale_axis_rangeframe() +
        no_x_text() +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y='# TFs')


    plot_df <- modules %>%
        distinct(tf, n_genes)

    p3 <- ggplot(plot_df, aes(1, n_genes)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        theme_rangeframe() + scale_axis_rangeframe() +
        no_x_text() +
        scale_y_continuous(
            trans=scales::pseudo_log_trans(base = 10)
        ) +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y=expression('# genes'))

    p_out <- p1 | p2 | p3 & no_margin()
    return(p_out)
}


#' Compute network graph embedding using UMAP.
#'
#' @importFrom tidyr pivot_longer
#' @import tidygraph
#'
#' @param object An object.
#' @param network Name of the network to use.
#' @param graph_name Name of the graph.
#' @param rna_assay Name of the RNA assay.
#' @param rna_slot Name of the RNA slot to use.
#' @param edge_method Method to compute edge weights:
#' * \code{'weighted'} - Correlation weighted by GRN coefficient.
#' * \code{'corr'} - Only correlation.
#' * \code{'coef'} - Only GRN coefficient.
#' @param ... Additional arguments for \code{\link[umap]{uwot}}.
#'
#' @return A SeuratPlus object.
#'
#' @rdname get_network_graph
#' @export
#' @method get_network_graph SeuratPlus
get_network_graph.SeuratPlus <- function(
    object,
    network = 'glm_network',
    graph_name = 'module_graph',
    rna_assay = 'RNA',
    rna_slot = 'data',
    edge_method = c('weighted', 'corr', 'coef'),
    features = NULL,
    random_seed = 111,
    ...
){
    edge_method <- match.arg(edge_method)
    modules <- NetworkModules(object, network=network)

    if (edge_method=='weighted'){
        rna_expr <- t(Seurat::GetAssayData(object, assay=rna_assay, slot=rna_slot))
        if (!is.null(features)){
            features <- intersect(features, colnames(rna_expr))
            rna_expr <- rna_expr[, features]
        }
        rna_expr <- rna_expr

        gene_cor <- sparse_cor(rna_expr)
        gene_cor_df <- gene_cor %>%
            as_tibble(rownames='source') %>%
            pivot_longer(!source, names_to='target', values_to='corr')

        modules_use <- modules@meta %>%
            filter(target%in%colnames(rna_expr), tf%in%colnames(rna_expr))

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target) %>%
            left_join(gene_cor_df, by=c('tf'='source', 'target')) %>%
            {.$corr[is.na(.$corr)] <- 0; .}

        reg_mat <- gene_net %>%
            select(target, tf, estimate) %>%
            pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>%
            column_to_rownames('target') %>% as.matrix() %>% Matrix::Matrix(sparse=T)

        # Layout with UMAP on adjacency matrix
        reg_factor_mat <- abs(reg_mat) + 1
        coex_mat <- gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)

    } else if (edge_method=='corr'){
        rna_expr <- t(GetAssayData(object, assay=rna_assay, slot=rna_slot))
        features <- intersect(features, colnames(rna_expr))
        rna_expr <- rna_expr[, features]

        coex_mat <- sparse_cor(rna_expr)
        gene_cor_df <- coex_mat %>%
            as_tibble(rownames='source') %>%
            pivot_longer(!source, names_to='target', values_to='corr')

        # Get adjacency df and matrix
        modules_use <- modules %>%
            filter(target%in%features, tf%in%features)

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target) %>%
            left_join(gene_cor_df, by=c('tf'='source', 'target')) %>%
            {.$corr[is.na(.$corr)] <- 0; .}

    } else if (edge_method=='coef'){
        modules_use <- modules %>%
            filter(target%in%features, tf%in%features)

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target)

        coex_mat <- gene_net %>%
            select(target, tf, estimate) %>%
            pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>%
            as_matrix() %>% Matrix::Matrix(sparse=T)
    }

    set.seed(random_seed)
    coex_umap <- get_umap(coex_mat, ...)

    gene_graph <- as_tbl_graph(gene_net) %>%
        activate(edges) %>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        mutate(dir=sign(estimate)) %>%
        activate(nodes) %>%
        mutate(centrality=centrality_pagerank()) %>%
        inner_join(coex_umap, by=c('name'='gene'))

    object@grn@networks[[network]]@graphs[[graph_name]] <- gene_graph
    return(object)
}


#' Plot network graph.
#'
#' @import tidygraph
#' @import ggraph
#'
#' @param object An object.
#' @param network Name of the network to use.
#' @param graph_name Name of the graph.
#' @param layout Layout for the graph. Can be 'umap' or any force-directed layout
#' implemented in \code{\link[ggraph]{ggraph}}
#' @param edge_width Edge width.
#' @param edge_color Edge color.
#' @param node_color Node color or color gradient.
#' @param node_size Edge size range.
#' @param text_size Text width.
#' @param color_nodes Logical, Whether to color nodes by centrality.
#' @param label_nodes Logical, Whether to label nodes with gene name.
#' @param color_edges Logical, Whether to color edges by direction.
#'
#' @return A SeuratPlus object.
#'
#' @rdname plot_network_graph
#' @export
#' @method plot_network_graph SeuratPlus
plot_network_graph.SeuratPlus <- function(
    object,
    network = 'glm_network',
    graph_name = 'module_graph',
    layout = 'umap',
    edge_width = 0.2,
    edge_color = c('-1'='darkgrey', '1'='orange'),
    node_color = pals::magma(100),
    node_size = c(1,5),
    text_size = 10,
    color_nodes = TRUE,
    label_nodes = TRUE,
    color_edges = TRUE
){
    gene_graph <- NetworkGraph(object)

    if (layout=='umap'){
        p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
    } else {
        p <- ggraph(gene_graph, layout=layout)
    }

    if (color_edges){
        p <- p + geom_edge_diagonal(aes(color=factor(dir)), width=edge_width) +
            scale_edge_color_manual(values=edge_color)
    } else {
        p <- p + geom_edge_diagonal(width=edge_width, color=edge_color[1])
    }

    if (color_nodes){
        p <- p + geom_node_point(aes(fill=centrality, size=centrality), color='darkgrey', shape=21) +
            scale_fill_gradientn(colors=node_color)
    } else {
        p <- p + geom_node_point(aes(size=centrality), color='darkgrey', shape=21, fill=node_color[1])
    }

    if (label_nodes){
        p <- p + geom_node_text(
            aes(label=name),
            repel=T, size=text_size/ggplot2::.pt, max.overlaps=99999
        )
    }
    p <- p + scale_size_continuous(range=node_size) +
        theme_void() + no_legend()

    return(p)
}

