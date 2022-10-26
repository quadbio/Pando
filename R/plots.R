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
    n_pcs = 30,
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
    network = DefaultNetwork(object),
    point_size = 0.5
){

    module_params <- NetworkModules(object, network=network)@params

    if (length(module_params)==0){
        stop('No modules found, please run `find_modules()` first.')
    }

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
    network = DefaultNetwork(object)
){
    modules <- NetworkModules(object, network=network)@meta

    if (nrow(modules)==0){
        stop('No modules found, please run `find_modules()` first.')
    }

    plot_df <- modules %>%
        distinct(target, n_regions)

    p1 <- ggplot(plot_df, aes(1, n_regions)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        no_x_text() +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y='# peaks') +
        ggtitle('# regions\nper target gene')


    plot_df <- modules %>%
        distinct(target, n_tfs)

    p2 <- ggplot(plot_df, aes(1, n_tfs)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        no_x_text() +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y='# TFs') +
        ggtitle('# TFs per\ntarget gene')


    plot_df <- modules %>%
        distinct(tf, n_genes)

    p3 <- ggplot(plot_df, aes(1, n_genes)) +
        geom_violin(size=0.2, fill='darkgrey', color='black') +
        theme_bw() +
        no_x_text() +
        scale_y_continuous(
            trans=scales::pseudo_log_trans(base = 10)
        ) +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y=expression('# genes')) +
        ggtitle('# target genes\nper TF')

    p_out <- p3 | p1 | p2 & no_margin()
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
#' @param umap_method Method to compute edge weights for UMAP:
#' * \code{'weighted'} - Correlation weighted by GRN coefficient.
#' * \code{'corr'} - Only correlation.
#' * \code{'coef'} - Only GRN coefficient.
#' * \code{'none'} - Don't compute UMAP and create the graph directly
#' from modules.
#' @param features Features to use to create the graph. If \code{NULL}
#' uses all features in the network.
#' @param random_seed Random seed for UMAP computation
#' @param ... Additional arguments for \code{\link[umap]{uwot}}.
#' @param verbose Print messages.
#'
#' @return A SeuratPlus object.
#'
#' @rdname get_network_graph
#' @export
#' @method get_network_graph SeuratPlus
get_network_graph.SeuratPlus <- function(
    object,
    network = DefaultNetwork(object),
    graph_name = 'module_graph',
    rna_assay = 'RNA',
    rna_slot = 'data',
    umap_method = c('weighted', 'corr', 'coef', 'none'),
    features = NULL,
    random_seed = 111,
    verbose = TRUE,
    ...
){
    umap_method <- match.arg(umap_method)
    modules <- NetworkModules(object, network=network)

    if (length(modules@params)==0){
        stop('No modules found, please run `find_modules()` first.')
    }

    if (is.null(features)){
        features <- NetworkFeatures(object, network=network)
    }

    if (umap_method=='weighted'){
        rna_expr <- t(Seurat::GetAssayData(object, assay=rna_assay, slot=rna_slot))
        features <- intersect(features, colnames(rna_expr))

        log_message('Computing gene-gene correlation', verbose=verbose)
        rna_expr <- rna_expr[, features]
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
            column_to_rownames('target') %>% as.matrix()

        log_message('Computing weighted regulatory factor', verbose=verbose)
        # Layout with UMAP on adjacency matrix
        reg_factor_mat <- abs(reg_mat) + 1
        coex_mat <- gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)

    } else if (umap_method=='corr'){
        net_features <- NetworkFeatures(object, network=network)
        rna_expr <- t(Seurat::GetAssayData(object, assay=rna_assay, slot=rna_slot))

        if (!is.null(features)){
            features <- intersect(intersect(features, colnames(rna_expr)), net_features)
        } else {
            features <- net_features
        }

        log_message('Computing gene-gene correlation', verbose=verbose)
        rna_expr <- rna_expr[, features]
        coex_mat <- sparse_cor(rna_expr)
        gene_cor_df <- coex_mat %>%
            as_tibble(rownames='source') %>%
            pivot_longer(!source, names_to='target', values_to='corr')

        # Get adjacency df and matrix
        modules_use <- modules@meta %>%
            filter(target%in%features, tf%in%features)

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target) %>%
            left_join(gene_cor_df, by=c('tf'='source', 'target')) %>%
            {.$corr[is.na(.$corr)] <- 0; .}

    } else if (umap_method=='coef'){
        modules_use <- modules@meta %>%
            filter(target%in%features, tf%in%features)

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target)

        coex_mat <- gene_net %>%
            select(target, tf, estimate) %>%
            pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>%
            column_to_rownames('target') %>% as.matrix()

    } else if (umap_method=='none' | is.null(umap_method)){
        modules_use <- modules@meta %>%
            filter(target%in%features, tf%in%features)

        gene_net <- modules_use %>%
            select(tf, target, everything()) %>%
            group_by(target)

        log_message('Getting network graph', verbose=verbose)
        gene_graph <- as_tbl_graph(gene_net) %>%
            activate(edges) %>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            mutate(dir=sign(estimate)) %>%
            activate(nodes) %>%
            mutate(centrality=centrality_pagerank())

        object@grn@networks[[network]]@graphs[[graph_name]] <- gene_graph
        return(object)
    }

    log_message('Computing UMAP embedding', verbose=verbose)
    set.seed(random_seed)
    coex_umap <- get_umap(as.matrix(coex_mat), ...)

    log_message('Getting network graph', verbose=verbose)
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
#' @param graph Name of the graph.
#' @param layout Layout for the graph. Can be 'umap' or any force-directed layout
#' implemented in \code{\link[ggraph]{ggraph}}
#' @param edge_width Edge width.
#' @param edge_color Edge color.
#' @param node_color Node color or color gradient.
#' @param node_size Node size range.
#' @param text_size Font size for labels.
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
    network = DefaultNetwork(object),
    graph = 'module_graph',
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
    gene_graph <- NetworkGraph(object, network=network, graph=graph)

    has_umap <- 'UMAP_1' %in% colnames(as_tibble(activate(gene_graph, 'nodes')))
    if (layout=='umap' & !has_umap){
        stop('No UMAP coordinates found, please run `get_network_graph()` first.')
    }

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
        p <- p + geom_node_point(
            color='darkgrey', shape=21, fill='lightgrey', size=node_size[1], stroke=0.5
        )
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


#' Get sub-network centered around one TF.
#'
#' @import tidygraph
#' @import ggraph
#'
#' @param object An object.
#' @param tf The transcription factor to center around.
#' @param network Name of the network to use.
#' @param graph Name of the graph.
#' @param features Features to use. If \code{NULL} uses all features in the graph.
#' @param order Integer indicating the maximal order of the graph.
#' @param keep_all_edges Logical, whether to maintain all edges to each leaf
#' or prune to the strongest overall connection.
#' @param verbose Logical. Whether to print messages.
#' @param parallel Logical. Whether to parallelize the computation with \code{\link[foreach]{foreach}}.
#'
#' @return A SeuratPlus object.
#'
#' @rdname get_tf_network
#' @export
#' @method get_tf_network SeuratPlus
get_tf_network.SeuratPlus <- function(
    object,
    tf,
    network = DefaultNetwork(object),
    graph = 'module_graph',
    features = NULL,
    order = 3,
    keep_all_edges = FALSE,
    verbose = TRUE,
    parallel = FALSE
){
    gene_graph <- NetworkGraph(object, network=network, graph=graph)
    gene_graph_nodes <- gene_graph %N>% as_tibble()

    if (is.null(features)){
        features <- NetworkFeatures(object, network=network)
    }

    features <- intersect(features, gene_graph_nodes$name)

    log_message('Getting shortest paths from TF', verbose=verbose)
    spaths <- igraph::all_shortest_paths(gene_graph, tf, features, mode='out')$res
    spath_list <- map_par(spaths, function(p){
        edg <- names(p)
        edg_graph <- gene_graph %>%
            filter(name%in%edg) %>%
            convert(to_shortest_path, from=which(.N()$name==edg[1]), to=which(.N()$name==edg[length(edg)])) %E>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            as_tibble()

        edg_dir <- edg_graph %>% pull(estimate) %>% sign() %>% prod()
        edg_est <- edg_graph %>% pull(estimate) %>% mean()
        path_df <- tibble(
            start_node = edg[1],
            end_node = edg[length(edg)],
            dir = edg_dir,
            path = paste(edg, collapse=';'),
            path_regions = paste(edg_graph$regions, collapse=';'),
            order = length(edg)-1,
            mean_estimate = edg_est
        )

        if ('padj' %in% colnames(edg_graph)){
            path_df$mean_log_padj <- edg_graph %>% pull(padj) %>% {-log10(.)} %>% mean()
        }

        return(
            list(
                path = path_df,
                graph = mutate(
                    edg_graph,
                    path=paste(edg, collapse=';'),
                    end_node=edg[length(edg)],
                    comb_dir=edg_dir
                )
            )
        )
    }, parallel=parallel)

    log_message('Pruning graph', verbose=verbose)
    spath_dir <- map_dfr(spath_list, function(x) x$path) %>%
        mutate(
            path_genes=str_split(path, ';'),
            path_regions=str_split(path_regions, ';')
        )
    spath_graph <- map_dfr(spath_list, function(x) x$graph)

    grn_pruned <- spath_dir %>%
        select(start_node, end_node, everything()) %>%
        group_by(end_node) %>% filter(order<=order)

    if (!keep_all_edges){
        if ('padj' %in% colnames(grn_pruned)){
            grn_pruned <- filter(grn_pruned, order==1 | mean_padj==max(mean_padj))
        } else {
            grn_pruned <- filter(grn_pruned, order==1 | mean_estimate==max(mean_estimate))
        }
    }

    spath_graph_pruned <- spath_graph %>%
        filter(path%in%grn_pruned$path) %>%
        select(from_node, to_node, end_node, comb_dir) %>% distinct()

    grn_graph_pruned <- gene_graph %E>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        as_tibble() %>% distinct()
    grn_graph_pruned <- suppressMessages(inner_join(grn_graph_pruned, spath_graph_pruned)) %>%
        select(from_node, to_node, everything(), -from, -to) %>%
        arrange(comb_dir) %>% as_tbl_graph()

    object@grn@networks[[network]]@graphs$tf_graphs[[tf]] <- grn_graph_pruned
    return(object)
}


#' Plot sub-network centered around one TF.
#'
#' @import tidygraph
#' @import ggraph
#'
#' @param object An object.
#' @param tf The transcription factor to center around.
#' @param network Name of the network to use.
#' @param graph Name of the graph.
#' @param circular Logical. Layout tree in circular layout.
#' @param edge_width Edge width.
#' @param edge_color Edge color.
#' @param node_size Node size.
#' @param text_size Font size for labels.
#' @param label_nodes String, indicating what to label.
#' * \code{'tfs'} - Label all TFs.
#' * \code{'all'} - Label all genes.
#' * \code{'none'} - Label nothing (except the root TF).
#' @param color_edges Logical, whether to color edges by direction.
#'
#' @return A SeuratPlus object.
#'
#' @rdname plot_tf_network
#' @export
#' @method plot_tf_network SeuratPlus
plot_tf_network.SeuratPlus <- function(
    object,
    tf,
    network = DefaultNetwork(object),
    graph = 'module_graph',
    circular = TRUE,
    edge_width = 0.2,
    edge_color = c('-1'='darkgrey', '1'='orange'),
    node_size = 3,
    text_size = 10,
    label_nodes = c('tfs', 'all', 'none'),
    color_edges = TRUE
){
    label_nodes <- match.arg(label_nodes)

    gene_graph <- NetworkGraph(object, network=network, graph='tf_graphs')
    gene_graph <- gene_graph[[tf]]

    p <- ggraph(gene_graph, layout='tree', circular=circular)

    if (color_edges){
        p <- p + geom_edge_diagonal(aes(color=factor(dir)), width=edge_width) +
            scale_edge_color_manual(values=edge_color)
    } else {
        p <- p + geom_edge_diagonal(width=edge_width, color=edge_color[1])
    }

    p <- p + geom_node_point(
        color='darkgrey', shape=21, fill='lightgrey', size=node_size, stroke=0.5
    )

    net_tfs <- colnames(NetworkTFs(object))
    if (label_nodes=='tfs'){
        p <- p + geom_node_label(
            aes(label=name, filter=name%in%net_tfs),
            size=text_size/ggplot2::.pt,
            label.padding=unit(0.1, 'line')
        )
    } else if (label_nodes=='all'){
        p <- p + geom_node_label(
            aes(label=name),
            size=text_size/ggplot2::.pt,
            label.padding=unit(0.1, 'line')
        )
    } else {
        p <- p + geom_node_label(
            aes(label=name, filter=name==tf),
            size=text_size/ggplot2::.pt,
            label.padding=unit(0.1, 'line')
        )
    }
    p <- p + scale_size_continuous(range=node_size) +
        theme_void() + no_legend()

    return(p)
}



