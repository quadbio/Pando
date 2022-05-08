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

#' Plot goodness-of-fit metrics.
#'
#' @import ggpointdensity patchwork
#' @return A ggplot2 object.
#'
#' @rdname plot_gof
#' @export
#' @method plot_gof SeuratPlus
plot_gof.SeuratPlus <- function(
    object,
    point_size = 0.5
){

    module_params <- NetworkModules(object)@params

    gof <- gof(object) %>%
        dplyr::filter(rsq<=1, rsq>=0) %>%
        mutate(nice=rsq>module_params$rsq_thresh&nvariables>module_params$nvar_thresh)

    models_use <- gof %>%
        filter(nice) %>%
        pull(target) %>%
        unique()

    coefs_use <- coef(object) %>%
        mutate(padj = p.adjust(pval, method = 'fdr')) %>%
        filter(target%in%models_use) %>%
        dplyr::filter(term!='(Intercept)') %>%
        mutate(
            tf_ = str_replace(term, '(.+):(.+)', '\\1'),
            peak_ = str_replace(term, '(.+):(.+)', '\\2')
        ) %>%
        mutate(
            tf = ifelse(str_detect(tf_, 'chr'), peak_, tf_),
            peak = ifelse(!str_detect(tf_, 'chr'), peak_, tf_)
        ) %>% dplyr::select(-peak_, -tf_)

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
#' @import patchwork
#' @return A ggplot2 object.
#'
#' @rdname plot_module_metrics
#' @export
#' @method plot_module_metrics SeuratPlus
plot_module_metrics.SeuratPlus <- function(
    object
){
    modules <- NetworkModules(object)@meta
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
        scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks=c(0, 1, 10, 100, 1000, 10000)) +
        theme(
            axis.line.x = element_blank(),
            axis.title.x = element_blank()
        ) +
        geom_boxplot(width=0.2, outlier.shape=NA, size=0.2) +
        labs(y=expression('# genes'))

    p_out <- p1 | p2 | p3 & no_margin()
    return(p_out)
}
