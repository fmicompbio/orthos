#' Visualize query results as a composite manhattan/density plot.
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#' @export
#'
#' @param queryResults A list containing the results of a query performed with
#'     \code{queryWithContrasts}
#' @param plot Logical scalar specifying if a plot should be generated.
#'
#' @return A composite manhattan/density plot for the scores of queries using
#' different contrast components against the respective contrast DBs.
#'
#' @examples
#' \dontrun{
#' qRES <- queryWithContrasts(MyDecomposedContrasts)
#' plotQueryResults.manh(qRES)
#' }
#'
#' @importFrom cowplot plot_grid
#'
plotQueryResultsManh <- function(queryResults, plot = TRUE) {
    .assertVector(x = queryResults, type = "list")
    stopifnot("'queryResults' must contain elements 'zscores' and 'TopHits'" = 
                  all(c("zscores", "TopHits") %in% names(queryResults)))
    .assertScalar(x = plot, type = "logical")
    CONTRASTS <- names(queryResults$TopHits)
    DATASETS <-  names(queryResults$TopHits[[1]])

    PLOTS <- list()
    for (dset in DATASETS) {
        DF <- data.frame(
            idx = seq_along(queryResults$zscores[[1]][1, ]),
            ACC = names(queryResults$zscores[[1]][1, ])
        )
        CONTR.PLOTS <- list()
        for (contrast in CONTRASTS) {
            DF <- cbind(DF, queryResults$zscores[[contrast]][dset, ])
            CONTR.PLOTS[[contrast]] <- .plotManhDens(
                queryResults$zscores[[contrast]][dset, ],
                queryResults$TopHits[[contrast]][[dset]],
                "z-score")
        }
        PLOTS[[dset]] <- cowplot::plot_grid(
            plotlist = CONTR.PLOTS, label_size = 8,
            labels = gsub("_CONTRASTS", "", CONTRASTS),
            label_x = c(0.35, 0.35, 0.35), vjust = 4,
            ncol = 3)
    }

    if (plot) {
        if (length(DATASETS) > 3) {
            warning("Too many datasets. Only the first three will be shown. ",
                    "The complete list of plots will however be returned by ",
                    "this function.", immediate. = TRUE)
        }
        combined_plot <- cowplot::plot_grid(plotlist = PLOTS[seq.int(min(length(DATASETS), 3))],
                                            label_size = 10,
                                            labels = DATASETS, ncol = 1,
                                            nrow = min(length(DATASETS), 3))
        print(combined_plot)
    }

    return(invisible(PLOTS))
}


#' Visualize query results as a composite Manhattan/Density plot
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#'
#' @param scores Numeric named vector (typically of z-scores) to use for
#'     visualization
#' @param annot Annotation dataframe for the topn results to highlight
#' @param scoreType Character scalar describing the type of \code{scores}
#'     (use to label the y axis).
#'
#' @return A composite manhattan/density plot
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#'  qRES <- queryWithContrasts(MyDecomposedContrasts)
#' .plotManhDens(qRES$zscores[[1]][1, ], qRES$TopHits[[contrast]][[dset]])
#' }
#'
#' @importFrom ggplot2 element_blank theme geom_point aes scale_fill_continuous
#'     theme_bw geom_bin2d ggplot labs
#' @importFrom ggpubr ggdensity rotate clean_theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom grid unit
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
.plotManhDens <- function(scores, annot = "", scoreType = "z-score") {
    # validate arguments
    .assertVector(x = scores, type = "numeric")
    .assertVector(x = annot, type = "DataFrame")

    # make sure scores have names
    if (is.null(names(scores))) {
        names(scores) <- as.character(seq_along(scores))
    }

    topn <- nrow(annot)
    min.score <- sort(scores, decreasing = TRUE)[topn]
    DF <- data.frame(
        idx = seq_along(scores),
        score = scores,
        ACC = names(scores)
    )

    blank.theme <-
        ggplot2::theme(
            axis.line = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            legend.position = "none",
            panel.background = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.background = ggplot2::element_blank()
        )

    expand_y <- 1.2
    dens.plot <-
        ggpubr::ggdensity(DF, "score", fill = "#33638DFF") +
        ggplot2::scale_x_continuous(limits = c(min(DF$score),
                                               max(DF$score) * expand_y)) +
        ggpubr::clean_theme() +
        ggplot2::geom_point(
            data = DF[annot$geo_accession, ],
            aes(x = .data$score, y = 0, color = annot$series_id),
            size = 1.5) +
        ggpubr::rotate() +
        ggplot2::theme(plot.margin = grid::unit(c(1, 0, 1, 0), "cm")) +
        blank.theme

    manh.plot <-
        ggplot2::ggplot(data = DF, aes(x = .data$idx, y = .data$score)) +
        ggplot2::scale_y_continuous(limits = c(min(DF$score),
                                               max(DF$score) * expand_y)) +
        ggplot2::geom_bin2d(bins = 200, na.rm = TRUE) +
        ggplot2::scale_fill_continuous(type = "viridis") +
        ggplot2::theme_bw() +
        ggplot2::geom_point(
            data = DF[annot$geo_accession, ],
            aes(color = annot$series_id),
            size = 1.5) +
        ggplot2::labs(x = "Contrast index",
                      y = paste0("Similarity (", scoreType, ")")) +
        ggrepel::geom_text_repel(
            data = DF[DF$ACC %in% annot$geo_accession, ],
            aes(x = .data$idx, y = .data$score, label = .data$ACC),
            size = 3, max.overlaps = 10000) +
        ggplot2::theme(legend.position = "none",
                       axis.ticks.x = ggplot2::element_blank(),
                       plot.margin = grid::unit(c(1, 0, 1, 0), "cm"),
                       panel.border = ggplot2::element_blank())

    P <- cowplot::plot_grid(
        manh.plot,
        NULL,
        dens.plot,
        ncol = 3,
        align = "hv",
        rel_widths = c(4, -0.2, 1))

    return(P)
}


#' Visualize query results as violin plots
#'
#' @author Panagiotis Papasaikas, Michael Stadler
#' @export
#'
#' @param queryResults A list containing the results of a query performed with
#'     \code{queryWithContrasts}
#' @param plot Logical scalar specifying if a plot should be generated.
#'
#'
#' @return A list of ggplot violin plots (one for each dataset) for the scores
#' of queries using different contrast components against the respective
#' contrast DBs.
#'
#' @examples
#' \dontrun{
#' qRES <- queryWithContrasts(MyDecomposedContrasts)
#' plotQueryResults.violin(qRES)
#' }
#'
#' @importFrom dplyr arrange group_by slice
#' @importFrom colorspace darken
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ggplot ggtitle geom_jitter geom_violin position_jitter
#'     scale_color_manual scale_fill_manual theme theme_minimal labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggsci pal_jco
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom plyr desc
#' @importFrom grDevices colorRampPalette
#'
plotQueryResultsViolin <- function(queryResults, plot = TRUE) {
    .assertVector(x = queryResults, type = "list")
    stopifnot("'queryResults' must contain elements 'zscores' and 'TopHits'" = 
                  all(c("zscores", "TopHits") %in% names(queryResults)))
    .assertScalar(x = plot, type = "logical")
    CONTRASTS <- names(queryResults$TopHits)
    DATASETS <-  names(queryResults$TopHits[[1]])
    topn <- nrow(queryResults$TopHits[[1]][[1]])
    TOP_META <- unique(as.data.frame(do.call(rbind,
                                             unlist(queryResults$TopHits))))

    PLOTS <- list()
    for (dset in DATASETS) {
        DF <- data.frame(
            idx = seq_along(queryResults$zscores[[1]][1, ]),
            ACC = names(queryResults$zscores[[1]][1, ])
        )
        for (contrast in CONTRASTS) {
            DF <- cbind(DF, queryResults$zscores[[contrast]][dset, ])
        }

        colnames(DF)[-seq_len(2)] <- gsub("_CONTRASTS", "", CONTRASTS)
        plot_df <- tidyr::pivot_longer(DF, cols = 3:ncol(DF),
                                       names_to = "COMPONENT",
                                       values_to = "score")
        plot_df$COMPONENT <- factor(plot_df$COMPONENT,
                                    levels = gsub("_CONTRASTS", "", CONTRASTS))
        # TopN highest values by COMPONENT
        plot_df2 <- plot_df |>
            dplyr::arrange(plyr::desc(.data$score)) |>
            dplyr::group_by(.data$COMPONENT) |>
            dplyr::slice(seq_len(topn))
        plot_df2$series <- TOP_META[plot_df2$ACC, "series_id"]

        mycolors <- grDevices::colorRampPalette(ggsci::pal_jco()(10))(length(
            unique(plot_df2$series)))
        mycolors <- colorspace::darken(mycolors, 0.25)
        pos <- ggplot2::position_jitter(width = 0.02, height = 0, seed = 2)

        PLOTS[[dset]] <-
            ggplot2::ggplot(plot_df, aes(.data$COMPONENT, .data$score,
                                         fill = .data$COMPONENT)) +
            ggplot2::geom_violin(trim = FALSE, size = 0.6) +
            ggplot2::scale_fill_manual(values = c("#8B451366", "#10701044",
                                                             "#FF550077")) +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none") +
            ggplot2::ggtitle(dset) +
            ggplot2::geom_jitter(
                data = plot_df2,
                aes(color = .data$series),
                size = 1.5, position = pos) +
            ggplot2::labs(x = "Contrast component",
                          y = "Similarity (z-score)") +
            ggrepel::geom_text_repel(
                data = plot_df2,
                aes(y = .data$score, label = .data$ACC, color = .data$series),
                size = 3, max.overlaps = 10000,
                position = pos) +
            ggplot2::scale_color_manual(values = mycolors)
    }

    if (plot) {
        if (length(PLOTS) > 4) {
            warning("Too many datasets. Only the first four will be shown. ",
                    "The complete list of plots will however be returned by ",
                    "this function.", immediate. = TRUE)
        }
        combined_plot <- cowplot::plot_grid(
            plotlist = PLOTS[seq_len(min(length(PLOTS), 4))])
        print(combined_plot)
    }

    return(invisible(PLOTS))
}

