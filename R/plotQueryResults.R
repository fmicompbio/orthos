#' Plot query results
#' 
#' @author Panagiotis Papasaikas
#' @export
#' 
#' @param scores Numeric named vector (typically of z-scores) to use for
#'  visualization
#' @param topn Highlight the topn results
#' @param annot Annotation vector used for coloring the topn results
#'
#' @return A plot
#' 
#' @examples 
#' x <- rnorm(100)
#' plotQueryResults(scores = x,
#'                  topn = 10,
#'                  annot = as.character(round(abs(x), 1)))
#' 
#' @importFrom ggplot2 element_blank theme geom_point aes scale_fill_continuous
#'     theme_bw geom_bin2d ggplot
#' @importFrom ggpubr ggdensity rotate clean_theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom grid unit
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' 
plotQueryResults <- function(scores, topn = 10, annot = "") {
    # validate arguments
    .assertVector(x = scores, type = "numeric")
    .assertScalar(x = topn, type = "numeric", rngIncl = c(1, length(scores)))
    .assertVector(x = annot)
    
    # make sure scores have names
    if (is.null(names(scores))) {
        names(scores) <- as.character(seq_along(scores))
    }
    
    min.score <- sort(scores, decreasing = TRUE)[topn]
    DF <- data.frame(
        idx = seq_along(scores),
        score = scores,
        ACC = names(scores),
        annot = annot
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
    
    dens.plot <-
        ggpubr::ggdensity(DF, "score", fill = "#33638DFF") + 
        ggpubr::clean_theme() +
        ggplot2::geom_point(
            data = DF |> dplyr::filter(.data$score >= min.score),
            aes(x = .data$score, y = 0, color = annot),
            size = 1.5) +
        ggpubr::rotate() +  
        ggplot2::theme(plot.margin = grid::unit(c(1, 0, 1, 0), "cm")) +
        blank.theme
    
    manh.plot <-
        ggplot2::ggplot(data = DF, aes(x = .data$idx, y = .data$score)) + 
        ggplot2::geom_bin2d(bins = 200) +
        ggplot2::scale_fill_continuous(type = "viridis") + 
        ggplot2::theme_bw() +
        ggplot2::geom_point(
            data = subset(DF, score >= min.score),
            aes(color = annot),
            size = 1.5) +
        ggrepel::geom_text_repel(
            data = subset(DF, score >= min.score),
            aes(x = .data$idx, y = .data$score , label = .data$ACC),
            size = 3) +
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
