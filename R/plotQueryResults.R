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
#' @importFrom ggplot2 element_blank theme geom_point aes scale_fill_continuous
#'     theme_bw geom_bin2d
#' @importFrom ggpubr ggdensity rotate rremove clean_theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom cowplot plot_grid
#' @importFrom grid unit
#' 
plotQueryResults <- function(scores, topn = 10, annot = "") {
  min.score <- sort(scores, decreasing = TRUE)[topn]
  DF <- data.frame(
    idx = 1:length(scores),
    score = scores,
    ACC = names(scores),
    annot = annot
  )
  
  blank.theme <-
    ggplot2::theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank()
    )
  
  dens.plot <-
    ggpubr::ggdensity(DF, "score", fill = "#33638DFF") + 
    ggpubr::clean_theme() +
    geom_point(
      data = DF %>% filter(score >= min.score),
      aes(x = score, y = 0, color = annot),
      size = 1.5,
    ) +
    ggpubr::rotate() +  
    ggplot2::theme(plot.margin = grid::unit(c(1, 0, 1, 0), "cm")) +
    blank.theme
  
  manh.plot <-
    ggplot(data = DF, aes(x = idx, y = score)) + 
    ggplot2::geom_bin2d(bins = 200) +
    ggplot2::scale_fill_continuous(type = "viridis") + 
    ggplot2::theme_bw() +
    ggplot2::geom_point(
      data = subset(DF, score >= min.score),
      aes(color = annot),
      size = 1.5
    ) +
    ggrepel::geom_text_repel(
      data = subset(DF, score >= min.score),
      aes(x = idx, y = score , label = ACC),
      size = 3
    ) +
    ggpubr::rremove("x.ticks") + 
    ggplot2::theme(legend.position = "none", 
                   plot.margin = grid::unit(c(1, 0, 1, 0), "cm"),
                   panel.border = ggplot2::element_blank())
  
  P <- cowplot::plot_grid(
    manh.plot,
    NULL,
    dens.plot,
    ncol = 3,
    align = "hv",
    rel_widths = c(4, -0.2, 1)
  )
  
  return(P)
}
