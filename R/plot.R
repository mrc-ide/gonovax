format.gonovax_grid <- function(grid) {
  data.frame(eff = grid$inputs$grid$eff * 100,
             dur = grid$inputs$grid$dur,
             a   = colMeans(grid$red_incid),
             b   = colMeans(grid$cum_vaccinated),
             c   = colMeans(grid$cum_red_incid),
             d   = colMeans(grid$cum_vaccinated / grid$cum_red_incid))
}

plot.gonovax_grid <- function(grid) {
  x <- format(grid)
  plot_heatmaps(x)
}

plot_heatmap <- function(x, what, title = "") {
  ggplot2::ggplot(x, ggplot2::aes(x = eff, y = dur, fill = .data[[what]])) +
    ggplot2::geom_tile() +
    ggplot2::xlab("Efficacy (%)") +
    ggplot2::ylab("Duration (years)") +
    ggplot2::ggtitle(title) +
    ggplot2::theme_grey() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::theme(plot.margin = ggplot2::unit(rep(0.5, 4), "cm"))
}

plot_heatmaps <- function(x) {
  p <- mapply(plot_heatmap,
              what = c("a", "b", "c", "d"),
              title = c("A) Reduction in incidence after 10 years",
                        "B) Courses of vaccine over 10 years",
                        "C) Infections averted over 10 years",
                        "D) Courses of vaccine per infection averted"),
              MoreArgs = list(x = x),
              SIMPLIFY = FALSE)
  gridExtra::grid.arrange(grobs = p)
}
