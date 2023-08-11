#' Plot the Euler pole solution
#'
#' @description Creates a simple ggplot to show the best fit small/great circle,
#' the data's misfit to the solution and some statistics of the solution.
#'
#' @inheritParams to_geomat
#' @param sm logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#' @param sigdig integer. Number of digits to show coordinates and apical half angle (default 4)
#' @param omerc logical. Whether the plot should be shown in the oblique
#' Mercator projection with the Euler pole at North (`TRUE`) or not (`FALSE`, the default).
#' @param expand numeric two element vector.  expand the map limits in latitude and longitude (`c(1, 1)` (degrees) by default)
#'
#' @returns ggplot
#'
#' @importFrom ggplot2 ggplot geom_sf scale_color_viridis_c labs aes
#' @importFrom sf st_cast st_bbox
#' @importFrom dplyr mutate
#' @importFrom tectonicr euler_pole geographical_to_PoR_sf
#'
#' @export
#'
#' @examples
#' data(tintina)
#' quick_plot(tintina)
#' quick_plot(tintina, omerc = TRUE)
#'
#' data(south_atlantic)
#' quick_plot(south_atlantic, sm = FALSE)
#' quick_plot(south_atlantic, sm = TRUE)
#' quick_plot(south_atlantic, sm = TRUE, omerc = TRUE)
quick_plot <- function(x, sm = TRUE, sigdig = 4, omerc = FALSE, expand = c(1, 1)) {
  res <- euler_solution(x, sm)
  deviation <- data_deviation(x, res)

  # pts <- to_geomat(x)
  suppressWarnings(
    x2 <- sf::st_cast(x, "POINT") |>
      dplyr::mutate(deviation = deviation)
  )

  circle <- smallcircle(res[1], res[2], res[3])
  suppressMessages(
    stats <- deviation_stats(x2$deviation)
  )

  if (omerc) {
    ep <- tectonicr::euler_pole(res[1], res[2])
    x <- tectonicr::geographical_to_PoR_sf(x, ep)
    x2 <- tectonicr::geographical_to_PoR_sf(x2, ep)
    circle <- tectonicr::geographical_to_PoR_sf(circle, ep)
  }

  box <- sf::st_bbox(x)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = circle, lty = 2, color = "darkgrey") +
    ggplot2::geom_sf(data = x) +
    ggplot2::geom_sf(data = x2, ggplot2::aes(color = abs(deviation))) +
    ggplot2::coord_sf(xlim = c(box[1] - expand[2], box[3] + expand[2]), ylim = c(box[2] - expand[1], box[4] + expand[1])) +
    ggplot2::scale_color_viridis_c("|Deviation| (\u00B0)") +
    ggplot2::labs(
      title = "Best fit Euler pole and cone",
      subtitle = paste0(
        "Pole: ",
        signif(res[1], digits = sigdig),
        "\u00B0 (lat), ",
        signif(res[2], digits = sigdig),
        "\u00B0 (lon) | Apical half angle of cone: ",
        signif(res[3], digits = sigdig),
        "\u00B"
      ),
      caption = paste0(
        "Misfit: ", signif(res[4], digits = 1),
        " | Dispersion of deviation from 0\u00B0: ", signif(stats$disp, digits = 1),
        " | Mean deviation: ",
        signif(stats$mean, digits = 1), "\u00B0 \u00B1 ",
        signif(stats$sd, digits = 2),
        "\u00B0 (sd.) | Rayleigh test: ",
        signif(stats$Rayleigh.test, digits = 3),
        " (p: ", signif(stats$p.value, digits = 2), ")"
      )
    )
}
