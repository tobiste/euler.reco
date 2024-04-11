#' Plot the Euler pole solution (ggplot2)
#'
#' @description Creates a simple ggplot to show the best fit small/great circle,
#' the data's misfit to the solution and some statistics of the solution.
#'
#' @inheritParams to_geomat
#' @param sc logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#' @param densify.x logical. Whether the points `x` along the lines structure
#' should be densified before analysis (`TRUE`) or not (`FALSE`, the default).
#' Densification would yield equally spaced points along the line.
#' @param proj character. Whether the plot should be shown in the geographic Mercator projection (`"geo"`, the default), the oblique
#' Mercator projection with the Euler pole at North (`"omerc"`), or the stereographic projection centered in the Euler pole (`"stereo"`).
#' @param expand numeric two element vector.  expand the map limits in latitude and longitude (`c(0, 0)` (degrees) by default)
#' @param ... optional arguments passed to [smoothr::densify()] (only if `densify.x = TRUE`).
#'
#' @importFrom sf st_cast st_bbox st_geometry sf.colors
#' @importFrom dplyr mutate
#' @importFrom tectonicr euler_pole geographical_to_PoR_sf
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' data(tintina)
#' quick_plot_gg(tintina)
#' quick_plot_gg(tintina, proj = "omerc")
#'
#' data(south_atlantic)
#' quick_plot_gg(south_atlantic)
#' quick_plot_gg(south_atlantic, proj = "omerc")
#' quick_plot_gg(south_atlantic, proj = "stereo")
quick_plot_gg <- function(x, sm = TRUE, densify.x = FALSE, ..., proj = c("geo", "omerc", "stereo"), expand = c(1, 1)) {
  proj <- match.arg(proj)

  if (densify.x) {
    is.sf <- inherits(x, "sf")
    is.df <- is.data.frame(x)
    is.mat <- is.matrix(x)
    stopifnot(any(is.sf, is.df, is.mat))
    if (!is.sf) {
      if (!is.df) {
        x <- as.data.frame(x)
        colnames(x) <- c("lat", "lon")
      } else {
        x <- sf::st_as_sf(x, coords = c("lon", "lat"), crs = "WGS84")
      }
    }
    x <- smoothr::densify(x, ...)
  }

  res <- euler_solution(x, sm)
  deviation <- data_deviation(x, res)

  suppressWarnings(
    x2 <- sf::st_cast(x, "POINT") |>
      dplyr::mutate(deviation = deviation)
  )

  circle <- smallcircle(res[1], res[2], res[3])
  suppressMessages(
    stats <- deviation_stats(x2$deviation)
  )

  ep <- tectonicr::euler_pole(res[1], res[2])
  if (proj == "omerc") {
    x <- tectonicr::geographical_to_PoR_sf(x, ep)
    x2 <- tectonicr::geographical_to_PoR_sf(x2, ep)
    #circle <- tectonicr::geographical_to_PoR_sf(circle, ep)
  } else if (proj == "stereo") {
    crs2 <- ep_stereo_crs(ep)
    x <- sf::st_transform(x, crs2)
    x2 <- sf::st_transform(x2, crs2)
    circle <- sf::st_transform(circle, crs2)
  }

  box <- sf::st_bbox(x)

  ggplot2::ggplot() +
    {
      if(proj == "omerc"){
        ggplot2::geom_hline(yintercept = 90 - res[3], lty = 2, color = "darkgrey")
      } else {
        ggplot2::geom_sf(data = circle, lty = 2, color = "darkgrey")
      }
      }+
    ggplot2::geom_sf(data = x) +
    ggplot2::geom_sf(data = x2, ggplot2::aes(color = abs(deviation))) +
    ggplot2::coord_sf(xlim = c(box[1] - expand[2], box[3] + expand[2]), ylim = c(box[2] - expand[1], box[4] + expand[1])) +
    ggplot2::scale_color_viridis_c("|Deviation| (\u00B0)") +
    ggplot2::labs(
      title = "Best fit Euler pole and cone",
      subtitle = paste0(
        "Lat./Lon.: ",
        round(res[1], digits = 1),
        "\u00B0 / ",
        round(res[2], digits = 1),
        "\u00B0 | Apical half angle of cone: ",
        round(res[3], digits = 1),
        "\u00B0"
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
