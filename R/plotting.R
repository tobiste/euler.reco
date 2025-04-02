#' Plot the Euler pole solution
#'
#' @description Creates a simple ggplot to show the best fit small/great circle,
#' the data misfit to the solution and some statistics of the solution.
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
#' @importFrom graphics abline lcm legend
#'
#' @export
#'
#' @examples
#' data(tintina)
#' quick_plot(tintina)
#' quick_plot(tintina, proj = "omerc")
#'
#' data(south_atlantic)
#' quick_plot(south_atlantic)
#' quick_plot(south_atlantic, proj = "omerc")
#' quick_plot(south_atlantic, proj = "stereo")
quick_plot <- function(x, sc = TRUE, densify.x = FALSE, ..., proj = c("geo", "omerc", "stereo"), expand = c(0, 0)) {
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

  res <- euler_solution(x, sc)
  deviation <- data_deviation(x, res)

  suppressWarnings(
    x2 <- sf::st_cast(x, "POINT") |>
      dplyr::mutate(deviation = deviation, abs_deviation = abs(deviation))
  )

  circle <- smallcircle(res[1], res[2], res[3])
  suppressMessages(
    stats <- deviation_stats(x2$deviation)
  )

  ep <- tectonicr::euler_pole(res[1], res[2])
  if (proj == "omerc") {
    x <- tectonicr::geographical_to_PoR(x, ep)
    x2 <- tectonicr::geographical_to_PoR(x2, ep)
    # circle <- tectonicr::geographical_to_PoR_sf(circle, ep)
  } else if (proj == "stereo") {
    crs2 <- ep_stereo_crs(ep)
    x <- sf::st_transform(x, crs2)
    x2 <- sf::st_transform(x2, crs2)
    circle <- sf::st_transform(circle, crs2)
  }

  box <- sf::st_bbox(x) + c(-expand[2], expand[2], -expand[1], expand[1])
  breaks <- pretty(x2$abs_deviation, n = 5)
  p.value <- ifelse(stats$p.value < 0.001, "<0.001", signif(stats$p.value, digits = 2))


  title <- paste0(
    "Best fit Euler pole and cone's apical half angle\n",
    "(Lat./Lon.: ",
    round(res[1], digits = 1),
    "\u00B0 / ",
    round(res[2], digits = 1),
    "\u00B0 | angle: ",
    round(res[3], digits = 1),
    "\u00B0)"
  )
  sub <- paste0(
    "Misfit: ", signif(res[4], digits = 1),
    " | Dispersion of deviation from 0\u00B0: ", signif(stats$disp, digits = 1),
    "\nMean deviation [95% CI]: ",
    signif(stats$mean, digits = 1), "\u00B0 [",
    signif(stats$mean - stats$CI95, digits = 2),
    "\u00B0, ", signif(stats$mean + stats$CI95, digits = 2), "\u00B0] | P value: ", p.value
  )

  plot(sf::st_geometry(x), extent = box, graticule = TRUE, axes = TRUE, main = title, sub = sub)
  if (proj == "omerc") {
    abline(h = 90 - res["angle"], lty = 2, col = "seagreen", lwd = 1.5)
  } else {
    plot(sf::st_geometry(circle), lty = 2, col = "seagreen", lwd = 1.5, reset = T, extent = box, add = TRUE)
  }
  plot(x2["abs_deviation"], fill = sf::sf.colors(length(breaks)), key.pos = 1, key.size = lcm(1.3), extent = box, add = TRUE)
  legend("topright", legend = breaks, fill = sf::sf.colors(length(breaks)), title = "|Deviation| (\u00B0)", bg = 'white')
}
