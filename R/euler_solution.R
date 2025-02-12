#' Helper functions for euler pole solution
#'
#' @param x matrix
#' @importFrom structr best_fit_plane
#' @name ep_help
NULL

#' @rdname ep_help
#' @importFrom stats setNames
ep_from_sc <- function(x) {
  res <- x |>
    geographical_to_cartesian2() |>
    structr::best_fit_plane()
  coords <- res$axis_c |>
    cartesian_to_geographical2()
  angle <- res$cone_angle * 180 / pi
  if (angle > 90) angle <- 180 - angle

  setNames(c(coords, angle,  res$r_s), c('lat', 'lon', 'angle', 'misfit'))
}

#' @rdname ep_help
#' @importFrom stats setNames
ep_from_gc <- function(x) {
  res <- x |>
    geographical_to_cartesian2() |>
    structr::best_fit_plane()
  coords <- res$axis_g |>
    cartesian_to_geographical2()
  misfit <- res$r_g
  setNames(c(coords, 90, misfit), c('lat', 'lon', 'angle', 'misfit'))
}

#' Euler pole solution for geological structures
#'
#' @description Finds the best Euler pole solution for a set of geographic locations that
#' are aligned along a geological structure.
#'
#' @inheritParams to_geomat
#'
#' @param sc logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#'
#' @param densify.x logical. Whether the points `x` along the lines structure
#' should be densified before analysis (`TRUE`) or not (`FALSE`, the default).
#'  Densification would yield equally spaced points along the line.
#'
#' @param ... optional arguments passed to [smoothr::densify()]
#'  (only if `densify.x = TRUE`).
#'
#' @importFrom sf st_coordinates st_as_sf
#' @importFrom smoothr densify
#'
#' @returns numeric vector given the latitude (in degrees), longitude (in degrees), the misfit (0 - low, 1 - high),
#' and (for small circle structures) the apical angle (in degrees) of the best fit Euler pole.
#'
#' @export
#'
#' @details
#' Based on Gray, N.H., Geiser, P.A., Geiser, J.R. (1980). On the least-square
#' fit of small and great circles to spherically projected data. Mathematical
#' Geology, Vol. 12, No. 3, 1980
#'
#'
#' @examples
#' test <- smallcircle(30, 150, 35)
#' euler_solution(test)
#'
#' # Example from Price and Carmicheal (1986), GEOLOGY:
#'
#' # from matrix
#' rmt_mat <- rbind(
#'   "yukon" = c(66.1, -147.8),
#'   "bigbend" = c(52.25, -122.65),
#'   "washington" = c(47.85, -121.85)
#' )
#' euler_solution(rmt_mat)
#'
#' # from data.frame
#' rmt_df <- data.frame("lat" = rmt_mat[, 1], "lon" = rmt_mat[, 2])
#' euler_solution(rmt_df)
#'
#' # from sf object
#' rmt_sf <- sf::st_as_sf(rmt_df, coords = c("lon", "lat"), crs = "WGS84")
#' euler_solution(rmt_sf)
#'
#' data(tintina)
#' euler_solution(tintina)
#' euler_solution(tintina, densify.x = TRUE)
#'
#' data(south_atlantic)
#' euler_solution(south_atlantic, sc = TRUE)
#' euler_solution(south_atlantic, densify.x = TRUE, sc = TRUE)
euler_solution <- function(x, sc = TRUE, densify.x = FALSE, ...) {
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
    x_coords <- sf::st_coordinates(x)
    x_coords <- cbind(x_coords[, 2], x_coords[, 1])
  } else {
    x_coords <- to_geomat(x)
  }

  if (sc) {
    ep_from_sc(x_coords)
  } else {
    ep_from_gc(x_coords)
  }
}
