#' Stress Dispersion for Given Euler Pole Coordinates
#'
#' @param x \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param grid grid given as `sfc` object
#' @param prd numeric. The theoretical orientation of \eqn{\sigma_{Hmax}}{SHmax} in PoR coordinates
#' @param prob Upper quantile threshold for circular dispersion to remove outliers. Default is 0.75
#'
#' @return SpatRaster
#' @export
#' @importFrom terra as.data.frame
#' @examples
#' data("san_andreas", package = "tectonicr")
#' euler_stress_dispersion(san_andreas, grid = latlon_grid(gridsize = 10), prd = 135)
euler_stress_dispersion <- function(x, grid, prd, prob = .75) {
  grid_df <- grid |>
    sf::st_coordinates() |>
    as.data.frame() |>
    dplyr::rename(lon = X, lat = Y)

  if (!"unc" %in% colnames(x)) {
    x$unc <- rep(1, nrow(x))
  }
  weights <- 1 / x$unc


  d <- numeric(nrow(grid_df))
  for (i in seq_along(grid_df$lon)) {
    azi <- tectonicr::PoR_shmax(x, grid_df[i, ])
    cd <- tectonicr::circular_distance(azi, prd)

    upp <- stats::quantile(cd, probs = prob)

    azi_in <- azi[cd >= upp]
    w_in <- weights[cd >= upp]
    d[i] <- tectonicr::circular_dispersion(azi_in, y = prd, w = w_in)
  }

  grid_df |>
    dplyr::bind_cols(dispersion = d) |>
    terra::rast(crs = terra::crs(sf::st_as_sf(grid)))
}

#' Extract Best Euler Pole from Dispersion Grid
#'
#' @param x SpatRaster
#' @param antipode logical. If `TRUE`, the antipode of the Euler pole is returned
#'
#' @return `sf` POINT object
#' @importFrom terra as.data.frame crs
#' @importFrom dplyr filter mutate
#' @importFrom sf st_as_sf
#' @importFrom tectonicr longitude_modulo
#' @export
#'
#' @examples
#' data("san_andreas", package = "tectonicr")
#' euler_stress_dispersion(san_andreas, grid = latlon_grid(gridsize = 10), prd = 135) |>
#'   extract_best_ep_from_grid()
extract_best_ep_from_grid <- function(x, antipode = FALSE) {
  ep <- terra::as.data.frame(x, xy = TRUE) |>
    dplyr::filter(dispersion == min(dispersion))
  if (antipode) {
    ep <- dplyr::mutate(ep,
      x = tectonicr::longitude_modulo(x + 180),
      y = -y
    )
  }
  sf::st_as_sf(ep, coords = c("x", "y"), crs = terra::crs(x))
}

#' Refine Euler Pole Solution by Densification
#'
#' @param x \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param grid grid given as `sfc` object
#' @param dispersion.threshold numeric. Threshold for dispersion
#' @param fact positive integer. Aggregation factor to increase number of cells in each direction (horizontally and vertically.
#' @param ... optional arguments passed to [euler_stress_dispersion()]
#'
#' @return SpatRast
#' @importFrom terra disagg crs
#' @importFrom sf st_as_sf st_as_sfc
refine_ep_dispersion <- function(x, grid, dispersion.threshold, fact = 2, ...) {
  grid2 <- grid |>
    tidyterra::filter(dispersion <= dispersion.threshold) |>
    terra::disagg(fact = fact)
  grid3 <- grid2 |>
    terra::as.data.frame(xy = TRUE) |>
    sf::st_as_sf(coords = c("x", "y"), crs = terra::crs(grid)) |>
    sf::st_as_sfc()
  euler_stress_dispersion(x, grid3, ...)
}

#' Best Euler Pole Solution for Shmax Data
#'
#' @description Iteratively refines the Euler pole solution for a given set of
#' \eqn{\sigma_{Hmax}}{SHmax} data by
#' densifying the grid and selecting the Euler pole with the lowest dispersion.
#'
#' @param x \code{data.frame} containing the coordinates of the point(s)
#' (\code{lat}, \code{lon}), the direction of
#' \eqn{\sigma_{Hmax}}{SHmax} \code{azi} and its standard deviation
#' \code{unc} (optional)
#' @param iter positive integer. Number of iterations to refine the Euler pole solution
#' @param fact positive integer. Aggregation factor to increase number of cells in each direction (horizontally and vertically, applied in each iteration.
#' @param ... optional arguments passed to [euler_stress_dispersion()]
#'
#' @return `sf` POINT object
#' @export
#'
#' @examples
#' \dontrun{
#' data("san_andreas", package = "tectonicr")
#' euler_from_stress(san_andreas, prd = 135)
#' }
euler_from_stress <- function(x, iter = 3, fact = 2, ...) {
  message(paste("Iteration 1"))
  start_grid <- latlon_grid(gridsize = 10)
  res <- euler_stress_dispersion(x, grid = start_grid, ...)
  i <- 2
  while (i <= (iter - 1)) {
    message(paste("Iteration", i))
    meddisp <- stats::median(terra::values(res), na.rm = TRUE)

    res <- refine_ep_dispersion(x, grid = res, dispersion.threshold = meddisp, fact = fact, ...)
    i <- i + 1
  }
  extract_best_ep_from_grid(res)
}
