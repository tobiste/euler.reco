#' Geographical coordinate matrix
#'
#' @description Creates the coordinates in matrix format
#'
#' @param x either a `matrix` or `data.frame` containing the geographical coordinates,
#' or a `sf` object.
#'
#' @returns `matrix` of transformed coordinates
#'
#' @importFrom sf st_transform st_coordinates
to_geomat <- function(x) {
  is.sf <- inherits(x, "sf")
  is.df <- is.data.frame(x)
  is.mat <- is.matrix(x)
  stopifnot(any(is.sf, is.df, is.mat))

  if (is.sf) {
    x <- sf::st_transform(x, crs = "WGS84")
    x_coords <- sf::st_coordinates(x)
    cbind(x_coords[, 2], x_coords[, 1]) # switch columns
  } else if (is.df) {
    cbind(x$lat, x$lon) # switch columns
  } else {
    x
  }
}


#' Cartesian and geographic coordinates
#'
#' @description Convert between the two coordinate systems
#'
#' @param x `matrix` containing the coordinates
#'
#' @returns `matrix` of the transformed coordinates
#'
#' @importFrom tectonicr deg2rad rad2deg
#'
#' @name coordinates
NULL

#' @rdname coordinates
geographical_to_cartesian2 <- function(x) {
  stopifnot(is.numeric(x))
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1]) * cos(x[, 2])
  cy <- cos(x[, 1]) * sin(x[, 2])
  cz <- sin(x[, 1])
  cbind(x = cx, y = cy, z = cz)
}

#' @rdname coordinates
cartesian_to_geographical2 <- function(x) {
  stopifnot(is.numeric(x))
  r <- sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2)
  lat <- asin(x[, 3] / r)
  lon <- atan2(x[, 2], x[, 1])
  tectonicr::rad2deg(cbind(lat, lon))
}

#' #' @rdname coordinates
#' geographical_to_acoscartesian <- function(x) {
#'   structr:::cartesian_to_acoscartesian(
#'     geographical_to_cartesian2(x)
#'   )
#' }
#'
#' #' @rdname coordinates
#' acoscartesian_to_geographical <- function(x) {
#'   cartesian_to_geographical2(
#'     structr:::acoscartesian_to_cartesian(x)
#'   )
#' }
#'

#' Small circle around a given pole
#'
#' @description Creates a small circle (or spherical cone) by specifying the axis and the half apical angle
#'
#' @param lat numeric. Latitude of the Euler pole (in degrees)
#' @param lon numeric. Longitude of the Euler pole  (in degrees)
#' @param angle numeric. Angle of small circle (half apical angle of the spherical cone).
#' If `angle = 90`, the great circle is returned.
#' @param n integer. resolution of the return. Amount of points along small circle
#' (1000 by default).
#'
#' @returns `sf` object
#'
#' @importFrom dplyr summarise
#'
#' @export
#'
#' @examples
#' smallcircle(45, 10, 12.3)
smallcircle <- function(lat, lon, angle = 90, n = 1000L) {
  stopifnot(is.numeric(lat), is.numeric(lon), is.numeric(angle), is.integer(n))
  pole <- tectonicr::euler_pole(lat, lon)

  sm_np <- data.frame(
    lon = seq(-180, 180, length.out = n),
    lat = rep(90 - angle, n)
  ) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = "WGS84") |>
    summarise(do_union = FALSE) |>
    sf::st_cast("MULTILINESTRING") |>
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )

  tectonicr::PoR_to_geographical_sf(sm_np, pole) |>
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )
}


#' Stereographic projection centered in Euler pole
#'
#' @param x two-element vector or data.frame specifying the latitude and longitude of the projections center.
#'
#' @return `sf::crs()` object
#'
#' @importFrom sf st_crs
#' @importFrom tectonicr longitude_modulo
#'
#' @export
#'
#' @examples
#' ep_stereo_crs(c(45, 10))
ep_stereo_crs <- function(x) {
  if (is.numeric(x) | is.matrix(x)) {
    x <- data.frame(lat = x[1], lon = x[2])
  }

  stopifnot(is.data.frame(x))
  if (x$lat < 0) {
    x$lat <- -x$lat
    x$lon <- tectonicr::longitude_modulo(x$lon + 180)
  }

  sf::st_crs(paste0("+proj=stere +lat_0=", x$lat, " +lat_0=", x$lon))
}




npts <- function(x) {
  sf::st_coordinates(x) |>
    nrow()
}

#' Generates a grid from max and min coordinates
#'
#' @param x_range,y_range two-column vector giving the minimum and maximum coordinates. Unit as implied by `crs`
#' @param gridsize size of grid points. Unit as implied by `crs`
#' @param crs Coordinate reference system. Should be parsed by [sf::st_crs()]. Default is `4326`
#' @param ... arguments passed to [sf::st_make_grid()]
#'
#' @return `sf` object
#' @export
#'
#' @examples
#' latlon_grid()
latlon_grid <- function(x_range = c(-180, 180), y_range = c(0, 90), gridsize = 10, crs = 4326, ...){
  sf::st_bbox(
    c(
      xmin = x_range[1],
      xmax = x_range[2],
      ymin = y_range[1],
      ymax = y_range[2]
    ),
    crs = sf::st_crs(crs)
  ) |>
    sf::st_make_grid(
      cellsize = gridsize,
      what = "centers",
      offset = c(x_range[1], y_range[1]), ...
    )
}
