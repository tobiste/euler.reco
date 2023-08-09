ep_from_sc <- function(x) {
  # x_cart <- geographical_to_acoscartesian(x) |>
  #   acoscartesian_to_cartesian()
  res <- x |>
    geographical_to_cartesian2() |>
    as.data.frame() |>
    structr::best_cone()
  coords <- res |>
    as.matrix() |>
    t() |>
    acoscartesian_to_geographical()
  names(coords) <- c("lat", "lon")
  angle <- res[4]
  names(angle) <- ("angle")
  c(coords, angle)
}

ep_from_gc <- function(x) {
  # x_cart <- geographical_to_acoscartesian(x) |>
  #   acoscartesian_to_cartesian()
  res <- x |>
    geographical_to_cartesian2()
  as.data.frame() |>
    structr::best_plane()
  coords <- res |>
    as.matrix() |>
    t() |>
    acoscartesian_to_geographical()
  names(coords) <- c("lat", "lon")
  R <- res[4]
  c(coords, R)
}

#' Euler pole solution for geological structures
#'
#' find the best Euler pole solution for a set of geographic locations that
#' are aligned along a geological structure.
#'
#' @inheritParams to_geomat
#' @param sm logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#' @param densify logical. Whether the points `x` along the lines structure
#' should be densified before analysis (`TRUE`) or not (`FALSE`, the default).
#'  Densification would yield equally spaced points along the line.
#' @param ... optional arguments passed to [smoothr::densify()]
#'  (only if `densify = TRUE`).
#' @importFrom sf st_coordinates st_as_sf
#' @importFrom smoothr densify
#' @returns three-element vector given the latitude, longitude, and apical
#' angle of the best fit Euler pole.
#' @export
#'
#' @examples
#' # Example from Price & Carmicheal (1986), GEOLOGY:
#'
#' # from matrix
#' rmt_mat <- rbind(
#'   "yukon" = c(66.1, -147.8),
#'   "bigbend" = c(52.25, -122.65),
#'   "washington" = c(47.85, -121.85)
#' )
#' # from data.frame
#' euler_solution(rmt_mat)
#'
#' rmt_df <- data.frame("lat" = rmt_mat[, 1], "lon" = rmt_mat[, 2])
#' euler_solution(rmt_df)
#'
#' # from sf object
#' rmt_sf <- sf::st_as_sf(rmt_df, coords = c("lon", "lat"), crs = "WGS84")
#' euler_solution(rmt_sf)
#'
#' data(tintina)
#' euler_solution(tintina)
#' euler_solution(tintina, densify = TRUE)
euler_solution <- function(x, sm = TRUE, densify = FALSE, ...) {
  if (densify) {
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



  if (sm) {
    ep_from_sc(x_coords)
  } else {
    ep_from_gc(x_coords)
  }
}


#' Position statistics
#'
#' Statistics on the distribution of geographic locations
#'
#' @inheritParams to_geomat
#' @returns three-element vector giving the mean geographic location, variance, and standard deviation
geo_distribution <- function(x) {
  x_coords <- to_geomat(x)
  x_cart <- geographical_to_cartesian2(x)
  xmean <- structr::to_vec(x_cart) |>
    structr::v_mean() |>
    cartesian_to_geographical()
  xsd <- structr::to_vec(x_cart) |>
    structr::v_delta() |>
    tectonicr::rad2deg()
  c(lat = xmean[1], lon = xmean[2], delta = xsd)
}



#' Convert object to geographic coordinate matrix
#'
#' @param x a `matrix` or `data.frame` containing the geographical coordinates,
#' or a `sf` object
#'
#' @returns `matrix`
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



#' Create a specified small circle around a given pole
#'
#' @param lat numeric. Latitude of the Euler pole (in degrees)
#' @param lon numeric. Longitude of the Euler pole  (in degrees)
#' @param angle numeric. Angle of small circle (if `angle = 90`, the great circle is returned)
#' @param n resolution of the return. Amount of points along small circle
#' (1000 by default).
#'
#' @returns `sf` object
#' @importFrom sf st_as_sf st_cast st_wrap_dateline
#' @importFrom dplyr summarise
#' @importFrom tectonicr euler_pole PoR_to_geographical_sf
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


# library(ggplot2)
#
# rmt.res <- euler_solution(rmt_mat)
# ep <- data.frame(lat = rmt.res[1], lon = rmt.res[2])
#
# ggplot()+
#   geom_sf(data = rmt_sf) +
#   geom_sf(data = tectonicr::eulerpole_smallcircles(ep, 50)) +
#   coord_sf(xlim = c(sf::st_bbox(rmt)[1], sf::st_bbox(rmt)[3]), ylim = c(sf::st_bbox(rmt)[2], sf::st_bbox(rmt)[4]))



#' Misfit of input data from Euler pole solution
#'
#' Calculates angular distance of the input points to the small or great circle
#' of the Euler pole solution
#'
#' @inheritParams to_geomat
#' @param ep three column vector. Euler pole solution
#' @param sm logical.
#'
#' @importFrom tectonicr geographical_to_cartesian
#' @returns numeric (angle in degrees)
#' @export
#'
#' @examples
#' data(tintina)
#' res <- euler_solution(tintina)
#' data_misfit(tintina, res)
data_misfit <- function(x, ep, sm = TRUE) {
  pts <- to_geomat(x) |>
    geographical_to_cartesian2()


  epc <- tectonicr::geographical_to_cartesian(c(ep[1], ep[2]))

  # distance to Euler pole
  epdist <- mapply(FUN = ep_pts_distance, x = pts[, 1], y = pts[, 2], z = pts[, 3], MoreArgs = list(ep = epc))
  ep[3] - epdist
}

ep_pts_distance <- function(x, y, z, ep) {
  tectonicr::angle_vectors(c(x, y, z), ep)
}

#' Distribution of misfit angles
#'
#' @param x Misfit angles in degrees
#'
#' @returns three column vector giving the mean misfit angle,
#' its standard deviation, the dispersion, the 95% confidence angle, and the minimum and maximum absolute angles
#'
#' @export
#'
#' @examples
#' data(tintina)
#' res <- euler_solution(tintina)
#' m <- data_misfit(tintina, res)
#' misfit_stats(m)
misfit_stats <- function(x) {
  #x <- tectonicr::deviation_norm(x)
  mean_misfit <- tectonicr::circular_mean(x)
  if (mean_misfit > 90) {
    mean_misfit <- 180 - mean_misfit
  }
  sd_misfit <- tectonicr::circular_sd(x)
  disp <- tectonicr::circular_dispersion(x)
  conf.angle <- tectonicr::confidence_angle(x)
  c(mean = mean_misfit, sd = sd_misfit, disp = disp, CI95 = conf.angle, min = min(abs(x)), max = max(abs(x)))
}


#' Plot the Euler pole solution
#'
#' @inheritParams to_geomat
#' @param sm logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#' @param sigdig integer indicating the number of decimal places (round) to be used.
#' @param ... optional arguments passed to [euler_solution()]
#' @returns ggplot
#' @export
#'
#' @examples
#' data(tintina)
#' quick_plot(tintina)
quick_plot <- function(x, sm = TRUE, sigdig = 1, ...) {
  res <- euler_solution(x, ...)
  misf <- data_misfit(x, res)

  pts <- to_geomat(x)
  x2 <- sf::st_cast(x, "POINT") |>
    dplyr::mutate(misfit = misf)

  circle <- smallcircle(res[1], res[2], res[3])

  stats <- misfit_stats(x2$misfit)

  ggplot() +
    geom_sf(data = circle, lty = 2, color = "grey") +
    geom_sf(data = x) +
    geom_sf(data = x2, aes(color = abs(misfit))) +
    coord_sf(xlim = c(sf::st_bbox(x)[1], sf::st_bbox(x)[3]), ylim = c(sf::st_bbox(x)[2], sf::st_bbox(x)[4])) +
    scale_color_viridis_c("|Misfit| (°)") +
    labs(
      title = "Best fit Euler pole and cone",
      subtitle = paste0(
        "Pole: ",
        round(res[1], sigdig),
        "° (lat), ",
        round(res[2], sigdig),
        "° (lon) | Apical half angle of cone: ",
        round(res[3], sigdig),
        "°"
      ),
      caption = paste0(
        "GOF: ", round(stats[3], 3),
        " | Mean misfit: ",
        round(stats[1], 2), "° ± ", round(stats[2], 2), "° (sd.)"
      )
    )
}


#' Coordinate conversion
#'
#' @param x `matrix`
#'
#' @returns `matrix`
#' @name coordinates
NULL

#' @rdname coordinates
geographical_to_cartesian2 <- function(x) {
  stopifnot(is.numeric(x))
  x <- tectonicr::deg2rad(x)
  cx <- cos(x[, 1]) * cos(x[, 2])
  cy <- cos(x[, 1]) * sin(x[, 2])
  cz <- sin(x[, 1])
  cbind(cx, cy, cz)
}

#' @rdname coordinates
cartesian_to_geographical2 <- function(x) {
  stopifnot(is.numeric(x))
  r <- sqrt(x[, 1]^2 + x[, 2]^2 + x[, 3]^2)
  lat <- asin(x[, 3] / r)
  lon <- atan2(x[, 2], x[, 1])
  tectonicr::rad2deg(cbind(lat, lon))
}

#' @rdname coordinates
geographical_to_acoscartesian <- function(x) {
  structr:::cartesian_to_acoscartesian(
    geographical_to_cartesian2(x)
  )
}

#' @rdname coordinates
acoscartesian_to_geographical <- function(x) {
  cartesian_to_geographical2(
    structr:::acoscartesian_to_cartesian(x)
  )
}
