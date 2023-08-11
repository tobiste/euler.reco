ep_from_sc <- function(x) {
  res <- x |>
    geographical_to_cartesian2() |>
    structr::best_fit_plane()
  coords <- res$axis_c |>
    cartesian_to_geographical2()
  names(coords) <- c("lat", "lon")
  angle <- res$cone_angle * 180 / pi
  if(angle > 90) angle = 180 - angle
  misfit <- res$r_s
  names(angle) <- "angle"
  names(misfit) <- "misfit"
  c(coords, angle, misfit)
}

ep_from_gc <- function(x) {
  res <- x |>
    geographical_to_cartesian2() |>
    structr::best_fit_plane()
  coords <- res$axis_g |>
    cartesian_to_geographical2()
  names(coords) <- c("lat", "lon")
  misfit <- res$r_g
  names(misfit) <- "misfit"
  c(coords, angle = 90, misfit)
}

#' Euler pole solution for geological structures
#'
#' Finds the best Euler pole solution for a set of geographic locations that
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
#' @returns numeric vector given the latitude (in degrees), longitude (in degrees), the misfit (0 - low, 1 - high),
#' and (for small circle structures) the apical angle (in degrees) of the best fit Euler pole.
#' @export
#'
#' @examples
#' test <- smallcircle(30, 150, 35)
#' euler_solution(test)
#'
#' # Example from Price & Carmicheal (1986), GEOLOGY:
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
#' euler_solution(tintina, densify = TRUE)
#'
#' data(south_atlantic)
#' euler_solution(south_atlantic, sm = TRUE)
#' euler_solution(south_atlantic, densify = TRUE, sm = TRUE)
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
  xmean <- x_cart |>
    structr::v_mean() |>
    tectonicr::cartesian_to_geographical()
  xsd <- x_cart |>
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



#' Deviation of input data from Euler pole solution
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
#' res1 <- euler_solution(tintina)
#' data_deviation(tintina, res1)
#'
#' data(south_atlantic)
#' res2 <- euler_solution(south_atlantic)
#' data_deviation(south_atlantic, res2)
data_deviation <- function(x, ep, sm = TRUE) {
  pts <- to_geomat(x) |>
    geographical_to_cartesian2()

  epc <- tectonicr::geographical_to_cartesian(c(ep[1], ep[2]))

  # distance to Euler pole
  epdist <- mapply(FUN = ep_pts_distance, x = pts[, 1], y = pts[, 2], z = pts[, 3], MoreArgs = list(ep = epc))
  ep[3] - epdist
}

ep_pts_distance <- function(x, y, z, ep) {
  d <- tectonicr::angle_vectors(c(x, y, z), ep)
  if(d > 90) {
    180-d
  } else {
      d
    }
}

#' Distribution of deviation angles
#'
#' returns the mean deviation, its standard deviation, the variance,
#' the dispersion from zero deviation, the 95% confidence angle, the minimum
#' and maximum absolute deviation angles, and the test results of the Rayleigh
#' test against a specified mean (i.e. zero deviation).
#'
#' @param x Misfit angles in degrees
#'
#' @returns `data.frame` giving the mean misfit angle,
#' its standard deviation, the dispersion from `0` degrees, the 95% confidence angle,
#' the minimum and maximum absolute angles,
#' and the V2 test results (alternative hypothesis is an unimodal distribution around `0` degrees)
#'
#' @export
#'
#' @examples
#' data(tintina)
#' res <- euler_solution(tintina)
#' m <- data_deviation(tintina, res)
#' deviation_stats(m)
deviation_stats <- function(x) {
  #x <- tectonicr::deviation_norm(x)
  mean_dev <- tectonicr::circular_mean(x)
  if (mean_dev > 90) {
    mean_dev <- 180 - mean_dev
  }
  R = tectonicr::rayleigh_test(x, mu = 0)

  data.frame(
    mean = mean_dev,
    sd =  tectonicr::circular_sd(x),
    var = tectonicr::circular_var(x),
    disp = tectonicr::circular_dispersion(x, 0),
    CI95 = tectonicr::confidence_angle(x),
    min = min(abs(x)), max = max(abs(x)),
    Rayleigh.test = R$statistic, p.value = R$p.value
    )
}


#' Plot the Euler pole solution
#'
#' @inheritParams to_geomat
#' @param sm logical. Whether the structure described by the points `x` is
#' expected to follow small (`TRUE`) or great circle (`FALSE`) arcs?
#' @param sigdig integer. Number of digits to show coordinates and apical half angle (default 4)
#' @param omerc logical. Whether the plot should be shown in the oblique
#' Mercator projection with the Euler pole at North (`TRUE`) or not (`FALSE`, the default).
#' @param expand two element vector expand the map limits in latitude and longitude (`c(1, 1)` (degrees) by default)
#' @returns ggplot
#' @import ggplot2
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

  pts <- to_geomat(x)
  suppressWarnings(
  x2 <- sf::st_cast(x, "POINT") |>
    dplyr::mutate(deviation = deviation)
  )

  circle <- smallcircle(res[1], res[2], res[3])
  suppressMessages(
  stats <- deviation_stats(x2$deviation)
  )

  if(omerc){
    ep <- tectonicr::euler_pole(res[1], res[2])
    x <- tectonicr::geographical_to_PoR_sf(x, ep)
    x2 <-  tectonicr::geographical_to_PoR_sf(x2, ep)
    circle <- tectonicr::geographical_to_PoR_sf(circle, ep)
  }

  box <- sf::st_bbox(x)

  ggplot2::ggplot() +
    ggplot2::geom_sf(data = circle, lty = 2, color = "darkgrey") +
    ggplot2::geom_sf(data = x) +
    ggplot2::geom_sf(data = x2, aes(color = abs(deviation))) +
    ggplot2::coord_sf(xlim = c(box[1] -expand[2], box[3] + expand[2]), ylim = c(box[2]-expand[1], box[4]+expand[1])) +
    ggplot2::scale_color_viridis_c("|Deviation| (°)") +
    ggplot2::labs(
      title = "Best fit Euler pole and cone",
      subtitle = paste0(
        "Pole: ",
        signif(res[1], digits = sigdig),
        "\u00b0 (lat), ",
        signif(res[2], digits = sigdig),
        "\u00b0 (lon) | Apical half angle of cone: ",
        signif(res[3], digits = sigdig),
        "°"
      ),
      caption = paste0(
        "Misfit: ", signif(res[4], digits = 1),
        " | Dispersion of deviation from 0\u00b0: ", signif(stats$disp, digits = 1),
        " | Mean deviation: ",
        signif(stats$mean, digits = 1), "\u00b0 \u00b1 ",
        signif(stats$sd, digits = 2),
        "\u00b0 (sd.) | Rayleigh test: ",
        signif(stats$Rayleigh.test, digits = 3),
        " (p: ", signif(stats$p.value, digits = 2), ")"
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
  cbind(x=cx, y=cy, z=cz)
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
  structr::cartesian_to_acoscartesian(
    geographical_to_cartesian2(x)
  )
}

#' @rdname coordinates
acoscartesian_to_geographical <- function(x) {
  cartesian_to_geographical2(
    structr::acoscartesian_to_cartesian(x)
  )
}
