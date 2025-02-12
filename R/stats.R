#' Distribution of Spherical Coordinates
#'
#' @description Statistics on the distribution of geographic locations
#'
#' @inheritParams to_geomat
#'
#' @returns three-element vector giving the mean geographic location, the standard deviation and the variance.
#'
#' @export
#'
#' @importFrom structr v_mean v_delta v_var
#' @importFrom tectonicr rad2deg cartesian_to_geographical
#'
#' @examples
#' xy <- cbind(lat = c(12, 14, 16), lon = c(-179, 179, 177))
#' geo_stats(xy)
geo_stats <- function(x) {
  x_cart <- to_geomat(x) |>
    geographical_to_cartesian2()
  xmean <- x_cart |>
    structr::v_mean() |>
    tectonicr::cartesian_to_geographical()
  xsd <- x_cart |>
    structr::v_delta() |>
    tectonicr::rad2deg()
  xvar <- structr::v_var(x_cart)
  c(lat = xmean[1], lon = xmean[2], delta = xsd, var = xvar)
}


#' Deviation of Input Data from Euler Pole Solution
#'
#' @description Calculates angular distance of the input points to the small or great circle
#' of the Euler pole solution
#'
#' @inheritParams to_geomat
#' @param ep three column vector. Euler pole solution
#' @param sm logical. Whether the solution is for a small circle or great circle geometries.
#'
#' @importFrom tectonicr geographical_to_cartesian
#' @returns numeric (angle in degrees). Angular difference of data points to the small or great circle.
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
  if (d > 90) {
    180 - d
  } else {
    d
  }
}

#' Distribution of deviation angles
#'
#' @description Mean deviation, standard deviation, variance, dispersion from zero deviation,
#' the 95% confidence angle, the minimum
#' and maximum absolute deviation angles, and the test results of the Rayleigh
#' test against a specified mean (i.e. zero deviation).
#'
#' @param x numeric. Misfit angles in degrees
#'
#' @returns `data.frame` giving the mean misfit angle,
#' its standard deviation, the dispersion from `0` degrees, the 95% confidence angle,
#' the minimum and maximum absolute angles,
#' and the V2 test results (alternative hypothesis is an unimodal distribution around `0` degrees)
#'
#' @export
#'
#' @importFrom tectonicr circular_mean rayleigh_test circular_sd circular_var circular_dispersion confidence_angle
#'
#' @examples
#' data(tintina)
#' res <- euler_solution(tintina)
#' m <- data_deviation(tintina, res)
#' deviation_stats(m)
deviation_stats <- function(x) {
  # x <- tectonicr::deviation_norm(x)
  mean_dev <- tectonicr::circular_mean(x)
  if (mean_dev > 90) {
    mean_dev <- 180 - mean_dev
  }
  R <- tectonicr::rayleigh_test(x, mu = 0)

  data.frame(
    mean = mean_dev,
    sd = tectonicr::circular_sd(x),
    var = tectonicr::circular_var(x),
    disp = tectonicr::circular_dispersion(x, 0),
    CI95 = tectonicr::confidence_angle(x),
    min = min(abs(x)), max = max(abs(x)),
    Rayleigh.test = R$statistic, p.value = R$p.value
  )
}
