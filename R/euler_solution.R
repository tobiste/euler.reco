ep_from_sc <- function(x) {
  # x_cart <- geographical_to_acoscartesian(x) |>
  #   acoscartesian_to_cartesian()
  res <- x |>
    structR::geographical_to_cartesian2()
  as.data.frame() |>
    structR::best_cone()
  coords <- res |>
    as.matrix() |>
    t() |>
    structR::acoscartesian_to_geographical()
  names(coords) <- c("lat", "lon")
  angle <- res[4]
  names(angle) <- ("angle")
  c(coords, angle)
}

ep_from_gc <- function(x) {
  # x_cart <- geographical_to_acoscartesian(x) |>
  #   acoscartesian_to_cartesian()
  res <- x |>
    structR::geographical_to_cartesian2()
    as.data.frame() |>
    structR::best_plane()
  coords <- res |>
    as.matrix() |>
    t() |>
    structR::acoscartesian_to_geographical()
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
#' @param sm logical. Are points `x` aligned on a small circle (`TRUE`) or great circle (`FALSE`)?
#' @importFrom sf st_coordinates
#' @return three-element vector given the latitude, longitude, and apical
#' angle of the best fit Euler pole.
#' @export
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
#' \dontrun{
#'  tf <- readRDS("C:/Users/tstephan/Documents/GitHub/cordillera-stress/data/faults.rds") %>%
#'  filter(Name %in% c("Rocky Mountain Trench", "Tintina Fault"))
#'  euler_solution(tf)
#'  }
euler_solution <- function(x, sm = TRUE) {
  x_coords <- to_geomat(x)

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
#' @return three-element vector giving the mean geographic location, variance, and standard deviation
geo_distribution <- function(x) {
  x_coords <- to_geomat(x)
  x_cart <- structR::geographical_to_cartesian2(x)
  xmean <- structR::to_vec(x_cart) |> structR::v_mean() |> structR::cartesian_to_geographical()
  xsd <- structR::to_vec(x_cart) |> structR::v_delta() |> tectonicr::rad2deg()
  c(lat = xmean[1], lon = xmean[2], delta = xsd)
}



#' Convert object to geographic coordinate matrix
#'
#' @param x a `matrix` or `data.frame` containing the geographical coordinates,
#' or a `sf` object
#'
#' @return matrix
to_geomat <- function(x){
  is.sf <- inherits(x, "sf")
  is.df <- is.data.frame(x)
  is.mat <- is.matrix(x)
  stopifnot(any(is.sf, is.df, is.mat))

  if(is.sf){
    x <- sf::st_transform(crs = "WGS84")
    x_coords <- sf::st_coordinates(x)
    cbind(x_coords[, 2], x_coords[, 1]) # switch columns
  } else if(is.df){
    cbind(x$lat, x$lon) # switch columns
  } else {
    x
  }
}

smallcircle <- function(lat, lon, angle){
  pole <- tectonicr::euler_pole(lat, lon)

  sm_np <- data.frame(
    lon = c(seq(-180, 180, 2), seq(-180, 180, 2)),
    lat = c(rep(90 - angle, 181), rep(-90 + angle, 181))
  ) %>%
    sf::st_as_sf(coords = c("lon", "lat")) %>%
    summarise(do_union = FALSE) %>%
    sf::st_cast("MULTILINESTRING") %>%
    sf::st_wrap_dateline(
      options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
      quiet = TRUE
    )

  tectonicr::PoR_to_geographical(x = sf::st_as_sf(sm_np), euler = pole) %>%
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

 tf <- readRDS("C:/Users/tstephan/Documents/GitHub/cordillera-stress/data/faults.rds") %>%
   filter(Name %in% c("Rocky Mountain Trench", "Tintina Fault"))
 euler_solution(tf)


data_misfit <- function(x, ep){
  x_coords <- to_geomat(x) |>
    as.data.frame() |>
    sf::st_as_sf(coords = c("V1", "V2"), crs = "WGS84")



}
