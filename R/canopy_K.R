# canopy_K.R
# Height-varying eddy diffusivity model
#   Above canopy (z > h): K(z) = k * u* * (z - d) / phi_h(z/L)   [MOST]
#   Within canopy (z <= h): K(z) = K(h) * exp(-alpha * (1 - z/h)) [Raupach attenuation]

#' MOST-based eddy diffusivity for heat/scalars
#' @param z Numeric. Height above ground (m)
#' @param u_star Numeric. Friction velocity (m/s)
#' @param L Numeric. Obukhov length (m)
#' @param disp_h Numeric. Displacement height (m)
#' @return K in m²/s
K_MOST <- function(z, u_star, L, disp_h) {
  k <- 0.4
  effective_z <- pmax(z - disp_h, 0.01)

  # Stability function for heat
  zL <- effective_z / L
  phih <- ifelse(
    is.na(L) | !is.finite(L), 1,                           # neutral
    ifelse(zL > 0, 1 + 5 * zL,                              # stable
           (1 - 16 * pmin(zL, 0))^(-0.5))                   # unstable
  )
  phih <- pmax(phih, 0.1)  # guard against blowup

  K <- k * u_star * effective_z / phih
  pmax(K, 1e-3)
}

#' Raupach (1989) canopy attenuation
#' @param z Numeric. Height above ground (m)
#' @param K_h Numeric. K at canopy top (m²/s)
#' @param canopy_h Numeric. Canopy height (m)
#' @param alpha Numeric. Decay parameter (~2-4 for forests, default 3)
#' @return K in m²/s
K_canopy_attenuate <- function(z, K_h, canopy_h, alpha = 3) {
  K_h * exp(-alpha * (1 - z / canopy_h))
}

#' Full K(z) profile combining MOST above and Raupach below canopy
#'
#' @param heights Numeric vector. Heights at which to evaluate K (m)
#' @param canopy_h Numeric. Canopy height (m). If NA or < 1, use MOST everywhere.
#' @param u_star Numeric. Friction velocity (m/s)
#' @param L Numeric. Obukhov length (m)
#' @param disp_h Numeric. Displacement height (m). Default: 0.66 * canopy_h.
#' @param alpha Numeric. Raupach decay parameter. Default: 3.
#' @return Numeric vector of K at each height (m²/s)
canopy_K_profile <- function(heights, canopy_h, u_star, L,
                              disp_h = NULL, alpha = 3) {

  if (is.na(canopy_h) || canopy_h < 1) {
    # No meaningful canopy — use MOST everywhere
    return(K_MOST(heights, u_star, L, 0))
  }

  if (is.null(disp_h) || is.na(disp_h)) disp_h <- 0.66 * canopy_h

  # K at canopy top (used as reference for attenuation below)
  K_h <- K_MOST(canopy_h, u_star, L, disp_h)

  sapply(heights, function(z) {
    if (is.na(z)) return(NA_real_)
    if (z >= canopy_h) {
      K_MOST(z, u_star, L, disp_h)
    } else {
      K_canopy_attenuate(z, K_h, canopy_h, alpha)
    }
  })
}
