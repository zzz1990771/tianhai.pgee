

cosf <- function(t, p) {
  (cos(t)) ^ (p - 2)
}

zeroDelt <- function(t, eps, pp) {
  stopifnot(t <= 1)
  given <- eps * integrate(cosf, 0, pi / 2, p = pp)$value
  integrate(cosf, 0, asin(t), p = pp)$value - given
}

computeDelta <- function(eps, p) {
  #compute delta value for Singular Value upper bound. Ref: Hochstenbach (2013)
  1 / uniroot(zeroDelt, interval = c(0, 0.5), eps = eps, pp = p)$root
}

