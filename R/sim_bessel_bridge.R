##########################
# these following functions simulate min and max points of Brownian bridges (BB)
# also includes functions to simulate Bessel bridges (BB conditioned on a minimum or maximum)
##########################

######################################## simulating minimum and maximum points of Brownian motion ########################################

#' Brownian bridge minimum point simulatiion
#'
#' Simulates the minimum value of a Brownian bridge
#'
#' @param x Start value of Brownian bridge
#' @param y End value of Brownian bridge
#' @param s Start time of Brownian bridge
#' @param t End time of Brownian bridge
#' @param a1 Lower bound of minimum point (a1 < a2 <= min(x,y))
#' @param a2 Upper bound of minimum point (a1 < a2 <= min(x,y))
#'
#' @return List with the simulated minimum, 'min', and time where minimum occurs, 'tau'
#'
#' @examples
#' # simulate a minimum of a Brownian bridge between 0 and 1 in time [0,1]
#' bb_min(x=0, y=0, s=0, t=1, a1=-1, a2 = 0)
#'
#' @export
bb_min <- function(x, y, s, t, a1, a2) {
  # returns a list
  # first element is the simulated minimum value of BB
  # second element is time at which it is attained
  # must make sure package 'statmod' is installed
  # a1 < a2 <= min(x,y)

  # Step 1: simulate two uniform random variables
  M_function <- function(a) exp(-2*(a-x)*(a-y)/(t-s))
  M_a1 <- M_function(a1)
  M_a2 <- M_function(a2)

  # note: if either a1 = a2 is equal to either x or y, then M_a1 = M_a2 = 1 and so u1 = 1
  # in this case, m = x or y
  u1 <- runif(1, min = M_a1, max = M_a2)
  u2 <- runif(1, min = 0, max = 1)

  # Step 2: set minimum value
  m <- x - ((sqrt((y-x)^(2)-2*(t-s)*log(u1)) - y + x) / 2)

  # Step 3: set V
  cond <- (x-m)/(x+y-2*m)
  if (u2 < cond) {
    mu <- (y-m)/(x-m)
    lambda <- ((y-m)^(2))/(t-s)
    V <- statmod::rinvgauss(1, mean = mu, shape = lambda)
  } else {
    mu <- (x-m)/(y-m)
    lambda <- ((x-m)^(2))/(t-s)
    V <- (1 / statmod::rinvgauss(1, mean = mu, shape = lambda))
  }

  # Step 4: set tau
  tau <- ((s*V)+t)/(1+V)

  return(list('min' = m, 'tau' = tau))
}

#' Brownian bridge maximum point simulatiion
#'
#' Simulates the maximum value of a Brownian bridge
#'
#' @param x Start value of Brownian bridge
#' @param y End value of Brownian bridge
#' @param s Start time of Brownian bridge
#' @param t End time of Brownian bridge
#' @param a1 Lower bound of mmaximum point (max(x,y) <= a1 < a2)
#' @param a2 Upper bound of maximum point (max(x,y) <= a1 < a2)
#'
#' @return List with the simulated maxmimum, 'max', and time where minimum occurs, 'tau'
#'
#' @examples
#' # simulate a maximum of a Brownian bridge between 0 and 1 in time [0,1]
#' bb_min(x=0, y=0, s=0, t=1, a1=0, a2 = 1)
#'
#' @export
bb_max <- function(x, y, s, t, a1, a2) {
  # returns a list
  # first element is the simulated minimum value of BB
  # second element is time at which it is attained
  # must make sure package 'statmod' is installed
  # max(x,y) <= a1 < a2

  # reflect the problem, and simulating a minimum
  sim_min <- bb_min(x = -x, y = -y, s = s, t = t, a1 = -a2, a2 = -a1)

  return(list('max' = - sim_min$min, 'tau' = sim_min$tau))
}



######################################## minimum Bessel Bridges simulation ########################################

min_bes_bridge_sim <- function(x, y, s, t, min, tau, q) {
  # returns a minimum Bessel Bridge simulation at time q given the Brownian bridge has minimum, min, at time tau

  # Step 1: set variable r
  if (q < tau) { r <- s; Wr <- x }
  else { r <- t; Wr <- y }

  # Step 2: simulate independent realisations of a Brownian bridge of unit length at time q'
  sd <- sqrt(abs(tau-q)*abs(q-r)) / abs(tau-r)
  b <- rnorm(n = 3, mean = 0, sd = sd)

  # Step 3: set simulated value
  first_term <- ((Wr-min)*abs(tau-q) / (abs(tau-r)^(3/2))) + b[1]
  W = min + sqrt(abs(tau-r) * (first_term^(2) + b[2]^(2) + b[3]^(2)))

  return(W)
}

min_bes_bridge_path <- function(x, y, s, t, min, tau, times, plot = FALSE) {
  # returns a matrix
  # first row is the values of the simulated Bessel bridge
  # fecond row is the corresponding time values

  # removing any time points where we already have values for
  times <- times[!(times %in% c(s, tau, t))]

  # splitting the time before and after the minimum value occured
  before_min_times <- times[times<=tau]
  after_min_times <- times[times>tau]

  # initialising matrices to store finite dim representation of path
  before_min <- matrix(c(x, s), nrow = 2, ncol = 1, byrow = T)
  # after_min starts at time 0, as we will flip time later on to end at t
  after_min <- matrix(c(y, 0), nrow = 2, ncol = 1, byrow = T)

  if (!(length(before_min_times)==0)) {
    for (i in 1:length(before_min_times)) {
      # simulating the first half of the Bessel Bridge
      l <- max(which(before_min[2,] <= before_min_times[i]))
      sim <- min_bes_bridge_sim(x = before_min[1,l], y,
                                s = before_min[2,l], t,
                                min, tau, q = before_min_times[i])
      before_min <- cbind(before_min, c(sim, before_min_times[i]))
    }
  }

  if (!(length(after_min_times)==0)) {
    for (i in rev(1:length(after_min_times))) {
      # simulating the second half of the Bessel bridge by reversing the time
      l <- max(which(after_min[2,] < abs(after_min_times[i]-t)))
      sim <- min_bes_bridge_sim(x = after_min[1,l], y = x,
                                s = after_min[2,l], t = abs(s-t),
                                min, tau = abs(tau-t), q = abs(after_min_times[i]-t))
      after_min <- cbind(after_min, c(sim, abs(after_min_times[i]-t)))
    }
  }

  # now need to transform the time back into the simulated time points we want and re-ordering
  after_min[2,] <- abs(after_min[2,]-t)
  after_min <- after_min[,order(after_min[2,])]

  # putting all the simulated points together to obtain simulated Bessel Bridge
  W <- cbind(before_min, c(min, tau), after_min)

  # removing any duplicates
  if (any(duplicated(W, MARGIN = 2))) {
    W <- W[, -which(duplicated(W, MARGIN = 2))]
  }

  # plot resulting simulated Bessel Bridge
  if (plot == TRUE) {
    plot(x = W[2,], y = W[1,], type = 'l', ylim = c(min(W[1,]), max(W[1,])), xlab = "Time", ylab = 'W')
    points(x = c(s, t, tau), y = c(x, y, min), pch = 4)
  }

  return(W)
}



######################################## maximum Bessel Bridges simulation ########################################

max_bes_bridge_sim <- function(x, y, s, t, max, tau, q) {
  # returns a maximum Bessel Bridge simulation at time q given the Brownian bridge has maximum at time tau

  # reflect the problem, and simulate Bessel bridge at time q, given a minimum point
  W <- min_bes_bridge_sim(x = -x, y = -y, s = s, t = t, min = -max, tau = tau, q = q)

  # reflect point on x-axis and return
  return(-W)
}

max_bes_bridge_path <- function(x, y, s, t, max, tau, times, plot = FALSE) {
  # returns a matrix
  # first row is the values of the simulated Bessel bridge
  # second row is the corresponding time values

  # reflect the problem, and simulate a Bessel bridge path for given times, given a minimum point
  W <- min_bes_bridge_path(x = -x, y = -y, s = s, t = t, min = -max, tau = tau, times = times, plot = FALSE)

  # reflect the simulated points on x-axis
  W[1,] <- - W[1,]

  # plot resulting simulated Bessel Bridge
  if (plot == TRUE) {
    plot(x = W[2,], y = W[1,], type = 'l', ylim = c(min(W[1,]), max(W[1,])), xlab = "Time", ylab = 'W')
    points(x = c(s, t, tau), y = c(x, y, max), pch = 4)
  }

  return(W)
}





######################################## alternative maximum Bessel bridge simulation ########################################

# max_bes_bridge_sim <- function(x, y, s, t, max, tau, q) {
#   # returns a maximum Bessel Bridge simulation at time q given the Brownian bridge has maximum at time tau
#
#   # Step 1: set variable r
#   if (q < tau) { r <- s; Wr <- -x }
#   else { r <- t; Wr <- -y }
#   min = -max
#
#   # Step 2: simulate independent realisations of a Brownian bridge of unit length at time q'
#   variance <- sqrt(abs(tau-q)*abs(q-r)) / (tau-r)
#   b <- rnorm(n = 3, mean = 0, sd = sqrt(variance))
#
#   # Step 3: set simulated value
#   first_term <- ((Wr-min)*abs(tau-q) / (abs(tau-r)^(3/2))) + b[1]
#   W = min + sqrt(abs(tau-r) * (first_term^(2) + b[2]^(2) + b[3]^(2)))
#
#   # Reflect point on x-axis
#   W = -W
#
#   return(W)
# }
#
# max_bes_bridge_path <- function(x, y, s, t, max, tau, times, plot = FALSE) {
#   # returns a matrix
#   # first row is the values of the simulated Bessel bridge
#   # second row is the corresponding time values
#
#   # removing any time points where we already have values for
#   remove <- c(s, tau, t)
#   times <- times[!(times %in% remove)]
#
#   # splitting the time before and after the maximum value occured
#   before_max_times <- times[times<tau]
#   after_max_times <- times[times>tau]
#
#   # initialising matrices to store finite dim representation of path
#   before_max <- matrix(c(x, s), nrow = 2, ncol = 1, byrow = T)
#   # after_max starts at time 0, as we will flip time later on to end at t
#   after_max <- matrix(c(y, 0), nrow = 2, ncol = 1, byrow = T)
#
#   if (!(length(before_max_times)==0)) {
#     for (i in 1:length(before_max_times)) {
#       # simulating the first half of the Bessel Bridge
#       l <- max(which(before_max[2,] < before_max_times[i]))
#       sim <- max_bes_bridge_sim(x = before_max[1,l], y,
#                                 s = before_max[2,l], t,
#                                 max, tau, q = before_max_times[i])
#       before_max <- cbind(before_max, c(sim, before_max_times[i]))
#     }
#   }
#
#   if (!(length(after_max_times)==0)) {
#     for (i in rev(1:length(after_max_times))) {
#       # simulating the second half of the Bessel bridge by reversing the time
#       l <- max(which(after_max[2,] < abs(after_max_times[i]-t)))
#       sim <- max_bes_bridge_sim(x = after_max[1,l], y = x,
#                                 s = after_max[2,l], t = abs(s-t),
#                                 max, tau = abs(tau-t), q = abs(after_max_times[i]-t))
#       after_max <- cbind(after_max, c(sim, abs(after_max_times[i]-t)))
#     }
#   }
#
#   # now need to transform the time back into the simulated time points we want and re-ordering
#   after_max[2,] <- abs(after_max[2,]-t)
#   after_max <- after_max[,order(after_max[2,])]
#
#   # putting all the simulated points together to obtain simulated Bessel Bridge
#   W <- cbind(before_max, c(max, tau), after_max)
#
#   # removing any duplicates
#   if (any(duplicated(W, MARGIN = 2))) {
#     W <- W[, -which(duplicated(W, MARGIN = 2))]
#   }
#
  # # plot resulting simulated Bessel Bridge
  # if (plot == TRUE) {
  #   plot(x = W[2,], y = W[1,], type = 'l', ylim = c(min(W[1,]), max(W[1,])), xlab = "Time", ylab = 'W')
  #   points(x = c(s, t, tau), y = c(x, y, max), pch = 4)
  # }
#
#   return(W)
# }
