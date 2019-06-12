##########################
# in the Unbounded Exact Algorithm, the 'phi function' does not have an upper bound
# in order to get an upper bound, we use 'layers' to find suitable information which establishes an interval so phi is bounded
# the following functions simulate the the 'Bessel layers' of a BB so that we can bound phi
# the first part of this script deals with the Cauchy sequences that are needed to simulate Bessel layers
# the function 'simulate_bessel_layer' simulates the Bessel layer
# the function 'simulate_layered_brownian_bridge_bessel' simulates Bessel layer intermediate points using the 'Bessel approach'
### note that there are other layering methods (eg. intersection layer approach), which requires more computation but aren't coded here
##########################

######################################## sums required for Cauchy sequences ########################################

#' Calculate interval: [S^{gamma}_{2k+1}, S^{gamma}_{2k}]
#'
#' This function calculates the interval [S^{gamma}_{2k+1}, S^{gamma}_{2k}] for given k
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start value of Brownian bridge
#' @param t end value of Brownian bridge
#' @param l lower bound of Brownian bridge
#' @param v upper bound of Brownian bridge
#' @param k integer value
#' @param previous_values defaults to NULL: but if k >= 2, then can pass the interval 
#'                        for (k-1), i.e. a vector of [S^{gamma}_{2k-1}, S^{gamma}_{2k-2}]
#'
#' @return vector of two values, S^{gamma}_{2k+1} and S^{gamma}_{2k}
#'
#' @examples
#' calc_SgammaK_intervals(x = 0, y = 0, s = 0, t = 1, l = -2, v = 1, k = 1)
#'
#' @export
calc_SgammaK_intervals <- function(x, y, s, t, l, v, k, previous_values = NULL) {
  # for given k, this function calculates and returns the interval (S_{2k+1}^{gamma}, S_{2k}^{gamma})
  # if previous values for (S_{2k+1}^{gamma}, S_{2k}^{gamma}) for (k-1) are given (VECTOR PASSED MUST BE IN THIS FORM)
  # then we don't re-calculate all sums and only add or subtract necessary functions
  # x,y,s,t,l,v,k are real numbers
  
  # defining functions sigma and phi
  sigma_bar <- function(j,l,v,x,y) exp(-(2/(t-s))*((abs(v-l)*j)+min(l,v)-x)*((abs(v-l)*j)+min(l,v)-y))
  phi_bar <- function(j,l,v,x,y) exp(-(2*j/(t-s))*(((abs(v-l)^(2))*j)+(abs(v-l)*(x-y))))
  sigma <- function(j) sigma_bar(j,l,v,x,y) + sigma_bar(j,-l,-v,-x,-y)
  phi <- function(j) phi_bar(j,l,v,x,y) + phi_bar(j,-l,-v,-x,-y)
  
  # if k = 0, then we need to calculate (S_{1}^{gamma}, S_{0}^{gamma}) which is equal to (1 - sigma(1), 1)
  if (k==0) {
    return(c(1 - sigma(1), 1))
  } 
  
  if (is.null(previous_values)) {
    # if no previous values are passed, then we calculate the sums needed
    S_2k <- 1 - sum(sapply(X = 1:k, FUN = function(j) (sigma(j) - phi(j))))
    S_2k_plus_1 <- S_2k - sigma(k+1)
    return(c(S_2k_plus_1, S_2k))
  } else if (is.vector(previous_values) & length(previous_values) == 2) {
    # first element of previous_values should be S_{2(k-1)+1}
    # second element of previous_values should be S_{2(k-1)}
    S_2k <- previous_values[2] - (sigma(k) - phi(k))
    S_2k_plus_1 <- S_2k - sigma(k+1)
    return(c(S_2k_plus_1, S_2k))
  } else {
    stop('previous_values should have the form c(S_{2k+1}^{gamma}, S_{2k}^{gamma}).')
  }
}

#' Calculate interval: [S^{delta,1}_{2k+1}, S^{delta,1}_{2k}]
#'
#' This function calculates the interval [S^{delta,1}_{2k+1}, S^{delta, 1}_{2k}] (case where min(x,y) > min)
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start value of Brownian bridge
#' @param t end value of Brownian bridge
#' @param min minimum of Brownian bridge
#' @param v upper bound of Brownian bridge
#' @param k integer value
#' @param previous_values defaults to NULL: but if k >= 2, then can pass the interval 
#'                        for (k-1), i.e. a vector of [S^{delta,1}_{2k-1}, S^{delta,1}_{2k-2}]
#'
#' @return vector of two values, S^{delta,1}_{2k+1} and S^{delta,1}_{2k}
#'
#' @examples
#' calc_SdeltaK_1_intervals(x = 0, y = 0, s = 0, t = 1, min = -2, v = 1, k = 1)
#'
#' @export
calc_SdeltaK_1_intervals <- function(x, y, s, t, min, v, k, previous_values = NULL) {
  # for given k, this function calculates and returns the interval (S_{2k+1}^{delta,1}, S_{2k}^{delta,1})
  # if previous values for (S_{2k+1}^{delta,1}, S_{2k}^{delta,1}) for (k-1) are given (VECTOR PASSED MUST BE IN THIS FORM)
  # then we don't re-calculate the sums and only add or subtract necessary functions
  # x,y,s,t,min,v,k are real numbers
  
  # defining functions sigma and phi (in the case we need to calculate from previous values)
  # we could alternatively mulitply the previous values by denom, pass it through to calc_SgammaK_intervals and convert,
  # but this requires more steps and calculations
  sigma_bar <- function(j,l,v,x,y) exp(-(2/(t-s))*((abs(v-l)*j)+min(l,v)-x)*((abs(v-l)*j)+min(l,v)-y))
  phi_bar <- function(j,l,v,x,y) exp(-(2*j/(t-s))*(((abs(v-l)^(2))*j)+(abs(v-l)*(x-y))))
  sigma <- function(j) sigma_bar(j,min,v,x,y) + sigma_bar(j,-min,-v,-x,-y)
  phi <- function(j) phi_bar(j,min,v,x,y) + phi_bar(j,-min,-v,-x,-y)
  
  # S_{k}^{delta,1} = S_{k}^{gamma} / denom
  denom <- 1 - (exp(-(2*(x-min)*(y-min))/(t-s)))
  
  # if k = 0, then we need to calculate (S_{1}^{gmma}, S_{0}^{gamma}) which is equal to (1 - sigma(1), 1) / denom
  if (k==0) {
    return(c(1 - sigma(1), 1) / denom)
  } 
  
  if (is.null(previous_values)) {
    # if no previous values are passed, then we calculate the sums needed by using calc_SgammaK_intervals
    return(calc_SgammaK_intervals(x, y, s, t, min, v, k) / denom)
  } else if (is.vector(previous_values) & length(previous_values) == 2) {
    # first element of previous_values should be S_{2(k-1)+1}
    # second element of previous_values should be S_{2(k-1)}
    S_2k <- previous_values[2] - ((sigma(k)- phi(k)) / denom)
    S_2k_plus_1 <- S_2k - (sigma(k+1) / denom)
    return(c(S_2k_plus_1, S_2k))
  } else {
    stop('previous_values should have the form c(S_{2k+1}^{delta,1}, S_{2k}^{delta,1}).')
  }
}

#' Calculate interval: [S^{delta,2}_{2k+1}, S^{delta,2}_{2k}]
#'
#' This function calculates the interval [S^{delta,2}_{2k+1}, S^{delta, 2}_{2k}] (case where min(x,y) == min)
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start value of Brownian bridge
#' @param t end value of Brownian bridge
#' @param min minimum of Brownian bridge
#' @param v upper bound of Brownian bridge
#' @param k integer value
#' @param previous_values defaults to NULL: but if k >= 2, then can pass the interval 
#'                        for (k-1), i.e. a vector of [S^{delta,2}_{2k-1}, S^{delta,2}_{2k-2}]
#'
#' @return vector of two values, S^{delta,2}_{2k+1} and S^{delta,2}_{2k}
#'
#' @examples
#' K = ceiling(sqrt((1)+(abs(1-(-2))*abs(1-(-2))))/(2*abs(1-(-2))))
#' calc_SdeltaK_2_intervals(x = -2, y = 0, s = 0, t = 0, min = -2, v = 1, k = K)
#'
#' @export
calc_SdeltaK_2_intervals <- function(x, y, s, t, min, v, k, previous_values = NULL) {
  # for given k, this function calculates and returns the interval (S_{2k+1}^{delta,2}, S_{2k}^{delta,2})
  # if previous values for (S_{2k+1}^{delta,2}, S_{2k}^{delta,2}) for (k-1) are given (VECTOR PASSED MUST BE IN THIS FORM)
  # then we don't re-calculate all sums and only add or subtract necessary functions
  # x,y,s,t,min,v,k are real numbers
  
  # defining functions psi and chi
  psi <- function(j) ((2*abs(v-min)*j)-(max(x,y)-min))*exp(-((2*abs(v-min)*j)*((abs(v-min)*j)-(max(x,y)-min)))/(t-s))
  chi <- function(j) ((2*abs(v-min)*j)+(max(x,y)-min))*exp(-((2*abs(v-min)*j)*((abs(v-min)*j)+(max(x,y)-min)))/(t-s))
  
  K <- sqrt((t-s)+(abs(v-min))^(2)) / (2*abs(v-min))
  if (k<K) stop('The given k is too small')
  
  if (is.null(previous_values)) {
    # if no previous values are passed, then we calculate the sums needed
    S_2k <- 1 - (sum(sapply(X = 1:k, FUN = function(j) (psi(j) - chi(j)))) / abs(x-y))
    S_2k_plus_1 <- S_2k - (psi(k+1)/abs(x-y))
    return(c(S_2k_plus_1, S_2k))
  } else if (is.vector(previous_values) & length(previous_values) == 2) {
    # first element of previous_values should be S_{2(k-1)+1}
    # second element of previous_values should be S_{2(k-1)}
    S_2k <- previous_values[2] - ((psi(j) - chi(j))/abs(x-y))
    S_2k_plus_1 <- S_2k - (psi(k+1)/(abs(x-y)))
    return(c(S_2k_plus_1, S_2k))
  } else {
    stop('previous_values should have the form c(S_{2k+1}^{delta,2}, S_{2k}^{delta,2}).')
  }
}

#' Calculate interval: [S^{delta}_{2k+1}, S^{delta}_{2k}]
#'
#' This function calculates the interval [S^{delta}_{2k+1}, S^{delta}_{2k}] (case where min(x,y) > min or where min(x,y) == min)
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start value of Brownian bridge
#' @param t end value of Brownian bridge
#' @param min minimum of Brownian bridge
#' @param v upper bound of Brownian bridge
#' @param k integer value
#' @param previous_values defaults to NULL: but if k >= 2, then can pass the interval 
#'                        for (k-1), i.e. a vector of [S^{delta}_{2k-1}, S^{delta}_{2k-2}]
#'
#' @return vector of two values, S^{delta}_{2k+1} and S^{delta}_{2k}
#'
#' @examples
#' # case where min(x,y ) > min
#' calc_SdeltaK_1_intervals(x = 0, y = 0, s = 0, t = 1, min = -2, v = 1, k = 1)
#'
#' # case where min(x,y) == min
#' K = ceiling(sqrt((1)+(abs(1-(-2))*abs(1-(-2))))/(2*abs(1-(-2))))
#' calc_SdeltaK_2_intervals(x = -2, y = 0, s = 0, t = 0, min = -2, v = 1, k = K)
#'
#' @export
calc_SdeltaK_intervals <- function(x, y, s, t, min, v, k, previous_values = NULL) {
  # for given k, this function calculates and returns the interval (S_{2k+1}^{delta}, S_{2k}^{delta})
  # if previous values for (S_{2k+1}^{delta}, S_{2k}^{delta}) for (k-1) are given (VECTOR PASSED MUST BE IN THIS FORM)
  # then we don't re-calculate all sums and only add or subtract necessary functions
  # x,y,s,t,min,v,k are real numbers
  
  if (min(x,y) > min) {
    return(calc_SdeltaK_1_intervals(x, y, s, t, min, v, k, previous_values = NULL))
  } else if (min(x,y) == min) {
    return(calc_SdeltaK_2_intervals(x, y, s, t, min, v, k, previous_values = NULL))
  } else {
    stop('min(x,y) < min: given minimum point is not the minimum of the BB.')
  }
}



######################################## Bessel layer simulation ########################################


#' Bessel Layer simulation
#'
#' Simulates a Bessel layer l for a given sequence a
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start time of Brownian bridge
#' @param t end time of Brownian bridge
#' @param a vector/sequence of numbers
#'
#' @examples
#' simulate_bessel_layer(x = 0, y = 0, s = 0, t = 1, a = seq(0.1, 0.5, 0.1))
#'
#' @export
simulate_bessel_layer <- function(x, y, s, t, a) {
  # given an increasing sequence of values, a, this function returns a Bessel layer,
  # for which a BB starting at x, ending at y, between [s,t] is contained
  
  # initialise values
  u <- runif(1,0,1); l <- 1; 
  
  repeat {
    k <- 2
    Sgamma_interval <- calc_SgammaK_intervals(x, y, s, t, l = min(x,y)-a[l], v = max(x,y)+a[l], k)
    while (Sgamma_interval[1] < u & u < Sgamma_interval[2]) {
      # if u is still in between (S_{2k+1}, S_{2k}), we keep looping until it is no longer in the interval
      k <- k+1
      Sgamma_interval <- calc_SgammaK_intervals(x, y, s, t, l = min(x,y)-a[l], v = max(x,y)+a[l], 
                                                k, previous_values = Sgamma_interval)
    }
    
    if (u <= Sgamma_interval[1]) {
      # i.e. u <= S_{2k+1}: return the sequence (as it might have changed) and the Bessel layer
      return(list('a' = a, 'l' = l))
    } else {
      # i.e. u >= S_{2k}: we try the next Bessel layer
      l <- l+1
    }
    
    # check that we have not reached the limit of our sequence of a: if so then we haven't accepted a layer yet
    # in this case, we extend the sequence of a, and carry on to try find a Bessel layer
    if (l > length(a)) {
      a <- c(a, a + a[length(a)])
    }
  }
}



######################################## intermediate point simulation for layered Brownian bridge ########################################

#' Layered Brownian Bridge sampler
#'
#' This function simulates a layered Brownian Bridge given a Bessel layer, at given times
#'
#' @param x start value of Brownian bridge
#' @param y end value of Brownian bridge
#' @param s start value of Brownian bridge
#' @param t end value of Brownian bridge
#' @param a vector/sequence of numbers
#' @param l integer number denoting Bessel layer, i.e. Brownian bridge is contained in [min(x,y)-a[l], max(x,y)+a[l]]
#' @param sim_times vector of real numbers to simulate Bessel bridge
#'
#' @return matrix of the simulated layered Brownian bridge path, first row is points X, second row are corresponding times
#'
#' @examples
#' # simulate Bessel layer
#' bes_layer <- simulate_bessel_layer(x = 0, y = 0, s = 0, t = 1, a = seq(0.1, 1.0, 0.1))
#' # simulate layered Brownian bridge
#' simulate_layered_brownian_bridge_bessel(x = 0, y = 0, s = 0, t = 1, a = bes_layer$a, l = bes_layer$l, sim_times = seq(0.2, 0.8, 0.2))
#'
#' @export
simulate_layered_brownian_bridge_bessel <- function(x, y, s, t, a, l, sim_times) {
  # x,y,s,t,l are real numbers
  # l is the Bessel layer that the simulated layered BB will be contained within
  # a is an increasing sequence - BB is contained within (min(x,y)-a[l], max(x,y)+a[l])
  # sim_times is a vector containing real numbers
  
  repeat {
    # initialise values
    u1 <- runif(1,0,1); u2 <- runif(1,0,1)
    
    if (u1 < 1/2) {
      # simulate minimum point and initialising values
      if (l==1 & a[l]==0) {
        # means that the simulated Bessel layer was [min(x,y), max(x,y)]
        # in this case, we cannot use bb_min, we just set the minimum to occur at min(x,y)
        # and when the minimum occured at s if x=min(x,y), or t if y=min(x,y)
        # note: there is no case that x=y, since if x=y and a[l]=0, this is a straight line, which won't be a BB
        # note: usually makes the algorithm very slow, since when simulating a Bessel bridge conditional on min point, 
        # it is not likely it will stay in this interval
        if (x < y) { m_hat = list('min' = x, 'tau' = s) } 
        else { m_hat = list('min' = y, 'tau' = t) }
        l1 <- m_hat$min; l2 <- m_hat$min;
        v1 <- max(x,y); v2 <- max(x,y)
      } else if (l==1 & a[l]!=0) {
        # means that simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        # but this is the first layer and we don't have a[l-1], so v1 = v2
        m_hat <- bb_min(x, y, s, t, a1 = min(x,y)-a[l], a2 = min(x,y))
        l1 <- m_hat$min; l2 <- m_hat$min
        v1 <- max(x,y) + a[l]; v2 <- max(x,y) + a[l]
      } else {
        # means that simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        m_hat <- bb_min(x, y, s, t, a1 = min(x,y)-a[l], a2 = min(x,y)-a[l-1])
        l1 <- m_hat$min; l2 <- m_hat$min
        v1 <- max(x,y) + a[l-1]; v2 <- max(x,y) + a[l]
      }
      
      # simulating Bessel bridge conditional on the minimum, m_hat$min, at intermediate time points
      S <- min_bes_bridge_path(x, y, s, t, min = m_hat$min, tau = m_hat$tau, times = sim_times, plot = FALSE)
      # checking that none of the simulated points are outside of the layer
      while (max(S[1,]) > (max(x,y) + a[l])) {
        S <- min_bes_bridge_path(x, y, s, t, min = m_hat$min, tau = m_hat$tau, times = sim_times, plot = FALSE)
      }
      
      # checking if BB remains in [l1,v1]
      j <- ceiling(sqrt((t-s)+(abs(v1-l1))^(2)) / (2*(abs(v1-l1))))
      # first row are the S_{2j+1}s // second row are the S_{2j}s
      Sdelta_interval <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = S[1,i], y = S[1,i+1], 
                                                                                  s = S[2,i], t = S[2,i+1], 
                                                                                  min = l1, v = v1, j))
      
      while (prod(Sdelta_interval[1,]) < u2 & u2 < prod(Sdelta_interval[2,])) {
        # if u is still in between (prod(S_{2j+1}), prod(S_{2j})), we keep looping unitl it is no longer in the interval
        j <- j + 1
        Sdelta_interval <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = S[1,i], y = S[1,i+1], 
                                                                                    s = S[2,i], t = S[2,i+1], 
                                                                                    min = l1, v = v1, j,
                                                                                    previous_values = Sdelta_interval[,i]))
      }
      
      if (u2 <= prod(Sdelta_interval[1,])) {
        # accept sample path
        return(S)
      } else if (u2 >= prod(Sdelta_interval[2,])) {
        # checking if BB remains in [l2,v2]
        k <- ceiling(sqrt((t-s)+(abs(v2-l2))^(2)) / (2*(abs(v2-l2))))
        # first row are the S_{2j+1}s // second row are the S_{2j}s
        Sdelta_interval2 <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = S[1,i], y = S[1,i+1], 
                                                                                     s = S[2,i], t = S[2,i+1], 
                                                                                     min = l2, v = v2, k))
        
        while (prod(Sdelta_interval2[1,]) < u2 & u2 < prod(Sdelta_interval2[2,])) {
          # if u is still in between (prod(S_{2j+1}), prod(S_{2j})), we keep looping unitl it is no longer in the interval
          j <- j + 1
          Sdelta_interval2 <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = S[1,i], y = S[1,i+1], 
                                                                                       s = S[2,i], t = S[2,i+1], 
                                                                                       min = l2, v = v2, j,
                                                                                       previous_values = Sdelta_interval2[,i]))
        }
        
        if (u2 <= prod(Sdelta_interval2[1,])) {
          if (runif(1,0,1) < 1/2) {
            # accept sample path
            return(S)
          } else {
            # start again
            next
          }
        } else if (u2 >= prod(Sdelta_interval2[2,])) {
          # start again
          next
        }
      }
    } 
    
    if (u1 > 1/2) {
      # simulate maximumt point and initialising values
      if (l==1 & a[l]==0) {
        # means that the simulated Bessel layer was [min(x,y), max(x,y)]
        # in this case, we cannot use bb_max, we just set the minimum to occur at min(x,y)
        # and when the minimum occured at s if x=max(x,y), or t if y=max(x,y)
        # note: there is no case that x=y, since if x=y and a[l]=0, this is a straight line, which won't be a BB
        # note: usually makes the algorithm very slow, since when simulating a Bessel bridge conditional on max point, 
        # it is not likely it will stay in this interval
        if (x < y) { m_check = list('max' = y, 'tau' = t) } 
        else { m_check = list('min' = x, 'tau' = s) }
        l1 <- min(x,y); l2 <- min(x,y);
        v1 <- m_check$max; v2 <- m_check$max
      } else if (l==1 & a[l]!=0) {
        # means that simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        # but this is the first layer and we don't have a[l-1], so l1 = l2
        m_check <- bb_max(x, y, s, t, a1 = max(x,y), a2 = max(x,y)+a[l])
        l1 <- min(x,y) - a[l]; l2 <- min(x,y) - a[l]
        v1 <- m_check$max; v2 <- m_check$max
      } else {
        # means that simulated Bessel layer was [min(x,y)-a[l], max(x,y)+a[l]]
        m_check <- bb_max(x, y, s, t, a1 = max(x,y)+a[l-1], a2 = max(x,y)+a[l])
        l1 <- min(x,y) - a[l-1]; l2 <- min(x,y) - a[l]
        v1 <- m_check$max; v2 <- m_check$max
      }
      
      # simulating Bessel bridge conditional on the maximum, m_check$max, at intermediate time points
      S <- max_bes_bridge_path(x, y, s, t, max = m_check$max, tau = m_check$tau, times = sim_times, plot = FALSE)
      # checking that none of the simulated points are outside of the layer
      while (min(S[1,]) < min(x,y) - a[l]) {
        S <- max_bes_bridge_path(x, y, s, t, max = m_check$max, tau = m_check$tau, times = sim_times, plot = FALSE)
      }

      # checking if BB remains in [l1,v1]
      j <- ceiling(sqrt((t-s)+(abs(v1-l1))^(2)) / (2*(abs(v1-l1))))
      # first row are the S_{2j+1}s // second row are the S_{2j}s
      Sdelta_interval <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = -S[1,i], y = -S[1,i+1], 
                                                                                  s = S[2,i], t = S[2,i+1], 
                                                                                  min = -v1, v = -l1, j))
      
      while (prod(Sdelta_interval[1,]) < u2 & u2 < prod(Sdelta_interval[2,])) {
        # if u is still in between (prod(S_{2j+1}), prod(S_{2j})), we keep looping unitl it is no longer in the interval
        j <- j + 1
        Sdelta_interval <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = -S[1,i], y = -S[1,i+1], 
                                                                                    s = S[2,i], t = S[2,i+1], 
                                                                                    min = -v1, v = -l1, j,
                                                                                    previous_values = Sdelta_interval[,i]))
      }
      
      if (u2 <= prod(Sdelta_interval[1,])) {
        # accept sample path
        return(S)
      } else if (u2 >= prod(Sdelta_interval[2,])) {
        # checking if BB remains in [l2,v2]
        k <- ceiling(sqrt((t-s)+(abs(v2-l2))^(2)) / (2*(abs(v2-l2))))
        # first row are the S_{2j+1}s // second row are the S_{2j}s
        Sdelta_interval2 <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = -S[1,i], y = -S[1,i+1], 
                                                                                     s = S[2,i], t = S[2,i+1], 
                                                                                     min = -v2, v = -l2, k))
        
        while (prod(Sdelta_interval2[1,]) < u2 & u2 < prod(Sdelta_interval2[2,])) {
          # if u is still in between (prod(S_{2j+1}), prod(S_{2j})), we keep looping unitl it is no longer in the interval
          j <- j + 1
          Sdelta_interval2 <- sapply(1:(ncol(S)-1), function(i) calc_SdeltaK_intervals(x = -S[1,i], y = -S[1,i+1], 
                                                                                       s = S[2,i], t = S[2,i+1], 
                                                                                       min = -v2, v = -l2, j,
                                                                                       previous_values = Sdelta_interval2[,i]))
        }
        
        if (u2 <= prod(Sdelta_interval2[1,])) {
          if (runif(1,0,1) < 1/2) {
            # accept sample path
            return(S)
          } else {
            # start again
            next
          }
        } else if (u2 >= prod(Sdelta_interval2[2,])) {
          # start again
          next
        }
      }
    }
  }
}
