% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_bessel_layers.R
\name{simulate_bessel_layer}
\alias{simulate_bessel_layer}
\title{Bessel Layer simulation}
\usage{
simulate_bessel_layer(x, y, s, t, a)
}
\arguments{
\item{x}{start value of Brownian bridge}

\item{y}{end value of Brownian bridge}

\item{s}{start time of Brownian bridge}

\item{t}{end time of Brownian bridge}

\item{a}{vector/sequence of numbers}
}
\description{
Simulates a Bessel layer l for a given sequence a
}
\examples{
simulate_bessel_layer(x = 0, y = 0, s = 0, t = 1, a = seq(0.1, 0.5, 0.1))

}
