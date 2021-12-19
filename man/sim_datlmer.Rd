\name{sim_datlmer}
\alias{sim_datlmer}
\docType{data}
\title{
Simulated Data Studying Social Mobility (Longitudinal)
}
\description{
This is a simulated data used to study social mobility under
longitudinal setting.
}
\usage{data("sim_datlmer")}
\format{
  A data frame with 50000 observations on the following 14 variables.
  \describe{
    \item{\code{id}}{an ordered factor with levels}
    \item{\code{obs}}{a numeric vector}
    \item{\code{eij}}{a numeric vector}
    \item{\code{origin}}{a factor with levels \code{1} \code{2} \code{3}}
    \item{\code{destination}}{a factor with levels \code{1} \code{2} \code{3}}
    \item{\code{origin_1}}{a numeric vector}
    \item{\code{origin_2}}{a numeric vector}
    \item{\code{origin_3}}{a numeric vector}
    \item{\code{destination_1}}{a numeric vector}
    \item{\code{destination_2}}{a numeric vector}
    \item{\code{destination_3}}{a numeric vector}
    \item{\code{yij}}{a numeric vector}
    \item{\code{age}}{a numeric vector}
    \item{\code{dir}}{a factor with levels \code{0} \code{1} \code{2}}
  }
}
\examples{
data(sim_datlmer)
## maybe str(sim_datlmer) ; plot(sim_datlmer) ...
}
\keyword{datasets}
