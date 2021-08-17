\name{sim_moderate_het}
\alias{sim_moderate_het}
\docType{data}
\title{
Simulated Data Studying Social Mobility (Cross-Sectional)
}
\description{
This is a simulated data used to study social mobility.
In this dataset, it is assumed that there exists a moderate social mobility.
}
\usage{data("sim_moderate_het")}
\format{
  A data frame with 30,000 observations on the following 6 variables.
  \describe{
    \item{\code{response}}{a numeric vector indicating the outcome variable}
    \item{\code{origin}}{a numeric vector indicating parents' socioeconomic status}
    \item{\code{destination}}{a numeric vector indicating child' socioeconomic status}
    \item{\code{mobility}}{a numeric vector indicating if child's socioeconomic
status is diffferent from that of parents'.}
  \item{\code{upmob}}{a numeric vector indicating child' socioeconomic status is higher than that of parents'.}
  \item{\code{dowmob}}{a numeric vector indicating child' socioeconomic status is lower than that of parents'.}

}
}
\examples{
data(sim_moderate_het)
head(sim_moderate_het)
}
\keyword{datasets}
