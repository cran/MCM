\name{sim_moderate_het}
\alias{sim_moderate_het}
\docType{data}
\title{
Simulated Cross-Sectional Data of Intergenerational Social Mobility
}
\description{
A simulated cross-sectional dataset of intergenerational social mobility. 
  The mobility effects are heterogeneous (e.g., with different directions and effect sizes) and 
  the mobility scale are moderate (e.g., the ratio of the number of nonmobile and the mobile is set at 3:2).
}
\usage{data("sim_moderate_het")}
\format{
  A data frame with 30,000 observations and 6 variables.
  \describe{
    \item{\code{response}}{a numeric vector indicating the outcome variable}
    \item{\code{origin}}{a numeric vector indicating parents' socioeconomic status group}
    \item{\code{destination}}{a numeric vector indicating child' socioeconomic status group}
    \item{\code{mobility}}{a numeric vector indicating if child's socioeconomic
status differs from that of their parents'.}
  \item{\code{upmob}}{a numeric vector indicating child' socioeconomic status is higher than that of parents'.}
  \item{\code{dowmob}}{a numeric vector indicating child' socioeconomic status is lower than that of parents'.}

}
}
\examples{
data(sim_moderate_het)
head(sim_moderate_het)
}
\keyword{datasets}
