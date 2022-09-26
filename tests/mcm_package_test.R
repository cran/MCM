
# MCM
rm(list = ls())
library(MCM)
data('sim_moderate_het')
mcm(response ~ origin * destination, data = sim_moderate_het,
    origin = "origin",destination="destination")

# # get the main effects
# fit <- mcm(response ~ origin * destination, data = sim_moderate_het,
#            origin = "origin",destination="destination",displayresult = T)
# fit$origin_main
# fit$destination_main

