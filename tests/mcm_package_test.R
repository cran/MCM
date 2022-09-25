
# MCM
rm(list = ls())
library(MCM)
data('sim_moderate_het')
mcm(response ~ origin * destination, data = sim_moderate_het,
    origin = "origin",destination="destination")

