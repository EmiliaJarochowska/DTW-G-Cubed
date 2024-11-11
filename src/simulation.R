# Evaluation of the performance of the DTW algorithm on simulated data
# Using the two scenarios presented in Hohmann et al. 2024 

# This code uses scenario A as distributed with the StratPal package
# and scenario B in the same format, obtained using https://github.com/MindTheGap-ERC/StratPal_data_prep

library(StratPal)
library(dtw)
library(admtools)
source("src/Custom_Step_Pattern.R")

data("scenarioA") # takes 2-3 min to read in
load("data/scenarioB.RData")

extract_and_align <- function(scenario, dist1, dist2) {
  #' align two records of water depth extracted from a forward model output
  #' 
  #' The water depth is transformed from the time domain into the strat
  #' domain using admtools pointing at specific strat positions
  #' the water depth curves are compared after removing wd = 0
  #' 
  #' @param scenario an object in the format obtained from https://github.com/MindTheGap-ERC/StratPal_data_prep
  #' @param dist1 can take values "2km", "4km" ... "12km"
  #' @param dist2 can take values "2km", "4km" ... "12km"
  #' @return an object of class dtw

  adm1 = admtools::tp_to_adm(t = scenario$t_myr,
                      h = scenario$h_m[,dist1],
                      T_unit = "Myr",
                      L_unit = "m")
  adm2 = admtools::tp_to_adm(t = scenario$t_myr,
                             h = scenario$h_m[,dist2],
                             T_unit = "Myr",
                             L_unit = "m")
  
  l1 = list(t = scenario$t_myr, val = scenario$wd[,dist1])
  l2 = list(t = scenario$t_myr, val = scenario$wd[,dist2])
  source <- l1 |> admtools::time_to_strat(adm1, destructive = T) 
  target <- l2 |> admtools::time_to_strat(adm2, destructive = T) 
  
  result <- dtw::dtw(source$val[(is.na(source$h) == F & source$val != 0)], 
                    target$val[(is.na(target$h) == F & target$val != 0)], 
                    keep.internals = T, 
                    step.pattern = asymmetricP1, 
                    open.begin = F, 
                    open.end = T)
  return(result)
}

#### Comparison without noise ####

# DTW alignment
extract_and_align(scenarioA, "6km", "12km") |> plot(type="threeway")
# perfect alignment
plot(scenarioA$h_m[,"6km"], scenarioA$h_m[,"8km"])

# DTW alignment
extract_and_align(scenarioB, "6km", "10km") |> plot(type="threeway")
# perfect alignment
plot(scenarioB$h_m[,"6km"], scenarioB$h_m[,"8km"])
