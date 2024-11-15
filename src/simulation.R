# Evaluation of the performance of the DTW algorithm on simulated data
# Using the two scenarios presented in Hohmann et al. 2024 

# This code uses scenario A as distributed with the StratPal package
# and scenario B in the same format, obtained using https://github.com/MindTheGap-ERC/StratPal_data_prep

library(StratPal)
library(dtw)
library(admtools)
library(astrochron)
library(DescTools)
source("src/Custom_Step_Pattern.R")

data("scenarioA") # takes 2-3 min to read in
load("data/scenarioB.RData")

standardize_geom <- function(df) {
  #' substract geometric mean and divide by geom. SD
  #' 
  #' @param df must be in the astrochron format, so second column is values
  #' @return a vector
  
  standardized <- (df[,2] - DescTools::Gmean(df[,2]))/DescTools::Gsd(df[,2])
  return(standardized)
}

extract_and_align <- function(scenario, dist1, dist2, dt1, dt2) {
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
  
  # Create a list object to project water depth values onto strat. domain
  # The list format is the best to pass to time_to_strat here
  l1 = list(t = scenario$t_myr, val = scenario$wd[,dist1])
  l2 = list(t = scenario$t_myr, val = scenario$wd[,dist2])
  source <- l1 |> admtools::time_to_strat(adm1, destructive = T) 
  target <- l2 |> admtools::time_to_strat(adm2, destructive = T) 
  
  # Make a df in the format require by linterp
  source <- do.call(cbind.data.frame, source)
  target <- do.call(cbind.data.frame, target)
  source <- source[,c(2,1)]
  target <- target[,c(2,1)]
  
  # Remove records that were emerged (no water, no deposition)
  source_val <- source[(is.na(source$h) == F & source$val != 0),]
  target_val <- target[(is.na(target$h) == F & target$val != 0),]
  
  # Interpolate and standardize
  # dt should be such that the lengths don't differ by a factor of >= 5
  source_val <- source_val %>%
    astrochron::linterp(dt = dt1,
                        verbose = F,
                        genplot = F) %>% standardize_geom()
  target_val <- target_val %>%
    astrochron::linterp(dt = dt2,
                        verbose = F,
                        genplot = F) %>% standardize_geom()
  
  result <- dtw::dtw(source_val, 
                    target_val, 
                    keep.internals = T, 
                    step.pattern = asymmetricP1, 
                    open.begin = F, 
                    open.end = T)
  return(result)
}

#### Comparison without noise ####

# DTW alignment
extract_and_align(scenarioB, 
                  "6km", 
                  "12km", 
                  dt1 = 0.4,
                  dt2 = 0.2) |> plot(type="threeway")
# perfect alignment
plot(scenarioB$h_m[,"6km"], scenarioB$h_m[,"12km"])

# DTW alignment
extract_and_align(scenarioA, 
                  "6km", 
                  "12km", 
                  dt1 = 0.4,
                  dt2 = 0.2) |> plot(type="threeway")
# perfect alignment
plot(scenarioA$h_m[,"6km"], scenarioA$h_m[,"8km"])


