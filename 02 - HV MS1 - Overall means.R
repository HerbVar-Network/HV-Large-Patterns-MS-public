# 02 - HV MS1 Overall means

# First run '01 - HV MS1 Load data and build tree.R' to load and prep data




# 1. Load packages ####
library(viridis)
library(patchwork)
library(ggExtra)
library(brms)
library(geonames)




# 2. Overall means for means and Gini ####

# Mean level of herbivory
averageMean = brm(plantMean ~ 1 + (1|species) + (1|gr(Genus_sp, cov = A)),
  data=data,
  data2=list(A = A),
  family='gaussian',
  cores=7, chains=7, iter=2000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(averageMean, 'averageMean.rds')

# Mean Gini
averageGiniAdj = brm(plantGiniAdj ~ 1 + (1|species) + (1|gr(Genus_sp, cov = A)),
  data=data,
  data2=list(A = A),
  family='beta',
  cores=16, chains=16, iter=1000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(averageGiniAdj, 'averageGiniAdj.rds')




# 3. Summarizing variability with proportions of plants ####

# What is the smallest proportion of plants with N percent of all percent herb?
mean(data$propPlantsHerb50) # only 11% of plants support 50% of all herbivory
sd(data$propPlantsHerb50) / sqrt(sum(!is.na(data$propPlantsHerb50)))

