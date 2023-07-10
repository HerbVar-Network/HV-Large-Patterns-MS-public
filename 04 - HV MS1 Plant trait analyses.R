# 04 - HV MS1 Plant trait analyses


# First run '01 - HV MS1 Load data and build tree.R' to load and prep data



# 1. Load packages #############################################################

# None needed here that aren't loaded in 01



# 2. Size models ###############################################################

## 2.1. Gini herbivory ####
sg = brm(plantGiniAdj ~
    sizeDiameterMeanLog +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=16, chains=16, iter=5000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))
saveRDS(sg, 'sg.rds')
r2_bayes(sg)
bayesfactor_parameters(sgb, null = 0, direction = "<")


## 2.2. Mean herbivory ####
sm = brm(plantMeanProp ~
    sizeDiameterMeanLog +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=16, chains=16, iter=5000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.985))
saveRDS(sm, 'sm.rds')
r2_bayes(sm)
bayesfactor_parameters(sm, null = 0)


## 2.3. Size with focal cover as a covariate ####
sg_density = brm(plantGiniAdj ~
    sizeDiameterMeanLog + focalPlantCoverMean +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=16, chains=16, iter=10000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))
saveRDS(sg_density, 'sg_density.rds')
r2_bayes(sg_density)
bayesfactor_parameters(sgb_density, null = 0, direction = "<")




# 3. Growth form models ############################################################

## 3.1. Gini herbivory ####
gg = brm(plantGiniAdj ~ growthForm_simp +
    (1|p|species) + (1|q|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=16, chains=16, iter=5000,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  control = list(adapt_delta = 0.98))
saveRDS(gg, 'gg.rds')
r2_bayes(gg)
bayesfactor_parameters(gg, null = 0, direction = "two-sided")

## 3.2. Mean herbivory ####
gm = brm(plantMeanProp ~ growthForm_simp +
    (1|p|species) + (1|q|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=16, chains=16, iter=5000,
  save_pars = save_pars(all = TRUE),
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99))
saveRDS(gm, 'gm.rds')
r2_bayes(gm)



