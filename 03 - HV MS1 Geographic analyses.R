# 03 - HV MS1 Geographic analyses
# Latitude and biome analyses


# First run '01 - HV MS1 Load data and build tree.R' to load and prep data



# 1. Load packages #############################################################

# None needed here that aren't loaded in 01




# 2. Latitude models ###########################################################

## 2.1. Latitude Gini models ####

### 2.1.1. Latitude Gini total effect ####
lg = brm(plantGiniAdj ~
    Lat_abs_scale +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  save_pars = save_pars(all = TRUE),
  cores=7, chains=7, iter=5750*2,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(lg, 'lg.rds')
r2_bayes(lg)
bayesfactor_parameters(lg, null = 0, direction = ">")


### 2.1.2. Latitude Gini direct effect w/ mean herbivory as covariate ####
lgm = brm(plantGiniAdj ~
    Lat_abs_scale + plantMeanProp +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior=beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=7, chains=7, iter=ceiling(40000 * 2 / 7),
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(lgm, 'lgm.rds')
r2_bayes(lgm)
bayesfactor_parameters(lgm, null = 0, direction = ">",
  parameters = 'Lat_abs_scale')


## 2.2. Latitude mean herbivory model ####

lm = brm(plantMeanProp ~
    Lat_abs_scale +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  cores=7, chains=7, iter=6000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.95))
saveRDS(lm, 'lm.rds')
r2_bayes(lm)
bayesfactor_parameters(lb_meantot, null = 0, direction = "<")


## 2.3. Latitude * hemisphere models ####

lgh = brm(plantGiniAdj ~
    Lat_abs_scale + Lat_abs_scale:hemi +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  save_pars = save_pars(all = TRUE),
  cores=7, chains=7, iter=5750,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(lgh, 'lgh.rds')
bayesfactor_parameters(lgh, null = 0, parameters = 'Lat_abs_scale:hemiS')

lmh = brm(plantMeanProp ~
    Lat_abs_scale + Lat_abs_scale:hemi +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior = beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  sample_prior = 'yes',
  save_pars = save_pars(all = TRUE),
  cores=7, chains=7, iter=5750,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(lmh, 'lmh.rds')


## 2.4. Latitude covariate models ####

### 2.4.1. Latitude-size model ####
lg_size= brm(plantGiniAdj ~
    Lat_abs_scale + sizeDiameterMeanLogScale +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  data=data,
  data2=list(A = A),
  family='beta',
  cores=7, chains=7, iter=1000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(lg_size, 'lg_size.rds')
r2_bayes(lg_size)

### 2.4.2. Latitude-size model ####
lg_density= brm(plantGiniAdj ~
    Lat_abs_scale + focalPlantCoverMean +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  data=data,
  data2=list(A = A),
  family='beta',
  prior = beta_priors,
  sample_prior = 'yes',
  save_pars = save_pars(all = TRUE),
  cores=16, chains=16, iter=5000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.85))
saveRDS(lg_density, 'lg_density.rds')
r2_bayes(lg_density)
bayesfactor_parameters(lg_density, null=0, direction = '>')




# 3. Biome model ##############################################################

b = brm(mvbind(plantGiniAdj, plantMeanProp) ~
    (1|r|Biome) +
    (1|p|species) + (1|q|gr(Genus_sp, cov = A)),
  data=data, data2=list(A=A),
  family='beta',
  sample_prior = 'yes',
  save_pars = save_pars(all = TRUE),
  cores=16, chains=16, iter=5000,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.9))
saveRDS(b, 'b.rds')
r2_bayes(bbre)



