# 06 - HV MS1 Sensitivity analyses

# First run '01 - HV MS1 Load data and build tree.R' to load and prep data



# 0. Load packages #############################################################

library(purrr)
library(broom)
library(broom.mixed)
library(lubridate)




# 1. Sample size sensitivity analyses ##########################################

dp = dp %>%
  group_by(surveyID) %>%
  mutate(
    totalPlantsPercHerbPlant = sum(!is.na(percHerbPlant))
  )


## Subsample dataset to N plants per survey ####

# SUBSET plant data dataframes to surveys with >= 30 plants
dp30 = dp %>%
  drop_na(percHerbPlant) %>%
  subset(totalPlantsPercHerbPlant >= 30) %>%
  group_by(surveyID)

# Summarize subsetted plant data dataframes to survey level
d30nonsub = dp30 %>%
  summarize(
    plantMean = mean(percHerbPlant),
    plantGini = Gini(percHerbPlant),
    Lat = first(Lat),
    Lat_abs = abs(first(Lat)),
    plants = n(),
    species = first(Genus_sp),
    Genus_sp = first(Genus_sp)
  ) %>%
  mutate(
    plantGiniAdj = pmax(pmin(plantGini, 0.99), 0.01),
    Lat_abs_scale = scale(Lat_abs, center=TRUE, scale=TRUE)
  ) %>%
  drop_na(plantGini)

# SUBSAMPLE plant data dataframes to 30 plants per survey and summarize them
d30l = list()
for(i in 1:100){
  dp30.temp =  slice_sample(dp30, n=30, replace=FALSE)
  d30l[[i]] = dp30.temp %>%
    summarize(
      plantMean = mean(percHerbPlant),
      plantGini = Gini(percHerbPlant),
      Lat = first(Lat),
      Lat_abs = abs(first(Lat)),
      plants = n(),
      species = first(Genus_sp),
      Genus_sp = first(Genus_sp)
    ) %>%
    mutate(
      bootRep = i,
      plantGiniAdj = pmax(pmin(plantGini, 0.99), 0.01),
      Lat_abs_scale = scale(Lat_abs, center=TRUE, scale=TRUE),
      sizeDiameterMeanLog =
        data$sizeDiameterMeanLog[match(surveyID, data$surveyID)]
    ) %>%
    drop_na(plantGini)
}
saveRDS(d30l, 'd30l.RDS')


# Bind dataframes in list into one big dataframe
d30 = d30l %>%
  bind_rows() %>%
  group_by(bootRep)



## Calculate summary stats ####
d30b = d30 %>%
  summarize(
    plantGiniMean = mean(plantGini, na.rm=T),
    plantMeanMean = mean(plantMean)
  )


### Mean herbivory ####
mean(data$plantMean, na.rm=TRUE)
sd(data$plantMean, na.rm=TRUE) / sqrt(sum(!is.na(data$plantMean)))

mean(d30b$plantMeanMean, na.rm=TRUE)
sd(d30b$plantMeanMean, na.rm=TRUE) / sqrt(sum(!is.na(d30b$plantMeanMean)))

averageMean = readRDS('averageMean.rds')


### Gini herbivory ####
averageGiniAdj = readRDS('averageGiniAdj.rds')

averageGiniAdj30 = list()
for(i in 1:100){
  averageGiniAdj30[[i]] = brm(plantGiniAdj ~ 1 +
      (1|species) + (1|gr(Genus_sp, cov = A)),
    prior=beta_priors_noslope,
    data=d30l[[i]],
    data2=list(A = A),
    family='beta',
    cores=16, chains=16, iter=800,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.85))
  print(i)
}
saveRDS(averageGiniAdj30, 'averageGiniAdj30.rds')
averageGiniAdj30 = readRDS('averageGiniAdj30.rds')

averageGiniAdj30tidy = averageGiniAdj30 %>%
  purrr::map_dfr(broom.mixed::tidy) %>%
  subset(term == '(Intercept)')

# Output summary
averageGiniAdj
plogis(.44)

mean(plogis(averageGiniAdj30tidy$estimate))
hist(plogis(averageGiniAdj30tidy$estimate))
quantile(plogis(averageGiniAdj30tidy$estimate), c(0.025, 0.975))




## Latitude models ####

lb_tot = readRDS('lg.rds')

lat30nonsub = brm(plantGiniAdj ~
    Lat_abs_scale +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior=beta_priors,
  data=d30nonsub,
  data2=list(A = A),
  family='beta',
  cores=12, chains=12, iter=800,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.85))
r2_bayes(lat30nonsub)

latd = brm(plantGiniAdj ~
    Lat_abs_scale +
    (1|species) + (1|gr(Genus_sp, cov = A)),
  prior=beta_priors,
  data=data,
  data2=list(A = A),
  family='beta',
  cores=12, chains=12, iter=800,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.85))
r2_bayes(latd)

lat30 = list()
for(i in 1:100){
  lat30[[i]] = brm(plantGiniAdj ~
      Lat_abs_scale +
      (1|species) + (1|gr(Genus_sp, cov = A)),
    prior=beta_priors,
    data=d30l[[i]],
    data2=list(A = A),
    family='beta',
    cores=12, chains=12, iter=600,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.85))
  print(i)
}
saveRDS(lat30, 'lat30.rds')
lat30 = readRDS('lat30.rds')

lat30tidy = lat30 %>%
  purrr::map_dfr(broom.mixed::tidy) %>%
  subset(term == 'Lat_abs_scale') %>%
  mutate(
    sig = ifelse(conf.low > 0, 1, 0), # is lat effect CI above zero?
  )
lat30tidy$R2_marg = sapply(lat30, function(x) r2_bayes(x)$R2_Bayes_marginal)

# Output summary
lg
r2_bayes(lg)

mean(lat30tidy$estimate)
hist(lat30tidy$estimate)
quantile(lat30tidy$estimate, c(0.025, 0.975))
mean(lat30tidy$R2_marg)
hist(lat30tidy$R2_marg)
quantile(lat30tidy$R2_marg, c(0.025, 0.975))
sum(lat30tidy$sig) / 100





## Plant size models ####

size30 = list()
for(i in 1:100){
  size30[[i]] = brm(plantGiniAdj ~
      sizeDiameterMeanLog +
      (1|species) + (1|gr(Genus_sp, cov = A)),
    prior=beta_priors,
    data=d30l[[i]],
    data2=list(A = A),
    family='beta',
    cores=12, chains=12, iter=600,
    backend = "cmdstanr",
    control = list(adapt_delta = 0.85))
  print(i)
}
saveRDS(size30, 'size30.rds')
size30 = readRDS('size30.rds')

size30tidy = size30 %>%
  purrr::map_dfr(broom.mixed::tidy) %>%
  subset(term == 'sizeDiameterMeanLog') %>%
  mutate(
    sig = ifelse(conf.high < 0, 1, 0), # is lat effect CI above zero?
  )
size30tidy$R2_marg = sapply(size30, function(x) r2_bayes(x)$R2_Bayes_marginal)

# size output summary
sg = readRDS('sg.RDS')
sg
r2_bayes(sg)

mean(size30tidy$estimate)
hist(size30tidy$estimate)
quantile(size30tidy$estimate, c(0.025, 0.975))
mean(size30tidy$R2_marg)
hist(size30tidy$R2_marg)
quantile(size30tidy$R2_marg, c(0.025, 0.975))
sum(size30tidy$sig) / 100



