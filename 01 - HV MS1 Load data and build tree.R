# 01 - HV MS1 Load data & build tree




# 1. Load packages ####

library(tidyverse)
library(V.PhyloMaker2)
library(phytools)
library(sensiPhy)
library(phyr)
library(DescTools)
library(rstan)
library(brms)
library(cmdstanr)
library(performance)
library(tidybayes)
library(bayestestR)




# 1. Set options ####

# Printing digits
options(scipen = 999)

# Stan options
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())




# 2. Load data ####

data = read.csv('data_survey_level_prep_2023-01-01 10-15-19.csv')
dsp = read.csv('data_species_level_prep_2023-01-01 10-15-19.csv')
dp = read.csv('data_plant_level_prep_2023-01-01 10-15-19.csv')

data$species = data$Genus_sp # 2nd species col for species level random effect


# 3. Build tree ####

# Extract species and rename columns to work with V.PhyloMaker2
spList0 <- subset(data, select=c("Genus_sp", "plantFamily"))
spList <- spList0[!duplicated(spList0),]
spList$genus <- gsub("(.+?)(\\_.*)", "\\1", spList$Genus_sp)
spList$species <- sub('_',' ',spList$Genus_sp)
spList$family <- spList$plantFamily
spList <- spList %>%
  select(species,genus,family)

# Build a tree using scenario S3 from V.PhyloMaker2
tree.a <- phylo.maker(sp.list=spList,
  tree=GBOTB.extended.TPL,
  nodes=nodes.info.1.TPL,
  scenarios="S3")

# Extract tree from the 'tree.a' list
hvtree <- tree.a$scenario.3

# Check number of spp matches
ifelse(length(hvtree$tip.label) == length(unique(data$Genus_sp)),
  "ONWARD", "PAUSE")

# Tally bound and pruned species, and where they are bound (for methods section)
table(tree.a$species.list$status)
bindSp <- subset(tree.a$species.list, status == "bind")

# Of bound spp, how many are bound at the genus vs family level?
# E.g., how many have a congener in the megatree?
megatreeGenera <- str_extract(GBOTB.extended.TPL$tip.label, "[^_]+")
# Number with a congener in the megatree
nrow(subset(bindSp, bindSp$genus %in% megatreeGenera))
# Number without a congener in the megatree
nrow(subset(bindSp, !bindSp$genus %in% megatreeGenera))




# 5. Build phylogenetic covariance matrix ####

A = ape::vcv.phylo(hvtree)




# 6. Set priors for brms models ####

beta_priors = c(
  set_prior("normal(0, 2)", class = "b"),
  set_prior("normal(0, 2)", class = 'Intercept'), # looks fairly flat on P scale
  set_prior('gamma(1, 0.05)', class = 'phi'),
  set_prior('cauchy(0,1)', class='sd')
)

beta_priors_noslope = beta_priors[-1,]
