# 05 - HV MS1 Phylogenetic signal


# First run '01 - HV MS1 Load data and build tree.R' to load and prep data



# 1. Load packages #############################################################

# None needed here that aren't loaded in 01




# 2. Prep data #################################################################
# Need named vectors

# Gini
spGini = dsp$plantGini
names(spGini) = dsp$Genus_sp

# Logit Gini
spGiniLogit = dsp$plantGiniAdjLogit
names(spGiniLogit) = dsp$Genus_sp

# Mean
spMean = dsp$plantMeanProp
names(spMean) = dsp$Genus_sp

# Logit mean
spMeanLogit = dsp$plantMeanPropLogit
names(spMeanLogit) = dsp$Genus_sp




# 3. Calculate lambda phylogenetic signal ######################################

# Gini logit
phylosig(hvtree, spGiniLogit, method='lambda', test=TRUE)

# Mean logit
phylosig(hvtree, spMeanLogit, method='lambda', test=TRUE)




# 4. Distribution of 1000 trees ################################################
# Uncertainty across binding sites

ranTrees0 <- phylo.maker(sp.list=spList,
  tree=GBOTB.extended.TPL,
  nodes=nodes.info.1.TPL,
  scenarios="S2", r=1000)
ranTrees <- ranTrees0$scenario.2

saveRDS(ranTrees, 'ranTrees.rds')

# Extract lambda for mean and Gini from each tree
nTrees = length(ranTrees)
out=list()

for(i in 1:nTrees){
  tempTree <- ranTrees[[i]]
  testGini <- phylosig(tempTree, spGiniLogit, method="lambda", test=TRUE)
  testMean <- phylosig(tempTree, spMeanLogit, method="lambda", test=TRUE)
  
  tree <- i
  scenario <- "S2"
  nTips <- length(tempTree$tip.label)
  
  lambda_gini <- testGini$lambda
  logL_gini <- testGini$logL
  logL0_gini <- testGini$logL0
  p_gini <- testGini$P
  
  lambda_mean <- testMean$lambda
  logL_mean <- testMean$logL
  logL0_mean <- testMean$logL0
  p_mean <- testMean$P
  
  out[[i]] <- data.frame(tree, scenario, nTips, lambda_gini, p_gini, logL_gini,
    logL0_gini, lambda_mean, p_mean, logL_mean, logL0_mean)
  rm(tempTree, testGini, testMean, tree, scenario, nTips, lambda_gini, p_gini,
    logL_gini, logL0_gini, lambda_mean, p_mean, logL_mean, logL0_mean)
}

treesOut <- do.call(rbind, out)

saveRDS(treesOut, 'treesOut.rds')

# Summary stats for Gini signal
mean(treesOut$lambda_gini)
ci(treesOut$lambda_gini, ci=0.95, method='ETI') # 0.51 (0.45, 0.51) P < 0.001
quantile(treesOut$lambda_gini, probs=c(0.025,0.975))
nrow(subset(treesOut, p_gini < 0.001)) # all runs (binding scenarios) have
# significant lambda @ p < 0.001

# Summary stats for mean signal
mean(treesOut$lambda_mean)
ci(treesOut$lambda_mean, ci=0.95, method='ETI') # 0.07 (0.06, 0.09) P = 1.0
nrow(subset(treesOut, p_mean < 0.05)) # no runs (binding scenarios) have
# significant signal in mean




# 5. Sensitivity to species sampling ###########################################
dsp = data.frame(dsp)
row.names(dsp) = dsp$Genus_sp

# Gini logit
sampGini <- samp_physig("plantGiniAdjLogit",
  data = dsp,
  n.sim = 1000,
  phy = hvtree,
  method="lambda")

saveRDS(sampGini, 'sampGiniPhyloSig.rds')

summary(sampGini)


# mean logit
sampMean <- samp_physig("plantMeanPropLogit",
  data = dsp,
  n.sim = 1000,
  phy = hvtree,
  method="lambda")

saveRDS(sampMean, 'sampMeanPhyloSig.rds')

summary(sampMean)

