###############################################################################
# 1. Do the SKAT analysis as in the doc,
#    to prove I can get it to work
###############################################################################
# Adapted from doc
library(SKAT)

# Use Tidyverse-proof variable names
data(SKAT.example)
d <- SKAT.example
rm(SKAT.example)

# Can use tibbles with named columns, so my mom can understand the data better
t_covariates <- tibble::as_tibble(d$X)
names(t_covariates) <- c("is_rocket", "speed")
d$X <- t_covariates

###############################################################################
# 2. Simulate data
###############################################################################

###############################################################################
# 2.1. Continuous phenotype is noise
###############################################################################

# Simulate the data
sim_data <- d

# No covariates
sim_data$X[, 1] <- 0.0
sim_data$X[, 2] <- 0.0

# No dichotomous phenotype
sim_data$y.b <- 0

# Continuous phenotype is all random noise
sim_data$y.c <- runif(n = length(sim_data$y.c))

# Linear null model based on continuous traits
linear_null_model_continuous <- SKAT_Null_Model(
  data = sim_data,
  formula = y.c ~ 1,
  out_type = "C" # continuous
)
high_p_value <- SKAT(Z = d$Z, linear_null_model_continuous)$p.value

###############################################################################
# 2.2. Continuous phenotype is the sum of the SNPs that are 1
#      i.e. all SNPs can contribute
###############################################################################
sim_data <- d

# No covariates
sim_data$X[, 1] <- 0.0
sim_data$X[, 2] <- 0.0

# No dichotomous phenotype
sim_data$y.b <- 0

# Continuous phenotype is the sum of your SNPs squared
# i.e. a non linear relationship
sim_data$y.c <- rowSums(sim_data$Z)

# Linear null model based on continuous traits
linear_null_model_continuous <- SKAT_Null_Model(
  data = sim_data,
  formula = y.c ~ 1,
  out_type = "C" # continuous
)
middle_p_value <- SKAT(Z = d$Z, linear_null_model_continuous)$p.value

testthat::expect_true(middle_p_value < high_p_value)

###############################################################################
# 2.3. Continuous phenotype is dependent on 1 SNP
###############################################################################
sim_data <- d

# No covariates
sim_data$X[, 1] <- 0.0
sim_data$X[, 2] <- 0.0

# No dichotomous phenotype
sim_data$y.b <- 0

# Find the first SNP that has a MAF between 5% and 20%
snp_index <- which(colSums(sim_data$Z) > 50 & colSums(sim_data$Z) < 200)[1]

# Continuous phenotype is zero for the standard,
# one for the minor allele
sim_data$y.c <- sim_data$Z[ , snp_index]

# Linear null model based on continuous traits
linear_null_model_continuous <- SKAT_Null_Model(
  data = sim_data,
  formula = y.c ~ 1,
  out_type = "C" # continuous
)
lowest_p_value <- SKAT(Z = d$Z, linear_null_model_continuous)$p.value

testthat::expect_true(lowest_p_value < middle_p_value)
