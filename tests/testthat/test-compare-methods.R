context("methods for the survival curve")
# set.seed(11)
source("./simulate_data.R")
# simulation
n_sim <- 2e2
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat
true_surv <- simulated$true_surv1

sl_lib_g <- c("SL.mean", "SL.glm")
sl_lib_censor <- c("SL.mean", "SL.glm")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.step.forward")
range(df$T.tilde)
df$T.tilde <- df$T.tilde + 1
k_grid <- 1:max(df$T.tilde)

sl_fit <- initial_sl_fit(
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  A = df$A,
  W = data.frame(df[, c("W", "W1")]),
  t_max = max(df$T.tilde),
  sl_treatment = sl_lib_g,
  sl_censoring = sl_lib_censor,
  sl_failure = sl_lib_failure
)
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

test_that("sl_1 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_1$survival, is.na)))
})
test_that("sl_0 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_0$survival, is.na)))
})

################################################################################
# moss hazard submodel
moss_hazard_l2 <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid,
  which_A_update = "obs"
)
moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
psi_moss_l2_1 <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)

moss_hazard_l2 <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid,
  which_A_update = "intervene"
)

moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
psi_moss_l2_1_alt <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)

library(data.table)
library(ggplot2)
dt <- data.table(psi_moss_l2_1,psi_moss_l2_1_alt)
dt[,t:=.I]
long <- melt(dt,id="t")
ggplot(long,aes(x=t, y=value, color=variable))+geom_line()+theme_bw()
