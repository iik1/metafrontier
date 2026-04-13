## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(metafrontier)

## ----simulate-panel-----------------------------------------------------------
set.seed(42)
panels <- lapply(1:4, function(t) {
  sim <- simulate_metafrontier(
    n_groups = 2,
    n_per_group = 50,
    beta_meta = c(1.0, 0.5, 0.3),
    tech_gap = c(0, 0.3 + 0.03 * t),  # G2 falls behind over time
    sigma_u = c(0.2, 0.3),
    sigma_v = 0.15,
    seed = 42 + t
  )
  sim$data$time <- t
  sim$data$id <- seq_len(nrow(sim$data))
  sim$data
})
panel_data <- do.call(rbind, panels)

table(panel_data$group, panel_data$time)

## ----malmquist----------------------------------------------------------------
malm <- malmquist_meta(
  log_y ~ log_x1 + log_x2,
  data = panel_data,
  group = "group",
  time = "time",
  orientation = "output",
  rts = "crs"
)

malm

## ----summary------------------------------------------------------------------
summary(malm)

## ----results-table------------------------------------------------------------
head(malm$malmquist, 10)

## ----verify-identity----------------------------------------------------------
m <- malm$malmquist
complete <- complete.cases(m[, c("MPI", "TEC", "TGC", "TC")])
all.equal(m$MPI[complete], m$TEC[complete] * m$TGC[complete] * m$TC[complete])

## ----group-vs-meta------------------------------------------------------------
# Within-group: MPI_group = EC_group x TC_group
head(malm$group_malmquist)

# Metafrontier: MPI_meta = EC_meta x TC_meta
head(malm$meta_malmquist)

## ----tgr-dynamics-------------------------------------------------------------
tgr_df <- malm$tgr

# Mean TGR by group and period
aggregate(cbind(TGR_from, TGR_to) ~ group, data = tgr_df, FUN = mean)

## ----tgc-by-group-------------------------------------------------------------
aggregate(TGC ~ group, data = tgr_df, FUN = mean)

## ----vrs-comparison-----------------------------------------------------------
malm_vrs <- malmquist_meta(
  log_y ~ log_x1 + log_x2,
  data = panel_data,
  group = "group",
  time = "time",
  rts = "vrs"
)

# Compare mean MPI under CRS vs VRS
data.frame(
  CRS = colMeans(malm$malmquist[, c("MPI", "TEC", "TGC", "TC")],
                 na.rm = TRUE),
  VRS = colMeans(malm_vrs$malmquist[, c("MPI", "TEC", "TGC", "TC")],
                 na.rm = TRUE)
)

## ----produc-example, eval = FALSE---------------------------------------------
# library(plm)
# data("Produc", package = "plm")
# 
# malm_us <- malmquist_meta(
#   gsp ~ pc + emp,
#   data = Produc,
#   group = "region",
#   time = "year",
#   rts = "crs"
# )
# summary(malm_us)

## ----utility-example, eval = FALSE--------------------------------------------
# library(sfaR)
# data("utility", package = "sfaR")
# 
# malm_util <- malmquist_meta(
#   y ~ k + labor + fuel,
#   data = utility,
#   group = "regu",
#   time = "year",
#   rts = "vrs"
# )
# summary(malm_util)

