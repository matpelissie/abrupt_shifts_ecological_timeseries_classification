###-###-###-###-###-###-###-###-###
#
# 04/07/2023 mathieu.pelissie@umontpellier.fr
#
# Run trajectory classification
#
# 02_run_classification.R
#
###-###-###-###-###-###-###-###-###

devtools::load_all()

library('tidyverse')
source("R/custom_functions.R")
source("R/functions_simu.R")
source("R/functions_trajclass.R")
source("R/functions_RAMLDB.R")


dir.create("analyses/classif", showWarnings = FALSE)
dir.create("analyses/classif/ts_example", showWarnings = FALSE)
dir.create("analyses/classif/library", showWarnings = FALSE)
dir.create("analyses/classif/data_transf", showWarnings = FALSE)

# Make sure you have enough core to run classification in parallel,
# otherwise reduce the number of cores just below:
ncores <- 12

# Noise levels:
noise_df <- get_data("data/00_simu/noise_levels.csv")
noise_names <- colnames(noise_df)

# Timeseries length:
name <- c("l20_1", "l25_1", "l33_1", "l50_1", "l100")

# Run Leave-One-Out procedure (step 3):
run_loo <- TRUE


# Full classification (steps 1, 2abc, 3) ----------------------------------

str <- "aic_asd"
asd_thr <- rep(0.15,5)

for (l in 1:length(name)){

  simu_list <- readRDS(paste0("data/03_simulations/all_simu_", name[l], ".rds"))

  # Run classification in parallel:
  classif <-
    parallel::mclapply(1:length(simu_list),
                       function (i) classif_noise_comb(simu_list, str,
                                                       run_loo=run_loo,
                                                       asd_thr=asd_thr[l], i),
                       mc.cores = ncores)

  # Store trajectories:
  traj_list <- lapply(classif, function(x) x$trajs)
  names(traj_list) <- names(simu_list)

  # Store matrices:
  mat_list <- lapply(classif, function(x) x$conf_mat)
  names(mat_list) <- names(simu_list)

  # Save all input and output:
  outlist_aicasd <- list("noise_df"=noise_df,
                  "simu_list"=simu_list,
                  "traj_list"=traj_list,
                  "mat_list"=mat_list)

  saveRDS(outlist_aicasd, paste0("analyses/classif/library/outlist_",
                                 name[l],"_",str,"_thr",asd_thr[l],".rds"))
}


# Classification based on AICc only (steps 1, 2bc, 3) ---------------------

str <- "aic"
asd_thr <- rep(0.15,5) # but won't be used

for (l in 1:length(name)){

  simu_list <- readRDS(paste0("data/03_simulations/all_simu_", name[l], ".rds"))

  # Run classification in parallel:
  classif <-
    parallel::mclapply(1:length(simu_list),
                       function (i) classif_noise_comb(simu_list, str,
                                                       run_loo=run_loo,
                                                       asd_thr=asd_thr[l], i),
                       mc.cores = ncores)

  # Store trajectories:
  traj_list <- lapply(classif, function(x) x$trajs)
  names(traj_list) <- names(simu_list)

  # Store matrices:
  mat_list <- lapply(classif, function(x) x$conf_mat)
  names(mat_list) <- names(simu_list)

  # Save all input and output
  outlist_aic <- list("noise_df"=noise_df,
                  "simu_list"=simu_list,
                  "traj_list"=traj_list,
                  "mat_list"=mat_list)

  saveRDS(outlist_aic, paste0("analyses/classif/library/outlist_",
                              name[l],"_",str,"_loo",run_loo,".rds"))
}



# Make timeseries transformations -----------------------------------------

# Noise levels:
noise_df <- get_data("data/00_simu/noise_levels.csv")
noise_names <- colnames(noise_df)

# Timeseries length:
name <- c("l20_1", "l25_1", "l33_1", "l50_1", "l100")

# Use step 2a (asdetect breakpoint validation):
str <- "aic_asd"

# Set detection threshold for validation:
asd_thr <- rep(0.15,5)

# Not necessary to run Leave-One-Out procedure in this case:
run_loo <- FALSE

# Transformation names:
transform <- c("raw","nrm","msz","std","log","sqr")

for (j in 1:length(transform)){

  for (l in 1:length(name)){

    simu_list <- readRDS(paste0("data/00_simu/all_simu_", name[l], ".rds"))

    if (transform[j] == "nrm"){
      simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              TB = (TB-min(TB))/(max(TB)-min(TB))
                            ) %>%
                            dplyr::ungroup()
      )
    } else if (transform[j] == "msz"){
      simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              TB = TB/mean(TB)
                            ) %>%
                            dplyr::ungroup()
      )
    } else if (transform[j] == "std"){
      simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              TB = (TB-mean(TB))/sd(TB)
                            ) %>%
                            dplyr::ungroup()
      )
    } else if (transform[j] == "log"){
      simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              TB = log(1+TB)
                            ) %>%
                            dplyr::ungroup()
      )
    } else if (transform[j] == "sqr"){
      simu_list <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::mutate(
                              TB = sqrt(TB)
                            ) %>%
                            dplyr::ungroup()
      )
    }

    # Run classification in parallel:
    classif <- parallel::mclapply(1:length(simu_list),
                                  function (i)
                                    classif_noise_comb(simu_list, str,
                                                       run_loo=run_loo,
                                                       asd_thr=asd_thr[l], i),
                                  mc.cores = ncores)

    # Store trajectories:
    traj_list <- lapply(classif, function(x) x$trajs)
    names(traj_list) <- names(simu_list)

    # Store matrices:
    mat_list <- lapply(classif, function(x) x$conf_mat)
    names(mat_list) <- names(simu_list)

    # Save all input and output
    outlist <- list(
      "noise_df"=noise_df,
      "simu_list"=simu_list, "traj_list"=traj_list,
      "mat_list"=mat_list)

    saveRDS(outlist, paste0("analyses/classif/data_transf/outlist_",
                            name[l],"_",str,"_",transform[j],
                            "_thr",asd_thr[l],".rds"))

  }
}
