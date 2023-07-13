###-###-###-###-###-###-###-###-###
#
# 04/07/2023 mathieu.pelissie@umontpellier.fr
#
# Generate library of time series
#
# 01_generate_timeseries.R
#
###-###-###-###-###-###-###-###-###

# devtools::load_all()

source("analyses/00_packages.R")
source("R/custom_functions.R")
source("R/functions_simu.R")
source("R/functions_trajclass.R")
source("R/functions_RAMLDB.R")

set.seed(2)

# Generate and store simulated timeseries ---------------------------------

# Noise levels
noise_df <- get_data("data/00_simu/noise_levels.csv")
noise_names <- colnames(noise_df)

param_df <- get_data("data/00_simu/param_df.csv")
name <- "l100"
### NB: For now values in iter column must be equal and all scenarios must be distinct

simu_list <- list()
scen_fct <- "make_scen"

for (n in 1:nrow(noise_df)){

  for (m in 1:length(noise_names)){

    assign(noise_names[m], noise_df %>%
             dplyr::pull(noise_names[m]) %>%
             `[`(n))
  }

  all_simu <- make_store_simu(param_df, scen_fct=scen_fct,
                              se=se, sr=sr, su=su, jfr=jfr, jsz=jsz)
  l <- nrow(all_simu %>% dplyr::filter(scen==unique(all_simu$scen)[1] & iter==1))
  noise_comb <- sprintf("l%g_se%g_sr%g_su%g_jfr%g_jsz%g",l,se,sr,su,jfr,jsz)

  simu_list[[ noise_comb ]] <- all_simu
  print(paste0(n,"/",nrow(noise_df)))
}

saveRDS(simu_list, paste0("data/00_simu/all_simu_", name, ".rds"))
print("Full length timeseries simulated and saved.")


# Make subsampling of stored time series ----------------------------------

simu_list <- readRDS(paste0("data/00_simu/all_simu_", name, ".rds"))

## Subsample 50 data points
simu_list_l50_1 <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::filter(year %% 2 == 1) %>%
                            dplyr::mutate(year=seq(1:50),
                                          scen = sub("l100", "l50", scen),
                                          scen = sub("brk50", "brk25", scen)) %>%
                            dplyr::ungroup()) %>%
  `names<-`(sub("l100", "l50", names(.)))
saveRDS(simu_list_l50_1, "data/00_simu/all_simu_l50_1.rds")


## Subsample 33 data points

simu_list_l33_1 <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::filter(year %% 3 == 0) %>%
                            dplyr::mutate(year=seq(1:33),
                                          scen = sub("l100", "l33", scen),
                                          scen = sub("brk50", "brk16", scen)) %>%
                            dplyr::ungroup()) %>%
  `names<-`(sub("l100", "l33", names(.)))

saveRDS(simu_list_l33_1, "data/00_simu/all_simu_l33_1.rds")


## Subsample 25 data points

simu_list_l25_1 <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::filter(year %% 4 == 1) %>%
                            dplyr::mutate(year=seq(1:25),
                                          scen = sub("l100", "l25", scen),
                                          scen = sub("brk50", "brk12", scen)) %>%
                            dplyr::ungroup()) %>%
  `names<-`(sub("l100", "l25", names(.)))

saveRDS(simu_list_l25_1, "data/00_simu/all_simu_l25_1.rds")


## Subsample 20 data points

simu_list_l20_1 <- lapply(simu_list,
                          function(x) x %>%
                            dplyr::group_by(scen, iter) %>%
                            dplyr::filter(year %% 5 == 1) %>%
                            dplyr::mutate(year=seq(1:20),
                                          scen = sub("l100", "l20", scen),
                                          scen = sub("brk50", "brk10", scen)) %>%
                            dplyr::ungroup()) %>%
  `names<-`(sub("l100", "l20", names(.)))

saveRDS(simu_list_l20_1, "data/00_simu/all_simu_l20_1.rds")

print("Subsampled timeseries generated and saved.")


# Generate 2 breakpoints timeseries ---------------------------------------

# Noise levels
set.seed(2)
noise_df <- get_data("data/00_simu/noise_levels.csv")
noise_names <- colnames(noise_df)

param_df <- get_data("data/00_simu/param_df_multibrk.csv")
name <- "l100_multibrk"

simu_list_mlt <- list()
scen_fct <- "make_scen_multibrk"

for (n in 1:nrow(noise_df)){

  for (m in 1:length(noise_names)){

    assign(noise_names[m], noise_df %>%
             dplyr::pull(noise_names[m]) %>%
             `[`(n))
  }

  all_simu <- make_store_simu(param_df, scen_fct=scen_fct,
                              se=se, sr=sr, su=su, jfr=jfr, jsz=jsz)
  l <- nrow(all_simu %>% dplyr::filter(scen==unique(all_simu$scen)[1] & iter==1))
  noise_comb <- sprintf("l%g_se%g_sr%g_su%g_jfr%g_jsz%g",l,se,sr,su,jfr,jsz)

  simu_list_mlt[[ noise_comb ]] <- all_simu
  print(paste0(n,"/",nrow(noise_df)))
}

saveRDS(simu_list_mlt, paste0("data/00_simu/all_simu_", name, ".rds"))
print("Multi-breaks timeseries simulated and saved.")
