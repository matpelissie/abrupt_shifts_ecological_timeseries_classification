###-###-###-###-###-###-###-###-###
#
# 17/03/2022 mathieu.pelissie@ens-lyon.fr
#
# Library of time series
#
# 03_ricker_library.R
#
###-###-###-###-###-###-###-###-###




# Run classification on simulated timeseries -------------------------------
# Maybe unnecessary because treated elsewhere?

simu_list <- readRDS(paste0("data/03_simulations/all_simu_", name, ".rds"))

mat_list <- list()
traj_list <- list()

for (i in 1:length(simu_list)){

  t <- Sys.time()
  print(paste0("noise combination: ",i,"/",length(simu_list)))
  noise_comb <- sub(".*?_", "", names(simu_list)[i])

  sets <- prep_data(df=simu_list[[i]], thr=-1, type="sim", apriori=TRUE)

  trajs <- traj_class(sets, str="aic_asd", abr_mtd=c("chg", "asd"), noise_comb=noise_comb, asd_chk=TRUE,
                      asd_thr=0.15, type="sim", showplots=FALSE, apriori=TRUE, run_loo=TRUE,
                      mad_thr=3, mad_cst=1.4826, edge_lim=5, congr_brk=5,
                      two_bkps = TRUE, smooth_signif = TRUE)

  traj_list[[noise_comb]] <- trajs

  conf_mat <- make_conf_mat(trajs$best_traj)
  mat_list[[noise_comb]] <- conf_mat

  print(paste0("completed in: ", round(Sys.time() - t, digits=2), " min/sec"))

}

outlist <- list("param_df"=param_df, "noise_df"=noise_df,
                "simu_list"=simu_list, "traj_list"=traj_list,
                "mat_list"=mat_list)

saveRDS(outlist, paste0("res/03_ricker_library/simu_library/outlist_",name,"_aic_scl_rep100.rds"))



# Individual classification -----------------------------------------------

name <- "l100"
simu_list <- readRDS(paste0("data/03_simulations/all_simu_", name, "_scl.rds"))


# CLASSIFICATION OUTPUT ---------------------------------------------------

# AIC only
outlist_aic_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_scl_loo.rds")
outlist_aic_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_scl_loo.rds")
outlist_aic_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_loo.rds")
outlist_aic <- list("l20" = outlist_aic_l20, "l50" = outlist_aic_l50, "l100" = outlist_aic_l100)

# AIC + step size >1
outlist_aicstz_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_stz_scl.rds")
outlist_aicstz_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_stz_scl.rds")
outlist_aicstz_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_stz_scl.rds")
outlist_aicstz <- list("l20" = outlist_aicstz_l20, "l50" = outlist_aicstz_l50, "l100" = outlist_aicstz_l100)

# as_detect single brk (sep=5)
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_loo_sep5.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_loo_sep5.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_loo_sep5.rds")
outlist_aicasd <- list("l20" = outlist_aicasd_l20, "l50" = outlist_aicasd_l50, "l100" = outlist_aicasd_l100)

# as_detect mlt brks (sep=5, no asd check)
outlist_asdmlt_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnck_scl_loo.rds")
outlist_asdmlt_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnck_scl_loo.rds")
outlist_asdmlt_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnck_scl_loo.rds")
outlist_aicasdmlt <- list("l20" = outlist_asdmlt_l20, "l50" = outlist_asdmlt_l50, "l100" = outlist_asdmlt_l100)

# as_detect mltbrks no check step size >1
outlist_asdmltz_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnckstz_scl.rds")
outlist_asdmltz_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnckstz_scl.rds")
outlist_asdmltz_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnckstz_scl.rds")
outlist_aicasdmltstz <- list("l20" = outlist_asdmltz_l20, "l50" = outlist_asdmltz_l50, "l100" = outlist_asdmltz_l100)

# as_detect mltbrks no check step size >1 MAD>2
outlist_asdmltmad2_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aicasdmltmad2_scl.rds")
outlist_asdmltmad2_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aicasdmltmad2_scl.rds")
outlist_asdmltmad2_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aicasdmltmad2_scl.rds")
outlist_asdmltmad2 <- list("l20" = outlist_asdmltmad2_l20, "l50" = outlist_asdmltmad2_l50, "l100" = outlist_asdmltmad2_l100)


# AIC vs AICc

# as_detect mlt brks (AIC)
outlist_AIC_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_AIC.rds")
outlist_AIC_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_AIC.rds")
outlist_AIC_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_AIC.rds")
outlist_AIC <- list("l20" = outlist_AIC_l20, "l50" = outlist_AIC_l50, "l100" = outlist_AIC_l100)

# as_detect mlt brks (AICc)
outlist_AICcasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_AICc.rds")
outlist_AICcasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_AICc.rds")
outlist_AICcasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_AICc.rds")
# outlist_AICc <- list("l20" = outlist_AICc_l20, "l50" = outlist_AICc_l50, "l100" = outlist_AICc_l100)

# outlist_AICc_l20_2 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_2_aic_asd_scl_AICc.rds")
# outlist_AICc_l50_2 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_2_aic_asd_scl_AICc.rds")
# outlist_AICc <- list("l20" = outlist_AICc_l20_2, "l50" = outlist_AICc_l50_2, "l100" = outlist_AICc_l100)

# as_detect mlt brks thr=0.15
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_thr0.15.rds")

outlist_aicasd <- list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)

# as_detect mlt brks thr=0.5
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_asd_scl_thr0.5.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_thr0.5.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_asd_scl_thr0.5.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_asd_scl_thr0.5.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_thr0.5.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_thr0.5.rds")

# as_detect mlt brks thr=0.15, improved where_as
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_asd_scl_thr0.15_21.03.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_thr0.15_21.03.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_asd_scl_thr0.15_21.03.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_asd_scl_thr0.15_21.03.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_thr0.15_21.03.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_thr0.15_21.03.rds")

outlist_aicasd_l25_0.15 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_2_aic_asd_scl_thr0.15.rds")
outlist_aicasd_l33_0.15 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_2_aic_asd_scl_thr0.15.rds")

outlist_aicasd_l25_0.14 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_2_aic_asd_scl_thr0.14.rds")
outlist_aicasd_l33_0.14 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_2_aic_asd_scl_thr0.14.rds")


outlist_aicasd_21.03 <- list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)
outlist_aicasd_21.03.2 <- list("l10"=outlist_aicasd_l10, "l25"=outlist_aicasd_l25, "l33"=outlist_aicasd_l33)

# as_detect mlt brks thr=0.2, improved where_as
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_asd_scl_thr0.2.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_thr0.2.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_asd_scl_thr0.2.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_asd_scl_thr0.2.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_thr0.2.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_thr0.2.rds")

outlist_aicasd_0.2 <- list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)

# as_detect mlt brks, improved where_as varying asd_thr
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_asd_scl_thr0.1_28.03.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_thr0.15_28.03.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_asd_scl_thr0.2_28.03.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_asd_scl_thr0.25_28.03.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_thr0.35_28.03.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_thr0.4_28.03.rds")

outlist_aicasd_28.03 <- list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)

# as_detect mlt brks, min thr variable
outlist_aicasd_l10 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l10_1_aic_asd_scl_thr0.99.rds")
outlist_aicasd_l20 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l20_1_aic_asd_scl_thr0.49.rds")
outlist_aicasd_l25 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l25_1_aic_asd_scl_thr0.24.rds")
outlist_aicasd_l33 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l33_1_aic_asd_scl_thr0.133.rds")
outlist_aicasd_l50 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l50_1_aic_asd_scl_thr0.073.rds")
outlist_aicasd_l100 <- readRDS("res/03_ricker_library/simu_library/out_minthrvar/outlist_l100_aic_asd_scl_thr0.025.rds")

outlist_aicasd_minthrvar <- list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)


# aic only
outlist_aic_l10 <- readRDS("res/03_ricker_library/simu_library/l10/outlist_l10_1_aic_scl.rds")
outlist_aic_l25 <- readRDS("res/03_ricker_library/simu_library/l25/outlist_l25_1_aic_scl.rds")
outlist_aic_l33 <- readRDS("res/03_ricker_library/simu_library/l33/outlist_l33_1_aic_scl.rds")
outlist_aic_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_scl_loo.rds")
outlist_aic_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_scl_loo.rds")
outlist_aic_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_loo.rds")



# Plot accuracy -----------------------------------------------------------

if(!exists("outlist")) outlist <- readRDS(paste0("res/03_ricker_library/simu_library/l100/outlist_",name,".rds"))
mat_list <- outlist$mat_list

accur <- sapply(1:length(mat_list),
                function(i) mat_list[[i]]$conf_mat_class$overall[["Accuracy"]],
                simplify="array") %>%
  matrix(nrow=4, byrow=TRUE) %>%
  `rownames<-`(c(".001", ".025", ".05", ".075")) %>%
  `colnames<-`(c(".001", ".025", ".05", ".075"))


accur_long <- accur %>%
  reshape2::melt(varnames=c("proc_error","env_noise"), value.name = "accur") %>%
  dplyr::mutate(env_noise = as.factor(env_noise),
                proc_error = as.factor(proc_error))


accur_plot <- ggplot(accur_long, aes(x=env_noise, y=forcats::fct_rev(proc_error), fill = accur)) +
  geom_tile()+
  geom_text(aes(label = round(accur, 3)), size=10)+
  colorspace::scale_fill_continuous_divergingx(palette = 'RdYlGn', mid = 0.75)+
  labs(title=names(accur), y="proc_error")+
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 20))+
  ggtitle("Overall accuracy")



# accurmat <- pheatmap::pheatmap(accur, cluster_rows = FALSE, cluster_cols = FALSE,
#                                legend = FALSE, display_numbers = TRUE,
#                                number_format = "%.2f", fontsize_number = 20)

# accurmat <- ComplexHeatmap::Heatmap(accur, column_title="Environmental noise",
#                                     row_title="Process error",
#                                     col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
#                                     cluster_rows = FALSE, cluster_columns = FALSE,
#                                     cell_fun = function(j, i, x, y, w, h, col){
#                                       grid::grid.text(round(accur[i,j],3), x, y)},
#                                     show_column_names = FALSE,
#                                     bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
#                                       text = ComplexHeatmap::anno_text(colnames(accur),
#                                                                        rot = 0, location = unit(1, "npc"), just = "right")),
#                                     show_heatmap_legend = FALSE)



# png(paste0("res/03_ricker_library/simu_library/noise_effect_accuracy_class_",name,"_aic_scl_rep100.png"), width = 600, height = 600)
# png(paste0("res/03_ricker_library/simu_library/noise_effect_accuracy_class_",name,"_asd_thr0.3_pval_scl.png"), width = 600, height = 600)
png(paste0("res/03_ricker_library/simu_library/noise_effect_accuracy_class_",name,"_aic_asd_thr0.15_scl.png"), width = 600, height = 600)
print(accur_plot)
dev.off()



# Plot the 4 scores for class: abrupt -------------------------------------

scores <- c("Balanced Accuracy", "Sensitivity", "Specificity", "Precision")

for (l in scores){

  score_abt <- sapply(1:length(mat_list),
                      function(i) mat_list[[i]]$conf_mat_class$byClass["Class: abrupt",][l],
                      simplify="array") %>% matrix(nrow=4, byrow=TRUE) %>%
    `rownames<-`(c(".001", ".025", ".05", ".075")) %>%
    `colnames<-`(c(".001", ".025", ".05", ".075"))

  gg_score_abt <- score_abt %>%
    reshape2::melt(x, varnames=c("proc_error","env_noise"), value.name = "score") %>%
    dplyr::mutate(env_noise = as.factor(env_noise),
                  proc_error = as.factor(proc_error))


  gl <- ggplot(gg_score_abt, aes(x=env_noise, y=forcats::fct_rev(proc_error), fill = score)) +
    geom_tile()+
    geom_text(aes(label = round(score, 2)), size=10)+
    colorspace::scale_fill_continuous_divergingx(palette = 'RdYlGn', mid = 0.75)+
    labs(title=names(score_abt), y="proc_error")+
    theme_minimal()+
    theme(legend.position="none", text = element_text(size = 20))+
    ggtitle(l)


  # Plots for all methods
  # png(paste0("res/03_ricker_library/simu_library/brkmtd_",l,"_class_",name,"_aic_scl_rep100.png"), width = 600, height = 600)
  # png(paste0("res/03_ricker_library/simu_library/brkmtd_",l,"_class_",name,"_asd_thr0.3_pval_scl.png"), width = 600, height = 600)
  png(paste0("res/03_ricker_library/simu_library/brkmtd_",l,"_class_",name,"_aic_asd_thr0.15_scl.png"), width = 600, height = 600)
  print(gl)
  dev.off()

}




# Figure SI confusion matrix noise combinations -----------------------------------------------------

conf_mat_noise(outlist_aic, "aic")
conf_mat_noise(outlist_aicstz, "aicstz")
conf_mat_noise(outlist_aicasd, "aicasd")
conf_mat_noise(outlist_aicasdmlt, "aicasdmlt")
conf_mat_noise(outlist_aicasdmltstz, "aicasdmltstz")
conf_mat_noise(outlist_asdmltmad2, "asdmltmad2")

conf_mat_noise(outlist_aicasd, "asdmltlowwl5")

conf_mat_noise(outlist_aicasd_21.03, "asdmlt_21.03(1)")
conf_mat_noise(outlist_aicasd_21.03.2, "asdmlt_21.03(2)")

conf_mat_noise(outlist_aicasd_0.2, "asdmlt_0.2")

conf_mat_noise(outlist_aicasd_28.03, "asdmlt_28.03")

conf_mat_noise(outlist_aicasd_minthrvar, "asdmlt_minthrvar")

conf_mat_noise(outlist_AIC, "AIC")
conf_mat_noise(outlist_AICc, "AICc")





# Confusion matrix + LOO distribution -------------------------------------

outlist <- readRDS(paste0("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_rep100_loo.rds"))
mat_list <- outlist$mat_list

gl <- lapply(seq_along(outlist$traj_list),
             function(x){
               ggplot(outlist$traj_list[[x]]$best_traj, aes(x=class_freq_loo))+
                 geom_histogram(binwidth=0.2)+
                 facet_grid(row = vars(class), col = vars(expected_class))+
                 theme_light()+
                 ggtitle(names(outlist$traj_list[x]))

             }
)

png(paste0("res/03_ricker_library/simu_library/noise_effect_confmat_class_",name,"_aic_asd_thr0.15_scl_loo.png"), width = 3000, height = 3000)

gridExtra::grid.arrange(grobs = gl, ncol = 4, clip = TRUE,
                        top = grid::textGrob("Environmental noise [0.001, 0.025, 0.05, 0.075]\n(each matrix: Reference)", gp=grid::gpar(fontsize=50,font=8)),
                        left = grid::textGrob("Process error [0.075, 0.05, 0.025, 0.001]\n(each matrix: Prediction)", rot=90, gp=grid::gpar(fontsize=50,font=8)))
dev.off()


# All noise combinations --------------------------------------------------

all_traj_loo <- do.call("rbind", (lapply(outlist$traj_list, function(x) x$best_traj)) )
ggplot(all_traj_loo, aes(x=class_freq_loo))+
  geom_histogram(binwidth=0.2)+
  facet_grid(row = vars(class), col = vars(expected_class))+
  theme_light()+
  ggtitle("All noise combinations")






# Figure 4_depr - Confusion matrix for varying noise, jump, length  ------------

outlist_l100 <- readRDS(paste0("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_loo.rds"))
mat_list_l100 <- outlist_l100$mat_list

outlist_l50 <- readRDS(paste0("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_scl_loo.rds"))
mat_list_l50 <- outlist_l50$mat_list

outlist_l20 <- readRDS(paste0("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_scl_loo.rds"))
mat_list_l20 <- outlist_l20$mat_list


# Extract class-specific indices into a data frame for all noise combinations
extract_index <- function(mat_list, index){

  df <- sapply(names(mat_list),
               function(i){
                 mat_list[[i]]$conf_mat_class$byClass[,index]
               }) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var="class") %>%
    dplyr::mutate(class = sub("Class: ", "", class)) %>%
    tidyr::pivot_longer(cols = !class,
                        names_to = "simu_batch",
                        values_to = index %>% sub(" ","_",.) %>%
                          tolower()) %>%
    dplyr::arrange(simu_batch)

  return(df)
}

# Combine all indices into a data frame for all noise combinations
extract_indices <- function(mat_list){

  # Class-specific indices
  indices <- c("Sensitivity","Specificity","Precision","Balanced Accuracy")
  indexlist <- lapply(indices, function(x)
    extract_index(mat_list, index=x))

  df_indices <- indexlist[[1]] %>%
    dplyr::left_join(indexlist[[2]],
                     by=c("class", "simu_batch")) %>%
    dplyr::left_join(indexlist[[3]],
                     by=c("class", "simu_batch")) %>%
    dplyr::left_join(indexlist[[4]],
                     by=c("class", "simu_batch")) %>%

    # Overall accuracy
    dplyr::left_join(
      sapply(names(mat_list),
             function(i) mat_list[[i]]$conf_mat_class$overall[["Accuracy"]]) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var="simu_batch") %>%
        `colnames<-`(c("simu_batch","accuracy")),
      by="simu_batch") %>%


    tidyr::separate(col=simu_batch, sep="_", remove=FALSE,
                    into=c("l","se","sr","su","jfr","jsz")) %>%
    dplyr::mutate(dplyr::across(!(simu_batch | class), gsub,
                                pattern = "[^0-9.-]",
                                replacement = ""),
                  dplyr::across(!(simu_batch | class), as.numeric))

  return(df_indices)
}


all_indices <- extract_indices(mat_list_l100) %>%
  rbind(extract_indices(mat_list_l50)) %>%
  rbind(extract_indices(mat_list_l20)) %>%
  dplyr::mutate(l = paste0("l",l) %>%
                  factor(levels = c("l100","l50","l20")))

sens_plot <- ggplot(all_indices, aes(x=sr, y=sensitivity))+
  geom_point(aes(shape=as.factor(jsz), color=class), alpha=0.5)+
  theme_light()+
  expand_limits(y=c(0,1))+
  scale_x_continuous(breaks = seq(0, 0.075, by = 0.025))+
  guides(shape="none", col="none")+
  # guides(shape=guide_legend(title="Jump size"))+
  facet_wrap(vars(l), ncol=1)+
  ggtitle("Sensitivity: TP / (TP + FN)")


spec_plot <- ggplot(all_indices, aes(x=sr, y=specificity))+
  geom_point(aes(shape=as.factor(jsz), color=class), alpha=0.5)+
  theme_light()+
  expand_limits(y=c(0,1))+
  scale_x_continuous(breaks = seq(0, 0.075, by = 0.025))+
  guides(shape="none", col="none")+
  # guides(shape=guide_legend(title="Jump size"))+
  facet_wrap(vars(l), ncol=1)+
  ggtitle("Specificity: TN / (TN + FP)")

prec_plot <- ggplot(all_indices, aes(x=sr, y=precision))+
  geom_point(aes(shape=as.factor(jsz), color=class), alpha=0.5)+
  theme_light()+
  expand_limits(y=c(0,1))+
  scale_x_continuous(breaks = seq(0, 0.075, by = 0.025))+
  guides(shape="none", col="none")+
  # guides(shape=guide_legend(title="Jump size"))+
  facet_wrap(vars(l), ncol=1)+
  ggtitle("Precision: TP / (TP + FP)")

balacc_plot <- ggplot(all_indices, aes(x=sr, y=balanced_accuracy))+
  geom_point(aes(shape=as.factor(jsz), color=class), alpha=0.5)+
  theme_light()+
  expand_limits(y=c(0,1))+
  scale_x_continuous(breaks = seq(0, 0.075, by = 0.025))+
  # guides(shape=FALSE, col=FALSE)+
  guides(shape=guide_legend(title="Jump size"))+
  facet_wrap(vars(l), ncol=1)+
  ggtitle("Balanced accuracy: (Sens + Spec)/2")

accu_plot <- ggplot(all_indices, aes(x=sr, y=accuracy))+
  geom_point(aes(shape=as.factor(jsz)), alpha=0.5)+
  theme_light()+
  expand_limits(y=c(0,1))+
  scale_x_continuous(breaks = seq(0, 0.075, by = 0.025))+
  guides(shape=guide_legend(title="Jump size"))+
  facet_wrap(vars(l), ncol=1)+
  ggtitle("Accuracy")


fig4 <- sens_plot + spec_plot + balacc_plot

cowplot::save_plot(filename = "res/03_ricker_library/simu_library/fig4.png",
                   plot=fig4, base_height = 8, base_width = 16)



accur <- sapply(1:length(mat_list_l100),
                function(i) mat_list_l100[[i]]$conf_mat_class$overall[["Accuracy"]],
                simplify="array") %>%
  matrix(nrow=4, byrow=TRUE) %>%
  `rownames<-`(c(".001", ".025", ".05", ".075")) %>%
  `colnames<-`(c(".001", ".025", ".05", ".075"))


accur_long <- accur %>%
  reshape2::melt(varnames=c("proc_error","env_noise"), value.name = "accur") %>%
  dplyr::mutate(env_noise = as.factor(env_noise),
                proc_error = as.factor(proc_error))


mat_list <- mat_list_l100
gl <- lapply(1:length(mat_list), function(i){

  byclass <- sum(mat_list[[i]]$conf_mat_class$table[,1])
  pheatmap::pheatmap(mat_list[[i]]$conf_mat_class$table/byclass,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     legend = FALSE, display_numbers = TRUE,
                     number_format = "%.2f", fontsize = 20,
                     number_color="black", fontsize_number = 20,
                     main = paste0("repl by class = ", byclass))

  # ComplexHeatmap::Heatmap(mat_list[[i]]$conf_mat_class$table, row_title="Prediction", column_title="Reference",
  #                         cluster_rows = FALSE, cluster_columns = FALSE,
  #                         cell_fun = function(j, i, x, y, w, h, col){
  #                           grid::grid.text(conf_mat$conf_mat_class$table[i,j], x, y)},
  #                         show_heatmap_legend = FALSE)
  grid::grid.grab()
})

gl[[14]] <- gl[[1]] ; gl[[1]] <- NULL # put no noise at the end




# Figure SI multibrks -----------------------------------------------------

name <- "l100_multibrk"
simu_list_mltbrk <- readRDS(paste0("data/03_simulations/all_simu_", name, ".rds"))

mat_list <- list()
traj_list <- list()
str_mbr <- "aic_asd"
abr_mtd_mbr <- c("chg", "asd")

for (i in 1:length(simu_list_mltbrk)){

  t <- Sys.time()
  print(paste0("noise combination: ",i,"/",length(simu_list_mltbrk)))
  noise_comb <- sub(".*?_", "", names(simu_list_mltbrk)[i])

  sets <- prep_data(df=simu_list_mltbrk[[i]], thr=-1, type="sim", apriori=TRUE)
  sets <- list("ts"=c(sets$ts[1], sets$ts[26], sets$ts[51], sets$ts[76],
                      sets$ts[101], sets$ts[126], sets$ts[151], sets$ts[176],
                      sets$ts[201], sets$ts[226], sets$ts[251], sets$ts[276]), "ts_type"=sets$ts_type)
  sets <- list("ts"=c(sets$ts[1], sets$ts[51], sets$ts[276]), "ts_type"=sets$ts_type)
  # sets <- list("ts"=c(sets$ts[201]), "ts_type"=sets$ts_type)

  ## workflow 3: AIC + asdetect:
  trajs <- traj_class(sets, str=str_mbr, abr_mtd=abr_mtd_mbr, noise_comb=noise_comb, asd_chk=TRUE,
                      asd_thr=0.15, type="sim", showplots=TRUE, apriori=TRUE, run_loo=TRUE,
                      two_bkps = TRUE, smooth_signif = TRUE, edge_lim=0, congr_brk=5, plot_one_in=1)
  # traj_list[[noise_comb]] <- trajs
  #
  # conf_mat <- make_conf_mat(trajs$best_traj)
  # mat_list[[noise_comb]] <- conf_mat

  print(paste0("completed in: ", round(Sys.time() - t, digits=2), " min/sec"))

}

outlist_mltbrk <- list("param_df"=param_df, "noise_df"=noise_df,
                "simu_list"=simu_list, "traj_list"=traj_list,
                "mat_list"=mat_list)


# Noise combination reproducible example ----------------------------------

noise_df <- get_data("data/03_simulations/noise_levels_df2.csv")
param_df <- get_data("data/03_simulations/param_df2.csv")

plot_ts_exp(param_df, noise_df, length="l100", seed=4)
plot_ts_exp(param_df, noise_df, length="l50", seed=4)
plot_ts_exp(param_df, noise_df, length="l20", seed=4)



# AIC distribution --------------------------------------------------------

res <- do.call("rbind", (lapply(outlist$traj_list, function(x) x$res)) )
# res <- outlist$traj_list[["l100_se0.025_sr0.025_su0"]]$res

res2 <- res %>%

  tibble::rownames_to_column(var="simu_id") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(min_aic = min(dplyr::c_across(dplyr::contains("aic")),
                              na.rm=TRUE),
                dplyr::across(dplyr::contains("aic") & !"min_aic", ~ .x - min_aic)) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_longer(cols = dplyr::contains("aic") & !("min_aic"|"aic_asd"),
                      names_to = "fit",
                      values_to = "delta_aic") %>%
  dplyr::mutate(expected_class = factor(expected_class, levels = c("no_change", "linear", "quadratic", "abrupt")),
                fit = sub("pol","qdr", fit), fit = sub("chg","abt", fit),
                fit = factor(fit, levels = c("aic_nch", "aic_lin", "aic_qdr", "aic_abt")),
                class = dplyr::case_when(
                  fit == "aic_nch" ~ "no_change",
                  fit == "aic_lin" ~ "linear",
                  fit == "aic_qdr" ~ "quadratic",
                  fit == "aic_abt" ~ "abrupt"),
                class = factor(class, levels = c("no_change", "linear", "quadratic", "abrupt"))
                )

misclassif <- res2 %>%
  dplyr::filter(delta_aic==0 & expected_class!=class) %>%
  dplyr::pull(simu_id) %>%
  unique()

res3 <- res2 %>%
  dplyr::filter(simu_id %in% misclassif)

res4 <- res3 %>%
  dplyr::group_by(simu_id) %>%
  dplyr::left_join(res3 %>%
                     dplyr::filter(delta_aic==0) %>%
                     dplyr::select(simu_id, class) %>%
                     dplyr::rename(misclassif_as = class),
                   by = "simu_id") %>%
  dplyr::filter(delta_aic!=0) %>%
  dplyr::slice(which.min(delta_aic)) %>%
  dplyr::ungroup()


aic_distrib <- ggplot(res4, aes(x=delta_aic, fill=fit))+
  geom_histogram(binwidth = 1)+
  theme_light()+
  coord_cartesian(x=c(0,20))+
  facet_wrap(vars(expected_class))

aic_distrib_by <- ggplot(res4, aes(x=delta_aic, fill=fit))+
  geom_histogram(binwidth = 1)+
  theme_light()+
  coord_cartesian(x=c(0,20))+
  facet_grid(row = vars(misclassif_as), col = vars(expected_class))+
  labs()


ggsave(aic_distrib_by, filename="res/03_ricker_library/simu_library/l100/aic_distrib_allnoisecomb_misclassif_zoom.png",
       width=10, height=7)




# No change performances and time series length ----------------

outlist_l20_1 <- readRDS("res/03_ricker_library/simu_library/l20_1/outlist_l20_1_aic_scl_noloo_nch.rds")
outlist_l50_1 <- readRDS("res/03_ricker_library/simu_library/l50_1/outlist_l50_1_aic_scl_noloo_nch.rds")
outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_noloo_nch.rds")

tot_mat_l20_1 <- outlist_l20_1$mat_list[[1]]$conf_mat_class$table
tot_mat_l50_1 <- outlist_l50_1$mat_list[[1]]$conf_mat_class$table
tot_mat_l100 <- outlist_l100$mat_list[[1]]$conf_mat_class$table
for(i in 2:length(outlist_l20_1$mat_list)){

  tot_mat_l20_1 <- tot_mat_l20_1 + outlist_l20_1$mat_list[[i]]$conf_mat_class$table
  tot_mat_l50_1 <- tot_mat_l50_1 + outlist_l50_1$mat_list[[i]]$conf_mat_class$table
  tot_mat_l100 <- tot_mat_l100 + outlist_l100$mat_list[[i]]$conf_mat_class$table

}

summ_df_nch <- as.data.frame(tot_mat_l20_1) %>%
  dplyr::mutate(length=20) %>%
  rbind(as.data.frame(tot_mat_l50_1) %>%
          dplyr::mutate(length=50)) %>%
  rbind(as.data.frame(tot_mat_l100) %>%
          dplyr::mutate(length=100)) %>%
  dplyr::filter(Reference=="no_change")

ggplot(summ_df_nch, aes(x="", y=Freq, fill=Prediction))+
  geom_bar(stat="identity", width=1)+
  coord_polar("y", start=0)+
  geom_text(aes(label = Freq), position = position_stack(vjust=0.5))+
  facet_wrap(~length, ncol = 3)+
  theme_void()+
  theme(legend.position = "bottom")


best_traj_l100 <- do.call("rbind", (lapply(outlist_l100$traj_list, function(x) x$best_traj))) %>%
  tibble::rowid_to_column("id") %>%
  dplyr::rename(class_l100=class)
best_traj_l50_1 <- do.call("rbind", (lapply(outlist_l50_1$traj_list, function(x) x$best_traj))) %>%
  tibble::rowid_to_column("id") %>%
  dplyr::rename(class_l50=class)
best_traj_l20_1 <- do.call("rbind", (lapply(outlist_l20_1$traj_list, function(x) x$best_traj))) %>%
  tibble::rowid_to_column("id") %>%
  dplyr::rename(class_l20=class)


best_traj_all <- best_traj_l100 %>%
  dplyr::select("id"|"class_l100") %>%
  dplyr::left_join(best_traj_l50_1 %>%
                     dplyr::select("id"|"class_l50"), by="id") %>%
  dplyr::left_join(best_traj_l20_1 %>%
                     dplyr::select("id"|"class_l20"), by="id")


library(ggalluvial)

ggplot(best_traj_all,
       aes(axis1 = class_l100, axis2 = class_l50, axis3 = class_l20)) +
  geom_alluvium(aes(fill=class_l100), width = 1/12, aes.bind = "alluvia") +
  geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("class_l100", "class_l50", "class_l20"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Effect of length on classification")




# LOO vs AIC weight -------------------------------------------------------





# [Temp] Combine magnitude output from different length:

# combine_mag <- function(outlist_all){
#
#   outlist <- lapply(outlist_all,
#                     function(z) mapply(function(x,y) x$res$res_abt$abt_res$chg %>%
#                                          tibble::rownames_to_column(var="simu_id") %>%
#                                          dplyr::mutate(length = stringr::str_extract(y, "[^_]+"),
#                                                        noise_lvl = sub(".*?\\_", "", y)),
#                                        z$traj_list, names(z$traj_list), SIMPLIFY=FALSE
#                     ) %>%
#                       # do.call(what="rbind") %>%
#                       dplyr::bind_rows() %>%
#                       tibble::remove_rownames()
#   ) %>%
#     # do.call(what="rbind") %>%
#     dplyr::bind_rows() %>%
#     dplyr::mutate(id = paste(simu_id, noise_lvl, sep="_"))
#
#   return(outlist)
# }


# mag_aic <- combine_mag(list(outlist_aicasd_l20,
#                             outlist_aicasd_l50,
#                             outlist_aicasd_l100))

outlist_aic <- combine_length(list(outlist_aic_l20,
                                   outlist_aic_l50,
                                   outlist_aic_l100))

# mag_aicasd <- combine_mag(list(outlist_aicasd_l20,
#                             outlist_aicasd_l50,
#                             outlist_aicasd_l100))

outlist_aicasd <- combine_length(list(outlist_aicasd_l20,
                                   outlist_aicasd_l50,
                                   outlist_aicasd_l100))

outlist_asdmlt <- combine_length(list(outlist_asdmlt_l20,
                                      outlist_asdmlt_l50,
                                      outlist_asdmlt_l100))


mag_aic <- outlist_aic %>%
  dplyr::left_join(mag_aic %>%
                     dplyr::select(id|contains("SD")|mag|step_size),
                   by="id")

mag_aicasd <- outlist_aicasd %>%
  dplyr::left_join(mag_aicasd %>%
                     dplyr::select(id|contains("SD")|mag|step_size),
                   by="id")



mag_aic %>%
  dplyr::filter(expected_class=="abrupt") %>%
  dplyr::slice(which(abs(step_size)<2)) %>% View()

mag_aic %>%
  dplyr::filter(expected_class=="no_change") %>%
  dplyr::slice(which(abs(step_size)>2)) %>% View()


mag_aic_plot <- mag_aic %>%
  ggplot(aes(x=abs(step_size), fill=expected_class))+
  # ggplot(aes(x=abs(step_size), fill=class))+
  geom_histogram(binwidth=0.1)+
  # facet_grid(rows = vars(length))+
  facet_grid(rows = vars(class))+
  scale_x_continuous(breaks = seq(0, 14, 2))+
  theme_light()
mag_aic_plot
ggsave(filename = paste0("res/03_ricker_library/simu_library/step_size_distrib_",
                         "expected_",
                         "facetclass_",
                         "class_simu_aiccoeff.png"),
       plot = mag_aic_plot, width=12, height=6)


mag_aicasd_plot <- mag_aicasd %>%
  ggplot(aes(x=abs(step_size), fill=expected_class))+
  # ggplot(aes(x=abs(step_size), fill=class))+
  geom_histogram(binwidth=0.1)+
  # facet_grid(rows = vars(length))+
  facet_grid(rows = vars(class))+
  scale_x_continuous(breaks = seq(0, 14, 2))+
  theme_light()
mag_aicasd_plot
ggsave(filename = paste0("res/03_ricker_library/simu_library/step_size_distrib_",
                         "expected_",
                         "facetclass_",
                         "class_simu_aicasdcoeff.png"),
       plot = mag_aicasd_plot, width=12, height=6)


mag_asdmlt_plot <- outlist_asdmlt %>%
  ggplot(aes(x=abs(step_size), fill=expected_class))+
  # ggplot(aes(x=abs(step_size), fill=class))+
  geom_histogram(binwidth=0.1)+
  # facet_grid(rows = vars(length))+
  facet_grid(rows = vars(class))+
  scale_x_continuous(breaks = seq(0, 14, 2))+
  theme_light()
mag_asdmlt_plot
ggsave(filename = paste0("res/03_ricker_library/simu_library/step_size_distrib_",
                         "expected_",
                         "facetclass_",
                         "class_simu_aicasdcoeff.png"),
       plot = mag_asdmlt_plot, width=12, height=6)



mag_aic_plot <- mag_aic %>%
  tidyr::separate(noise_lvl, sep="_",
                  into=c("se", "sr","su","jfr","jsz"),
                  remove=FALSE) %>%
  # dplyr::filter(class!=expected_class) %>%
  # dplyr::filter(class=="abrupt") %>%
  # dplyr::filter(expected_class=="abrupt") %>%
  ggplot(aes(x=abs(step_size), fill=expected_class))+
  geom_histogram(binwidth=0.1)+
  # facet_grid(rows = vars(length))+
  # facet_grid(rows = vars(class))+
  facet_grid(cols = vars(sr), rows = vars(jsz), scales="fixed")+
  scale_x_continuous(breaks = seq(0, 14, 2))+
  theme_light()
mag_aic_plot



# Detection rate with as_detect only -----

outlist_aicasd <- combine_length(list(outlist_aicasd_l20,
                                      outlist_aicasd_l50,
                                      outlist_aicasd_l100))

asd_rate <- outlist_aicasd %>%
  dplyr::mutate(asd_brk = ifelse(!is.na(loc_brk_asd),"abrupt","smooth") %>%
                  as.factor(),
                expected_class = ifelse(expected_class=="abrupt","abrupt","smooth") %>%
                  as.factor())

asd_rate_l20 <- asd_rate %>% dplyr::filter(length=="l20")
asd_rate_l50 <- asd_rate %>% dplyr::filter(length=="l50")
asd_rate_l100 <- asd_rate %>% dplyr::filter(length=="l100")
asd_mat <- caret::confusionMatrix(data=asd_rate$asd_brk, reference=asd_rate$expected_class)
asd_mat_l20 <- caret::confusionMatrix(data=asd_rate_l20$asd_brk, reference=asd_rate_l20$expected_class)
asd_mat_l50 <- caret::confusionMatrix(data=asd_rate_l50$asd_brk, reference=asd_rate_l50$expected_class)
asd_mat_l100 <- caret::confusionMatrix(data=asd_rate_l100$asd_brk, reference=asd_rate_l100$expected_class)
asd_mat$table
asd_mat_l20$table
asd_mat_l50$table
asd_mat_l100$table


# Combine LOO and AIC weight info in a long format:

loo_weight_long <- function(outlist){

  loo_weight <- outlist %>%
    dplyr::select(id | dplyr::contains("loo_")) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "class_tested",
                        values_to = "loo") %>%
    dplyr::mutate(class_tested = sub("loo_","",class_tested)) %>%

    left_join(
      outlist %>%
        dplyr::select(id | dplyr::contains("weight_aic_")) %>%
        tidyr::pivot_longer(cols = -id,
                            names_to = "class_tested",
                            values_to = "weight_aic") %>%
        dplyr::mutate(class_tested = sub("weight_aic_","",class_tested)),
      by=c("id","class_tested")
    ) %>%
    dplyr::mutate(class_tested = factor(class_tested,
                                        levels = c("no_change","linear","quadratic","abrupt"))) %>%

    dplyr::left_join(outlist %>%
                       dplyr::select(dplyr::contains(c("id","class")) | noise_lvl | length) %>%
                       dplyr::mutate(correct_classif = ifelse(class==expected_class, TRUE, FALSE)),
                     by = "id") %>%
    dplyr::mutate(length = factor(length, levels = c("l20","l50","l100")))
}


# * * Comparison of workflows --------

outlist_aic <- combine_length(list(outlist_aic_l20,
                                   outlist_aic_l50,
                                   outlist_aic_l100))
outlist_aicasd <- combine_length(list(outlist_aicasd_l20,
                                      outlist_aicasd_l50,
                                      outlist_aicasd_l100))


workflow_comp_simu <- outlist_aicasd %>%
  dplyr::select(id, class, expected_class) %>%
  dplyr::rename(class_aicasd = class) %>%
  dplyr::left_join(outlist_aic %>%
                     dplyr::select(id, class) %>%
                     dplyr::rename(class_aic = class),
                   by = "id") %>%
  dplyr::group_by(class_aicasd, class_aic) %>%
  dplyr::summarise(freq=n())

library(ggalluvial)

ggplot(workflow_comp_simu,
       aes(y = freq, axis1 = class_aic, axis2 = class_aicasd)) +
  ggalluvial::geom_alluvium(aes(fill=class_aicasd), width = 1/12, aes.bind = "alluvia") +
  ggalluvial::geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("AIC", "AIC+asdetect"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Classification output from both workflows simulated time series")


workflow_comp_simu_exp <- outlist_aicasd %>%
  dplyr::select(id, class, expected_class) %>%
  dplyr::rename(class_aicasd = class) %>%
  dplyr::left_join(outlist_aic %>%
                     dplyr::select(id, class) %>%
                     dplyr::rename(class_aic = class),
                   by = "id") %>%
  dplyr::group_by(class_aicasd, class_aic, expected_class) %>%
  dplyr::summarise(freq=n())

ggplot(workflow_comp_simu_exp,
       aes(y = freq, axis1 = class_aic, axis2 = expected_class, axis3 = class_aicasd)) +
  ggalluvial::geom_alluvium(aes(fill=expected_class), width = 1/12, aes.bind = "alluvia") +
  ggalluvial::geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("AIC", "expectation", "AIC+asdetect"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Classification output from both workflows on simulated time series and expectation")


ggplot(workflow_comp_simu_exp,
       aes(y = freq, axis1 = class_aic, axis2 = expected_class)) +
  ggalluvial::geom_alluvium(aes(fill=expected_class), width = 1/12, aes.bind = "alluvia") +
  ggalluvial::geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("AIC", "expectation"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Classification output from AIC workflow on simulated time series and expectation")


ggplot(workflow_comp_simu_exp,
       aes(y = freq, axis1 = class_aicasd, axis2 = expected_class)) +
  ggalluvial::geom_alluvium(aes(fill=expected_class), width = 1/12, aes.bind = "alluvia") +
  ggalluvial::geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("AIC+asdetect", "expectation"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Classification output from AIC + asdetect workflow on simulated time series and expectation")



# Graph sensitivity-specificity --------------

sens_spec_bylength <- function(outlist, workflow, class){

  pool_mat <- outlist$mat_list[[1]]$conf_mat_class$table
  for(i in 2:length(outlist$mat_list)){
    pool_mat <- pool_mat + outlist$mat_list[[i]]$conf_mat_class$table
  }

  TP <- pool_mat[class,class]
  FP <- sum(pool_mat[class,]) - TP
  FN <- sum(pool_mat[,class]) - TP
  TN <- sum(pool_mat) - (TP+FP+FN)

  sens <- TP/(FN+TP)
  spec <- TN/(FP+TN)
  # spec <- lapply(outlist$mat_list, # Equivalent since the sum of all columns is equal among matrices
  #        function(x) x$conf_mat_class$byClass[paste0("Class: ",class),"Specificity"]) %>%
  #   unlist() %>% mean()

  df <- data.frame(
    workflow = workflow,
    length = stringr::str_extract(names(outlist$mat_list)[1], "[^_]+"),
    class = class,
    sens = sens,
    spec = spec)

  return(df)
  }


sens_spec_df <-
  sapply(list("no_change", "linear", "quadratic", "abrupt"),
         function(i)
           mapply(
             function(j, k)
               sens_spec_bylength(
                 outlist = j,
                 workflow = k,
                 class = i
               ),
             list(
               # outlist_aic_l20,
               # outlist_aic_l50,
               # outlist_aic_l100,
               # outlist_aicasd_l20,
               # outlist_aicasd_l50,
               # outlist_aicasd_l100,
               # outlist_asdmlt_l20,
               # outlist_asdmlt_l50,
               # outlist_asdmlt_l100
               # # outlist_asdmltz_l20,
               # outlist_asdmltz_l50,
               # outlist_asdmltz_l100

               outlist_AICcasd_l10,
               outlist_AICcasd_l20,
               outlist_AICcasd_l25,
               outlist_AICcasd_l33,
               outlist_AICcasd_l50,
               outlist_AICcasd_l100,



               outlist_AICc_l10,
               outlist_AICc_l20,
               outlist_AICc_l25,
               outlist_AICc_l33,
               outlist_AICc_l50,
               outlist_AICc_l100


             ),
             c(rep("aic_asd_mltpl", 6), rep("aic", 6)),
             # c(rep("aic", 3), rep("aic_asd", 3)),
             # c(rep("aic", 3), rep("aic_asd", 3), rep("aic_asd_mltpl", 3), rep("aic_asd_mltpl_stpsz", 3)),
             SIMPLIFY = FALSE
           ) %>% do.call(rbind, .),
         simplify = FALSE) %>% do.call(rbind, .)


sens_spec_plot <- sens_spec_df %>%
  # dplyr::mutate(length = factor(length, levels=c("l100","l50","l20"))) %>%
  dplyr::mutate(length = factor(length, levels=c("l100","l50","l33","l25","l20","l10"))) %>%
  # dplyr::filter(class=="abrupt") %>%
  ggplot(aes(x=1-sens, y=1-spec))+
  geom_point(aes(shape=workflow, color=length), size=3, alpha=0.7)+
  # geom_point(aes(color=workflow, shape=length), size=3)+
  expand_limits(x=c(0,1), y=c(0,1))+
  facet_wrap(vars(class))+
  labs(x="FNR (1-sensitivity)", y="FPR (1-specificity)")+
  geom_abline(slope=-1, intercept=1, lty=2)+
  theme_light()
sens_spec_plot
ggsave(filename="res/03_ricker_library/simu_library/sens_spec_plot_workflows_mlt_AICc_only.png",
       plot=sens_spec_plot,
       width=8, height=7)


# By noise combination
sens_spec_bynoise <- function(outlist, workflow, class){

  df <- data.frame()
  for(i in 1:length(outlist$mat_list)){

    mat <- outlist$mat_list[[i]]$conf_mat_class$table

    TP <- mat[class,class]
    FP <- sum(mat[class,]) - TP
    FN <- sum(mat[,class]) - TP
    TN <- sum(mat) - (TP+FP+FN)

    sens <- TP/(FN+TP)
    spec <- TN/(FP+TN)

    dline <- data.frame(
      workflow = workflow,
      noise = names(outlist$mat_list)[i],
      class = class,
      sens = sens,
      spec = spec) %>%
      tidyr::separate(col=noise, sep="_", remove=FALSE,
                      into=c("length","se","sr","su","jfr","jsz")) %>%
      dplyr::mutate(dplyr::across(!(workflow | noise | class | length), gsub,
                                  pattern = "[^0-9.-]",
                                  replacement = ""),
                    dplyr::across(!(workflow | noise | class | length), as.numeric))

    df <- dplyr::bind_rows(df, dline)

  }

  return(df)
}


sens_spec_noise_df <-
  sapply(list("no_change", "linear", "quadratic", "abrupt"),
         function(i)
           mapply(
             function(j, k)
               sens_spec_bynoise(
                 outlist = j,
                 workflow = k,
                 class = i
               ),
             list(
               # outlist_aic_l20,
               # outlist_aic_l50,
               # outlist_aic_l100
               outlist_asdmlt_l20,
               outlist_asdmlt_l50,
               outlist_asdmlt_l100
             ),
             c(rep("aic_asd_mltpl", 3)),
             # c(rep("aic", 3), rep("aic_asd", 3)),
             # c(rep("aic", 3), rep("aic_asd", 3), rep("aic_asd_mltpl", 3), rep("aic_asd_mltpl_stpsz", 3)),
             SIMPLIFY = FALSE
           ) %>% do.call(rbind, .),
         simplify = FALSE) %>% do.call(rbind, .)


sens_spec_noise_plot <- sens_spec_noise_df %>%
  dplyr::mutate(length = factor(length, levels=c("l100","l50","l20")),
                # class = factor(class, levels=c("abrupt","no_change","linear","quadratic")),
                sr = factor(sr),
                jsz = factor(jsz)) %>%
  # dplyr::filter(class=="abrupt") %>%
  ggplot(aes(x=1-sens, y=1-spec))+
  geom_point(aes(shape=jsz, color=sr), size=3, alpha=0.7)+
  # geom_point(aes(color=workflow, shape=length), size=3)+
  expand_limits(x=c(0,1), y=c(0,1))+
  facet_grid(cols=vars(class), rows=vars(length))+
  labs(x="FNR (1-sensitivity)", y="FPR (1-specificity)")+
  geom_abline(slope=-1, intercept=1, lty=2)+
  theme_light()
sens_spec_noise_plot
ggsave(filename="res/03_ricker_library/simu_library/sens_spec_noise_plot_asd.png",
       plot=sens_spec_noise_plot,
       width=10, height=7)


sens_spec_noise_choice_plot <- sens_spec_noise_df %>%
  dplyr::filter((sr==0.001 & se==0.025) | (sr==0.075 & se==0.025)) %>%
  dplyr::mutate(length = factor(length, levels=c("l100","l50","l20")),
                # class = factor(class, levels=c("abrupt","no_change","linear","quadratic")),
                noise_lvl = ifelse(sr==0.001 & se==0.025, "low noise", "high noise") %>% as.factor()) %>%
  # dplyr::filter(class=="abrupt") %>%
  ggplot(aes(x=1-sens, y=1-spec))+
  geom_point(aes(shape=noise_lvl, color=length), size=3, alpha=0.7)+
  # geom_point(aes(color=workflow, shape=length), size=3)+
  expand_limits(x=c(0,1), y=c(0,1))+
  facet_wrap(vars(class))+
  labs(x="FNR (1-sensitivity)", y="FPR (1-specificity)")+
  geom_abline(slope=-1, intercept=1, lty=2)+
  theme_light()+
  # ggtitle("AIC only")+
  ggtitle("AIC + as_detect multiple")


sens_spec_noise_choice_plot
ggsave(filename="res/03_ricker_library/simu_library/sens_spec_noise_choice_plot_aicasd.png",
       plot=sens_spec_noise_choice_plot,
       width=8, height=7)


# With LOO
library(ggridges)

loo_weight_aic <- loo_weight_long(outlist_aicasd)

loo_weight_aic %>%
  dplyr::filter(class_tested == class) %>%
  ggplot()+
  geom_histogram(aes(x=loo))


loo <- loo_weight_aic %>%
  dplyr::filter(class_tested == class) %>%
  dplyr::mutate(loo_lvl = dplyr::case_when(
    loo < 0.5 ~ "loo_0-.5",
    loo < 0.75 ~ "loo_.5-.75",
    loo < 0.9 ~ "loo_.75-.9",
    loo <= 1 ~ "loo_.9-1"),
    loo_lvl = factor(loo_lvl, levels = c("loo_0-.5","loo_.5-.75",
                                         "loo_.75-.9","loo_.9-1"))
  )

loo %>%
  ggplot(aes(x=loo, y=class, fill=class), alpha=0.5)+
  geom_density_ridges(stat = "binline",
                      binwidth = 0.05, draw_baseline = FALSE)+
  coord_flip()+
  facet_wrap(vars(expected_class), ncol=4)+
  expand_limits(y=c(0,1))+
  theme_light()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))


loo %>%
  ggplot()+
  geom_histogram(aes(x=loo, fill=expected_class),
                 position="stack",binwidth = 0.05, alpha=0.5)+
  coord_flip()+
  facet_wrap(vars(class), ncol=4)+
  expand_limits(y=c(0,1))+
  theme_light()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  labs(y="Output classes (frequency)")



workflow_comp_simu_loo <- loo %>%
  dplyr::select(id, class, expected_class, loo_lvl) %>%
  dplyr::rename(class_aic = class) %>%
  # dplyr::left_join(outlist_aic %>%
  #                    dplyr::select(id, class) %>%
  #                    dplyr::rename(class_aic = class),
  #                  by = "id") %>%
  dplyr::group_by(class_aic, expected_class, loo_lvl) %>%
  dplyr::summarise(freq=n())

ggplot(workflow_comp_simu_loo,
       aes(y = freq, axis1 = expected_class, axis2 = class_aic, axis3 = loo_lvl)) +
  ggalluvial::geom_alluvium(aes(fill=expected_class), width = 1/12, aes.bind = "alluvia") +
  ggalluvial::geom_stratum(fill="grey30", width = 1/12, color = "white") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("expected", "output", "LOO"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme_minimal()+
  theme(legend.position = "none")+
  ggtitle("Classification output from both workflows simulated time series")





loo_weight_aicasd <- loo_weight_long(outlist_aicasd)


# [different] diff_loc in abrupt trajectories

outlist %>%
  dplyr::filter(expected_class=="abrupt") %>%
  pull(diff_loc) %>%
  # mean(na.rm = TRUE)
  hist(breaks=seq(0,100,by=1))


# Check maximal AIC weight corresponds to best class (for aic without asd)

test <- outlist %>%
  # dplyr::filter(length=="l100" & simu_id %in% paste0("l100_cst_F3_r1.0_H0.75_iter",1:4,"0")) %>%
  dplyr::select(id | contains("weight") | class | best_model) %>%
  tidyr::pivot_longer(cols = contains("weight"),
                      names_to = "class_weight",
                      values_to = "weight") %>%
  mutate(class_weight = sub("weight_aic_","",class_weight)) %>%
  dplyr::group_by(id) %>%
  dplyr::slice(which.max(weight)) %>%
  dplyr::group_by(class, class_weight) %>%
  dplyr::summarise(n())

# The "exceptions" are due to non-unique AIC weights due to rounding to 3 digits
expt <- outlist %>%
  # dplyr::filter(length=="l100" & simu_id %in% paste0("l100_cst_F3_r1.0_H0.75_iter",1:4,"0")) %>%
  dplyr::select(id | contains("weight") | class | best_model) %>%
  tidyr::pivot_longer(cols = contains("weight"),
                      names_to = "class_weight",
                      values_to = "weight") %>%
  mutate(class_weight = sub("weight_aic_","",class_weight)) %>%
  dplyr::group_by(id) %>%
  dplyr::slice(which.max(weight)) %>%
  dplyr::filter((class != class_weight) & class != "no_change") %>%
  dplyr::pull(id)

outlist %>%
  filter(id %in% expt)




# Compare confusion matrix of following two (and/or filter for LOO/AIC weight plot)

loo_weight %>%
  dplyr::group_by(id) %>%
  dplyr::filter(weight_aic==max(weight_aic)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(correct_classif = ifelse(class_tested==expected_class, TRUE, FALSE)) %>%
  dplyr::count(correct_classif)

#     correct_classif     n
#     <lgl>           <int>
#   1 FALSE           11174
#   2 TRUE            20035


loo_weight %>%
  dplyr::filter(class_tested==class) %>%
  dplyr::count(correct_classif)

#     correct_classif     n
#     <lgl>           <int>
#   1 FALSE            7804
#   2 TRUE            23396


# For confusion matrix
loo_weight %>%
  dplyr::group_by(id) %>%
  dplyr::filter(weight_aic==max(weight_aic)) %>%
  dplyr::ungroup() %>%
  dplyr::count(class, expected_class, .drop = FALSE)


loo_weight %>%
  dplyr::filter(class_tested==class) %>%
  dplyr::count(class, expected_class, .drop = FALSE)




loo_aicw_plot <- ggplot(loo_weight
                        # %>% dplyr::filter(class_tested=="abrupt")
                        # %>% dplyr::filter(class_tested==class) # specificity
                        # %>% dplyr::filter(class_tested==expected_class) # sensitivity
                        # %>% dplyr::filter(class==expected_class)
                        %>% dplyr::group_by(id) %>%
                          dplyr::filter(weight_aic==max(weight_aic)) %>%
                          dplyr::ungroup()
                        ,
                        aes(x=weight_aic, y=loo, shape=expected_class))+

  geom_point(alpha=0.3, aes(color=correct_classif))

loo_aicw_plot



loo_aicw_plot <- ggplot(loo_weight
       # %>% dplyr::filter(class_tested=="abrupt")
       %>% dplyr::filter(class_tested==class) # specificity
       # %>% dplyr::filter(class_tested==expected_class) # sensitivity
       # %>% dplyr::filter(class==expected_class)
       ,
       aes(x=weight_aic, y=loo, shape=class))+

  geom_point(alpha=0.3, aes(color=correct_classif))+
  # geom_point(alpha=0.3)+

  # Sensitivity:
  # labs(title = "Sensitivity (scores of time series fits of expected classes)", x = "AIC weight", y = "Leave-One-Out frequency")+
  # guides(shape=guide_legend(title="classification output"),
         # color=guide_legend(title="Fit the expected class"))+

  # Specificity:
  labs(title = "Specificity (scores of time series fits actually obtained)", x = "AIC weight", y = "Leave-One-Out frequency")+
  guides(shape=guide_legend(title="classification output"),
         color=guide_legend(title="Fit the expected class"))+

  theme_light()+
  facet_grid(rows=vars(length), cols=vars(class_tested))
  # facet_wrap(facets=vars(class_tested, length), ncol=4, as.table=FALSE)

loo_aicw_plot

ggsave(
  "res/03_ricker_library/simu_library/loo_aicw_sensitivity.png"
  # "res/03_ricker_library/simu_library/loo_aicw_specificity.png"
  , plot = loo_aicw_plot,
       width = 4000, height = 2000, units = "px")



### Specificity

loo_aicw_plot <- ggplot(loo_weight
                        # %>% dplyr::filter(noise_lvl=="se0.025_sr0.075_su0_jfr0.1_jsz5")
                        %>% dplyr::filter(class_tested==class), # specificity
                        aes(x=weight_aic, y=loo, shape=expected_class))+

  geom_point(alpha=0.3, aes(color=correct_classif))+
  # geom_point(alpha=0.3)+

  # Specificity:
  labs(title = "Specificity (scores of time series fits actually obtained)", x = "AIC weight", y = "Leave-One-Out frequency")+
  guides(shape=guide_legend(title="expected classification"),
         color=guide_legend(title="Fit the expected class"))+


  theme_light()+
  facet_grid(rows=vars(length), cols=vars(class))

loo_aicw_plot

ggsave(
  "res/03_ricker_library/simu_library/loo_aicw_output.png"
  , plot = loo_aicw_plot,
  width = 4000, height = 2000, units = "px")



# LOO by type -------------------------------------------------------------

outlist_aic <- combine_length(list(outlist_aic_l20,
                                   outlist_aic_l50,
                                   outlist_aic_l100))

loo_boxplot_fun(outlist, save=FALSE)

loo_l20 <- loo_boxplot_fun(outlist %>%
                  dplyr::filter(length=="l20"), save=FALSE)
loo_l50 <- loo_boxplot_fun(outlist %>%
                  dplyr::filter(length=="l50"), save=FALSE)
loo_l100 <- loo_boxplot_fun(outlist %>%
                  dplyr::filter(length=="l100"), save=FALSE)

# png("res/03_ricker_library/simu_library/loo_type_aic_asd_sep3_thr0.15_scl_length.png", width = 1500, height = 500)
# png("res/03_ricker_library/simu_library/loo_type_aic_asd_sep5_thr0.15_scl_length.png", width = 1500, height = 500)
png("res/03_ricker_library/simu_library/loo_type_aic_scl_length.png", width = 1500, height = 500)
cowplot::plot_grid(loo_l100, loo_l50, loo_l20,
                   ncol= 3, labels=c("full time series","subsample 1/2", "subsample 1/5"))
dev.off()


# LOO by class -------------------------------------------------------------

outlist_aic <- combine_length(list(outlist_aic_l20,
                                   outlist_aic_l50,
                                   outlist_aic_l100))

outlist_aicasdmlt <- combine_length(list(outlist_asdmlt_l20,
                                   outlist_asdmlt_l50,
                                   outlist_asdmlt_l100))

loo_aicasdmlt <- loo_boxplot_class_fun(outlist_aicasdmlt, save=TRUE,
                                     workflow="aicasdmlt", class_type="expected_class",
                                     title="AIC + as_detect multiple")
loo_aic <- loo_boxplot_class_fun(outlist_aic, save=TRUE, class_type="expected_class",
                               workflow="aic", title = "AIC only")


# AIC weight by class -----------

wAIC_aicasdmlt <- weight_boxplot_class_fun(outlist_aicasdmlt, save=TRUE,
                                     workflow="aicasdmlt", title="AIC + as_detect multiple")
wAIC_aic <- weight_boxplot_class_fun(outlist_aic, save=TRUE,
                               workflow="aic", title = "AIC only")


# Conf mat in one --------------------------------------

# aic_asd sep = 3
# outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_loo_sep3.rds")
# outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_loo_sep3.rds")
# outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_loo_sep3.rds")

# aic_asd sep = 5
outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_loo_sep5.rds")
outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_loo_sep5.rds")
outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_loo_sep5.rds")

# aic
# outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_scl_loo.rds")
# outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_scl_loo.rds")
# outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_loo.rds")

# aic_asd_mlt
# outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_asdmlt.rds")
# outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_asdmlt.rds")
# outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_asdmlt.rds")

# aic_asd_mlt_nck
outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_asdmltnck_scl.rds")
outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_asdmltnck_scl.rds")
outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_asdmltnck_scl.rds")

# aic_asd_mlt_nck_stz<1
outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnckstz_scl.rds")
outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnckstz_scl.rds")
outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnckstz_scl.rds")



mat_list <- c(outlist_l20$mat_list,
              outlist_l50$mat_list,
              outlist_l100$mat_list)

tot_mat <- mat_list[[1]]$conf_mat_class$table


for(i in 2:length(mat_list)){
  tot_mat <- tot_mat + mat_list[[i]]$conf_mat_class$table
}
byclass <- sum(tot_mat[,1])

# png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_sep5_aic_asd_thr0.15_scl.png"), width = 3000, height = 3000)
# png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_sep3_aic_asd_thr0.15_scl.png"), width = 3000, height = 3000)
# png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_aic_scl.png"), width = 3000, height = 3000)
# png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_sep5_aic_asdmlt_thr0.15_scl.png"), width = 3000, height = 3000)
# png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_sep5_aic_asdmltnck_thr0.15_scl.png"), width = 3000, height = 3000)
png(paste0("res/03_ricker_library/simu_library/confmat_class_pooled_sep5_aic_asdmltnckstz_thr0.15_scl.png"), width = 3000, height = 3000)


pheatmap::pheatmap(tot_mat/byclass,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   legend = FALSE, display_numbers = TRUE,
                   number_format = "%.2f", fontsize = 170,
                   number_color="black", fontsize_number = 170)
dev.off()



# All classes -------------------------------------------------------------




# aic_asd sep = 3
p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_loo_sep3.rds")
p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_loo_sep3.rds")
p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_loo_sep3.rds")

# aic_asd sep = 5
# p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_loo_sep5.rds")
# p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_loo_sep5.rds")
# p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_loo_sep5.rds")

# aic
# p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_scl_loo.rds")
# p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_scl_loo.rds")
# p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_scl_loo.rds")

# aic_asd_mlt sep = 5
# p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asd_scl_asdmlt.rds")
# p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asd_scl_asdmlt.rds")
# p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asd_scl_asdmlt.rds")

# aic_asd_mlt_nocheck sep = 5 loo
p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnck_scl_loo.rds")
p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnck_scl_loo.rds")
p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnck_scl_loo.rds")

# aic_asd_mlt_nocheck sep = 5 step size <1
p20 <- conf_mat_pool("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnckstz_scl.rds")
p50 <- conf_mat_pool("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnckstz_scl.rds")
p100 <- conf_mat_pool("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnckstz_scl.rds")



# png("res/03_ricker_library/simu_library/confmat_class_pool_aic_asd_sep3_thr0.15_scl_length.png", width = 9000, height = 3000)
# png("res/03_ricker_library/simu_library/confmat_class_pool_aic_asd_sep5_thr0.15_scl_length.png", width = 9000, height = 3000)
# png("res/03_ricker_library/simu_library/confmat_class_pool_aic_scl_length.png", width = 9000, height = 3000)
# png("res/03_ricker_library/simu_library/confmat_class_pool_aic_asdmlt_sep5_thr0.15_scl_length.png", width = 9000, height = 3000)
png("res/03_ricker_library/simu_library/confmat_class_pool_aic_asdmltnckstz_sep5_thr0.15_loo_scl_length.png", width = 9000, height = 3000)


cowplot::plot_grid(p20$gtable, p50$gtable, p100$gtable,
                   ncol= 3, labels=LETTERS[1:3])
dev.off()



# Why not abrupt ----------------------------------------------------------

outlist <- list("l10" = outlist_aicasd_l10, "l20" = outlist_aicasd_l20,
                "l25" = outlist_aicasd_l25, "l33" = outlist_aicasd_l33,
                "l50" = outlist_aicasd_l50, "l100" = outlist_aicasd_l100)

outlist <- list("l25_thr0.15" = outlist_aicasd_l25_0.15, "l33_thr0.15" = outlist_aicasd_l33_0.15,
                "l25_thr0.14" = outlist_aicasd_l25_0.14, "l33_thr0.14" = outlist_aicasd_l33_0.14)

outlist <- list("l25_thr0.14" = outlist_aicasd_l25_0.14, "l33_thr0.14" = outlist_aicasd_l33_0.14)

out <- combine_length(outlist)

out %>%
  dplyr::filter(expected_class=="abrupt") %>%
  dplyr::group_by(length, asd_thr, not_abr) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::mutate(length=factor(length, levels=c("l10","l20","l25","l33","l50","l100")),
                not_abr=factor(not_abr, levels=c("no_congr_brk","no_asd_brk","not_lowest_aic"))) %>%
  tidyr::drop_na() %>%
  ggplot(aes(x=length, y=n, fill=not_abr))+
  geom_col()+
  # facet_wrap(vars(asd_thr))+
  theme_light()



# Classif according asdetect only -----------------------------------------

outlist <- list("l10" = outlist_aicasd_l10, "l20" = outlist_aicasd_l20,
                "l25" = outlist_aicasd_l25, "l33" = outlist_aicasd_l33,
                "l50" = outlist_aicasd_l50, "l100" = outlist_aicasd_l100)

out <- combine_length(outlist)

out_asd <- out %>%
  dplyr::mutate(asd_class = ifelse(!is.na(loc_brk_asd), "abrupt","not_abrupt"),
                aicasd_class = ifelse(class!="abrupt", "not_abrupt", "abrupt"),
                abr_expected_class = ifelse(expected_class!="abrupt", "not_abrupt", "abrupt"),
                asd_class = factor(asd_class, levels=c("abrupt","not_abrupt")),
                aicasd_class = factor(aicasd_class, levels=c("abrupt","not_abrupt")),
                abr_expected_class = factor(abr_expected_class, levels=c("abrupt","not_abrupt")),
                length = factor(length, levels=c("l10","l20","l25","l33","l50","l100")))

outlist_aic <- list("l10" = outlist_aic_l10, "l20" = outlist_aic_l20,
                "l25" = outlist_aic_l25, "l33" = outlist_aic_l33,
                "l50" = outlist_aic_l50, "l100" = outlist_aic_l100)

out_aic <- combine_length(outlist_aic)

out_aic <- out_aic %>%
  dplyr::mutate(aic_class = ifelse(class!="abrupt", "not_abrupt", "abrupt"),
                abr_expected_class = ifelse(expected_class!="abrupt", "not_abrupt", "abrupt"),
                aic_class = factor(aic_class, levels=c("abrupt","not_abrupt")),
                abr_expected_class = factor(abr_expected_class, levels=c("abrupt","not_abrupt")),
                length = factor(length, levels=c("l10","l20","l25","l33","l50","l100")))

conf_mat_aic <- caret::confusionMatrix(data=out_aic$aic_class, reference=out_aic$abr_expected_class)
conf_mat_asd <- caret::confusionMatrix(data=out_asd$asd_class, reference=out_asd$abr_expected_class)
conf_mat_aicasd <- caret::confusionMatrix(data=out_asd$aicasd_class, reference=out_asd$abr_expected_class)


conf_mat_aic$table[,1] <- conf_mat_aic$table[,1]/colSums(conf_mat_aic$table)[1]
conf_mat_aic$table[,2] <- conf_mat_aic$table[,2]/colSums(conf_mat_aic$table)[2]

p_aic <- pheatmap::pheatmap(conf_mat_aic$table,
                            breaks=seq(0,1, by=0.01),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            legend = FALSE, display_numbers = TRUE,
                            number_format = "%.2f", fontsize = 170,
                            number_color="black", fontsize_number = 170,
                            main="AICc only")


conf_mat_asd$table[,1] <- conf_mat_asd$table[,1]/colSums(conf_mat_asd$table)[1]
conf_mat_asd$table[,2] <- conf_mat_asd$table[,2]/colSums(conf_mat_asd$table)[2]

p_asd <- pheatmap::pheatmap(conf_mat_asd$table,
                            breaks=seq(0,1, by=0.01),
                        cluster_rows = FALSE, cluster_cols = FALSE,
                        legend = FALSE, display_numbers = TRUE,
                        number_format = "%.2f", fontsize = 170,
                        number_color="black", fontsize_number = 170,
                        main="asdetect only")


conf_mat_aicasd$table[,1] <- conf_mat_aicasd$table[,1]/colSums(conf_mat_aicasd$table)[1]
conf_mat_aicasd$table[,2] <- conf_mat_aicasd$table[,2]/colSums(conf_mat_aicasd$table)[2]

p_aicasd <- pheatmap::pheatmap(conf_mat_aicasd$table,
                               breaks=seq(0,1, by=0.01),
                            cluster_rows = FALSE, cluster_cols = FALSE,
                            legend = FALSE, display_numbers = TRUE,
                            number_format = "%.2f", fontsize = 170,
                            number_color="black", fontsize_number = 170,
                            main="AICc + asdetect")

png("res/03_ricker_library/simu_library/abr_detection_propbis.png", width = 9000, height = 3000)

cowplot::plot_grid(p_aic$gtable, p_aicasd$gtable, p_asd$gtable,
                   ncol= 3, label_size=170)
dev.off()



## Broken down by length ------------------------------------------------

out_aic_length <- out_aic %>%
  split(out_aic$length)

out_asd_length <- out_asd %>%
  split(out_asd$length)


conf_mat_length <- function(list_length, colname){


  mats <- lapply(1:length(list_length),
                 function(i){
                   m <- caret::confusionMatrix(data=list_length[[i]][[colname]],
                                          reference=list_length[[i]]$abr_expected_class) %>%
                   `[`("table") %>%
                   `[[`(1)
                   m[,1] <- m[,1]/colSums(m)[1]
                   m[,2] <- m[,2]/colSums(m)[2]
                   m
                   }) %>%
    `names<-`(names(list_length))

  plots <- lapply(1:length(mats),
                  function(i){
                    pheatmap::pheatmap(mats[[i]],
                                       cluster_rows = FALSE, cluster_cols = FALSE,
                                       legend = FALSE, display_numbers = TRUE,
                                       number_format = "%.2f", fontsize = 10,
                                       main = names(mats)[i],
                                       number_color="black", fontsize_number = 15)
                  grid::grid.grab()}
                  )

  pl <- cowplot::plot_grid(plotlist = plots, nrow=1)

  return(pl)
}

p_aic <- conf_mat_length(out_aic_length, "aic_class")
p_asd <- conf_mat_length(out_asd_length, "asd_class")
p_aicasd <- conf_mat_length(out_asd_length, "aicasd_class")

by_length <- cowplot::plot_grid(plotlist = list(p_aic, p_asd, p_aicasd), nrow=3)

ggsave(by_length, filename="res/03_ricker_library/simu_library/abr_detection_by_length_thr0.2_impr_prop.pdf",
       width=14, height=7)

# Max asdetect line -------------------------------------------------------


# All time series
distr_asd <- max_asdetect(list("l20"=outlist_aicasd_l20, "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100))

distr_asd <- max_asdetect(list("l10"=outlist_aicasd_l10, "l20"=outlist_aicasd_l20,
                               "l25"=outlist_aicasd_l25,"l33"=outlist_aicasd_l33,
                               "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100))


outlist_aic_all <- list("l10"=outlist_aic_l10, "l20"=outlist_aic_l20,
                        "l25"=outlist_aic_l25,"l33"=outlist_aic_l33,
                        "l50"=outlist_aic_l50, "l100"=outlist_aic_l100)

outlist_aicasd_all <- list("l10"=outlist_aicasd_l10, "l20"=outlist_aicasd_l20,
                           "l25"=outlist_aicasd_l25,"l33"=outlist_aicasd_l33,
                           "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100)

outlist_aicasd_abrupt <- keep_aic(outlist_aic_all, outlist_aicasd_all, class_kept="abrupt")

distr_asd <- max_asdetect(outlist_aicasd_abrupt)



p_distr_asd <- distr_asd %>%
  dplyr::mutate(length = factor(length, levels = c("l10","l20","l25","l33","l50","l100")),
                expected_class = dplyr::case_when(
                  grepl("cst",simu_id) ~ "no_change",
                  grepl("lin_pos_F1-5_r1.6_H0.75",simu_id) |
                    grepl("lin_neg_F5-1_r1.6_H0.75",simu_id) ~ "linear",
                  grepl("lin_pos_F2.5-4_r0.75_H2",simu_id) |
                    grepl("lin_pos_F0-5_r1.0_H2",simu_id) ~ "quadratic",
                  grepl("abt",simu_id) |
                    grepl("lin_pos_F3-4_r1.0_H0.75",simu_id) ~ "abrupt"
                )) %>%
  # expected class
  # dplyr::group_by(length) %>%
  ggplot(aes(x=max_asd, fill=length))+
  # ggplot(aes(x=max_asd, fill=noise_lvl))+
  geom_histogram()+
  theme_light()+
  facet_grid(cols=vars(length), rows=vars(expected_class))+
  guides(fill="none")+
  ggtitle("Maximal (absolute) detection value from asdetect time series")

ggsave(p_distr_asd, filename="res/03_ricker_library/simu_library/max_asdetect_distrib_bestaic_madcstdefault.png",
       width=10, height=7)


# ROC curve ---------------------------------------------------------------

max_asd <- max_asdetect(list("l10"=outlist_aicasd_l10, "l20"=outlist_aicasd_l20,
                               "l25"=outlist_aicasd_l25,"l33"=outlist_aicasd_l33,
                               "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100))

roc_df <- max_asd %>%
  dplyr::mutate(length = factor(length, levels = c("l10","l20","l25","l33","l50","l100")),
                expected_class = dplyr::case_when(
                  grepl("cst",simu_id) ~ "no_change",
                  grepl("lin_pos_F1-5_r1.6_H0.75",simu_id) |
                    grepl("lin_neg_F5-1_r1.6_H0.75",simu_id) ~ "linear",
                  grepl("lin_pos_F2.5-4_r0.75_H2",simu_id) |
                    grepl("lin_pos_F0-5_r1.0_H2",simu_id) ~ "quadratic",
                  grepl("abt",simu_id) |
                    grepl("lin_pos_F3-4_r1.0_H0.75",simu_id) ~ "abrupt"
                ),
                expected_class = ifelse(expected_class=="abrupt", "abrupt", "not abrupt"),
                expected_class = factor(expected_class, levels=c("abrupt","not abrupt")))

outlist_aic <- combine_length(list(outlist_aic_l10,
                                   outlist_aic_l20,
                                   outlist_aic_l25,
                                   outlist_aic_l33,
                                   outlist_aic_l50,
                                   outlist_aic_l100))

roc_df2 <- roc_df %>%
  dplyr::left_join(outlist_aic %>%
                     dplyr::select(id, class),
                   by="id")

# For equal plan

# roc_df2 <- roc_df2 %>%
#   dplyr::filter(expected_class=="not abrupt") %>%
#   dplyr::group_by(length) %>%
#   dplyr::slice_sample(prop=1/3) %>%
#   dplyr::bind_rows(roc_df2 %>%
#                      dplyr::filter(expected_class=="abrupt"))

roc <- data.frame()

for (i in 0:101){

  # thres_df <- roc_df2 %>%
  #   dplyr::mutate(class_thr = ifelse(class=="abrupt" & max_asd >= i/100, "abrupt","not abrupt"),
  #                 class_thr = factor(class_thr, levels=c("abrupt","not abrupt")))

  thres_df <- roc_df2 %>%
    dplyr::mutate(class_thr = ifelse(max_asd >= i/100, "abrupt","not abrupt"),
                  class_thr = factor(class_thr, levels=c("abrupt","not abrupt")))

  # thres_df <- thres_df %>%
  #   split(.$length) %>%
  #   lapply(function(x)
  #     caret::confusionMatrix(data=x$class_thr, reference=x$expected_class))

  roc_line <- thres_df %>%
    split(.$length) %>%
    lapply(function(x){
      mat <- caret::confusionMatrix(data=x$class_thr, reference=x$expected_class)
      df <- data.frame(thres = i/100,
                       TP = mat$byClass[["Sensitivity"]],
                       FP = 1 - mat$byClass[["Specificity"]])}) %>%
    do.call(what="rbind") %>% tibble::rownames_to_column(var="length")



  # roc_line <- sapply(thres_df,
  #                    function(x){
  #                      data.frame(thres=i/100, TP=x$byClass[["Sensitivity"]],
  #                                           FP=1-x$byClass[["Specificity"]])},
  #                    simplify = TRUE
  #                    ) %>%
  #   t() %>% as.data.frame() %>% tibble::rownames_to_column(var="length")

  # perf <- caret::confusionMatrix(data=thres_df$class_thr, reference=thres_df$expected_class)

  # roc_line <- data.frame(thres=i/100,
  #                        TP=perf$byClass[["Sensitivity"]],
  #                        FP=1-perf$byClass[["Specificity"]])

  roc <- dplyr::bind_rows(roc, roc_line)

}

p_roc <- roc %>%
  dplyr::mutate(length = factor(length, levels=c("l10","l20","l25","l33","l50","l100"))) %>%
  ggplot(aes(x=FP, y=TP, color=length))+
  geom_line()+
  geom_point(alpha=0.2)+
  geom_abline(slope = 1, lty=2)+
  expand_limits(x=c(0,1), y=c(0,1))+
  ggrepel::geom_label_repel(aes(label = ifelse(thres%%0.2==0,
                                               as.character(thres), "")),
                   box.padding = 0.3,
                   segment.color = "grey50",
                   max.overlaps=100)+
  theme_light()

ggsave(plot=p_roc, filename = "res/03_ricker_library/simu_library/ROC_curve_labels.png",
       width=16, height=14, units="cm")

# library(plotly)

pl_roc <- roc %>%
  dplyr::mutate(length = factor(length, levels=c("l10","l20","l25","l33","l50","l100"))) %>%
  plotly::plot_ly(x=~FP, y=~TP, color=~length, type="scatter", text=~paste("thresh: ",thres))

htmltools::save_html(pl_roc, file = "res/03_ricker_library/simu_library/ROC_curve.html")


# ROC curve

roc_df <- readRDS("res/03_ricker_library/simu_library/roc_thr/roc_df_thr.rds") %>%
  dplyr::mutate(class = ifelse(class=="abrupt", "abrupt", "not abrupt"),
                expected_class = ifelse(expected_class=="abrupt", "abrupt", "not abrupt"),
                class = factor(class, levels=c("abrupt","not abrupt")),
                expected_class = factor(expected_class, levels=c("abrupt","not abrupt")))


roc_thr <- data.frame()

for (thr in seq(0,1,0.05)){

  roc_line <- roc_df %>%
    dplyr::filter(asd_thr==thr) %>%
    split(.$length) %>%
    lapply(function(x){
      mat <- caret::confusionMatrix(data=x$class, reference=x$expected_class)
      df <- data.frame(thres = thr,
                       TP = mat$byClass[["Sensitivity"]],
                       FP = 1 - mat$byClass[["Specificity"]])}) %>%
    do.call(what="rbind") %>% tibble::rownames_to_column(var="length")

  roc_thr <- dplyr::bind_rows(roc_thr, roc_line)

}


roc_aic_only <- combine_length(list(outlist_aic_l10,
                                    outlist_aic_l20,
                                    outlist_aic_l25,
                                    outlist_aic_l33,
                                    outlist_aic_l50,
                                    outlist_aic_l100)) %>%
  dplyr::mutate(class = ifelse(class=="abrupt", "abrupt", "not abrupt"),
                expected_class = ifelse(expected_class=="abrupt", "abrupt", "not abrupt"),
                class = factor(class, levels=c("abrupt","not abrupt")),
                expected_class = factor(expected_class, levels=c("abrupt","not abrupt"))) %>%
  split(.$length) %>%
  lapply(function(x){
    mat <- caret::confusionMatrix(data=x$class, reference=x$expected_class)
    df <- data.frame(thres = thr,
                     TP = mat$byClass[["Sensitivity"]],
                     FP = 1 - mat$byClass[["Specificity"]])}) %>%
  do.call(what="rbind") %>% tibble::rownames_to_column(var="length")





p_roc_thr <- roc_thr %>%
  dplyr::mutate(length = factor(length, levels=c("l10","l20","l25","l33","l50","l100"))) %>%
  ggplot(aes(x=FP, y=TP, color=length))+
  geom_line()+
  geom_point(alpha=0.2)+
  geom_abline(slope = 1, lty=2)+
  expand_limits(x=c(0,1), y=c(0,1))+
  ggrepel::geom_label_repel(aes(label = ifelse(thres%%0.2==0,
                                               as.character(thres), "")),
                            box.padding = 0.3,
                            segment.color = "grey50",
                            max.overlaps=100)+
  geom_point(data=roc_aic_only, aes(x=FP, y=TP, color=length), shape=2, size=4)+
  theme_light()
p_roc_thr


# Find tangent
roc_thr %>%
  split(.$length) %>%
  lapply(function(x) x %>% mutate(slope = (lag(TP)-TP)/(lag(FP)-FP)))




# Classification congruence across lengths --------------------------------

# Whole classification

outlist_aicasd <- combine_length(list(outlist_aicasd_l10,
                                         outlist_aicasd_l20,
                                         outlist_aicasd_l25,
                                         outlist_aicasd_l33,
                                         outlist_aicasd_l50,
                                         outlist_aicasd_l100))

out_length_wide <- outlist_aicasd %>%
  dplyr::select(id, length, class, expected_class) %>%
  dplyr::mutate(id = sub(".*?_", "", id),
                id = sub("_brk(\\d+)_", "_", id),
                length = factor(length, levels=c("l10","l20","l25","l33","l50","l100"))) %>%
  # dplyr::filter(expected_class=="abrupt") %>%
  dplyr::group_by(id) %>%
  tidyr::pivot_wider(names_from=length, values_from=class)


out_length_long <- outlist_aicasd %>%
  dplyr::select(id, length, class, expected_class, noise_lvl) %>%
  dplyr::mutate(id = sub(".*?_", "", id),
                id = sub("_brk(\\d+)_", "_", id),
                length = factor(length, levels=c("l10","l20","l25","l33","l50","l100"))) %>%
  # dplyr::filter(expected_class=="abrupt") %>%
  dplyr::mutate(id_short = paste(expected_class, as.integer(factor(id)))) %>%
  dplyr::arrange(expected_class, noise_lvl, id)


among_length <- catmaply::catmaply(
  out_length_long,
  x = id_short,
  # x_order = id_short,
  y = length,
  z = class
  )

htmltools::save_html(among_length, file = "res/03_ricker_library/simu_library/classif_among_length_test.html")

# Proportion of misclassified whatever the length
out_length_long %>%
  dplyr::filter(length != "l10") %>%
  dplyr::mutate(correct = expected_class==class) %>%
  dplyr::filter(correct==FALSE) %>%
  dplyr::group_by(id, expected_class) %>%
  dplyr::summarise(n=n()) %>%
  dplyr::filter(n==5) %>%
  dplyr::group_by(expected_class) %>%
  dplyr::summarise(n=n()/2600)

#   expected_class      n
#   <fct>           <dbl>
# 1 no_change      0.0273
# 2 linear         0.0769
# 3 quadratic      0.0204
# 4 abrupt         0.0385

# detection value
distr_asd <- max_asdetect(list("l10"=outlist_aicasd_l10, "l20"=outlist_aicasd_l20,
                               "l25"=outlist_aicasd_l25,"l33"=outlist_aicasd_l33,
                               "l50"=outlist_aicasd_l50, "l100"=outlist_aicasd_l100))

asd_thr <- 0.15

detect_length_long <- distr_asd %>%
  dplyr::mutate(expected_class = dplyr::case_when(
    grepl("cst",id) ~ "no_change",
    grepl("lin_pos_F1-5_r1.6_H0.75",id) |
      grepl("lin_neg_F5-1_r1.6_H0.75",id) ~ "linear",
    grepl("lin_pos_F2.5-4_r0.75_H2",id) |
      grepl("lin_pos_F0-5_r1.0_H2",id) ~ "quadratic",
    grepl("abt",id) |
      grepl("lin_pos_F3-4_r1.0_H0.75",id) ~ "abrupt"),
    id = sub(".*?_", "", id),
    id = sub("_brk(\\d+)_", "_", id),
    length = factor(length, levels=c("l10","l20","l25","l33","l50","l100")),
    max_asd = ifelse(max_asd>asd_thr, "abrupt", "not abrupt"),
    max_asd = factor(max_asd, levels = c("not abrupt", "abrupt"))
    ) %>%
  # dplyr::select(-id) %>%
  dplyr::filter(expected_class=="abrupt") %>%
  dplyr::arrange(expected_class, noise_lvl, id)

maxasd_among_length <- catmaply::catmaply(
  detect_length_long,
  x = id,
  x_order = id,
  y = length,
  y_order = length,
  z = max_asd
)

maxasd_among_length <- plotly::plot_ly(
  detect_length_long,
  x = ~id,
  # x_order = id,
  y = ~length,
  z = ~max_asd,
  type="heatmap"
)

htmltools::save_html(maxasd_among_length, file = "res/03_ricker_library/simu_library/maxasd_among_length.html")



# Figure 4 - Performances -------------------------------------------------

# aic_asd_mlt_nck_stz<1
outlist_l20 <- readRDS("res/03_ricker_library/simu_library/l20/outlist_l20_1_aic_asdmltnckstz_scl.rds")
outlist_l50 <- readRDS("res/03_ricker_library/simu_library/l50/outlist_l50_1_aic_asdmltnckstz_scl.rds")
outlist_l100 <- readRDS("res/03_ricker_library/simu_library/l100/outlist_l100_aic_asdmltnckstz_scl.rds")

mat_list <- c(outlist_l20$mat_list,
              outlist_l50$mat_list,
              outlist_l100$mat_list)

tot_mat <- mat_list[[1]]$conf_mat_class$table

for(i in 2:length(mat_list)){
  tot_mat <- tot_mat + mat_list[[i]]$conf_mat_class$table
}
byclass <- sum(tot_mat[,1])

mat <- pheatmap::pheatmap(tot_mat/byclass,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   legend = FALSE, display_numbers = TRUE,
                   number_format = "%.2f", fontsize = 170,
                   number_color="black", fontsize_number = 170)


sens_spec_df <-
  sapply(list("no_change", "linear", "quadratic", "abrupt"),
         function(i)
           mapply(
             function(j, k)
               sens_spec(
                 outlist = j,
                 workflow = k,
                 class = i
               ),
             list(
               outlist_l20,
               outlist_l50,
               outlist_l100
             ),
             c(rep("aic_asd_mltpl_stpsz", 3)),
             SIMPLIFY = FALSE
           ) %>% do.call(rbind, .),
         simplify = FALSE) %>% do.call(rbind, .)

sens_spec_plot <- sens_spec_df %>%
  dplyr::mutate(length = factor(length, levels=c("l100","l50","l20"))) %>%
  # dplyr::filter(class=="abrupt") %>%
  ggplot(aes(x=1-sens, y=1-spec))+
  geom_point(aes(shape=class, color=length), size=3, alpha=0.7)+
  # geom_point(aes(color=workflow, shape=length), size=3)+
  expand_limits(x=c(0,1), y=c(0,1))+
  # facet_wrap(vars(class))+
  labs(x="FNR (1-sensitivity)", y="FPR (1-specificity)")+
  geom_abline(slope=-1, intercept=1, lty=2)+
  theme_light()
sens_spec_plot


# Vert-Pre et al. 2012 "confusion matrix" ---------------------------------

vertpre_mat <- matrix(c(0.54,0.04,0.05,0.12,
                        0.14,0.81,0.33,0.13,
                        0.08,0.11,0.57,0.04,
                        0.24,0.04,0.04,0.71),
                      ncol=4, nrow=4, byrow=TRUE) %>%
  `rownames<-`(c("abundance", "regimes", "mixed", "random")) %>%
  `colnames<-`(c("abundance", "regimes", "mixed", "random"))

png("res/03_ricker_library/simu_library/vertpre_conf_mat.png", width = 3000, height = 3000)
pheatmap::pheatmap(vertpre_mat,
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   legend = FALSE, display_numbers = TRUE,
                   number_format = "%.2f", fontsize = 120,
                   number_color="black", fontsize_number = 120,
                   main="Vert-pre et al. confusion matrix")
dev.off()


# Time series parameters --------------------------------------------------

simu_list <- readRDS("data/03_simulations/all_simu_l100.rds")


df_list <- list()
for (i in 1:13){

  df_list <- c(df_list,
               simu_list[[i]] %>%
                 split(list(.$iter, .$scen))
  )
}

vals_simu <- summary_param_bis(df_list)

vals_simu_long <- vals_simu %>%
  tidyr::pivot_longer(cols=-ts, names_to="param", values_to = "val") %>%
  mutate(type = "simu")

ggplot(vals_long, aes(x=val))+
  geom_histogram()+
  theme_light()+
  facet_wrap(vars(param), scales = "free")
