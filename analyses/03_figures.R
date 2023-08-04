###-###-###-###-###-###-###-###-###
#
# 08/03/2023 mathieu.pelissie@ens-lyon.fr
#
# Make figures for paper 1
#
# fig_art1.R
#
###-###-###-###-###-###-###-###-###


library('tidyverse')
source("R/functions_simu.R")
source("R/functions_trajclass.R")
source("R/functions_output.R")


dir.create("analyses/figs", showWarnings = FALSE)
dir.create("analyses/figs/supp_mat", showWarnings = FALSE)


# Main text figures -------------------------------------------------------


# Figure 1 - Trajectory classification directly made on Inkscape ----------



# Figure 2 - Overview of simulated trajectories ---------------------------

# Initialization:
noise_df <- get_data("data/00_simu/noise_levels.csv")
param_df_fig2 <- get_data("data/00_simu/param_df_fig2.csv")
scen_fct <-  "make_scen"
titles <- c("No change","Linear decrease","Linear increase",
            "Convex","Accelerated decrease","Accelerated increase",
            "Concave","Decelerated decrease","Decelerated increase",
            "Abrupt decrease","Abrupt increase","Abrupt decrease")
gl <- list()
set.seed(1)

# Import noise values:
nz <- sapply(colnames(noise_df),
             function(i) noise_df %>%
               dplyr::pull(i) %>% `[`(5),
             simplify = FALSE,
             USE.NAMES = TRUE)

for (j in 1:nrow(param_df_fig2)){ # For each scenario j

  # Import parameter values:
  prm <- sapply(colnames(param_df_fig2),
                function(i) param_df_fig2 %>%
                  dplyr::pull(i) %>% `[`(j),
                simplify = FALSE,
                USE.NAMES = TRUE)

  # Generate mortality scenario timeseries:
  Fseq <- make_scen(shape=prm$shape, trend=prm$trend, min=prm$min,
                    max=prm$max, velocity=prm$velocity, l=prm$l,
                    breaktime = prm$breaktime)

  # Generate population timeseries:
  sim <- run_simu(name=Fseq$name, r=prm$r, F=Fseq$scen, Ts=length(Fseq$scen),
                  K=prm$K, P=prm$P, H=prm$H, iter=prm$iter, thr=prm$thr,
                  init=prm$init, expected_class=prm$expected_class,
                  sr=nz$sr, se=nz$se, su=nz$su, jfr=nz$jfr, jsz=nz$jsz)

  # Plot biomass timeseries:
  pbiom <- plot_simu_simple(sim$df, var="TB", tbcolor="black",
                            alpha=0.5, xname="Year", yname="Biomass")+
    theme(axis.title = element_text(size=18),
          axis.text = element_text(size=20),
          axis.line = element_line(colour = "black", linewidth = 1))

  # Plot exploitation scenario timeseries:
  pmort <- plot_simu_simple(sim$df, var="F", tbcolor="blue",
                            alpha=0.5, xname="Year", yname="Mortality",
                            ylim=8)+
    theme(axis.title = element_text(size=10),
          axis.text = element_text(size=7))

  # Combine both plots:
  p <- cowplot::ggdraw() +
    cowplot::draw_plot(pbiom, x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_plot(pmort, x = 0.62, y = .67, width = .35, height = .3) +
    cowplot::draw_text(titles[j], x = 0.55, y = 0.25,
                       hjust = 0.5, size = 20, fontface="bold")

  gl <- c(gl, list(p))

  # Save one example for figure 3:
  if(j==12){
    sim$df <- sim$df %>%
      dplyr::mutate(expected_traj=prm$expected_traj,
                    expected_class=prm$expected_class)
    sim_fig3 <- sim
    saveRDS(sim_fig3, "analyses/classif/ts_example/ts_example_abr.rds")
  }
}

# Keep four of the plots:
p <- cowplot::plot_grid(plotlist = gl[c(1,2,5,12)], nrow=2, ncol=2,
                        labels=LETTERS[1:4], label_size=25)

# Save the figure:
cowplot::save_plot(filename = "analyses/figs/fig2.pdf",
                   plot=p, base_height = 8, base_width = 12)

## Then arranged on Inkscape



# Figure 3 - Example of classification output -----------------------------

# Load an example timeseries:
sim_fig3 <- readRDS("analyses/classif/ts_example/ts_example_abr.rds")

# Reshape the data:
sets <- prep_data(df=sim_fig3$df, thr=-1, type="sim", apriori=TRUE)

# Run the classification:
trajs <- traj_class(sets, str="aic_asd", abr_mtd=c("chg", "asd"),
                    noise_comb=noise_comb, asd_chk=TRUE,
                    asd_thr=0.15, type="sim", showplots=TRUE,
                    apriori=TRUE, run_loo=TRUE, save_plot=FALSE,
                    two_bkps = TRUE, smooth_signif = TRUE,
                    outplot=TRUE, mad_thr=3, mad_cst=1.4826)

# Save the figure:
png(filename = "analyses/figs/fig3.png", width=8, height=6, units="in", res=300)
print(trajs$class_plot)
dev.off()

# Letters added with Inkscape



# Figure 4 - Confusion matrix ---------------------------------------------

lib_dir <- "analyses/classif/library/"

# Make confusion matrices for classification outputs with AICc only:
p100_aic <- conf_mat_pool(paste0(lib_dir, "outlist_l100_aic.rds"),
                          "100 timepoints")
p50_aic <- conf_mat_pool(paste0(lib_dir, "outlist_l50_1_aic.rds"),
                         "50 timepoints")
p20_aic <- conf_mat_pool(paste0(lib_dir, "outlist_l20_1_aic.rds"),
                         "20 timepoints")

# Make confusion matrices for classification outputs with AICc + asdetect:
p100_aicasd <- conf_mat_pool(paste0(lib_dir,
                                    "outlist_l100_aic_asd_thr0.15_looTRUE.rds"),
                             "")
p50_aicasd <- conf_mat_pool(paste0(lib_dir,
                                   "outlist_l50_1_aic_asd_thr0.15_looTRUE.rds"),
                            "")
p20_aicasd <- conf_mat_pool(paste0(lib_dir,
                                   "outlist_l20_1_aic_asd_thr0.15_looTRUE.rds"),
                            "")

# Save the figure:
pdf("analyses/figs/fig4.pdf", width = 180, height = 120)
cowplot::plot_grid(p100_aic$gtable, p50_aic$gtable, p20_aic$gtable,
                   p100_aicasd$gtable, p50_aicasd$gtable, p20_aicasd$gtable,
                   ncol = 3, labels=LETTERS[1:6], label_size=250)
dev.off()

# Make confusion matrices with colour scale lengend:
p100_aic_leg <- conf_mat_pool(paste0(lib_dir, "outlist_l100_aic.rds"),
                              "100 timepoints", show_legend = TRUE)
legend <- p100_aic_leg$gtable$grobs[[5]]

# Save the colour scale:
pdf("analyses/figs/fig4_leg.pdf", width = 2, height = 5)
cowplot::plot_grid("", "",
                   "", legend,
                   "", "", nrow=3, ncol=2) %>% suppressWarnings()
dev.off()

# Legend added with Inkscape and other minor changes made



# Figure 5 - Real world examples ------------------------------------------

# Load fisheries data:
path <- "data/01_RAMLDB/RAMLDB v4.61/R Data/"
load(paste0(path,"DBdata[asmt][v4.61].RData"))

# Select and reshape the data:
# ATBTUNAEATL: Atlantic bluefin tuna Eastern Atlantic
ts_type <- "TCbest"
thr <- 0
stk <- "ATBTUNAEATL"

set <- extract_RAM(stk, ts_type) %>%
  prep_data(thr=thr, type="RAM", apriori=FALSE)

# Run the classification:
output_RAM <- traj_class(sets=set, str="aic_asd", abr_mtd=c("chg", "asd"),
                         asd_thr=0.15, asd_chk=FALSE, type="RAM",
                         showplots=TRUE, apriori=FALSE, run_loo=TRUE,
                         save_plot=FALSE, two_bkps=TRUE, smooth_signif=TRUE,
                         outplot=TRUE, ind_plot="abt", dirname="analyses/")

# Plot the best trajectory shape:
p_RAM <- output_RAM$class_plot$plots[[1]] +
  theme(plot.background = element_rect(fill = "white"))


# Load Bird PECBMS index:
df_bird <- readr::read_csv(
  "data/02_PECBMS/europe-indicesandtrends-till2021_mod.csv")

# Select and reshape the data:
# Ortolan bunting (Emberiza hortulana) in Europe:
species <- unique(df_bird$PECBMS_species_name) %>% sort()

df_list_bird <- df_bird %>%
  dplyr::group_by(PECBMS_species_name) %>%
  dplyr::group_split() %>%
  stats::setNames(species)

df_list_hort <- df_list_bird["Emberiza_hortulana"]

# Run the classification:
output_bird <- run_classif_data(df_list_hort, min_len=20, str="aic_asd",
                           normalize=FALSE, showplots=TRUE, save_plot=FALSE,
                           run_loo=TRUE, two_bkps=TRUE, smooth_signif=TRUE,
                           group="PECBMS_species_name", time="Year",
                           variable="Index", outplot=TRUE, ind_plot="abt",
                           dirname="analyses/")

# Plot the best trajectory shape:
p_bird <- output_bird$outlist$Emberiza_hortulana$class_plot$plots[[1]] +
  theme(plot.background = element_rect(fill = "white"))


# Load Odonate index from Termaat et al. 2019:
df_odo <- readr::read_csv("data/03_Termaat2019/ddi12913-sup-0003-datas1.csv")%>%
  tidyr::pivot_longer(cols = (!region & !species & !nplot &
                         !trend & !trend_se & !trend_category),
               names_to = "year",
               values_to = "index") %>%
  dplyr::mutate(year = as.numeric(year)) %>%
  tidyr::drop_na()

# Select and reshape the data
# Small red-eyed damselfly (Erythromma viridulum) in Britain:
df_odo_brit <- df_odo %>%
  dplyr::filter(region %in% c("Britain"))

species <- unique(df_odo_brit$species) %>% sort()

df_list_odo <- df_odo_brit %>%
  dplyr::group_by(species) %>%
  dplyr::group_split() %>%
  stats::setNames(species)

df_list_erth <- df_list_odo["Erythromma viridulum"]

# Run the classification:
output_odo <- run_classif_data(df_list_erth, min_len=20, str="aic_asd",
                           normalize=FALSE, showplots=TRUE, save_plot=FALSE,
                           run_loo=TRUE, two_bkps=TRUE, smooth_signif=TRUE,
                           group="species", time="year", variable="index",
                           outplot=TRUE, ind_plot="abt", dirname="analyses/")

# Plot the best trajectory shape:
p_odo <- output_odo$outlist$`Erythromma viridulum`$class_plot$plots[[1]] +
  theme(plot.background = element_rect(fill = "white"))

# Combine plots and save the figure:
fig5 <- cowplot::plot_grid(p_RAM, p_bird, p_odo, nrow=1, ncol=3,
                        labels=c("A","B","C"), label_size=15, label_x = -0.0)
cowplot::save_plot(filename = "analyses/figs/fig5.pdf",
                   plot=fig5, base_height = 3, base_width = 12)

## Then completed on Inkscape with silhouettes



# Supplementary figures ---------------------------------------------------

# Figure S1 - Confidence matrices by noise AICc + asdetect ------------------

lib_dir <- "analyses/classif/library/"

# Load classification output:
outlist_aicasd_l100 <- readRDS(
  paste0(lib_dir, "outlist_l100_aic_asd_thr0.15_looTRUE.rds"))
outlist_aicasd_l50 <- readRDS(
  paste0(lib_dir,"outlist_l50_1_aic_asd_thr0.15_looTRUE.rds"))
outlist_aicasd_l20 <- readRDS(
  paste0(lib_dir,"outlist_l20_1_aic_asd_thr0.15_looTRUE.rds"))

outlist_aicasd_0.15 <- list("l100"=outlist_aicasd_l100,
                            "l50"=outlist_aicasd_l50,
                            "l20"=outlist_aicasd_l20)

# Make matrices:
all_aicasd <- conf_mat_noise(outlist_aicasd_0.15, "fig_S1_aicasd", save=FALSE)

# Save figure:
pdf(paste0("analyses/figs/supp_mat/fig_s1_aicasd.pdf"), width=40, height=100)
print(all_aicasd)
dev.off()

# Then adjusted with Inkscape



# Figure S2 - Confidence matrices by noise AICc only ------------------------

lib_dir <- "analyses/classif/library/"

# Load classification output:
outlist_aic_l100 <- readRDS(paste0(lib_dir,"outlist_l100_aic.rds"))
outlist_aic_l50 <- readRDS(paste0(lib_dir,"outlist_l50_1_aic.rds"))
outlist_aic_l20 <- readRDS(paste0(lib_dir,"outlist_l20_1_aic.rds"))

outlist_aic_0.15 <- list("l100"=outlist_aic_l100,
                         "l50"=outlist_aic_l50,
                         "l20"=outlist_aic_l20)

# Make matrices:
all_aic <- conf_mat_noise(outlist_aic_0.15, "fig_S2_aic", save=FALSE)

# Save figure:
pdf(paste0("analyses/figs/supp_mat/fig_s2_aic.pdf"), width = 40, height = 100)
print(all_aic)
dev.off()

# Then adjusted with Inkscape



# Figure S3 - Example classification by noise AICc + asdetect -------------

# Load simulated timeseries
simu_list <- readRDS("data/00_simu/all_simu_l100.rds")

class_list <- c("nch","lin","pol")
id_list <- c(1, 201, 401)

# Classify examples for non-abrupt simulated timeseries:
no_abt <-
  mapply(function(class, id){

    p <- lapply(1:length(simu_list),
                function(i){
                  noise_comb <- sub(".*?_", "", names(simu_list)[i])

                  sets <- prep_data(df=simu_list[[i]], thr=-1,
                                    type="sim", apriori=TRUE)
                  sets <- list("ts"=c(sets$ts[id]), "ts_type"=sets$ts_type)
                  trajs <- traj_class(sets, str="aic_asd",
                                      abr_mtd=c("chg", "asd"),
                                      noise_comb=noise_comb, asd_chk=TRUE,
                                      asd_thr=0.15, type="sim", showplots=TRUE,
                                      apriori=TRUE, run_loo=FALSE,
                                      two_bkps = TRUE, smooth_signif = TRUE,
                                      outplot=TRUE, ind_plot=class)
                  trajs$class_plot[[1]]
                }
    )

    pl <- (cowplot::plot_grid(plotlist = p[c(1,4,7,10)], nrow=1) /
             cowplot::plot_grid(plotlist = p[c(2,5,8,11)], nrow=1) /
             cowplot::plot_grid(plotlist = p[c(3,6,9,12)], nrow=1))

    title <- cowplot::ggdraw() + cowplot::draw_label(class, size=30)

    pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))

  }, class_list, id_list, SIMPLIFY=FALSE)


# Classify examples for abrupt simulated timeseries:
abt <- lapply(1:length(simu_list),
              function(i){
                noise_comb <- sub(".*?_", "", names(simu_list)[i])

                sets <- prep_data(df=simu_list[[i]], thr=-1,
                                  type="sim", apriori=TRUE)
                sets <- list("ts"=c(sets$ts[601]), "ts_type"=sets$ts_type)
                trajs <- traj_class(sets, str="aic_asd",
                                    abr_mtd=c("chg", "asd"),
                                    noise_comb=noise_comb, asd_chk=TRUE,
                                    asd_thr=0.15, type="sim", showplots=TRUE,
                                    apriori=TRUE, run_loo=FALSE,
                                    two_bkps = TRUE, smooth_signif = TRUE,
                                    outplot=TRUE, ind_plot="abt")
                trajs$class_plot$plots[[1]]
              }
)

abt_pl <- (cowplot::plot_grid(plotlist = abt[c(1,4,7,10)], nrow=1) /
             cowplot::plot_grid(plotlist = abt[c(2,5,8,11)], nrow=1) /
             cowplot::plot_grid(plotlist = abt[c(3,6,9,12)], nrow=1))

title <- cowplot::ggdraw() + cowplot::draw_label("abt", size=30)

abt_pl2 <- cowplot::plot_grid(title, abt_pl, ncol=1, rel_heights=c(0.1, 1))

# Combine non-abrupt and abrupt plots:
fig_s3 <- cowplot::plot_grid(no_abt[[1]]) /
  cowplot::plot_grid(no_abt[[2]]) /
  cowplot::plot_grid(no_abt[[3]]) /
  cowplot::plot_grid(abt_pl2)

# Save the figure:
pdf(paste0("analyses/figs/supp_mat/fig_s3.pdf"), width = 14, height = 28)
print(fig_s3)
dev.off()

# Then adjusted with Inkscape



# Figure S4 - ROC curve to adjust asdetect threshold -----------------------

dir.create("analyses/figs/supp_mat/fig_s4", showWarnings = FALSE)

# Load classification for different threshold values:
roc_df <- readRDS("analyses/classif/roc/roc_df_thr.rds") %>%
  # Simplify classification as abrupt/not abrupt:
  dplyr::mutate(class = ifelse(class=="abrupt", "abrupt", "not abrupt"),
                expected_class = ifelse(expected_class=="abrupt",
                                        "abrupt", "not abrupt"),
                class = factor(class, levels=c("abrupt","not abrupt")),
                expected_class = factor(expected_class,
                                        levels=c("abrupt","not abrupt")))

# Make confusion matrices for each length and threshold value:
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


# Load classification output for AICc only:
outlist_aic_l20 <- readRDS("analyses/classif/library/outlist_l20_1_aic.rds")
outlist_aic_l25 <- readRDS("analyses/classif/library/outlist_l25_1_aic.rds")
outlist_aic_l33 <- readRDS("analyses/classif/library/outlist_l33_1_aic.rds")
outlist_aic_l50 <- readRDS("analyses/classif/library/outlist_l50_1_aic.rds")
outlist_aic_l100 <- readRDS("analyses/classif/library/outlist_l100_aic.rds")

# Make confusion matrices for each length:
roc_aic_only <- combine_length(list(outlist_aic_l20,
                                    outlist_aic_l25,
                                    outlist_aic_l33,
                                    outlist_aic_l50,
                                    outlist_aic_l100)) %>%
  dplyr::mutate(class = ifelse(class=="abrupt", "abr.", "not abr."),
                expected_class = ifelse(expected_class=="abrupt",
                                        "abr.", "not abr."),
                class = factor(class, levels=c("abr.","not abr.")),
                expected_class = factor(expected_class,
                                        levels=c("abr.","not abr.")),
                length = factor(length,
                                levels=c("l20","l25","l33","l50","l100"))) %>%
  split(.$length) %>%
  lapply(function(x){
    mat <- caret::confusionMatrix(data=x$class, reference=x$expected_class)
    df <- data.frame(TP = mat$byClass[["Sensitivity"]],
                     FP = 1 - mat$byClass[["Specificity"]])

    return(list(mat, df))
  })

# Summarise TPR and FPR into a data frame:
roc_aic_only_df <- roc_aic_only %>%
  lapply(function(x) x %>% `[[`(2)) %>%
  do.call(what="rbind") %>%
  tibble::rownames_to_column(var="length")

# Plot confusion matrices for each length:
roc_aic_only_conf_mat <- roc_aic_only %>%
  lapply(function(x){
    mat <- x %>%
      `[[`(1) %>%
      `[[`("table")

    mat[,1] <- mat[,1]/2400
    mat[,2] <- mat[,2]/7200

    # Reverse for consistency with full confusion matrix:
    mat <- mat[nrow(mat):1, nrow(mat):1]

    mat %>%
      pheatmap::pheatmap(
        color = grDevices::colorRampPalette(
          rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
        breaks = seq(0, 1, .01),
        cluster_rows = FALSE, cluster_cols = FALSE,
        legend = FALSE, display_numbers = TRUE,
        number_format = "%.2f", fontsize = 40,
        number_color="black", fontsize_number = 50)

  }
  )

# Load classification output for AICc + asdetect:

lib_dir <- "analyses/classif/library/"

outlist_aicasd_l20 <- readRDS(
  "analyses/classif/library/outlist_l20_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l25 <- readRDS(
  "analyses/classif/library/outlist_l25_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l33 <- readRDS(
  "analyses/classif/library/outlist_l33_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l50 <- readRDS(
  "analyses/classif/library/outlist_l50_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l100 <- readRDS(
  "analyses/classif/library/outlist_l100_aic_asd_thr0.15_looTRUE.rds")

# Make confusion matrices for each length:
roc_aicasd <- combine_length(list(outlist_aicasd_l20,
                                  outlist_aicasd_l25,
                                  outlist_aicasd_l33,
                                  outlist_aicasd_l50,
                                  outlist_aicasd_l100)) %>%
  dplyr::mutate(class = ifelse(class=="abrupt", "abr.", "not abr."),
                expected_class = ifelse(expected_class=="abrupt",
                                        "abr.", "not abr."),
                class = factor(class, levels=c("abr.","not abr.")),
                expected_class = factor(expected_class,
                                        levels=c("abr.","not abr.")),
                length = factor(length,
                                levels=c("l20","l25","l33","l50","l100"))) %>%
  split(.$length) %>%
  lapply(function(x){
    mat <- caret::confusionMatrix(data=x$class, reference=x$expected_class)
    df <- data.frame(TP = mat$byClass[["Sensitivity"]],
                     FP = 1 - mat$byClass[["Specificity"]])

    return(list(mat, df))
  })

# Plot confusion matrices for each length:
roc_aicasd_conf_mat <- roc_aicasd %>%
  lapply(function(x){
    mat <- x %>%
      `[[`(1) %>%
      `[[`("table")

    mat[,1] <- mat[,1]/2400
    mat[,2] <- mat[,2]/7200

    # Reverse for consistency with full confusion matrix:
    mat <- mat[nrow(mat):1, nrow(mat):1]

    mat %>%
      pheatmap::pheatmap(
        color = grDevices::colorRampPalette(
          rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
        breaks = seq(0, 1, .01),
        cluster_rows = FALSE, cluster_cols = FALSE,
        legend = FALSE, display_numbers = TRUE,
        number_format = "%.2f", fontsize = 40,
        number_color="black", fontsize_number = 50)

  }
  )

# Save confusion matrices:
pdf("analyses/figs/supp_mat/fig_s4/fig_s4cd.pdf", width = 30, height = 12)

cowplot::plot_grid(roc_aic_only_conf_mat$l20$gtable,
                   roc_aic_only_conf_mat$l25$gtable,
                   roc_aic_only_conf_mat$l33$gtable,
                   roc_aic_only_conf_mat$l50$gtable,
                   roc_aic_only_conf_mat$l100$gtable,
                   roc_aicasd_conf_mat$l20$gtable,
                   roc_aicasd_conf_mat$l25$gtable,
                   roc_aicasd_conf_mat$l33$gtable,
                   roc_aicasd_conf_mat$l50$gtable,
                   roc_aicasd_conf_mat$l100$gtable,
                   ncol = 5, byrow=TRUE)
dev.off()


# Remove failing limit cases (thr=0 and thr=0.05 for lengths 100 and 50):
roc_thr <- roc_thr[-c(1,6,5,10),]

# Plot and save complete ROC curve:
p_roc_thr <- roc_thr %>%
  dplyr::mutate(length = factor(length,
                                levels=c("l20","l25","l33","l50","l100"))) %>%
  ggplot(aes(x=FP, y=TP, color=length))+
  geom_line()+
  geom_point(alpha=0.2)+
  geom_abline(slope = 1, lty=2)+
  expand_limits(x=c(0,1), y=c(0,1))+
  geom_point(data=roc_aic_only_df, aes(x=FP, y=TP, color=length),
             shape=2, size=4)+
  theme_light()+
  theme(legend.position = "none")+
  labs(x="False Positive rate", y="True Positive rate")

ggsave(filename = "analyses/figs/supp_mat/fig_s4/fig_s4a.pdf",
       width=10, height=10,
       plot = p_roc_thr)

# Plot and save a close-up of the ROC curve with threshold values annotated:
p_roc_thr_zoom <- p_roc_thr +
  coord_cartesian(xlim=c(0,0.25), ylim=c(0.75, 1)) +
  theme(legend.position = "none", panel.grid.minor = element_blank())+
  ggrepel::geom_label_repel(aes(label = ifelse(thres%%0.1==0,
                                               as.character(thres), "")),
                            box.padding = 0.3,
                            segment.color = "grey50",
                            max.overlaps=100)

ggsave(filename = "analyses/figs/supp_mat/fig_s4/fig_s4b.pdf",
       width=4, height=4, plot = p_roc_thr_zoom)

# All combined with Inkscape


# Figure S5 - Distribution of classification scores -----------------------

# Load classification output:
outlist_aicasd_l100 <- readRDS(
  "analyses/classif/library/outlist_l100_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l50 <- readRDS(
  "analyses/classif/library/outlist_l50_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l33 <- readRDS(
  "analyses/classif/library/outlist_l33_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l25 <- readRDS(
  "analyses/classif/library/outlist_l25_1_aic_asd_thr0.15_looTRUE.rds")
outlist_aicasd_l20 <- readRDS(
  "analyses/classif/library/outlist_l20_1_aic_asd_thr0.15_looTRUE.rds")

# Combine classification output into a single dataframe:
outlist_aicasd <- combine_length(list(outlist_aicasd_l20,
                                      outlist_aicasd_l25,
                                      outlist_aicasd_l33,
                                      outlist_aicasd_l50,
                                      outlist_aicasd_l100))

# Reshape scores:
loo_wAICc_nmrse <- loo_wAICc_nmrse_long(outlist_aicasd)

loo_wAICc_nmrse_obs <- loo_wAICc_nmrse %>%
  dplyr::filter(class_tested==class) %>%
  dplyr::mutate(correct = ifelse(class==expected_class, TRUE, FALSE))

loo_wAICc_nmrse_obs_long <- loo_wAICc_nmrse_obs %>%
  dplyr::rename(wAICc = weight_aic, LOO = loo, NRMSE = nrmse) %>%
  tidyr::pivot_longer(cols = c("wAICc", "LOO", "NRMSE"),
                      names_to = "metric", values_to = "value") %>%
  dplyr::mutate(metric = factor(metric, levels=c("wAICc", "LOO", "NRMSE")))

# Compute deciles:
deciles <- loo_wAICc_nmrse_obs_long %>%
  dplyr::filter(class_tested==class) %>%
  dplyr::mutate(correct = ifelse(class==expected_class, TRUE, FALSE),
                correct = factor(correct, levels = c("TRUE", "FALSE"))) %>%
  dplyr::group_by(metric, class_tested, correct) %>%
  dplyr::summarise(lower = quantile(value, probs = .1),
                   upper = quantile(value, probs = .9))

# Plot score distributions:
fig_s5 <- ggplot(loo_wAICc_nmrse_obs_long %>%
                   dplyr::filter(class_tested==class) %>%
                   dplyr::mutate(correct = ifelse(class==expected_class,
                                                  TRUE, FALSE)),
                 aes(x = value, fill = correct)) +
  geom_histogram(alpha = 0.4, binwidth = 0.1, position="identity") +
  theme_classic() +
  facet_grid(rows=vars(metric), cols=vars(class_tested), scales="free_x")+
  theme(legend.position = "none")+
  geom_vline(data = deciles, aes(xintercept=lower, col=correct, lty=correct)) +
  geom_vline(data = deciles, aes(xintercept=upper, col=correct, lty=correct))+
  scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
  scale_colour_manual(values = c("TRUE" = "#00bfc4", "FALSE" = "#f8766d"))

# Save the figure:
ggsave(filename = "analyses/figs/supp_mat/fig_s5.png",
       plot = fig_s5, width=8, height=6)



# Figure S6 - Impact of data transformation -------------------------------

dir.create("analyses/figs/supp_mat/fig_s6", showWarnings = FALSE)

# Make confusion matrices:

# For untransformed 'raw' timeseries:
outlist_aic_raw_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_raw_thr0.15.rds")
outlist_aic_raw_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_raw_thr0.15.rds")
outlist_aic_raw_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_raw_thr0.15.rds")
outlist_aic_raw_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_raw_thr0.15.rds")
outlist_aic_raw_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_raw_thr0.15.rds")

outlist_aic_raw <- combine_length(list(outlist_aic_raw_l20, outlist_aic_raw_l25,
                                       outlist_aic_raw_l33, outlist_aic_raw_l50,
                                       outlist_aic_raw_l100))
conf_mat_raw <- conf_mat(outlist_aic_raw, "no_transfo",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)

# For timeseries scaled by mean value:
outlist_aic_msz_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_msz_thr0.15.rds")
outlist_aic_msz_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_msz_thr0.15.rds")
outlist_aic_msz_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_msz_thr0.15.rds")
outlist_aic_msz_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_msz_thr0.15.rds")
outlist_aic_msz_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_msz_thr0.15.rds")

outlist_aic_msz <- combine_length(list(outlist_aic_msz_l20, outlist_aic_msz_l25,
                                       outlist_aic_msz_l33, outlist_aic_msz_l50,
                                       outlist_aic_msz_l100))
conf_mat_msz <- conf_mat(outlist_aic_msz, "divide_by_mean",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)

# For rescaled timeseries:
outlist_aic_nrm_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_nrm_thr0.15.rds")
outlist_aic_nrm_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_nrm_thr0.15.rds")
outlist_aic_nrm_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_nrm_thr0.15.rds")
outlist_aic_nrm_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_nrm_thr0.15.rds")
outlist_aic_nrm_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_nrm_thr0.15.rds")

outlist_aic_nrm <- combine_length(list(outlist_aic_nrm_l20, outlist_aic_nrm_l25,
                                       outlist_aic_nrm_l33, outlist_aic_nrm_l50,
                                       outlist_aic_nrm_l100))
conf_mat_nrm <- conf_mat(outlist_aic_nrm, "normalize_btwn0and1",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)

# For standardized timeseries:
outlist_aic_std_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_std_thr0.15.rds")
outlist_aic_std_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_std_thr0.15.rds")
outlist_aic_std_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_std_thr0.15.rds")
outlist_aic_std_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_std_thr0.15.rds")
outlist_aic_std_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_std_thr0.15.rds")

outlist_aic_std <- combine_length(list(outlist_aic_std_l20, outlist_aic_std_l25,
                                       outlist_aic_std_l33, outlist_aic_std_l50,
                                       outlist_aic_std_l100))
conf_mat_std <- conf_mat(outlist_aic_std, "standardize",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)

# For timeseries after square root transformation:
outlist_aic_sqr_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_sqr_thr0.15.rds")
outlist_aic_sqr_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_sqr_thr0.15.rds")
outlist_aic_sqr_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_sqr_thr0.15.rds")
outlist_aic_sqr_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_sqr_thr0.15.rds")
outlist_aic_sqr_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_sqr_thr0.15.rds")

outlist_aic_sqr <- combine_length(list(outlist_aic_sqr_l20, outlist_aic_sqr_l25,
                                       outlist_aic_sqr_l33, outlist_aic_sqr_l50,
                                       outlist_aic_sqr_l100))
conf_mat_sqr <- conf_mat(outlist_aic_sqr, "squareroot",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)

# For timeseries after logarithm transformation:
outlist_aic_log_l20 <- readRDS(
  "analyses/classif/data_transf/outlist_l20_1_aic_asd_log_thr0.15.rds")
outlist_aic_log_l25 <- readRDS(
  "analyses/classif/data_transf/outlist_l25_1_aic_asd_log_thr0.15.rds")
outlist_aic_log_l33 <- readRDS(
  "analyses/classif/data_transf/outlist_l33_1_aic_asd_log_thr0.15.rds")
outlist_aic_log_l50 <- readRDS(
  "analyses/classif/data_transf/outlist_l50_1_aic_asd_log_thr0.15.rds")
outlist_aic_log_l100 <- readRDS(
  "analyses/classif/data_transf/outlist_l100_aic_asd_log_thr0.15.rds")

outlist_aic_log <- combine_length(list(outlist_aic_log_l20, outlist_aic_log_l25,
                                       outlist_aic_log_l33, outlist_aic_log_l50,
                                       outlist_aic_log_l100))
conf_mat_log <- conf_mat(outlist_aic_log, "logarithm",
                         dirname="analyses/figs/supp_mat/fig_s6", save=TRUE)


# Plot examples of transformed timeseries:
transform <- c("raw","nrm","msz","std","log","sqr")

for (j in 1:length(transform)){ # for all 5 transformations + 'raw'

  # Load the same timeseries and noise level combination:
  sim_fig3 <- readRDS("analyses/classif/ts_example/ts_example_abr.rds")
  noise_comb <- "sr0.025_se0.025_su0_jsz0_jfr0.2"

  if (transform[j] == "nrm"){ # for rescaling

    sim_fig3$df <- sim_fig3$df %>%
      dplyr::group_by(scen, iter) %>%
      dplyr::mutate(
        TB = (TB-min(TB))/(max(TB)-min(TB))
      ) %>%
      dplyr::ungroup()

  } else if (transform[j] == "msz"){ # for scaling by mean value

    sim_fig3$df <- sim_fig3$df %>%
      dplyr::group_by(scen, iter) %>%
      dplyr::mutate(
        TB = TB/mean(TB)
      ) %>%
      dplyr::ungroup()

  } else if (transform[j] == "std"){ # for standadization

    sim_fig3$df <- sim_fig3$df %>%
      dplyr::group_by(scen, iter) %>%
      dplyr::mutate(
        TB = (TB-mean(TB))/sd(TB)
      ) %>%
      dplyr::ungroup()

  } else if (transform[j] == "log"){ # for log transformation

    sim_fig3$df <- sim_fig3$df %>%
      dplyr::group_by(scen, iter) %>%
      dplyr::mutate(
        TB = log(1+TB)
      ) %>%
      dplyr::ungroup()

  } else if (transform[j] == "sqr"){ # for square root transformation

    sim_fig3$df <- sim_fig3$df %>%
      dplyr::group_by(scen, iter) %>%
      dplyr::mutate(
        TB = sqrt(TB)
      ) %>%
      dplyr::ungroup()
  }

  # Reshape the timeseries:
  sets <- prep_data(df=sim_fig3$df, thr=NULL, type="sim", apriori=TRUE)

  # Run the classification:
  trajs <- traj_class(sets, str="aic_asd", abr_mtd=c("chg", "asd"),
                      asd_thr=0.15, noise_comb=noise_comb, asd_chk=TRUE,
                      type="sim", showplots=TRUE, apriori=TRUE, run_loo=TRUE,
                      two_bkps = TRUE, smooth_signif = TRUE, outplot=TRUE,
                      plot_one_in=1, mad_thr=3, mad_cst=1.4826, save_plot=FALSE)

  # Save plots:
  png(filename = paste0("analyses/figs/supp_mat/fig_s6/ts_example_",
                        transform[j], "_aicasd_thr0.15.png"),
      width=8, height=6, units="in", res=300)
  print(trajs$class_plot)
  dev.off()

}

# Then combined with Inkscape



# Figure S7 - Multi break point examples -----------------------------------

dir.create("analyses/figs/supp_mat/fig_s7", showWarnings = FALSE)
dir.create("analyses/figs/supp_mat/fig_s7/abrupt", showWarnings = FALSE)

# Load timeseries with two breakpoints:
simu_list_mltbrk <- readRDS("data/00_simu/all_simu_l100_multibrk.rds")

# Focus on one noise level combination:
i <- 4
noise_comb <- sub(".*?_", "", names(simu_list_mltbrk)[i])

# Reshape the timeseries:
sets <- prep_data(df=simu_list_mltbrk[[i]], thr=-1, type="sim", apriori=TRUE)

# Select only 12 timeseries with different trajectories:
subsets <- list("ts"=c(sets$ts[1], sets$ts[26], sets$ts[51], sets$ts[76],
                    sets$ts[101], sets$ts[126], sets$ts[151], sets$ts[176],
                    sets$ts[201], sets$ts[226], sets$ts[251], sets$ts[276]),
             "ts_type"=sets$ts_type)


# Run the classification:
trajs <- traj_class(subsets, str="aic_asd", abr_mtd=c("chg", "asd"),
                    noise_comb=noise_comb, asd_chk=FALSE, asd_thr=0.15,
                    type="sim", showplots=TRUE, apriori=TRUE, run_loo=TRUE,
                    two_bkps = TRUE, smooth_signif = TRUE, edge_lim=5,
                    congr_brk=5, outplot=TRUE, save_plot=TRUE, plot_one_in=1,
                    save_plot_bis=FALSE,
                    dirname="analyses/figs/supp_mat/fig_s7/")

# Then combined with Inkscape



# Figure S8 - Additional empirical examples -------------------------------

# Load fisheries data:
path <- "data/01_RAMLDB/RAMLDB v4.61/R Data/"
load(paste0(path,"DBdata[asmt][v4.61].RData"))

# Select and reshape the data:
# COD3NO - Atlantic cod (Gadus morhua) - Southern Grand Banks:
ts_type <- "TCbest"
thr <- 0
stk <- "COD3NO"
str <- "aic_asd"

# Run the classification:
set <- extract_RAM(stk, ts_type) %>%
  prep_data(thr=thr, type="RAM", apriori=FALSE)

# Run the classification:
output_RAM <- traj_class(sets=set, str=str, abr_mtd=c("chg", "asd"),
                         asd_thr=0.15, asd_chk=FALSE, type="RAM",
                         showplots=TRUE, apriori=FALSE, run_loo=TRUE,
                         mad_cst=1, save_plot = FALSE, two_bkps=TRUE,
                         smooth_signif=TRUE, outplot=TRUE,ind_plot=NULL,
                         dirname="analyses/")

# Plot the trajectory shapes:
p_RAM <- output_RAM$class_plot+
  theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("Atlantic cod (Gadus morhua) catch in Southern Grand Banks (Canada)")


# Load Bird PECBMS index
df_bird <- readr::read_csv(
  "data/02_PECBMS/europe-indicesandtrends-till2021_mod.csv")

# Select and reshape the data:
# European greenfinch (Chloris chloris) in Europe:
species <- unique(df_bird$PECBMS_species_name) %>% sort()

df_list_bird <- df_bird %>%
  dplyr::group_by(PECBMS_species_name) %>%
  dplyr::group_split() %>%
  stats::setNames(species)

df_list_grf <- df_list_bird["Chloris_chloris"]

# Run the classification:
output_bird <- run_classif_data(df_list_grf, min_len=20, str="aic_asd",
                                normalize=FALSE, showplots=TRUE,
                                run_loo=TRUE, two_bkps=TRUE, smooth_signif=TRUE,
                                group="PECBMS_species_name", time="Year",
                                variable="Index", outplot=TRUE, ind_plot=NULL,
                                save_plot = FALSE, dirname = "analyses/")

# Plot the trajectory shapes:
p_bird <- output_bird$outlist$Chloris_chloris$class_plot+
  theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("European greenfinch (Chloris chloris) in Europe")


# Load Odonate index from Termaat et al. 2019:
df_odo <- readr::read_csv("data/03_Termaat2019/ddi12913-sup-0003-datas1.csv")%>%
  tidyr::pivot_longer(cols = (!region & !species & !nplot &
                         !trend & !trend_se & !trend_category),
               names_to = "year",
               values_to = "index") %>%
  dplyr::mutate(year = as.numeric(year)) %>%
  tidyr::drop_na()

# Select and reshape the data:
# White-megged damselfly (Platycnemis pennipes) in Britain
df_odo_brit <- df_odo %>%
  dplyr::filter(region %in% c("Britain"))

species <- unique(df_odo_brit$species) %>% sort()

df_list_odo <- df_odo_brit %>%
  dplyr::group_by(species) %>%
  dplyr::group_split() %>%
  stats::setNames(species)

df_list_plty <- df_list_odo["Platycnemis pennipes"]

# Run the classification:
output_odo <- run_classif_data(df_list_plty, min_len=20, str="aic_asd",
                               normalize=FALSE, showplots=TRUE,
                               run_loo=TRUE, two_bkps=TRUE, smooth_signif=TRUE,
                               group="species", time="year", variable="index",
                               outplot=TRUE, ind_plot=NULL,
                               save_plot = FALSE, dirname = "analyses/")

# Plot the trajectory shapes:
p_odo <- output_odo$outlist$`Platycnemis pennipes`$class_plot+
  theme_update(plot.title = element_text(hjust = 0.5))+
  ggtitle("White-megged damselfly (Platycnemis pennipes) in Britain")


# Combine plots and save the figure:
fig_s8 <- cowplot::plot_grid(p_RAM, p_bird, p_odo, nrow=3, ncol=1,
                             labels=c("A","B","C"), label_size=15, label_x=-0.0)
cowplot::save_plot(filename = "analyses/figs/supp_mat/fig_s8.pdf",
                   plot=fig_s8, base_height = 20, base_width = 10)

## Then completed on Inkscape with silhouettes



# Figure S9 - Scenario space ----------------------------------------------

library(patchwork)

# Bifurcation diagram for h=0.75:
bifurc_2D_h0.75_p2 <- bifurc_analysis(Ts=100, K=10, H=0.75, P=2,
                                      rmax=2.0, rstep=0.1,
                                      Fmax=7.5, Fstep=0.1,
                                      thr=0, init=0.9)

# Bifurcation diagram for h=2:
bifurc_2D_h2_p2 <- bifurc_analysis(Ts=100, K=10, H=2, P=2,
                                   rmax=2.0, rstep=0.1,
                                   Fmax=7.5, Fstep=0.1,
                                   thr=0, init=0.9)

# Combine and save plot:
fig_s9 <- cowplot::plot_grid(bifurc_2D_h0.75_p2$plot_inc,
                             bifurc_2D_h2_p2$plot_inc, nrow=1, ncol=2,
                            labels=c("A","B"), label_size=15, label_x = -0.0)

cowplot::save_plot(filename = "analyses/figs/supp_mat/fig_s9.pdf",
                   plot=fig_s9, base_height = 4, base_width = 10)

# Then edited with Inkscape
