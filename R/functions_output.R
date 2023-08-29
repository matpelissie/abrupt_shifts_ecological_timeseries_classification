###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# Functions to analyse outputs
#
# functions_output.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###


#' Combine classification output from different length
#'
#' @param outlist_all list of "outlists" corresponding to different
#' timeseries lengths
#'
#' @return data frame combining classification output info from timeseries
#' of all lengths
#'
#' @export

combine_length <- function(outlist_all){

  outlist <- lapply(outlist_all, # for different timeseries lengths
                    function(z) mapply( # for different noise combinations
                      function(x,y) x$best_traj %>%
                        dplyr::mutate(length = stringr::str_extract(y, "[^_]+"),
                                      noise_lvl = sub(".*?\\_", "", y)),
                      z$traj_list, names(z$traj_list), SIMPLIFY=FALSE
                    ) %>%
                      dplyr::bind_rows() %>%
                      tibble::remove_rownames()
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(id = paste(simu_id, noise_lvl, sep="_"))

  return(outlist)
}



#' Make a confusion matrix
#'
#' @param combine_length_out "outlist" from the traj_class function
#'  corresponding to a given timeseries lengths or combine_length function
#' @param save logical, should the figure be saved automatically?
#' @param dirname directory name where to save plots
#' @param filename character to name the file saved
#'
#' @return confusion matrix as heatmap
#'
#' @export

conf_mat <- function(combine_length_out, save, dirname, filename){

  # Make the confusion matrix from the data:
  mat <- caret::confusionMatrix(data=combine_length_out$class,
                                reference=combine_length_out$expected_class) %>%
    `[`("table") %>%
    `[[`(1) %>%
    `colnames<-`(sub("_"," ", colnames(.))) %>%
    `rownames<-`(sub("_"," ", rownames(.)))

  # Normalize to get proportions:
  mat <- mat/(nrow(combine_length_out)/4)

  # Make the heatmap:
  plot <- pheatmap::pheatmap(mat,
                             color = grDevices::colorRampPalette(
                               rev(RColorBrewer::brewer.pal(
                                 n = 7, name ="RdYlBu")))(100),
                             breaks = seq(0, 1, .01),
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             legend = FALSE, display_numbers = TRUE,
                             number_format = "%.2f", fontsize = 30,
                             number_color="black", fontsize_number = 30)

  # Save the heatmap:
  if (save){
    png(paste0(dirname, "/confmat_", filename ,".png"),
        width = 800, height = 800)
    print(plot)
    dev.off()
    }

  return(plot)
}



#' Make confusion matrices by timeseries length and noise combination
#'
#' @param outlist list of "outlist" from the traj_class function
#'  corresponding to  timeseries lengths or combine_length function
#' @param save logical, should the figure be saved automatically?
#' @param dirname directory name where to save plots
#' @param filename character to name the file saved
#'
#' @return patchwork of cowplots
#'
#' @export


conf_mat_noise <- function(outlist, save, dirname, filename){

  # Make the confusion matrices from the data:
  mat <- lapply(outlist,
                function(x) # for different timeseries lengths

                  lapply(1:length(x$mat_list),
                         function(i){ # for different noise combinations

                    byclass <- sum(x$mat_list[[i]]$conf_mat_class$table[,1])
                    title_list <- c("sr = 0.001","","",
                                    "sr = 0.025","","",
                                    "sr = 0.050","","",
                                    "sr = 0.075","","")

                    mat2 <- x$mat_list[[i]]$conf_mat_class$table %>%
                      `colnames<-`(sub("_"," ", colnames(.))) %>%
                      `rownames<-`(sub("_"," ", rownames(.)))

                    # Make heatmaps:
                    pheatmap::pheatmap(mat2/byclass,
                                       color = grDevices::colorRampPalette(
                                         rev(RColorBrewer::brewer.pal(
                                           n = 7, name ="RdYlBu")))(100),
                                       breaks = seq(0, 1, .01),
                                       cluster_rows = FALSE,
                                       cluster_cols = FALSE,
                                       legend = FALSE, display_numbers = TRUE,
                                       number_format = "%.2f", fontsize = 50,
                                       fontsize_row = 30, fontsize_col = 30,
                                       fontsize_number = 35,
                                       main = title_list[i],
                                       number_color="black")

                    grid::grid.grab()
                  })
  )

  titles <- c("100 timepoints", "50 timepoints", "20 timepoints")

  # Arrange plots:
  p <- mapply(function(gl, tl){

    pl <- (cowplot::plot_grid(plotlist = gl[c(1,4,7,10)], nrow=1) /
             cowplot::plot_grid(plotlist = c(gl[c(2,5,8,11)]), nrow=1) /
             cowplot::plot_grid(plotlist = c(gl[c(3,6,9,12)]), nrow=1)) %>%
      suppressWarnings()

    title <- cowplot::ggdraw() + cowplot::draw_label(tl, size=80)

    # Add titles by length:
    pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }, mat, titles, SIMPLIFY=FALSE
  )

  all <- cowplot::plot_grid(p$l100) /
    cowplot::plot_grid(p$l50) /
    cowplot::plot_grid(p$l20)

  # Save plots:
  if (save){
    pdf(paste0(dirname, "/confmat_noise_", filename, ".pdf"),
        width = 50, height = 100)
    print(all)
    dev.off()
  }

  return(all)
}



#' Combine classification output from different length
#'
#' @param outlist data frame combining classification output info from
#' all timeseries of all lengths
#'
#' @return data frame with wAICc, LOO, and NRMSE in a long format
#'
#' @export

loo_wAICc_nmrse_long <- function(outlist){

  # Make long format for LOO:
  df <- loo_wAICc_nmrse <- outlist %>%
    dplyr::select(id | dplyr::contains("loo_")) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "class_tested",
                        values_to = "loo") %>%
    dplyr::mutate(class_tested = sub("loo_","",class_tested)) %>%

    # Add long format for AICc weight:
    dplyr::left_join(
      outlist %>%
        dplyr::select(id | dplyr::contains("weight_aic_")) %>%
        tidyr::pivot_longer(cols = -id,
                            names_to = "class_tested",
                            values_to = "weight_aic") %>%
        dplyr::mutate(class_tested = sub("weight_aic_","",class_tested)),
      by=c("id","class_tested")
    ) %>%

    # Add long format for NRMSE:
    dplyr::left_join(
      outlist %>%
        dplyr::select(id | dplyr::contains("nrmse_")) %>%
        tidyr::pivot_longer(cols = -id,
                            names_to = "class_tested",
                            values_to = "nrmse") %>%
        dplyr::mutate(class_tested = sub("nrmse_","",class_tested)),
      by=c("id","class_tested")
    ) %>%

    dplyr::mutate(class_tested = factor(class_tested,
                                        levels = c("no_change", "linear",
                                                   "quadratic", "abrupt"))) %>%

    # Add whether the classification fitted expected class:
    dplyr::left_join(outlist %>%
                       dplyr::select(dplyr::contains(c("id","class")) |
                                     noise_lvl | length) %>%
                       dplyr::mutate(correct_classif =
                                       ifelse(class == expected_class,
                                              TRUE, FALSE)),
                     by = "id") %>%

    dplyr::mutate(length = factor(length, levels = c("l20","l25","l33","l50","l100")))

  return(df)
}



#' Extract a stock RAMLBD timeseries
#' (the RAMLBD should be loaded in the environment)
#'
#' @param id short stock ID character
#' @param ts_type type of timeseries (TBbest, ERbest...) character
#'
#' @return data frame of the selected timeseries
#'
#' @export

extract_RAM <- function(id, ts_type){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(year, tidyselect::all_of(ts_type)) %>%
    dplyr::mutate(scen = paste(ts_type, id, sep="_")) %>%
    dplyr::relocate(scen) %>%
    na.omit()

  return(ts)

}
