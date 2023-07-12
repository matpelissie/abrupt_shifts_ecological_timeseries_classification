###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# 11/04/2023 mathieu.pelissie@ens-lyon.fr
#
# Functions to analyse outputs
#
# functions_output.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###


#' Combine classification output from different length
#'
#' @param outlist_all a list of "outlists" corresponding to different time series lengths
#'
#' @return a data frame combining classification output info from all time series of all lengths
#' @export

combine_length <- function(outlist_all){

  outlist <- lapply(outlist_all,
                    function(z) mapply(function(x,y) x$best_traj %>%
                                         dplyr::mutate(length = stringr::str_extract(y, "[^_]+"),
                                                       noise_lvl = sub(".*?\\_", "", y)),
                                       z$traj_list, names(z$traj_list), SIMPLIFY=FALSE
                    ) %>%
                      # do.call(what="rbind") %>%
                      dplyr::bind_rows() %>%
                      tibble::remove_rownames()
  ) %>%
    # do.call(what="rbind") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(id = paste(simu_id, noise_lvl, sep="_"))

  return(outlist)
}


#' Get the maximum absolute value of the detection line for each time series
#'
#' @param outlist_all a list of "outlists" corresponding to different time series lengths
#'
#' @return a data frame combining the maximum absolute values of the detection line of all lengths
#' @export

max_asdetect <- function(outlist_all){

  df <- lapply(outlist_all,
               function(z) mapply(
                 function(x,y) sapply(x$res_detail$res_abt$shifts_res,
                                      function(t) t$asd_detect %>% abs() %>% max(), simplify=TRUE) %>%
                   data.frame() %>%
                   tibble::rownames_to_column(var="simu_id") %>%
                   dplyr::as_tibble() %>%
                   # dplyr::rename(max_asd=".") %>%
                   dplyr::mutate(length = stringr::str_extract(y, "[^_]+"),
                                 noise_lvl = sub(".*?\\_", "", y)),
                 z$traj_list, names(z$traj_list), SIMPLIFY=FALSE
               ) %>%
                 # do.call(what="rbind") %>%
                 dplyr::bind_rows() %>%
                 tibble::remove_rownames()
  ) %>%
    # do.call(what="rbind") %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(id = paste(simu_id, noise_lvl, sep="_")) %>%
    dplyr::rename(max_asd=".")

  return(df)
}

#' Make confusion matrices by noise combination
#'
#' @param list a "outlist" from the traj_class function corresponding to a given time series lengths or combine_length function
#' @param name a character to name the file saved
#'
#' @return No return value
#' @export

conf_mat <- function(combine_length_out, name, dirname, save){

  mat <- caret::confusionMatrix(data=combine_length_out$class,
                                reference=combine_length_out$expected_class) %>%
    `[`("table") %>%
    `[[`(1) %>%
    `colnames<-`(sub("_"," ", colnames(.))) %>%
    `rownames<-`(sub("_"," ", rownames(.)))

  mat <- mat/12000

  plot <- pheatmap::pheatmap(mat,
                             color = grDevices::colorRampPalette(
                               rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
                             breaks = seq(0, 1, .01),
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             legend = FALSE, display_numbers = TRUE,
                             number_format = "%.2f", fontsize = 30,
                             number_color="black", fontsize_number = 30)

  if (save){
    # pdf(paste0(dirname, "/confmat_",name,".pdf"), width = 10, height = 10)
    png(paste0(dirname,"/confmat_",name,".png"), width = 800, height = 800)
    print(plot)
    dev.off()
    }

  return(plot)
}




#' Make confusion matrices by noise combination
#'
#' @param outlist a "outlist" from the traj_class function corresponding to a given time series lengths
#' @param name a character to name the file saved
#' @param save boolean, should the figure be saved automatically?
#'
#' @return No return value
#' @export

conf_mat_noise <- function(outlist, name, save=TRUE){

  mat <- lapply(outlist,
                function(x)
                  lapply(1:length(x$mat_list), function(i){

                    byclass <- sum(x$mat_list[[i]]$conf_mat_class$table[,1])
                    title_list <- c("sr=.001 se=.001","sr=.001 se=.025","","",
                                    "sr=.025 se=.025","","","sr=.050 se=.025","","",
                                    "sr=.075 se=.025","","")

                    mat2 <- x$mat_list[[i]]$conf_mat_class$table %>%
                      `colnames<-`(sub("_"," ", colnames(.))) %>%
                      `rownames<-`(sub("_"," ", rownames(.)))

                    pheatmap::pheatmap(mat2/byclass,
                                       color = grDevices::colorRampPalette(
                                         rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
                                       breaks = seq(0, 1, .01),
                                       cluster_rows = FALSE, cluster_cols = FALSE,
                                       legend = FALSE, display_numbers = TRUE,
                                       number_format = "%.2f", fontsize = 30,
                                       main = title_list[i],
                                       number_color="black", fontsize_number = 35)

                    grid::grid.grab()
                  })
  )

  # p <- mat %>%
  #   lapply(function(gl){
  #
  #     pl <- (cowplot::plot_grid(plotlist = gl[c(1,2,5,8,11)], nrow=1) /
  #       cowplot::plot_grid(plotlist = c("",gl[c(3,6,9,12)]), nrow=1) /
  #       cowplot::plot_grid(plotlist = c("",gl[c(4,7,10,13)]), nrow=1)) %>%
  #       suppressWarnings()
  #
  #     title <- cowplot::ggdraw() + cowplot::draw_label("length 20")
  #
  #     pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  #     }
  #   )

  titles <- c("length 20", "length 50", "length 100")
  # titles <- c("length 10", "length 25", "length 33")

  p <- mapply(function(gl, tl){

    pl <- (cowplot::plot_grid(plotlist = gl[c(1,2,5,8,11)], nrow=1) /
             cowplot::plot_grid(plotlist = c("",gl[c(3,6,9,12)]), nrow=1) /
             cowplot::plot_grid(plotlist = c("",gl[c(4,7,10,13)]), nrow=1)) %>%
      suppressWarnings()

    title <- cowplot::ggdraw() + cowplot::draw_label(tl)

    pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }, mat, titles, SIMPLIFY=FALSE
  )


  all <- cowplot::plot_grid(p$l20) /
    cowplot::plot_grid(p$l50) /
    cowplot::plot_grid(p$l100)


  # all <- cowplot::plot_grid(p$l10) /
  #   cowplot::plot_grid(p$l25) /
  #   cowplot::plot_grid(p$l33)

  if (save){
    # png(paste0("res/03_ricker_library/simu_library/confmat_noise_",name,".png"), width = 6000, height = 10000)
    pdf(paste0("res/03_ricker_library/simu_library/confmat_noise_",name,".pdf"), width = 50, height = 100)
    print(all)
    dev.off()
    }

  return(all)

}


conf_mat_noise_final <- function(outlist, name, save=TRUE){

  mat <- lapply(outlist,
                function(x)
                  lapply(1:length(x$mat_list), function(i){

                    byclass <- sum(x$mat_list[[i]]$conf_mat_class$table[,1])
                    title_list <- c("sr = 0.001","","",
                                    "sr = 0.025","","",
                                    "sr = 0.050","","",
                                    "sr = 0.075","","")

                    mat2 <- x$mat_list[[i]]$conf_mat_class$table %>%
                      `colnames<-`(sub("_"," ", colnames(.))) %>%
                      `rownames<-`(sub("_"," ", rownames(.)))

                    pheatmap::pheatmap(mat2/byclass,
                                       color = grDevices::colorRampPalette(
                                         rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100),
                                       breaks = seq(0, 1, .01),
                                       cluster_rows = FALSE, cluster_cols = FALSE,
                                       legend = FALSE, display_numbers = TRUE,
                                       number_format = "%.2f", fontsize = 50,
                                       fontsize_row = 30, fontsize_col = 30,
                                       main = title_list[i],
                                       number_color="black", fontsize_number = 35)

                    grid::grid.grab()
                  })
  )

  # p <- mat %>%
  #   lapply(function(gl){
  #
  #     pl <- (cowplot::plot_grid(plotlist = gl[c(1,2,5,8,11)], nrow=1) /
  #       cowplot::plot_grid(plotlist = c("",gl[c(3,6,9,12)]), nrow=1) /
  #       cowplot::plot_grid(plotlist = c("",gl[c(4,7,10,13)]), nrow=1)) %>%
  #       suppressWarnings()
  #
  #     title <- cowplot::ggdraw() + cowplot::draw_label("length 20")
  #
  #     pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  #     }
  #   )

  titles <- c("100 time points", "50 time points", "20 time points")
  # titles <- c("length 10", "length 25", "length 33")

  p <- mapply(function(gl, tl){

    pl <- (cowplot::plot_grid(plotlist = gl[c(1,4,7,10)], nrow=1) /
             cowplot::plot_grid(plotlist = c(gl[c(2,5,8,11)]), nrow=1) /
             cowplot::plot_grid(plotlist = c(gl[c(3,6,9,12)]), nrow=1)) %>%
      suppressWarnings()

    title <- cowplot::ggdraw() + cowplot::draw_label(tl, size=80)

    pl2 <- cowplot::plot_grid(title, pl, ncol=1, rel_heights=c(0.1, 1))
  }, mat, titles, SIMPLIFY=FALSE
  )


  all <- cowplot::plot_grid(p$l100) /
    cowplot::plot_grid(p$l50) /
    cowplot::plot_grid(p$l20)


  # all <- cowplot::plot_grid(p$l10) /
  #   cowplot::plot_grid(p$l25) /
  #   cowplot::plot_grid(p$l33)

  if (save){
    # png(paste0("res/03_ricker_library/simu_library/confmat_noise_",name,".png"), width = 6000, height = 10000)
    pdf(paste0("res/03_ricker_library/simu_library/confmat_noise_",name,".pdf"), width = 50, height = 100)
    print(all)
    dev.off()
  }

  return(all)

}




#' Keep only time series with a given class as best AICc
#'
#' @param outlist_aic_all a list of "outlists" classified with AICc only to different time series lengths
#' @param outlist_aicasd_all a list of "outlists" classified with AICc + as.detect to different time series lengths
#' @param class_kept a character for the type of class kept
#'
#' @return A list of "outlists" with class kept only corresponding to different time series lengths
#' @export

keep_aic <- function(outlist_aic_all, outlist_aicasd_all, class_kept="abrupt"){

  kept_res <- mapply(function(z,t){

    id <- mapply(function(x,y){

      x$best_traj %>%
        dplyr::select(simu_id, class) %>%
        dplyr::rename(class_aic = class) %>%
        dplyr::left_join(y$best_traj, by="simu_id") %>%
        dplyr::filter(class_aic==class_kept) %>%
        dplyr::pull(simu_id)
    }
    ,
    z$traj_list, t$traj_list, SIMPLIFY = FALSE)

  }
  ,
  outlist_aic_all, outlist_aicasd_all, SIMPLIFY = FALSE
  )

  temp <- outlist_aicasd_all
  for (i in 1:length(temp)){
    for (j in 1:length(temp[[i]]$traj_list)){

      temp[[i]]$traj_list[[j]]$res_detail$res_abt$shifts_res <-
        names(temp[[i]]$traj_list[[j]]$res_detail$res_abt$shifts_res) %in% kept_res[[i]][[j]] %>%
        keep(temp[[i]]$traj_list[[j]]$res_detail$res_abt$shifts_res, .)
    }
  }

  return(temp)
}



# Compute time series parameters
summary_param <- function(df_list, ts_type="TB"){

  summ <- lapply(df_list,
                 function(x){

                   ar <- stats::acf(x[[ts_type]], lag=2, plot=FALSE)

                   data.frame(ts = paste0(x$scen[1],"_iter",x$iter[1]),
                              cov = sd(x[[ts_type]])/mean(x[[ts_type]])*100,
                              skw = moments::skewness(x[[ts_type]]),
                              kts = moments::kurtosis(x[[ts_type]]),
                              ar1 = ar$acf[2],
                              ar2 = ar$acf[3],
                              d = D(x[[ts_type]])
                              )
                 }) %>%
    do.call('bind_rows', args=.)

  return(summ)
}

summary_param_bis <- function(df_list, ts_type="TB"){

  summ <- mapply(function(x, name){

                   ar <- stats::acf(x[[ts_type]], lag=2, plot=FALSE)

                   data.frame(ts = name,
                              par_cov = sd(x[[ts_type]])/mean(x[[ts_type]])*100,
                              par_skw = moments::skewness(x[[ts_type]]),
                              par_kts = moments::kurtosis(x[[ts_type]]),
                              par_ar1 = ar$acf[2],
                              par_ar2 = ar$acf[3],
                              par_d = D(x[[ts_type]])
                   )
                 }, df_list, names(df_list), SIMPLIFY=FALSE) %>%
    do.call('bind_rows', args=.)

  return(summ)
}


#' Combine classification output from different length
#'
#' @param outlist a data frame combining classification output info from all time series of all lengths
#'
#' @return a data frame with wAICc, LOO, and NRMSE in a long format
#' @export

loo_wAICc_nmrse_long <- function(outlist){

  loo_wAICc_nmrse <- outlist %>%
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

    left_join(
      outlist %>%
        dplyr::select(id | dplyr::contains("nrmse_")) %>%
        tidyr::pivot_longer(cols = -id,
                            names_to = "class_tested",
                            values_to = "nrmse") %>%
        dplyr::mutate(class_tested = sub("nrmse_","",class_tested)),
      by=c("id","class_tested")
    ) %>%
    dplyr::mutate(class_tested = factor(class_tested,
                                        levels = c("no_change","linear","quadratic","abrupt"))) %>%

    dplyr::left_join(outlist %>%
                       dplyr::select(dplyr::contains(c("id","class")) | noise_lvl | length) %>%
                       dplyr::mutate(correct_classif = ifelse(class==expected_class, TRUE, FALSE)),
                     by = "id") %>%
    dplyr::mutate(length = factor(length, levels = c("l20","l25","l33","l50","l100")))
}

