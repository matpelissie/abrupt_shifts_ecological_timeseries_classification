###-###-###-###-###-###-###-###-###
#
# 17/03/2022 mathieu.pelissie@ens-lyon.fr
#
# Library of time series
#
# 03_ricker_library.R
#
###-###-###-###-###-###-###-###-###


# ROC curve ---------------------------------------------------------------

outlist_aic_l100 <- readRDS("analyses/classif/library/outlist_l100_aic.rds")
outlist_aic_l50 <- readRDS("analyses/classif/library/outlist_l50_1_aic.rds")
outlist_aic_l33 <- readRDS("analyses/classif/library/outlist_l33_1_aic.rds")
outlist_aic_l25 <- readRDS("analyses/classif/library/outlist_l25_1_aic.rds")
outlist_aic_l20 <- readRDS("analyses/classif/library/outlist_l20_1_aic.rds")

max_asd <- max_asdetect(list("l20"=outlist_aicasd_l20, "l25"=outlist_aicasd_l25,
                             "l33"=outlist_aicasd_l33,"l50"=outlist_aicasd_l50,
                             "l100"=outlist_aicasd_l100))

roc_df <- max_asd %>%
  dplyr::mutate(length = factor(length, levels = c("l20","l25","l33","l50","l100")),
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

outlist_aic <- combine_length(list(outlist_aic_l20,
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


roc_aic_only <- combine_length(list(outlist_aic_l20,
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
  dplyr::mutate(length = factor(length, levels=c("l20","l25","l33","l50","l100"))) %>%
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
