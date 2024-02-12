###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# Functions for trajectory classification
#
# functions_trajclass.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###

# I) Data preparation --------------------------------------------------------

## I-A) Simulated data -----------------------------------------------------

#' Get csv data
#'
#' @param file Relative path and file name with noise scenarios.
#'
#' @return Tibble with noise or simulation paramaters for scenarios in row.
#'
#' @export

get_data <- function(file){
  readr::read_csv(file)
}



#' Make mortality scenario
#'
#' @param shape Shape of the timeseries to generate, either "cst" (constant),
#' "lin" (linear), "qdr" (quadratic monotonous), "ccv" (concave),
#' "cvx" (convex), or "abt" (abrupt).
#' @param l Length of the timeseries.
#' @param trend (lin, qdr, abt) Character "pos" (positive) or "neg" (negative).
#' @param min (lin, qdr, abt) Minimal value.
#' @param max (lin, qdr, abt) Maximal value.
#' @param velocity (qdr) "acc" (accelerating) or "dec" (decelerating).
#' @param breaktime (abt) Timestep at which the break occurred.
#'
#' @return Vector with the desired features and string with the parameters.
#'
#' @export

make_scen <- function(shape, l=100, trend=NULL,
                      min=NULL, max=NULL, velocity=NULL,
                      breaktime=floor(l/2)){

  if(shape=="cst") {

    scen <- rep(min, l)
    name <- sprintf("l%g_%s_X%g", l, shape, min)

  } else if(shape=="lin") { # Shape linear

    scen <- seq(min, max, length.out = l)
    name <- sprintf("l%g_%s_%s_X%g-%g", l, shape, trend, min, max)
    if(trend=="neg"){
      scen <- rev(scen)
      name <- sprintf("l%g_%s_%s_X%g-%g", l, shape, trend, max, min)
    }

  } else if(shape=="ccv") { # Shape concave

    scen <- c(seq(max, min, length.out = breaktime),
              seq(min, max, length.out = l-breaktime))
    name <- sprintf("l%g_%s_X%g-%g-%g", l, shape, max, min, max)

  } else if(shape=="cvx") { # Shape convex

    scen <- c(seq(min, max, length.out = breaktime),
              seq(max, min, length.out = l-breaktime))
    name <- sprintf("l%g_%s_X%g-%g-%g", l, shape, min, max, min)

  } else if(shape=="qdr") { # Shape quadratic

    if(velocity=="dec"){

      if(trend=="pos"){
        a <- (min-max)/l^2
        b <- 2*(max-min)/l
        c <- min
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (max-min)/l^2
        b <- 2*(min-max)/l
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, max, min)
      }
    }

    if(velocity=="acc"){

      if(trend=="pos"){
        a <- (max-min)/l^2
        b <- 0
        c <- min
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, min, max)
      }

      if(trend=="neg"){
        a <- (min-max)/l^2
        b <- 0
        c <- max
        name <- sprintf("l%g_%s_%s_%s_X%g-%g",
                        l, shape, trend, velocity, max, min)
      }

    }

    pol2deg <- function(x) a*x^2+b*x+c
    scen <- pol2deg(1:l)


  } else if(shape=="abt") { # Shape abrupt

    if(breaktime>l) stop("Break point should occur within the timeseries")

    if(trend=="pos"){
      scen <- c(rep(min, breaktime-1), rep(max, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g",
                      l, shape, trend, breaktime, min, max)
    }

    if(trend=="neg"){
      scen <- c(rep(max, breaktime-1), rep(min, l-breaktime+1))
      name <- sprintf("l%g_%s_%s_brk%g_X%g-%g",
                      l, shape, trend, breaktime, max, min)
    }
  }

  return(list("scen"=scen, "name"=name))

}



#' Make mortality scenario with multiple breakpoints
#'
#' @param l Length of the timeseries.
#' @param trend Character "ccv" (concave), "cvx" (convex),
#' "pos" (positive), or "neg" (negative).
#' @param val* Value of mortality for the first (val1), middle (val2),
#' and last (val3) steps.
#' @param brk* Timestep at which the first (brk1) and second (brk2) breaks
#' occurred.
#'
#' @return Vector with the desired features and string with the parameters.
#'
#' @export

make_scen_multibrk <- function(l=100, trend=NULL,
                               val1=NULL, val2=NULL, val3=NULL,
                               brk1=NULL, brk2=NULL){

  if(brk1>l | brk2>l) stop("Breakpoints should occur within the timeseries")

  shape <- "abt"
  scen <- c(rep(val1, brk1), rep(val2, brk2-brk1), rep(val3, l-brk2))
  name <- sprintf("l%g_%s_%s_brk%g-%g_X%g-%g-%g",
                  l, shape, trend, brk1, brk2, val1, val2, val3)

  return(list("scen"=scen, "name"=name))
}



#' Make simulations and plots
#'
#' @param name Character for scenario name.
#' @param Ts Length of timeseries.
#' @param r Vector of growth rates (same length as timeseries).
#' @param F Vector of mortality rates (same length as timeseries).
#' @param K Carrying capacity.
#' @param P Exponent for sigmoid functional response.
#' @param H Half-saturation rate.
#' @param sr Demographic stochasticity.
#' @param se Environmental stochasticity.
#' @param su Observational noise (i.e. measurement error), not used.
#' @param iter Number of repetitions of timeseries.
#' @param thr Unfeasibility threshold.
#' @param init Percentage of carrying capacity as initial size.
#' @param trs Length of the transient period to be removed (default = 50 steps).
#' @param expected_class Expected trajectory class.
#' @param jfr Frequency of jump events.
#' @param jsz Base size of the jump event (buffered by demographic
#' stochasticity and population size).
#'
#' @return List of two elements:
#' - list of matrix with time varying inputs and outputs for all iterations
#' - data frame with time varying inputs and outputs for all iterations
#'
#' @export

run_simu <- function(name, Ts, r, F, K, P, H, sr, se, su,
                     iter, thr, init, trs=50, expected_class, jfr, jsz){

  # if and only if F is variable:
  full_name <- sprintf("%s_r%#.2g_H%g", sub("X", "F", name), r, H)

  if (length(r)==1) r <- rep(r,Ts)
  if (length(H)==1) H <- rep(H,Ts)

  if (length(r)!=length(F) | length(F)!=length(H) | length(F)!= Ts){
    stop("Parameters should be vectors the same length as timeseries length Ts")
  }

  TBts <- matrix(NA, nrow=Ts, ncol=iter)
  TCts <- matrix(NA, nrow=Ts, ncol=iter)
  rts <- matrix(NA, nrow=Ts, ncol=iter)
  Fts <- matrix(NA, nrow=Ts, ncol=iter)
  Hts <- matrix(NA, nrow=Ts, ncol=iter)
  srts <- matrix(NA, nrow=Ts, ncol=iter)
  sets <- matrix(NA, nrow=Ts, ncol=iter)
  suts <- matrix(NA, nrow=Ts, ncol=iter)

  b <- r/K

  i <- 1
  extra <- 0
  while (i <= iter){

    # Begin at a given % of carrying capacity:
    x <- K*init

    # Let the equilibrium be reached:
    for (t in 1:trs){

      N <- x*exp(r[1]-b[1]*x)
      TC <- F[1]*x^P/(H[1]^P+x^P)
      x <- N - TC

      if (x<thr){
        if (x<(-10^1)) x <- -10^1 # to avoid reaching -Inf
        print("Warning: the systems fails to stabilize, check the parameters")
      }

    }

    # First timepoint with stochasticity:
    TBts[1, i] <- x*exp(su*rnorm(1))

    for (t in 2:Ts){

      r1 <- r[t]*exp(rnorm(n=1, sd=sr)) # demographic stochasticity
      N <- x*exp(r1-b[t]*x+rnorm(n=1, sd=se))
      TC <- F[t]*x^P/(H[t]^P+x^P)

      if (runif(1) < jfr){ # In case of jump is present

        jump_ok <- FALSE
        jump_count <- 0

        while(!jump_ok & jump_count<50){ # Avoid jumping to extinction

          J <- rnorm(n=1, sd=se)*jsz*x
          if (N - TC + J > 0) jump_ok <- TRUE
          if (N - TC + J < 0) jump_count <- jump_count + 1
        }

      } else {
        J <- 0
      }

      x <- N - TC + J


      if (x<thr){

        extra <- extra+1
        break

      } else {

        TBts[t, i] <- x*exp(rnorm(n=1, sd=su)) # observation error (not used)
        TCts[t, i] <- TC
        rts[t, i] <- r1
        Fts[t, i] <- F[t]
        Hts[t, i] <- H[t]
        srts[t, i] <- sr
        sets[t, i] <- se
        suts[t, i] <- su

      }

    }

    # Is constant timeseries not collapsed:
    if (!expected_class %in% c("abrupt") &
        x>0.1*TBts[1, i] & TBts[90, i]>0.1*TBts[1, i]) i <- i+1

    # Is timeseries complete:
    else if (expected_class %in% c("abrupt") & x>thr) i <- i+1

    # To escape if not possible to complete the timeseries:
    if (extra>2*max(iter,50)){
      return(warning("Warning: too many simulations fail to avoid extinction,
                     check parameters and noise"))
    }
  }

  lib_ts <- data.frame(scen = as.vector(full_name),
                       iter = rep(1:iter, each=Ts),
                       year = rep(1:Ts, iter),
                       TB = as.vector(TBts),
                       r = as.vector(rts),
                       F = as.vector(Fts),
                       H = as.vector(Hts),
                       P = rep(P, Ts*iter),
                       TC = as.vector(TCts),
                       sr = as.vector(srts),
                       se = as.vector(sets),
                       su = as.vector(suts)
  )

  scen_mat <- list("TB"=TBts, "TC"=TCts, "r"=rts, "F"=Fts)

  return(list(
    "scen_mat"=scen_mat,
    "df"=lib_ts
  ))

}



#' Make and save simulated data
#'
#' @param param_df Data frame with scenarios parameters.
#' @param scen_fct Character to specify the scenario-making function to use.
#' @param sr Demographic stochasticity.
#' @param se Environmental stochasticity.
#' @param su Observational noise (i.e. measurement error), not used.
#' @param jfr Frequency of jump events.
#' @param jsz Base size of the jump event (buffered by demographic
#' stochasticity and population size).
#'
#' @return Data frame with simulated data for all scenarios.
#'
#' @export

make_store_simu <- function(param_df, scen_fct, sr, se, su, jfr, jsz){

  param_names <- colnames(param_df)
  all_simu <- data.frame(df=matrix(ncol = 17))[-1,]

  for (j in 1:nrow(param_df)){ # For each scenario

    for (i in 1:length(param_names)){ # import parameter values

      assign(param_names[i], param_df %>%
               dplyr::pull(param_names[i]) %>%
               `[`(j))
    }

    if (scen_fct == "make_scen"){
      Fseq <- make_scen(shape=shape, trend=trend, min=min, max=max,
                        velocity=velocity, l=l, breaktime = breaktime)

    } else if (scen_fct == "make_scen_multibrk"){
      Fseq <- make_scen_multibrk(trend=trend, l=l,
                                 val1=val1, val2=val2, val3=val3,
                                 brk1=brk1, brk2=brk2)
    }

    sim <- run_simu(name=Fseq$name, r=r, F=Fseq$scen, Ts=length(Fseq$scen),
                    K=10, P=2, H=H, sr=sr, se=se, su=su, iter=iter,
                    thr=0, init=0.9, expected_class=expected_class,
                    jfr=jfr, jsz=jsz)

    if(class(sim) != "character"){ # if no warnings

      sim$df <- sim$df %>%
        dplyr::mutate(expected_traj=expected_traj,
                      expected_class=expected_class)
      all_simu <- all_simu %>% rbind(sim$df)
    }

  }

  return(all_simu)

}



### Bifurcation plot --------------------------------------------------------

#' Give the state at the equilibrium for a given timeseries
#'
#' @param lib_ts Data frame with inputs and outputs for all iterations.
#' @param thr Unfeasibility threshold (default 0).
#'
#' @return Character indicating state at equilibrium (unfeasible, stable,
#' cycles, chaos).
#'
#' @export

eq_state <- function(lib_ts, thr=0){

  sample_length <- nrow(lib_ts)/10

  # Unique lag-1 difference of the last 10% timepoints:
  unique_diff <- lib_ts %>%
    tail(n=sample_length) %>%
    dplyr::mutate(diff_lag1=TB-dplyr::lag(TB)) %>%
    tidyr::drop_na(diff_lag1) %>%
    dplyr::distinct(diff_lag1) %>%
    dplyr::pull()

  # Determine equilibrium state based on a series of criteria:
  state <- NA

  if (lib_ts %>%
      dplyr::pull(TB) %>%
      tail(1) < thr) { # last value below unfeasibility threshold
    state <- "unfeasible"

  } else if (unique_diff %>% length() == 1 | # fully stabilized
             max(abs(unique_diff)) < 10^(-5) | # if oscillating very little
             unique_diff %>% sign() %>% unique() %>%
             length() == 1) { # if not yet stabilized
    state <- "stable"

  } else if (unique_diff %>% length() < (sample_length-1) | # regular cycles
             sign(abs(unique_diff[-1]) -
                  abs(stats::lag(unique_diff, k=2)[-1])) %>% unique() %>%
             length() == 1) { # if lag-2 cycles of decreasing amplitude
    state <- "cycle"

  } else if (unique_diff %>%
             length() == sample_length-1) { # all lag-1 diff values different
    state <- "chaos"
  }

  return(state)

}



#' Plot stable states and evidence potential alternative stable states
#'
#' @param Ts Length of timeseries.
#' @param rval Value of growth rate.
#' @param P Exponent for sigmoid functional response.
#' @param H Half-saturation rate.
#' @param K Carrying capacity.
#' @param Fmax Maximal mortality rate.
#' @param Fstep Interval between two values of mortality rate.
#' @param TBhigh Starting value for high biomass.
#'
#' @return Table and ggplot with values of stable states along a range of F.
#'
#' @export

ASS_row <- function(Ts, rval, P, H, K, Fmax, Fstep, TBhigh){

  Fseq <- seq(0, Fmax, Fstep)

  row_inc <- row_dec <- states_inc <- states_dec <- stab_inc <- stab_dec <-
    data.frame(matrix(NA, length(rval), length(Fseq))) %>%
    `rownames<-`(paste0("r", rval)) %>%
    `colnames<-`(paste0("F", Fseq))

  # Find equilibrium for increasing F & high starting biomass:
  init <- TBhigh
  for (j in Fseq){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name="test", Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1,
                     thr=0, init,expected_class="abrupt", jfr=0, jsz=0)[[2]] %>%
      try() # expected_class abrupt to avoid removing collapsed timeseries

    if(class(simu) !="try-error"){

      states_inc[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)

      # take mean of last 10% timepoints:
      init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K
      init <- max(init, 10^(-5)) # to avoid reaching zero

      row_inc[paste0("r", rval), paste0("F", j)] <- init*K

    } else {
      states_inc[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_inc[paste0("r", rval), paste0("F", j)] <- NA
    }

  }

  states_inc_long <- states_inc %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "states_inc") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  row_inc_long <- row_inc %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "value_inc") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  # Find equilibrium for decreasing F & low starting biomass:
  for (j in rev(Fseq)){

    F <- rep(j, Ts)
    r <- rep(rval, Ts)
    simu <- run_simu(name=NA, Ts, r, F, K, P, H=H, sr=0, se=0, su=0, iter=1,
                     thr=0, init,expected_class="abrupt", jfr=0, jsz=0)[[2]] %>%
      try()

    if(class(simu) !="try-error"){

      states_dec[paste0("r",rval), paste0("F",j)] <- eq_state(simu, thr=0)

      # take mean of last 10% timepoints:
      init <- mean(simu[(nrow(simu)*0.9):nrow(simu), "TB"])/K
      init <- max(init, 10^(-5)) # to avoid reaching zero

      row_dec[paste0("r", rval), paste0("F", j)] <- init*K

    } else {
      states_dec[paste0("r",rval), paste0("F",j)] <- "unfeasible"
      row_dec[paste0("r", rval), paste0("F", j)] <- NA
    }

  }

  states_dec_long <- states_dec %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "states_dec") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  row_dec_long <- row_dec %>%
    tibble::rownames_to_column("r") %>%
    tidyr::pivot_longer(cols = !r,
                        names_to = "F",
                        names_prefix = "F",
                        values_to = "value_dec") %>%
    dplyr::mutate(r = sub("r","",r) %>% as.numeric(),
                  F = as.numeric(F))

  # Compare equilibria:
  bif_wide <- row_inc_long %>%
    dplyr::left_join(row_dec_long, by=c("r", "F")) %>%
    dplyr::left_join(states_inc_long, by=c("r", "F")) %>%
    dplyr::left_join(states_dec_long, by=c("r", "F")) %>%
    # arbitrary threshold of 10^(-1) for equal equilibrium value:
    dplyr::mutate(ASS = ifelse(abs(value_inc - value_dec)>10^(-1),
                               TRUE, FALSE)) %>%
    dplyr::mutate(diffstate = ifelse(states_inc != states_dec,
                                     TRUE, FALSE))


  bif_long <- bif_wide %>%
    tidyr::pivot_longer(cols = contains("_"),
                        names_to = c(".value","sense"),
                        names_pattern = "(\\w+)_(\\w+)")


  bif_plot <- ggplot(bif_long, aes(x=F, y=value, color=as.factor(sense)))+
    geom_line()+
    geom_point(size=0.5)+
    labs(y="TB at equilibrium", color="Sense")+
    theme_light()+
    expand_limits(y = 0)+
    theme(legend.position = c(0.9,0.8))


  stab_plot <- ggplot(bif_long, aes(x=F, y=stab, color=as.factor(sense)))+
    geom_hline(yintercept=0, lty=3)+
    geom_line(lty=1)+
    labs(y="Stability")+
    theme_light()+
    expand_limits(y = 0)+
    theme(legend.position = "none")+
    labs(caption = paste0("r = ", rval, " ; H = ", H, " ; P = ", P))

  bif_stab <- bif_plot / stab_plot +
    patchwork::plot_layout(heights = unit(c(4, 1), c('null')))

  return(list("table" = bif_wide, "plot" = bif_stab))

}



#' Plot system stable states under a range of constant r and F
#'
#' @param Ts Length of timeseries.
#' @param K Carrying capacity.
#' @param rmax Maximal growth rate value (default 3).
#' @param rstep Interval between two values of growth rate (default 0.1).
#' @param Fmax Maximal mortality rate value (default 3).
#' @param Fstep Interval between two values of mortality rate (default 0.1).
#' @param P Exponent for sigmoid functional response.
#' @param H Half-saturation rate.
#' @param thr Extinction threshold.
#' @param init Percentage of carrying capacity as initial size.
#'
#' @return List of plots with stable states and along ranges of r and F and
#' the data frame used for the plots.
#'
#' @export

bifurc_analysis <- function(Ts, K, rmax=3, rstep=0.1,
                            Fmax=3, Fstep=0.1, P, H, thr, init){

  rseq <- seq(0, rmax, rstep)
  Fseq <- seq(0, Fmax, Fstep)
  ASS <- data.frame()

  # Find equilibrium states for every growth rate value:
  for (i in rseq){

    ASS <- ASS %>%
      rbind(ASS_row(Ts=Ts, rval = i, Fmax=Fmax, K=K,
                    Fstep=Fstep, H=H, P=P, TBhigh=init)$table)

    print(paste0(i, "/", rmax))
  }

  # Plot types and values of equilibrium state:
  bif_ana <- ASS %>%
    dplyr::mutate(ASS = ifelse(ASS, TRUE, NA))

  bif_ana$states_inc <- factor(bif_ana$states_inc,
                               levels = c("unfeasible","stable","cycle","chaos"))

  bif_ana$states_dec <- factor(bif_ana$states_dec,
                               levels=c("unfeasible","stable","cycle","chaos"))
  state.col <- c(unfeasible="#c44601", stable="#5ba300",
                 cycle="#0073e3", chaos="#e79b01")

  # Types of equilibrium and ASS when increasing F:
  plot_inc <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_inc, alpha=value_inc))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10), na.value=1)+
    geom_point(data=subset(bif_ana, !is.na(ASS) &
                             !states_inc %in% c("unfeasible", "chaos")),
               # display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ",
                          init*100, "% K ; increasing F"))


  # Types of equilibrium and ASS when decreasing F:
  plot_dec <- ggplot(bif_ana, aes(x=F, y=r))+
    geom_tile(aes(fill=states_dec, alpha=value_dec))+
    scale_fill_manual(values = state.col)+
    scale_alpha(range=c(0.3,1), limits=c(0,10))+
    geom_point(data=subset(bif_ana, !is.na(ASS) &
                             !states_dec %in% c("unfeasible", "chaos")),
               # display ASS only for stable equilibrium
               aes(x=F, y=r), size=1)+
    scale_color_manual(values = c("TRUE"="black"))+
    theme_minimal()+
    labs(caption = paste0("h = ", H, " ;  init TB inc = ",
                          init*100, "% K ; decreasing F"))


  return(list("plot_inc" = plot_inc,
              "plot_dec" = plot_dec,
              "bif_ana" = bif_ana))
}



## I-B) Empirical data -----------------------------------------------------



#' Extract stocks RAMLBD timeseries
#' (the RAMLBD should be loaded in the environment)
#'
#' @param id Short stock ID character.
#' @param ts_type Type of timeseries (TBbest, ERbest...) character.
#' @param drop_na Keep only years with all types of timeseries available.
#'
#' @return Data frame of the selected timeseries.
#'
#' @export

extract_RAM <- function(id, ts_type, drop_na=TRUE){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid %in% id) %>%
    dplyr::mutate(scen = paste(ts_type, stockid, sep="_")) %>%
    dplyr::select(scen, year, tidyselect::all_of(ts_type)) %>%
    dplyr::relocate(scen)

  if (drop_na){
    ts <- ts %>% na.omit()
  }

  return(ts)

}


## I-C) Data standardization -----------------------------------------------

#' Process timeseries data prior classification
#'
#' @param df Data frame with timeseries to classify.
#' @param thr Unfeasibility threshold (default 0).
#' @param type Character specifying whether the series come from simulations
#' ("sim") or from empirical data ("RAM" or "data").
#' @param apriori Logical, whether expected trajectories are indicated.
#'
#' @return List of data frames ready for analyses and the type of timeseries.
#'
#' @export

prep_data <- function(df, thr=0, type="sim", apriori){

  if (type=="sim"){ # For simulated timeseries

    # Reshape the data frame:
    if ("iter" %in% names(df)){

      ts_type <- "TB"
      time_type <- names(df)[3]
      dataset <- df %>%
        {if (apriori==TRUE) dplyr::select(., scen, iter, year, TB,
                                          expected_traj,
                                          expected_class) else .} %>%
        {if (apriori==FALSE) dplyr::select(., scen, iter, year, TB) else .} %>%
        dplyr::rename(Y=TB)

      iter <- seq_len(max(df$iter))
      scen_list <- df %>% dplyr::distinct(scen) %>% dplyr::pull()
    }

  } else if (type=="RAM" | type=="data"){ # For empirical data

    df <- df %>% dplyr::rename_with(.cols=1, ~"scen")
    ts_type <- names(df)[3]
    time_type <- names(df)[2]
    iter <- 1
    dataset <- df %>%
      dplyr::rename(Y=dplyr::all_of(ts_type),
                    year=dplyr::all_of(time_type)) %>%
      dplyr::mutate(iter = 1)

    scen_list <- df %>% dplyr::distinct(scen) %>% dplyr::pull()

  }

  # Break down the data frame into a list of timeseries:
  sets <- list()
  for (j in 1:length(scen_list)){

    by_scen <- dataset %>%
      dplyr::filter(scen==scen_list[j]) %>%
      split(.$iter)

    dats <- lapply(by_scen, function(dat){

      X <- dat$year
      Y <- dat$Y

      if (apriori==TRUE){

        expected_traj <- dat$expected_traj
        expected_class <- dat$expected_class

        set <- data.frame(X=X, Y=Y,
                          expected_traj=expected_traj,
                          expected_class=expected_class)

      } else if (apriori==FALSE){

        set <- data.frame(X=X, Y=Y)
      }

    })

    sets <- c(sets, dats)
  }

  # Rename each timeseries:
  if (type=="sim"){
    names(sets) <- paste0(rep(scen_list, each=max(df$iter)),"_iter",
                          sprintf("%02d", 1:max(df$iter)))
  } else {
    names(sets) <- scen_list
  }

  if (!is.null(thr)){
    sets <- sets %>%
      lapply(function(x) x %>% dplyr::filter(Y>thr))
  }

  return(list("ts"=sets, "ts_type"=ts_type, "time_type"=time_type))

}



# III) Trajectory classification -----------------------------------------------

## III-A) Smooth trajectory models ------------------------------------------------

#' Trajectory classification (from Rigal et al. 2020)
#'
#' @param Y Vector of the timeseries values.
#' @param X Vector of the timesteps.
#' @param dataset Two-column dataframe with timesteps and timeseries values.
#' @param interval_size Set to 0.5 for timesteps regularly spaced by 1 unit.
#'
#' @return Three-row data frame with infos about the trajectory
#' (no change, linear, and polynomial).
#'
#' @export

class_trajectory_mod <- function (Y = NULL, X = NULL, dataset = NULL,
                                   interval_size = 0.5){

  ## Initiation and data check
  if (is.null(X) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }

  if (is.null(X) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,2]
    X <- dataset[,1]

  } else {
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      } else {
        Y <- dataset[, Y]
        X <- dataset[, X]
      }

    } else {
      if (!(class(Y) %in% c("numeric","integer")) == TRUE &
          !(class(X) %in% c("numeric","integer")) == TRUE) {
        stop("'Y' and 'X' must be either characters or vector but 'class'
             must be similar")}
    }
  }

  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),] # ordering the X values

  if (length(X)<4){
    stop("timeseries length must be at least 4")
  }

  Y <- data$Y
  X <- data$X


  ## Trajectory fitting
  # No change model:
  null.model <- lm(Y~1)
  nrmse_nch <- sqrt(sum(summary(null.model)$residuals^2)/length(Y))/sd(Y)
  aic_nch <- MuMIn::AICc(null.model)

  # Linear model:
  linear.model <- lm(Y~X)
  nrmse_lin <- sqrt(sum(summary(linear.model)$residuals^2)/length(Y))/sd(Y)
  aic_lin <- MuMIn::AICc(linear.model)

  # Quadratic model:
  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))
  nrmse_pol <- sqrt(sum(summary(orthogonal_polynomial)$residuals^2)/
                      length(Y))/sd(Y)
  aic_pol <- MuMIn::AICc(orthogonal_polynomial)

  ## Get relevant values for quadratic output
  # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values
  # in the X interval for first_order_coefficient, second_order_coefficient
  # and intercept, knowing that X'= alpha*X + beta and chi = eta*X'^2 + theta

  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]

  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1] ~ X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1] ~ X)$coef[1]

  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2 ~
                  orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2 ~
                    orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta

  Y2 <- Y*(max(X)-min(X))/(max(Y)-min(Y))
  # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3
  polynomial_orthonormal_basis <- lm(Y2~poly(X,2, raw=T))$coefficients

  # Quadratic model output:
  classification <-
    data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
               first_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[2, 4],
               second_order_coefficient = (alpha^2)*gammab*eta,
               second_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[3, 4],
               strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
               intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
               x_m = (X[length(X)]-X[1])/2+X[1],
               # points of interest:
               p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),
               p2 = (-polynomial_orthonormal_basis[2]+1)/
                 (2*polynomial_orthonormal_basis[3]),
               p3 = (-polynomial_orthonormal_basis[2]-1)/
                 (2*polynomial_orthonormal_basis[3]),
               aic = aic_pol,
               nrmse = nrmse_pol)

  # Linear model output:
  classification[2,] <-
    data.frame(first_order_coefficient = delta*alpha,
               first_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[3, 4],
               strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
               intercept = epsilon+delta*beta,
               x_m = (X[length(X)]-X[1])/2+X[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_lin,
               nrmse = nrmse_lin)

  # No change model output:
  classification[3,] <-
    data.frame(first_order_coefficient = 0,
               first_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[2, 4],
               second_order_coefficient = 0,
               second_order_pvalue =
                 summary(orthogonal_polynomial)$coefficients[3, 4],
               strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
               intercept = null.model$coefficients,
               x_m = (X[length(X)]-X[1])/2+X[1],
               p1 = NA,
               p2 = NA,
               p3 = NA,
               aic = aic_nch,
               nrmse = nrmse_nch)


  # Significance of the coefficients:
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){
    classification$best_model <- c("best", NA, NA) # quadratic model
  } else if(summary(orthogonal_polynomial)$coefficients[2, 4]<=0.05){
    classification$best_model <- c(NA, "best", NA) # linear model
  } else {
    classification$best_model <- c(NA, NA, "best") # no change model
  }

  # retrieve the adjusted coefficient of determination
  # classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared

  # Give classification for each model fitted:
  for (i in 1:3){

    # Compute the derivative at xm-delta and at xm + delta with delta being
    # half of the input interval size
    derivative <-
      2*(classification$x_m[i] - (X[length(X)]-X[1])*(interval_size/2))*
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]
    derivative2 <-
      2*(classification$x_m[i] + (X[length(X)]-X[1])*(interval_size/2))*
      classification$second_order_coefficient[i] +
      classification$first_order_coefficient[i]


    if(sign(derivative) != sign(derivative2) | i==3){
      # non consistent direction around x_m
      classification$derivative[i]  <-  NA
      classification$intercept_derivative[i]  <-  NA

    } else {
      # consistent direction around x_m
      classification$derivative[i] <- mean(c(derivative, derivative2))
      classification$intercept_derivative[i] <-
        (classification$second_order_coefficient[i]*classification$x_m[i]^2+
           classification$first_order_coefficient[i]*classification$x_m[i]+
           classification$intercept[i]) -
        classification$x_m[i]*classification$derivative[i]
    }

    # Compute the derivative of the curvature function:
    classification$derivated_curvature[i] <-
      -12*(classification$second_order_coefficient[i]^2)*
      (2*classification$second_order_coefficient[i]*classification$x_m[i]+
         classification$first_order_coefficient[i])*
      (classification$second_order_coefficient[i]/
         abs(classification$second_order_coefficient[i]))/
      ((1+(2*classification$second_order_coefficient[i]*classification$x_m[i]+
             classification$first_order_coefficient[i])^2)^(2.5))

    # Keep derivated curvature even if not significant for polynomial fit:
    if(classification$second_order_pvalue[i]>0.05 & i != 1){
      classification$derivated_curvature[i] <- NA
    }

    # Classify the direction:
    classification$direction[i] <- NA
    classification$direction[i][which(
      classification$derivative[i] > 0)] <- "increase"
    classification$direction[i][which(
      classification$derivative[i] < 0)] <- "decrease"
    classification$direction[i][which(
      is.na(classification$derivative[i]))] <- "stable"
    classification$direction[i][which(
      as.numeric(classification$first_order_pvalue[i])>0.05 &
        as.numeric(classification$second_order_pvalue[i])>0.05)] <- "stable"

    # Classify the acceleration:
    classification$acceleration[i] <- NA
    classification$acceleration[i][which(
      classification$derivated_curvature[i] < 0)] <- "accelerated"
    classification$acceleration[i][which(
      classification$derivated_curvature[i] > 0)] <- "decelerated"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] < 0)] <- "concave"
    classification$acceleration[i][which(
      classification$direction[i] == "stable" &
        classification$second_order_coefficient[i] > 0)] <- "convex"
    classification$acceleration[i][which(
      is.na(classification$derivated_curvature[i]))] <- "constant"

    # Give the final classification combining direction and acceleration:
    classification$shape_class[i] <- paste(classification$direction[i],
                                           classification$acceleration[i],
                                           sep="_")
  }

  # To remove what corresponds to one useless parameter:
  if(classification$shape_class[2] == "stable_constant"){
    aic_lin <- classification[2,]$aic <- MuMIn::AICc(linear.model)-2
  }

  # Provide the linear approach results for comparison:
  linear.model.summary <- summary(linear.model)

  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]

  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]

  row.names(classification) <- c("Y_pol","Y_lin", "Y_nch")

  return(classification)

}



#' Markov-chain simulations of classification
#'
#' @param dataset Data frame (timeseries) ready for analyses.
#' @param niter Number of MC simulations (option to run more than one disabled).
#' @param ref_year Reference year (default the middle year of the interval).
#' @param correction Logical, to the reference value to 100 and correct values
#' below 0 before logtransformation.
#' @param fit Character to specify the type of fit
#' (either "nch", "lin", or "pol", default is best fit).
#'
#' @return Data frame with as many rows as iterations.
#'
#' @export

mc_trend <- function(dataset,
                     niter,
                     ref_year=NULL,
                     correction=FALSE,
                     fit=NULL){

  # Initiate output data frame:
  b <- data.frame(t(rep(NA, 15)))
  attributes(b)$names <- c("second_order_coefficient",
                           "first_order_coefficient",
                           "strd_error",
                           "shape_class",
                           "intercept",
                           "p_1",
                           "p_2",
                           "p_3",
                           "slope_p_value",
                           "slope",
                           "ref_year",
                           "aic",
                           "best_model",
                           "second_order_pvalue",
                           "nrmse"
  )

  # Define reference year and correct accordingly if requested:
  if(is.null(ref_year)){
    ref_year <- dataset$X[round(nrow(dataset)/2)+1]
  }

  if(correction == TRUE){
    ref_value <- dataset$Y[dataset$X == ref_year]
    if(ref_value == 1){
      dataset$Y <- 100*dataset$Y # set reference year value to 100
      dataset$Y[which(dataset$Y <= 1)] <- 1 # set values < 1 to 1
    }

    if(ref_value!=1 & ref_value!=100){
      if(ref_value>1){
        dataset$Y <- dataset$Y/ref_value
      }

      if(ref_value<1){
        stop("use 'correction = FALSE' when value of the
             reference year is strictly below 1")
      }

      dataset$Y <- 100*dataset$Y
      if(length(which(dataset$Y <= 1))>0){
        print("caution, low values corrected, if strongly decreasing or
              increasing trajectory, use respectively  first or last year
              as referential")}
      dataset$Y[which(dataset$Y <= 1)] <- 1
    }

    if(ref_value == 100){
      dataset$Y[which(dataset$Y <= 1)] <- 1
    }
    dataset$log <- log(dataset$Y) # log transform Y

  }

  if(correction == FALSE){
    if(min(dataset$Y)<0){
      min_value <- abs(min(dataset$Y))
      dataset$Y <- dataset$Y+min_value
    } else {min_value <- 0}
    dataset$log <- dataset$Y
  }

  # Run classification (option to make several normal sampling of Y disabled):
  for(i in 1:niter){

    if (niter==1){ # no resampling if only one iteration
      a <- dataset$Y
    }

    if(correction == TRUE){
      # set reference year value to 100 and retransform values if logtranformed
      a <- exp(a)/exp(a[which(dataset$X == ref_year)])*100
    } else {a <- a-min_value}

    # Run classification:
    a <- class_trajectory_mod(a, dataset$X)

    best_model <- a %>%
      dplyr::filter(best_model=="best") %>%
      rownames() %>%
      sub("Y_","",.)

    # Keep the fitting output chosen:
    if(is.null(fit)){
      a <- a %>% dplyr::filter(best_model=="best")
    }else if(fit == "nch"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_nch")
    }else if(fit == "lin"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_lin")
    }else if(fit == "pol"){
      a <- a %>% dplyr::filter(rownames(.)=="Y_pol")
    }


    # Store classification output from each iteration (only one here):
    b[i, 1] <- a$second_order_coefficient
    b[i, 2] <- a$first_order_coefficient
    b[i, 3] <- a$strd_error
    b[i, 4] <- a$shape_class
    b[i, 5] <- a$intercept

    # Record changing point inside timeseries:
    if(a$second_order_coefficient!=0){
      if(findInterval(a$p1,  c(min(dataset$X), max(dataset$X))) == 1){
        b[i, 6] <- a$p1} else {b[i, 6] <- NA}
      if(findInterval(a$p2,  c(min(dataset$X), max(dataset$X))) == 1){
        b[i, 7] <- a$p2} else {b[i, 7] <- NA}
      if(findInterval(a$p3,  c(min(dataset$X), max(dataset$X))) == 1){
        b[i, 8] <- a$p3} else {b[i, 8] <- NA}

    } else {
      b[i, 6] <- NA
      b[i, 7] <- NA
      b[i, 8] <- NA
    }

    b[i, 9] <- a$linear_slope_pvalue
    b[i, 10] <- a$linear_slope
    b[i, 12] <- a$aic
    b[i, 13] <- best_model
    b[i, 14] <- a$second_order_pvalue
    b[i, 15] <- a$nrmse
  }

  b[, 4] <- as.factor(b[, 4])
  b[, 11] <- rep(ref_year, nrow(b))

  return(b)
}



#' Summary of MC simulations of classification
#'
#' @param sets List of data frame (timeseries) ready for analyses.
#' @param niter Number of MC simulations (option to run more than one disabled).
#' @param ref_year Reference year (default the middle year of the interval).
#' @param correction Logical, to the reference value to 100 and correct values
#' below 0 before logtransformation.
#' @param fit Character to specify the type of fit
#' (either "nch", "lin", or "pol", default is best fit).
#'
#' @return Data frame with as many rows as time series.
#'
#' @export

res_trend <- function(sets,
                      niter,
                      ref_year=NULL,
                      correction=FALSE,
                      fit=NULL){

  res <- data.frame()

  # Run classification for each timeseries in the list:
  for (i in seq_len(length(sets))){

    if(nrow(sets[[i]])>3){ # Errors occur if timeseries length is below 4

      simulated <- mc_trend(sets[[i]], niter, ref_year=NULL,
                            correction=correction, fit) %>%
        dplyr::mutate(best_model=as.factor(best_model))

      # Retrieve best shape and model (only one iteration):
      if(length(levels(simulated$shape_class)) == 1){
        max_shape <- levels(simulated$shape_class)
        best_model <- levels(simulated$best_model)}

      # Summarize classification outputs (useful if multiple iterations):
      alpha2 <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 1]))
      sd_alpha2 <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 1]))
      alpha1 <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 2]))
      sd_alpha1 <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 2]))
      inter <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 5]))
      strd <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 3]))
      p_1 <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 6]), na.rm=T)
      sd_p_1 <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 6]), na.rm=T)
      if(!is.na(p_1) &&
          findInterval(p_1, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){
        p_1 <- sd_p_1 <- as.numeric(NA)}
      p_2 <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 7]), na.rm=T)
      sd_p_2 <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 7]), na.rm=T)
      if(!is.na(p_2) &&
          findInterval(p_2, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){
        p_2 <- sd_p_2 <- as.numeric(NA)}
      p_3 <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 8]), na.rm=T)
      sd_p_3 <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 8]), na.rm=T)
      if(!is.na(p_3) &&
          findInterval(p_3, c(min(sets[[i]]$X), max(sets[[i]]$X))) != 1){
        p_3 <- sd_p_3 <- as.numeric(NA)}
      slope_p_value <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 9]), na.rm=T)
      slope <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 10]), na.rm=T)
      slope_sd <- sd(
        as.numeric(simulated[simulated$shape_class == max_shape, 10]), na.rm=T)
      aic <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 12]))
      second_order_pvalue <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 14]), na.rm=T)
      nrmse <- mean(
        as.numeric(simulated[simulated$shape_class == max_shape, 15]))

    } else { # If the timeseries is below 4 timesteps:
      alpha2 <- alpha1 <- sd_alpha1 <- inter <- strd <- p_1 <- sd_p_1 <- p_2 <-
        sd_p_2 <- p_3 <- sd_p_3 <- max_shape <- slope_p_value <- slope <-
        slope_sd <- aic <- best_model <- second_order_pvalue <- nrmse <- NA}

    ref_year <- simulated[1,11]


    res <- res %>%
      rbind(data.frame(alpha2, alpha1, sd_alpha1, inter, strd, p_1, sd_p_1,
                       p_2, sd_p_2, p_3, sd_p_3, second_order_pvalue,
                       slope_p_value, slope, slope_sd, nrmse, ref_year,
                       max_shape=as.factor(max_shape), aic=aic,
                       best_model=best_model, best_class=NA, trend=NA))

  }

  # Define best class and trend from best shape:
  res <- res %>%
    dplyr::mutate(
      best_class = ifelse(!is.na(max_shape),
                          dplyr::case_when(
                            stringr::str_detect(
                              max_shape,"stable_constant") ~ "no_change",
                            stringr::str_detect(
                              max_shape, "constant") ~ "linear",
                            stringr::str_detect(
                              max_shape, "acc|dec|conc|conv") ~ "quadratic"
                          ), NA),
      trend = ifelse(!is.na(max_shape),
                     sub("_.*", "", max_shape), NA)
    )

  return(res)
}



## III-B) Abrupt models ----------------------------------------------------

#' Breakpoints analysis using 'asdetect' package
#'
#' @param ts Timeseries to analyse as 'ts' object.
#' @param asd_thr Numeric threshold in detection timeseries.
#' @param check_true_shift Logical parameter to perform or not the check that
#' the detected shift(s) do not correspond to a flat section.
#' @param lowwl Lowest window length used in algorithm (default 5).
#' @param highwl Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of timeseries.
#' @param mad_thr Threshold of anomalous change in number
#' of median absolute deviations (default 3).
#' @param mad_cst Correction factor for asymptotic normal consistency.
#'
#' @return List of two objects:
#' - one-row data frame with info about potentially detected breakpoints
#' - vector of detection values from as_detect method
#' - list of runs corresponding to local extrema in detection timeseries
#'
#' @references Boulton C. A. & Lenton T. M. (2019) 'A new method for detecting
#' abrupt shifts in time series'.
#'
#' @export

asd_fct <- function(ts, asd_thr, check_true_shift,
                    lowwl, highwl, mad_thr, mad_cst){

  time <- ts %>%
    dplyr::pull(1) %>%
    as.numeric()

  vals <- ts %>%
    dplyr::pull(ncol(ts))

  # Make the detection timeseries:
  len <- length(vals)
  detect <- as_detect_mad(vals, lowwl=lowwl, highwl=highwl,
                          mad_thr=mad_thr, mad_cst=mad_cst)

  # Return position(s) of shift(s) (or of the maximum detected value):
  where <- try(
    where_as_quiet(detect, thresh = asd_thr, quiet=TRUE), silent=TRUE)

  # In case of error due to a threshold too low,
  # a threshold of 1 returns the maximum (absolute) value:
  if(class(where)[1]=="try-error"){

    where_low <- try(where_as_quiet(detect, thresh = 1, quiet=TRUE))
    where <- where_low

  }


  # Summarize as_detect outputs in one line (with potentially multiple shifts):
  asd_out <- data.frame(abbr="asd",
                        mtd="asdetect",
                        n_brk=0,
                        loc_brk=NA,
                        aic=NA,
                        trend=NA,
                        mag=NA,
                        rel_chg=NA,
                        SDbef=NA,
                        SDaft=NA,
                        nrmse=NA,
                        abruptness=NA)
  run <- list()

  if (length(where$as_pos)>0){

    for (i in 1:length(where$as_pos)){ # If multiple shifts detected

      if(check_true_shift){ # If the additional check is performed

        if(abs(where$dt_val[i])>asd_thr){ # If detected shift is above threshold

          # Returns 1 if true abrupt shift detected,
          # using fixed section length for short timeseries
          if (len<20) true_shift <-
              custom_shift_type(vals, where$as_pos[i], width = 2)

          else true_shift <-
              custom_shift_type(vals, where$as_pos[i], width = "tenth")

        } else {true_shift <- 0}

        # Adds confirmed shift(s) to the summary:
        if (abs(where$dt_val[i])>asd_thr & true_shift == 1){

          asd_out["loc_brk"] <- paste0(asd_out["loc_brk"],";",
                                       c(time[where$as_pos][i], NA)[1])
          asd_out["n_brk"] <- asd_out["n_brk"] + 1
          run <- c(run, list(time[where$as_run[[i]]]))
        }

      } else { # Without additional check

        # Adds shift(s) to the summary:
        if (abs(where$dt_val[i])>asd_thr){

          asd_out["loc_brk"] <- paste0(asd_out["loc_brk"],";",
                                       c(time[where$as_pos][i], NA)[1])
          asd_out["n_brk"] <- asd_out["n_brk"] + 1
          run <- c(run, list(time[where$as_run[[i]]]))
        }
      }
    }
  }

  # To keep length non null:
  if(length(run)==0) run <- NA

  asd_out["loc_brk"] <- sub("NA;", "", asd_out["loc_brk"])

  return(list("df" = asd_out,
              "detect" = detect,
              "run" = run))

}



#' Breakpoints analysis using 'chngpt' package
#'
#' @param ts Timeseries to analyse as 'ts' object.
#'
#' @return List of two objects:
#' - one-row data frame with info about potentially detected breakpoints
#' - vector of model fitted values
#'
#' @references Fong Y., Huang Y., Gilbert P. B. & Permar S. R. (2017) 'chngpt:
#' threshold regression model estimation and inference' BMC Bioinformatics.
#'
#' @export

chg_fct <- function(ts){

  # Fit breakpoint model on the timeseries:
  Y <- tail(names(ts),1)
  chg <- chngpt::chngptm(formula.1 = as.formula(paste(Y,"~1")),
                         formula.2 = ~year,
                         type="step", family="gaussian",
                         data=ts %>%
                           dplyr::mutate(year=as.numeric(year)) # for dates
                         )

  # Get predicted values from the model:
  pred_chg <- data.frame(year = ts$year,
                         bp = chg$best.fit$fitted.values)

  # Compute root mean square error normalized by standard deviation of the fit:
  nrmse <- sqrt(sum(residuals(chg)^2)/length(ts$Y))/sd(ts$Y)

  # Summarize chngt outputs in one line:
  chg_out <- data.frame(abbr = "chg",
                        mtd = "chgnpt",
                        n_brk = 1,
                        loc_brk = chg$chngpt,
                        aic = MuMIn::AICc(chg),
                        trend = ifelse(pred_chg$bp[1] >
                                         pred_chg$bp[length(pred_chg$bp)],
                                       "decrease", "increase"),
                        mag = pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1],
                        rel_chg =
                          (pred_chg$bp[length(pred_chg$bp)] - pred_chg$bp[1]) /
                          max(abs(pred_chg$bp[length(pred_chg$bp)]),
                                    abs(pred_chg$bp[1])),
                        SDbef = ts %>%
                          dplyr::filter(year<=chg$chngpt) %>%
                          dplyr::pull(Y) %>%
                          sd(),
                        SDaft = ts %>%
                          dplyr::filter(year>=chg$chngpt) %>%
                          dplyr::pull(Y) %>%
                          sd(),
                        nrmse = nrmse
  )
  chg_out <- chg_out %>% mutate(abruptness = mag/((SDbef+SDaft)/2))

  return(list("chg_out" = chg_out, "pred_chg" = pred_chg))

}



#' Multiple breakpoints analysis
#'
#' @param ts Timeseries to analyses as data frame (2 columns: year, Y).
#' @param abr_mtd Vector with abbreviation(s) corresponding to the
#' breakpoints method(s) to use ("asd" and/or "chg").
#' @param asd_thr (asdetect) Numeric threshold in detection timeseries.
#' @param asd_chk (asdetect) Logical parameter for check_true_shift in asd_fct.
#' @param lowwl (asdetect) Lowest window length used in algorithm (default 5).
#' @param highwl (asdetect) Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of timeseries.
#' @param mad_thr (asdetect) Threshold of anomalous change in number
#' of median absolute deviations (default 3).
#' @param mad_cst (asdetect) Correction factor for asymptotic
#' normal consistency.
#'
#' @return List of four objects (if both chg and asd methods used):
#' - data frame with info about potentially detected breakpoints
#' - list of chg method output with model fitted values
#' - vector of detection values from asd method
#' - list of runs of uncertainty for asdetect breakpoints
#'
#' @export

shifts <- function(ts, abr_mtd, asd_thr, asd_chk,
                   lowwl, highwl, mad_thr, mad_cst){

  res_table <- data.frame()

  # Run asdetect method:
  if ("asd" %in% abr_mtd){

    asd_out <- asd_fct(ts, asd_thr, check_true_shift=asd_chk,
                       lowwl, highwl, mad_thr, mad_cst)
    res_table <- rbind(res_table, asd_out$df)
  }

  # Run chngpt method:
  if ("chg" %in% abr_mtd){

    chg_outlist <- chg_fct(ts)
    chg_out <- chg_outlist$chg_out

    res_table <- rbind(res_table, chg_out)
  }

  # Return the output depending on the methods used:
  if ("asd" %in% abr_mtd & "chg" %in% abr_mtd) {
    return(list("res_table" = res_table, "chg_outlist" = chg_outlist,
                "asd_detect" = asd_out$detect, "asd_run" = asd_out$run))

  } else if ("chg" %in% abr_mtd) {
    return(list("res_table" = res_table, "chg_outlist" = chg_outlist))

  } else {
    return(list("res_table" = res_table))
  }

}



#' Summary of breakpoints analyses
#'
#' @param sets List of data frame ready for analyses.
#' @param abr_mtd Vector with abbreviation(s) corresponding to the
#' breakpoints method(s) to use ("asd" and/or "chg").
#' @param asd_thr (asdetect) Numeric threshold in detection timeseries.
#' @param asd_chk (asdetect) Logical parameter for check_true_shift in asd_fct.
#' @param lowwl (asdetect) Lowest window length used in algorithm (default 5).
#' @param highwl (asdetect) Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of timeseries.
#' @param mad_thr (asdetect) Threshold of anomalous change in number
#' of median absolute deviations (default 3).
#' @param mad_cst (asdetect) Correction factor for asymptotic
#' normal consistency.
#'
#' @return List of two objects:
#' - list of results from breakpoints analyses for each timeseries
#' - list of main output from breakpoints analyses by breakpoint method
#'
#' @export

abrupt_classif <- function(sets, abr_mtd,
                           asd_thr, asd_chk, lowwl, highwl, mad_thr, mad_cst){

  ts <- sets %>%
    lapply(function(x) dplyr::select(x, c(X,Y)) %>%
             dplyr::rename(year=X))

  # Perform breakpoint analyses on each timeseries:
  shifts_res <- ts %>%
    lapply(function(x) shifts(x, abr_mtd, asd_thr, asd_chk,
                              lowwl, highwl, mad_thr, mad_cst))

  # Make list of analyses main outputs by breakpoint method:
  abt_res <- list()
  for(k in abr_mtd){

    abt_res_mtd <- shifts_res %>%
      lapply(function(x) x$res_table %>%
               dplyr::filter(abbr==k)) %>%
      do.call(what=rbind.data.frame)
    abt_res[[k]] <- abt_res_mtd
  }


  return(list("shifts_res" = shifts_res, "abt_res" = abt_res))
}



### Custom functions from asdetect -----------------------------------

#' Create detection time series
#'
#' @description Adapted from 'as_detect' from the asdetect package to control
#' the MAD (median absolute deviation) threshold
#'
#' @param ts Time series to detect shifts in.
#' @param dt Time step of the time series. Default is NA.
#' @param lowwl Lowest window length used in algorithm. Default is 5.
#' @param highwl Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of ts.
#' @param mad_thr Threshold of anomalous change in number
#' of median absolute deviations. Default is 3.
#' @param mad_cst Correction factor for asymptotic normal consistency (1.4826).
#'
#' @return Returns a numeric array which is the detection time series,
#' the same length as the inputted time series. If dt is supplied rather
#' than kept as the default NA then a data frame is returned including:
#'
#' t  Time values
#' detect  Detection time series
#'
#' @export

as_detect_mad <- function (ts, dt = NA, lowwl = 5, highwl = "default",
                           mad_thr=3, mad_cst=1.4826){

  l <- length(ts)
  if (highwl == "default") {
    higwl <- floor(l/3)
  } else {
    higwl <- highwl
  }

  # tip_distrib: counts how often timepoints are involved in windows
  # exceeding mad_thr*MAD:
  tip_distrib <- rep(0, l)
  wls <- lowwl:higwl # set of window lengths

  for (j in 1:length(wls)) { # For the different window length
    breaks <- l/wls[j]
    if (breaks != floor(breaks)) {
      remainder <- l - floor(breaks) * wls[j]
      ts_temp <- ts[(floor(remainder/2) + 1):(floor(remainder/2) +
                                                floor(breaks) * wls[j])]
      # ts_temp: longest timeseries divisible by window length
    }
    if (breaks == floor(breaks)) {
      remainder <- NA
      ts_temp <- ts
    }

    lm_coeffs <- array(NA, dim = c(breaks, 2))
    # lm_coeff: array to store intercept and slope of lm from each window
    # of a given length

    for (i in 1:breaks) { # For the different windows
      lm_mod <- lm(ts_temp[1:wls[j] + ((i - 1) * wls[j])] ~
                     c(1:wls[j]))
      # lm_mod: result of the linear regression made on a given window
      lm_coeffs[i, ] <- lm_mod$coefficients
    }

    outlier_ind_up <-
      which((lm_coeffs[, 2] - median(lm_coeffs[, 2]))/
              mad(lm_coeffs[, 2], constant=mad_cst) > mad_thr)
    outlier_ind_down <-
      which((lm_coeffs[, 2] - median(lm_coeffs[, 2]))/
              mad(lm_coeffs[, 2], constant=mad_cst) < -mad_thr)
    # outlier_ind_up: list of windows with a slope above 3*MAD
    # outlier_ind_down: list of windows with a slope below -3*MAD

    if (is.na(remainder)) {
      for (k in outlier_ind_up) {
        tip_distrib[1:wls[j] + ((k - 1) * wls[j])] <-
          tip_distrib[1:wls[j] + ((k - 1) * wls[j])] + 1
      }
      for (k in outlier_ind_down) {
        tip_distrib[1:wls[j] + ((k - 1) * wls[j])] <-
          tip_distrib[1:wls[j] + ((k - 1) * wls[j])] - 1
      }
    }
    if (!is.na(remainder)) {
      for (k in outlier_ind_up) {
        tip_distrib[(floor(remainder/2) + 1) +
                      1:wls[j] + ((k - 1) * wls[j])] <-
          tip_distrib[(floor(remainder/2) + 1) +
                        1:wls[j] + ((k - 1) * wls[j])] + 1
      }
      for (k in outlier_ind_down) {
        tip_distrib[(floor(remainder/2) + 1) +
                      1:wls[j] + ((k - 1) * wls[j])] <-
          tip_distrib[(floor(remainder/2) + 1) +
                        1:wls[j] + ((k - 1) * wls[j])] - 1
      }
    }
  }

  # dt: timestep of the timeseries, if equals one leave default dt=NA
  if (is.na(dt)) {
    result <- tip_distrib/length(wls)
    # result: cumulative contribution from all window lengths divided by
    # the number of window lengths
  } else {
    t <- rep(NA, length(tip_distrib))
    t[1] <- 0
    for (i in 2:length(t)) {
      t[i] <- t[i - 1] + dt
    }
    result <- data.frame(t = t, detect = tip_distrib/length(wls))
  }
  return(result)
}

# To call hidden functions from asdetect:
environment(as_detect_mad) <- asNamespace('asdetect')



#' Determine where abrupt shift occurs
#'
#' @description Adapted from 'where_as' from the asdetect package to
#' retrieve the extent of flat extrema and to make it silent
#'
#' @param ts Time series, most commonly a detection time series created
#' from as_detect(). Any numerical vector will work.
#' @param dt Time step of the time series. Default is NA.
#' @param thresh The threshold value to look for maximum values in. It will
#' also look for negative values such that their absolute value is above
#' the threshold.
#' @param quiet Logical to allow to print if no shift detected above threshold.
#'
#' @return A list is returned comprising of the components:
#'
#' as_pos	 The position(s) of the detected abrupt shifts. If dt is supplied,
#' it will list the time rather than position.
#' dt_val  The detection value at the position or time of detection.
#' as_run The position(s) of extrema in detection timeseries (as uncertainty).
#'
#' If the maximum value of the (absolute) time series is less than the
#' threshold searched for, the maximum detection value and corresponding
#' position or time is returned, along with the message
#' 'Threshold not detected, maximum returned instead' is printed.
#'
#' @export

where_as_quiet <- function (ts, dt = NA, thresh = 0.7, quiet = FALSE){

  inds <- which(abs(ts) > thresh)
  # inds: position of timepoints with detection line above threshold

  if (length(inds) > 0) { # if any point above threshold

    encoding <- rle(diff(inds))
    # encoding: lengths and values of runs of lag-1 difference

    numtip <- length(which(encoding$value == 1))
    # numtip: number of runs above threshold

    whichtip <- which(encoding$value == 1)
    # whichtip: position of runs in encoding

    tip_run <- list()
    # tip_run: list of run(s) with min or max detection score
    tip_pos <- rep(NA, numtip)
    # tip_pos: list of position(s) with min or max detection score (median of runs)
    tip_prob <- rep(NA, numtip)
    # tip_prob: list of detection score(s) of min or max

    for (k in 1:numtip) {
      if (k == 1) {
        if (ts[inds[1]] > 0) {

          tip_run[k] <- list(
            inds[which(ts[ inds[1:(encoding$length[1] + 1)] ] ==
                         max( ts[ inds[1:(encoding$length[1] + 1)] ]))])
          # if the maximum is spread over several timepoint,
          # it returns the median position:
          tip_pos[k] <- floor(median(tip_run[[k]]))

        }
        if (ts[inds[1]] < 0) {

          tip_run[k] <- list(
            inds[which(ts[ inds[1:(encoding$length[1] + 1)] ] ==
                         min( ts[ inds[1:(encoding$length[1] + 1)] ]))])
          tip_pos[k] <- floor(median(tip_run[[k]]))

        }
        tip_prob[k] <- ts[tip_pos[k]]
      } else {

        inds_temp <-
          inds[sum(encoding$length[1:(whichtip[k] - 1)]):
                 sum(encoding$length[1:whichtip[k]]) + 1]
        # inds_temp: for run k, positions of the kth run

        if (ts[inds_temp[1]] > 0) {

          tip_run[k] <-
            list(inds_temp[which(ts[inds_temp] == max(ts[inds_temp]))])
          tip_pos[k] <- floor(median(tip_run[[k]]))

        }
        if (ts[inds_temp[1]] < 0) {

          tip_run[k] <-
            list(inds_temp[which(ts[inds_temp] == min(ts[inds_temp]))])
          tip_pos[k] <- floor(median(tip_run[[k]]))

        }
        tip_prob[k] <- ts[tip_pos[k]]
      }
    }

    if (!is.na(dt)) {
      t <- rep(NA, length(ts))
      t[1] <- 0
      for (i in 2:length(t)) {
        t[i] <- t[i - 1] + dt
      }
      tip_pos <- t[tip_pos]
    }
    results <- list()
    results$as_pos <- tip_pos
    results$dt_val <- tip_prob
    results$as_run <- tip_run
    return(results)
  }

  # If no point above threshold, still return the maximum with a warning:
  if (length(inds) == 0) {
    results <- list()
    results$as_pos <- which.max(abs(ts))
    if (max(abs(ts)) == max(ts)) {
      results$dt_val <- max(ts)
    }
    if (-max(abs(ts)) == min(ts)) {
      results$dt_val <- min(ts)
    }
    results$as_run <- which(ts == results$dt_val)
    if (!quiet) print("Threshold not detected, maximum returned instead")

    return(results)
  }
}

# To call hidden functions from asdetect:
environment(where_as_quiet) <- asNamespace('asdetect')



#' Determine abrupt shift type
#'
#' @description Adapted from 'shift_type' from the asdetect package to solution
#' a limit case
#'
#' @param ts Time series to detect shifts in.
#' @param where_as_pos Output from the where_as() function.
#' @param dt Time step of the time series ts. As default is FALSE but is
#' used in calculations if supplied as a numerical value into the function.
#' It must be supplied if it was used in the where_as() function.
#' @param width The width of the section used to measure the gradient around
#' the detected shift. Default is 'tenth' which will use 10 percent of the
#' length of the time series. A numerical value can be supplied, which will
#' be either in points, or time if dt is also supplied. Note user supplied
#' values are not in percentages as the default is. The value of width is
#' divided in half and the section used to measure the gradient around the
#' detected tipping point is spread evenly across the detected position.
#'
#' @return A value of 1 is returned if the function determines that a
#' true abrupt shift has occurred, and 0 otherwise.
#'
#' @export

custom_shift_type <- function (ts, where_as_pos, dt = FALSE, width = "tenth")
{
  if (width == "tenth") {
    w <- floor(length(ts)/20) # shift_type needs timeseries of length 20 or more
  }
  else {
    w <- floor(width/2)
  }
  if (dt == FALSE) {
    pos <- where_as_pos
  }
  else {
    pos <- where_as_pos/dt
    if (width != "tenth") {
      w <- w/dt
    }
  }

  if (pos <= w) {
    # Include case when pos == w otherwise (pos - w) will equal 0 and will cause
    # length problem below
    poss <- 1:(pos + w)
  }
  else if (pos > (length(ts) - w)) {
    poss <- (pos - w):length(ts)
  }
  else {
    poss <- (pos - w):(pos + w)
  }
  as_mod <- lm(ts[poss] ~ c(1:length(poss)))
  full_mod <- lm(ts ~ c(1:length(ts)))
  as_grad <- as_mod$coefficients[2]
  full_grad <- full_mod$coefficients[2]
  if (abs(as_grad) < abs(full_grad)) {
    marker <- 0
  }
  else {
    marker <- 1
  }
  return(marker)
}

# To call hidden functions from asdetect:
environment(custom_shift_type) <- asNamespace('asdetect')



## III-C) Classification ---------------------------------------------------

#' Make trajectory fitting
#'
#' @param sets List of data frame ready for analyses.
#' @param abr_mtd Vector of abbreviation with breakpoints
#' methods to use ("asd" and/or "chg").
#' @param type Character specifying whether the series come
#' from simulations ("sim"), RAMLDB ("RAM"), or other empirical data ("data").
#' @param asd_thr (asdetect) Numeric threshold in detection timeseries.
#' @param asd_chk (asdetect) Logical parameter for check_true_shift in asd_fct.
#' @param lowwl (asdetect) Lowest window length used in algorithm (default 5).
#' @param highwl (asdetect) Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of timeseries.
#' @param mad_thr (asdetect) Threshold of anomalous change in number
#' of median absolute deviations (default 3).
#' @param mad_cst (asdetect) Correction factor for asymptotic
#' normal consistency.
#' @param apriori Logical to state whether expected trajectories is indicated.
#'
#' @return List of two objects:
#' - data frame with the results of different trajectory fitting
#' - list by trajectory type of the main outputs
#'
#' @export

fit_models <- function(sets, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl,
                       mad_thr, mad_cst, apriori){

  # Running the different model fitting:
  res_nch <- res_trend(sets, niter=1, correction=FALSE, fit="nch")
  res_lin <- res_trend(sets, niter=1, correction=FALSE, fit="lin")
  res_pol <- res_trend(sets, niter=1, correction=FALSE, fit="pol")
  res_abt <- abrupt_classif(sets, abr_mtd, asd_thr, asd_chk,
                            lowwl, highwl, mad_thr, mad_cst)
  lengths <- data.frame(first = sapply(sets, function(x) min(x$X),
                                       simplify="array"),
                        last = sapply(sets, function(x) max(x$X),
                                      simplify="array"),
                        length = sapply(sets, function(x) nrow(x),
                                        simplify="array"))

  # Combine all results:
  summ_res <- cbind(
    res_nch %>% dplyr::select(contains("max_shape")|contains("aic")|
                                contains("trend")|contains("nrmse")) %>%
      dplyr::rename_with(~str_c(., "_nch")),
    res_lin %>% dplyr::select(contains("max_shape")|contains("aic")|
                                contains("trend")|contains("nrmse")|slope) %>%
      dplyr::rename_with(~str_c(., "_lin")),
    res_pol %>% dplyr::select(contains("max_shape")|contains("aic")|
                                contains("best_model")|contains("trend")|
                                contains("nrmse")) %>%
      dplyr::rename_with(~str_c(., "_pol")),
    do.call("cbind",
            lapply(abr_mtd, function(x) res_abt$abt_res[[x]] %>%
                     dplyr::select(contains("brk")|aic|trend|nrmse) %>%
                     dplyr::mutate(max_shape=paste0(n_brk,"_breakpoint")) %>%
                     dplyr::rename_with(~str_c(., paste0("_", x)))
                   )
            ),
    res_abt$abt_res[["chg"]] %>%
      dplyr::select(mag, rel_chg, SDbef, SDaft, abruptness),
    lengths
  ) %>% dplyr::rename(signif_model = best_model_pol,
                      slope = slope_lin)

  # Add expectations if known:
  if(type=="sim" & apriori){

    expected <- do.call(rbind, (lapply(sets, function(x) x[1,]))) %>%
      dplyr::select(expected_traj, expected_class)
    summ_res <- cbind(summ_res, expected)

  }

  return(list("summ_res"=summ_res,
              "res_fit"=list("res_nch"=res_nch,
                             "res_lin"=res_lin,
                             "res_pol"=res_pol,
                             "res_abt"=res_abt)
              )
         )
}



#' Run timeseries classification
#'
#' @param sets List of data frame ready for analyses.
#' @param str Character indicating whether to use 'asdetect'
#' breakpoint validation ("aic" or "aic_asd").
#' @param abr_mtd Vector of abbreviation with breakpoints
#' methods to use ("asd" and/or "chg").
#' @param type Character specifying whether the series come
#' from simulations ("sim"), RAMLDB ("RAM"), or other empirical data ("data").
#' @param noise_comb For simulated timeseries, a character indicating
#' the noise combination used to generate the timeseries.
#' @param apriori Logical, whether expected trajectories are indicated.
#' @param run_loo Logical, whether to perform leave-one-out process.
#' @param two_bkps Logical, if true, looks for more than one breakpoints.
#' @param smooth_signif Logical, should significance of the coefficient be taken
#'  into account in the selection of best model.
#' @param asd_thr (asdetect) Numeric threshold in detection timeseries.
#' @param asd_chk (asdetect) Logical parameter for check_true_shift in asd_fct.
#' @param lowwl (asdetect) Lowest window length used in algorithm (default 5).
#' @param highwl (asdetect) Highest window length used in algorithm.
#' If 'default' then highwl is set to 1/3 of the length of timeseries.
#' @param mad_thr (asdetect) Threshold of anomalous change in number
#' of median absolute deviations (default 3).
#' @param mad_cst (asdetect) Correction factor for asymptotic
#' normal consistency.
#' @param edge_lim Numeric, minimal breakpoint distance to start or end dates
#' (default 5 timesteps).
#' @param congr_brk Numeric, maximal acceptable distance between chngpt and
#' as_detect breaks (default 5 timesteps).
#' @param showplots Logical that indicate whether to show plots
#' (slows down computation if true).
#' @param save_plot Logical, if plots are shown, whether to save plots
#' (either the four trajectory fits or one of a given fit).
#' @param save_plot_bis Logical, if plots are shown,whether to save plot_bis
#' (best fit with asdetect detection curve).
#' @param outplot Logical, last plot in return value.
#' @param ind_plot Character to return the last plot for a given trajectory
#' (either "nch", "lin", "pol", "abt", or "best").
#' @param detection_plot Logical to plot detection plot when abrupt.
#' @param plot_one_in For simulations, to limit the number of plots saved,
#' the number of timeseries out of which to save the plot (default 10).
#' @param dirname Directory name where to save plots.
#'
#' @return List of three to four objects:
#' - data frame with the results of different trajectory fitting,
#' - data frame with the output of the best fitting model,
#' - data frame with detailed results of different trajectory fitting,
#' - (optional) plot with either the four trajectory fits or one of a given fit
#'
#' @export

traj_class <- function(sets, str, abr_mtd, type="sim", noise_comb=NULL,
                       apriori=FALSE, run_loo, two_bkps, smooth_signif,
                       asd_thr, asd_chk, lowwl=5, highwl="default",
                       mad_thr=3, mad_cst=1.4826, edge_lim=5, congr_brk=5,
                       showplots=FALSE, save_plot=TRUE, save_plot_bis=TRUE,
                       outplot=FALSE, ind_plot=NULL, detection_plot=TRUE,
                       dirname=NULL, plot_one_in=10){

  ### Original timeseries

  # Make trajectory fitting:
  res <- fit_models(sets$ts, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl,
                    mad_thr, mad_cst, apriori)

  # Define best trajectory according to decision tree:
  best_traj <- best_traj_aic(class_res=res, type=type, apriori=apriori,
                             aic_selec=str, smooth_signif=smooth_signif,
                             edge_lim=edge_lim, congr_brk=congr_brk)

  # If looking for two breakpoints:
  if(two_bkps){

    # Keep timeseries with breakpoints:
    list_brk_ts <- best_traj %>%
      dplyr::filter(class=="abrupt") %>%
      dplyr::select(simu_id, loc_brk_chg)

    # If some timeseries have breakpoints:
    if (nrow(list_brk_ts) != 0){

      # Split timeseries (keeping the breakpoint in each):
      sets2 <- list()
      for (i in 1:nrow(list_brk_ts)){

        sets2[[paste0(list_brk_ts$simu_id[i],"_1st")]] <-
          sets$ts[[list_brk_ts$simu_id[i]]] %>%
          dplyr::filter(X <= list_brk_ts$loc_brk_chg[i])

        sets2[[paste0(list_brk_ts$simu_id[i],"_2nd")]] <-
          sets$ts[[list_brk_ts$simu_id[i]]] %>%
          dplyr::filter(X >= list_brk_ts$loc_brk_chg[i])
      }

      # Remove too short sub timeseries:
      sets2 <- sets2[lapply(sets2, nrow)>=15]

      # If any timeseries eligible:
      if(length(sets2)>0){

        # Fit models of splitted timeseries:
        res2 <- fit_models(sets2, abr_mtd, type, asd_thr, asd_chk, lowwl, highwl,
                           mad_thr, mad_cst, apriori)
        best_traj_split <- best_traj_aic(class_res=res2, type=type,
                                         apriori=apriori, aic_selec=str,
                                         smooth_signif=smooth_signif,
                                         edge_lim=edge_lim, congr_brk=congr_brk)

        # Locate auxiliary breakpoints:
        list_brk_ts2 <- best_traj_split %>%
          dplyr::mutate(part = factor(gsub("^.*_", "", simu_id),
                                      levels=c("1st", "2nd")),
                        simu_id = gsub("_1st|_2nd", "", simu_id)) %>%
          dplyr::filter(traj == "abrupt") %>%
          dplyr::select(simu_id, loc_brk_chg, part) %>%
          # tidyr::pivot_wider(names_from = part, values_from = loc_brk_chg) %>%
          tidyr::spread(key = part, value = loc_brk_chg, drop=FALSE) %>%
          dplyr::rename(loc_aux1_chg = "1st",
                        loc_aux2_chg = "2nd")

        best_traj <- dplyr::left_join(best_traj, list_brk_ts2, by="simu_id") %>%
          dplyr::mutate(nb_brk = dplyr::case_when(
            class != "abrupt" ~ 0,
            is.na(loc_aux1_chg) & is.na(loc_aux2_chg) ~ 1,
            is.na(loc_aux1_chg) | is.na(loc_aux2_chg) ~ 2,
            !(is.na(loc_aux1_chg) & is.na(loc_aux2_chg)) ~ 3
          ))

        # Specify trajectory of splitted timeseries

        best_traj <- best_traj %>%
          dplyr::left_join(
            best_traj_split %>%
              dplyr::mutate(part = factor(gsub("^.*_", "", simu_id),
                                          levels=c("1st", "2nd")),
                            simu_id = gsub("_1st|_2nd", "", simu_id)) %>%
              dplyr::select(simu_id, part, traj, trend, class) %>%
              tidyr::pivot_wider(names_from = part, values_from = -c(simu_id, part)),
            by="simu_id")


        # tidyr::pivot_wider(names_from = part, values_from = loc_brk_chg) %>%

      } else {

        best_traj <- best_traj %>%
          dplyr::mutate(loc_aux1_chg = NA,
                        loc_aux2_chg = NA,
                        nb_brk = ifelse(class == "abrupt", 1, 0),
                        traj_1st = NA, traj_2nd = NA,
                        trend_1st = NA, trend_2nd = NA,
                        class_1st = NA, class_2nd = NA
          )
      }

    } else {

      best_traj <- best_traj %>%
        dplyr::mutate(loc_aux1_chg = NA,
                      loc_aux2_chg = NA,
                      nb_brk = ifelse(class == "abrupt", 1, 0),
                      traj_1st = NA, traj_2nd = NA,
                      trend_1st = NA, trend_2nd = NA,
                      class_1st = NA, class_2nd = NA
        )
    }
  }


  ### Leave-One-Out timeseries
  if (run_loo){

    # Make all leave-one-out timeseries and name them accordingly:
    loo <- lapply(sets$ts, function(x) lapply(1:nrow(x), function(k) x[-k,]))
    for (i in 1:length(loo)) names(loo[[i]]) <-
        paste0(names(loo)[i], "_", paste0("loo", 1:length(loo[[i]])))

    # Make trajectory fitting for loo:
    res_loo <- lapply(loo, function(x) fit_models(x, abr_mtd, type, asd_thr,
                                                  asd_chk, lowwl, highwl,
                                                  mad_thr, mad_cst, apriori))

    # Define best trajectory according to decision tree for loo:
    best_traj_loo <- lapply(res_loo, function(x)
      best_traj_aic(class_res=x, type=type, apriori=apriori,
                    aic_selec=str,smooth_signif=smooth_signif,
                    edge_lim=edge_lim, congr_brk=congr_brk))

    best_traj_loo <- lapply(1:length(sets$ts), function(i)
      best_traj_loo[[i]] %>%
        dplyr::mutate(X=sets$ts[[i]]$X))

    # LOO frequencies for class:
    loo_class <- do.call("rbind",  lapply(1:length(sets$ts), function(i) {

      class_freq <- best_traj_loo[[i]] %>%
        dplyr::group_by(class) %>%
        dplyr::summarise(class_freq_loo = n()/nrow(best_traj_loo[[i]])) %>%
        tidyr::pivot_wider(names_from = class, values_from = class_freq_loo)

      missing <- setdiff(levels(best_traj_loo[[i]]$class), names(class_freq))

      class_freq[missing] <- 0

      class_freq <- class_freq[levels(best_traj_loo[[i]]$class)] %>%
        dplyr::select_all(list(~ paste0("loo_", .))) %>%
        dplyr::mutate(simu_id = names(sets$ts)[i])


    }))

    best_traj <- best_traj %>%
      dplyr::left_join(loo_class, by="simu_id")

  }

  # store asd parameters:
  if (str == "aic_asd") best_traj <- best_traj %>%
      dplyr::mutate(asd_thr=asd_thr,
                    lowwl=lowwl,
                    highwl=highwl,
                    mad_thr=mad_thr)

  # Plot the 4 different fits with best one colored:
  if(showplots == TRUE){

    if(!run_loo) best_traj_loo <- NULL
    if(sapply(sets$ts[[1]], lubridate::is.Date)[[1]]==TRUE){ is_date <- TRUE
    } else { is_date <- FALSE }

    ## Make plots:
    plots_traj_nch <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_nch,
                                          best_traj, plot_class="no_change",
                                          best_traj_loo=best_traj_loo,
                                          ind_plot=ind_plot)
    plots_traj_lin <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_lin,
                                          best_traj, plot_class="linear",
                                          best_traj_loo=best_traj_loo,
                                          ind_plot=ind_plot)
    plots_traj_pol <- plot_traj_multi_abt(sets, rslt=res$res_fit$res_pol,
                                          best_traj, plot_class="quadratic",
                                          best_traj_loo=best_traj_loo,
                                          ind_plot=ind_plot)
    plots_traj_chg <- plots_traj_abt <-
      plot_traj_multi_abt(sets, rslt=res$res_fit$res_abt, best_traj,
                          plot_class="abrupt", asd_thr,
                          best_traj_loo=best_traj_loo,
                          detection_plot=detection_plot,
                          ind_plot=ind_plot,
                          is_date=is_date)

    if(!is.null(ind_plot)){
      if(ind_plot=="best"){

        ind_plots <- get(paste0("plots_traj_",
                                sub("max_shape_","",best_traj$best_model)))
      } else {
        ind_plots <- get(paste0("plots_traj_",ind_plot))
      }
    }


    for (i in 1:nrow(res$summ_res)){


      if(str=="aic_asd") class_plot <-
          cowplot::plot_grid(plots_traj_nch$plots[[i]], plots_traj_lin$plots[[i]],
                             plots_traj_pol$plots[[i]], plots_traj_abt$plots[[i]],
                             nrow=2, ncol=2, align = 'h')
      if(str=="aic") class_plot <-
          cowplot::plot_grid(plots_traj_nch$plots[[i]], plots_traj_lin$plots[[i]],
                             plots_traj_pol$plots[[i]], plots_traj_abt$plots[[i]],
                             nrow=2, ncol=2, align = 'h')

      if(str=="aic_asd") endname <- paste0("_", str, "_coef", smooth_signif,
                                           "asd_thr", asd_thr, "l",lowwl, "h",
                                           highwl, "mad", mad_thr,".png")
      if(str=="aic") endname <- paste0("_",str,"_coef",smooth_signif,".png")

      if(type=="sim" & save_plot & i%%plot_one_in==0){

        png(filename =
              paste0(dirname,
                     sets$ts[[i]]$expected_class[1],"/",
                     paste(head(stringr::str_split(
                       names(sets$ts)[i],"_")[[1]],-1), collapse="_"),"_",
                     noise_comb,"_", gsub("^.*_", "", names(sets$ts)[i]),"__",
                     endname), width=8, height=6, units="in", res=300)
        print(class_plot)
        dev.off()


        if(str=="aic_asd" & save_plot_bis){

          png(filename =
                paste0(dirname,
                       sets$ts[[i]]$expected_class[1],"abr_",
                       paste(head(stringr::str_split(
                         names(sets$ts)[i],"_")[[1]],-1), collapse="_"),"_",
                       noise_comb,"_", gsub("^.*_", "", names(sets$ts)[i]),
                       endname), width=6, height=6, units="in", res=300)
          print(plots_traj_abt$plot_bis[[i]])
          dev.off()
        }

      }

      if ((type=="RAM" | type=="data") & save_plot){

        class_plot <- class_plot +
          patchwork::plot_annotation(sub("_iter01", "", names(sets$ts)[i]) %>%
                                       sub("_"," ",.),
                                     theme = theme(plot.title =
                                                     element_text(hjust = 0.5)))

        png(filename = paste0(dirname,
                              sub("_iter01","",names(sets$ts)[i]), endname),
            width=8, height=6, units="in", res=300)
        print(class_plot)
        dev.off()


        if(str=="aic_asd" & save_plot_bis){

          png(filename = paste0(dirname,"abr_",
                                sub("_iter01","",names(sets$ts)[i]), endname),
              width=6, height=6, units="in", res=300)
          print(plots_traj_abt$plot_bis[[i]])
          dev.off()
          }
        }
    }

  }

  if (outplot){ # For last plot only

    if(!is.null(ind_plot)){
      return(list("res"=res$summ_res, "best_traj"=best_traj,
                  "res_detail"=res$res_fit, "class_plot"=ind_plots))

    } else {
      return(list("res"=res$summ_res, "best_traj"=best_traj,
                  "res_detail"=res$res_fit, "class_plot"=class_plot))
    }
  } # If no plot to include in the output (lighter output):

  return(list("res"=res$summ_res, "best_traj"=best_traj,
              "res_detail"=res$res_fit))
}



#' Determine best trajectory based on AIC comparisons
#'
#' @param class_res Data frame with results from the different classifications.
#' @param type Character specifying whether the series come
#' from simulations ("sim"), RAMLDB ("RAM"), or other empirical data ("data").
#' @param apriori Logical, whether expected trajectories are indicated.
#' @param aic_selec Character specifying how the best trajectory is chosen,
#' based on lowest aic only ("aic") or also the confirmation
#' with asdetect ("aic_asd").
#' @param smooth_signif Logical, should significance of the coefficient be taken
#'  into account in the selection of best model.
#' @param edge_lim Numeric, minimal breakpoint distance to start or end dates.
#' @param congr_brk Numeric, maximal acceptable distance between chngpt and
#' as_detect breaks.
#'
#' @return Data frame with the output of the best fitting model,
#' trajectory, and class.
#'
#' @export

best_traj_aic <- function(class_res, type, apriori, aic_selec,
                          smooth_signif, edge_lim, congr_brk){

  traj_lvl <- c("stable_constant","increase_constant","decrease_constant",
                "stable_concave","stable_convex","decrease_accelerated",
                "decrease_decelerated","increase_accelerated",
                "increase_decelerated","1_breakpoint")

  class_lvl <- c("no_change","linear","quadratic","abrupt")

  best_traj <- class_res$summ_res %>%

  # Find minimal AIC and subtract to all AICc:
  tibble::rownames_to_column(var="simu_id") %>%
  dplyr::rowwise() %>%

  dplyr::mutate(min_aic = min(dplyr::c_across(dplyr::contains("aic")),
                              na.rm=TRUE),
                dplyr::across(dplyr::contains("aic") & !"min_aic",
                              ~ .x - min_aic),
                dplyr::across(.names = "weight_{.col}",
                              .cols = dplyr::contains("aic") & !"min_aic",
                              .fns = ~ exp(-0.5*.x)),
                dplyr::across(.cols = dplyr::contains("weight_aic"),
                              .fns = ~ round(.x / sum(dplyr::c_across(
                                dplyr::contains("weight_aic")), na.rm=TRUE),
                                digits=2))) %>%
  dplyr::ungroup()

  # (if as_detect involved) Divide multiple asd breakpoints into several columns
  # and compute distance between chg and asd breakpoints:
  if (aic_selec=="aic_asd"){

    best_traj <- best_traj %>%
      dplyr::mutate(loc_brk_asd_sep = strsplit(loc_brk_asd, ";")) %>%
      tidyr::unnest(loc_brk_asd_sep) %>%
      dplyr::mutate(loc_brk_asd_sep = as.numeric(loc_brk_asd_sep),
                    loc_brk_chg = as.numeric(loc_brk_chg)) %>%
      dplyr::group_by(simu_id) %>%
      dplyr::mutate(row = row_number()) %>%
      tidyr::pivot_wider(names_from=row, names_prefix="loc_brk_asd_",
                         values_from=loc_brk_asd_sep) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::contains('loc_brk_asd_'),
                                  .names = "dup_{.col}")) %>%
      dplyr::rename_with(~sub('dup_loc_brk_asd_',"diff_loc_",.),
                         dplyr::contains('dup_loc_brk_asd_')) %>%
      dplyr::mutate(dplyr::across(.cols = dplyr::contains("diff_loc_"),
                                  ~abs(loc_brk_chg-.x))) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(loc_brk_asd = as.numeric(loc_brk_asd),
                    diff_loc_min = min(dplyr::c_across(
                      dplyr::contains("diff_loc_")), na.rm=TRUE),
                    loc_brk_asd_min = min(dplyr::c_across(
                      dplyr::contains("loc_brk_asd_")), na.rm=TRUE),
                    loc_brk_asd_max = max(dplyr::c_across(
                      dplyr::contains("loc_brk_asd_")), na.rm=TRUE),
                    dplyr::across(.cols = dplyr::contains("loc"),
                                  ~ dplyr::na_if(., Inf)),
                    loc_brk_asd = as.character(loc_brk_asd)) %>%
    suppressWarnings()
  }

  # Expand in length to select by AICc below:
  best_traj <- best_traj %>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(contains("aic") & !contains("weight_aic"),
                        names_to = "best_aic", values_to = "vals")


  # Keep model with lowest AICc:
  if (aic_selec=="aic"){

    # Arrange models according to their AICc:
    rank_aic <- best_traj %>%
      dplyr::group_by(simu_id) %>%
      dplyr::filter(!best_aic %in% c("aic_asd", "min_aic")) %>%
      dplyr::mutate(aic_rank = dplyr::dense_rank(vals)) %>%
      dplyr::select(contains("aic") & !contains("weight_aic") | "simu_id") %>%
      tidyr::pivot_wider(names_from = aic_rank, names_prefix = "aic_rank",
                         values_from = best_aic)

    # Keep best model:
    best_traj <- best_traj %>%
      dplyr::filter(vals == 0) %>%
      dplyr::select(-vals) %>%
      dplyr::left_join(rank_aic, by="simu_id")

    # If best AICc model is abrupt but do not fill in the edge limit criterion
    # take 2nd best model in the ranking:
    best_traj <- best_traj %>%
      dplyr::mutate(best_aic = ifelse(best_aic=="aic_chg" &
                                        (loc_brk_chg-first<edge_lim |
                                           last-loc_brk_chg<edge_lim),
                                      aic_rank2, best_aic)) %>%
      dplyr::mutate(best_aic = sub("aic_","max_shape_", best_aic))

  }


  # Keep lowest AICc but confirm with asdetect if abrupt:
  if (aic_selec=="aic_asd"){

    # Arrange models according to their AICc:
    rank_aic <- best_traj %>%
      dplyr::group_by(simu_id) %>%
      dplyr::filter(!best_aic %in% c("aic_asd", "min_aic")) %>%
      dplyr::mutate(aic_rank = dplyr::dense_rank(vals)) %>%
      dplyr::select(contains("aic") & !contains("weight_aic") | "simu_id") %>%
      tidyr::pivot_wider(names_from = aic_rank, names_prefix = "aic_rank",
                         values_from = best_aic)

    # Keep the model with the lowest AIC:
    best_traj <- best_traj %>%
      dplyr::filter(vals == 0) %>%
      dplyr::select(-vals)

    # If best AICc model is abrupt but do not fill in the following criteria
    # take 2nd best model in the ranking:
    best_traj <- best_traj %>%
      dplyr::left_join(rank_aic, by="simu_id") %>%

      dplyr::mutate(
        not_abr = "not_lowest_aic",

        # If asdetect doesn't find any breakpoints:
        not_abr = ifelse(best_aic=="aic_chg" &
                           max_shape_asd=="0_breakpoint",
                         "no_asd_brk", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" &
                            max_shape_asd=="0_breakpoint",
                          aic_rank2, best_aic),

        # If the breakpoints found are not congruent:
        not_abr = ifelse(best_aic=="aic_chg" & diff_loc_min > congr_brk,
                         "no_congr_brk", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" & diff_loc_min > congr_brk,
                                        aic_rank2, best_aic),

        # If the breakpoints are too close to the edges (for chg breakpoint):
        not_abr = ifelse(best_aic=="aic_chg" &
                           (loc_brk_chg-first<edge_lim |
                              last-loc_brk_chg<edge_lim),
                         "edge_lim", not_abr),

        best_aic = ifelse(best_aic=="aic_chg" &
                            (loc_brk_chg-first<edge_lim |
                               last-loc_brk_chg<edge_lim),
                          aic_rank2, best_aic),

        not_abr = ifelse(best_aic=="aic_chg", NA, not_abr)
      ) %>%

      dplyr::mutate(best_aic = sub("aic_","max_shape_", best_aic))

  }


  # Display best model and trajectory for each simulation:
  best_traj <- best_traj %>%

    dplyr::select(simu_id|contains("max_shape")|contains("expected")|
                    best_aic|signif_model|contains("loc")|contains("weight")|
                    mag|rel_chg|contains("SD")|contains("nrmse")|slope|
                    abruptness|contains("not_abr")) %>%

    tidyr::pivot_longer(contains("max_shape"),
                        names_to = "best_model",
                        values_to = "traj") %>%

    # Keep the trends for each fit (and then only keep the best):
    dplyr::left_join(best_traj %>%
                       dplyr::select(simu_id | contains("trend")) %>%
                       {if(aic_selec == "aic_asd")
                         dplyr::select(., -"trend_asd") else .} %>%
                       tidyr::pivot_longer(contains("trend"),
                                           names_to = "trend_model",
                                           values_to = "trend") %>%
                       dplyr::rename(best_model = trend_model) %>%
                       dplyr::mutate(best_model = sub("trend_", "max_shape_",
                                                      best_model)),
                     by = c("simu_id", "best_model"))


  # Check significance of coefficients for smooth (non-abrupt) trajectories:
  if(smooth_signif){

    best_traj <- best_traj %>%
      dplyr::mutate(signif_model = paste0("max_shape_", signif_model),
                    best_aic = ifelse(best_aic != "max_shape_chg" &
                                        signif_model != best_aic,
                                      signif_model, best_aic))
  }

  # Final reshaping of the results:
  best_traj <- best_traj %>%
    dplyr::filter(best_aic==best_model) %>%
    {if(aic_selec == "aic_asd") dplyr::select(., -best_aic & -weight_aic_asd)
      else dplyr::select(., -best_aic)} %>%
    dplyr::mutate(simu_id = factor(simu_id),
                  traj = traj %>% factor(levels = traj_lvl),
                  class = dplyr::case_when(
                    stringr::str_detect(traj,"stable_constant") ~ "no_change",
                    # acts as smooth significant for 1st order
                    stringr::str_detect(traj, "constant") ~ "linear",
                    stringr::str_detect(traj, "breakpoint") ~ "abrupt",
                    stringr::str_detect(traj, "acc|dec|conc|conv") ~ "quadratic"
                  ) %>% factor(levels = class_lvl),
                  trend = trend %>% factor(levels = c("stable", "decrease",
                                                      "increase"))) %>%
    dplyr::rename(weight_aic_no_change = weight_aic_nch,
                  weight_aic_linear = weight_aic_lin,
                  weight_aic_quadratic = weight_aic_pol,
                  weight_aic_abrupt = weight_aic_chg,
                  nrmse_no_change = nrmse_nch,
                  nrmse_linear = nrmse_lin,
                  nrmse_quadratic = nrmse_pol,
                  nrmse_abrupt = nrmse_chg) %>%
    dplyr::relocate(where(is.factor)) %>%
    dplyr::relocate(c(simu_id, best_model))

  # Make expected trajectories factors:
  if(type=="sim" & apriori){
    best_traj <- best_traj %>%
      dplyr::mutate(expected_traj = expected_traj %>%
                      factor(levels = traj_lvl),
                    expected_class = expected_class %>%
                      factor(levels = class_lvl))
  }

  return(best_traj)

}



#' Run timeseries classification for a given noise combination (for simulations)
#'
#' @param simu_list List of simulation data frames for each noise combination.
#' @param str Character specifying the type of structure to use
#' for classification ("aic","aic_asd").
#' @param asd_thr Numeric threshold for as_detect method.
#' @param run_loo Logical to state whether to perform leave-one-out process.
#' @param edge_lim Numeric, minimal breakpoint distance to start or end dates.
#' @param i Integer to indicate the noise combination (to run in parallel).
#'
#' @return List two objects:
#' - list with results of different trajectory fitting and best model
#' - list of confusion matrices considering trajectories and classes
#'
#' @export

classif_noise_comb <- function(simu_list, str, run_loo, asd_thr, edge_lim, i){

  # Define the noise combination:
  print(paste0("noise combination: ",i,"/",length(simu_list)))
  noise_comb <- names(simu_list)[i]

  # Reshape the data:
  sets <- prep_data(df=simu_list[[i]], thr=NULL, type="sim", apriori=TRUE)

  # Define the breakpoint algorithms used:
  if(str == "aic") abr_mtd <- c("chg")
  if(str == "aic_asd") abr_mtd <- c("chg", "asd")

  # Run the classification:
  trajs <- traj_class(sets, str=str, abr_mtd=abr_mtd, asd_thr=asd_thr,
                      asd_chk=FALSE, type="sim", noise_comb=noise_comb,
                      smooth_signif=TRUE, showplots=FALSE, apriori=TRUE,
                      run_loo=run_loo, two_bkps = TRUE, lowwl=5,
                      highwl="default", mad_thr=3, mad_cst=1.4826,
                      edge_lim=5, congr_brk=5, plot_one_in=10)

  # Make the confusion matrices:
  conf_mat <- make_conf_mat(trajs$best_traj)

  return(list("trajs" = trajs,
              "conf_mat" = conf_mat))
}



#' Run timeseries classification (for empirical data)
#'
#' @param df_list List of data frames with timeseries to classify.
#' @param min_len Minimal timeseries length to consider.
#' @param str Character specifying the type of structure to use
#' for classification ("aic","aic_asd").
#' @param asd_thr Numeric threshold for as_detect method.
#' @param showplots Logical that indicate whether to show plots
#' (slows down computation if true).
#' @param run_loo Logical, whether to perform leave-one-out process.
#' @param two_bkps Logical, if true, looks for more than one breakpoints.
#' @param smooth_signif Logical, should significance of the coefficient be taken
#'  into account in the selection of best model.
#' @param group Column name containing the timeseries identification name (e.g.
#' species or population name).
#' @param time Column name containing timeseries time units.
#' @param variable Column name containing timeseries values.
#' @param outplot Logical, last plot in return value.
#' @param ind_plot Character to return the last plot for a given trajectory
#' (either "nch", "lin", "pol", or "abt").
#' @param dirname Directory name where to save plots.
#' @param save_plot Logical, if plots are shown, whether to save plots
#' (either the four trajectory fits or one of a given fit).
#'
#' @return List of two objects:
#' - Data frame with classification output with one row for each timeseries
#' - List of detailed classification output for each timeseries
#'
#' @export


run_classif_data <- function(df_list, min_len=20, str, asd_thr,
                             showplots=TRUE, run_loo, two_bkps, smooth_signif,
                             group, time, variable, outplot=FALSE,
                             ind_plot=NULL, dirname=NULL, save_plot=TRUE){

  # List of timeseries meeting the length criterion:
  df_list <- df_list[df_list %>% lapply(nrow)>=min_len]

  # Classify for all timeseries of this type:
  traj_ts <- data.frame()
  outlist <- list()
  thr <- NULL

  for (i in 1:length(df_list)){

    set <- df_list[[i]] %>%
      dplyr::rename(scen = .data[[group]],
                    year = .data[[time]]) %>%
      dplyr::select(scen, year, all_of(variable)) %>%
      prep_data(thr=thr, type="data", apriori=FALSE)

    set_length <- nrow(set$ts[[1]])

    if(str == "aic") abr_mtd <- c("chg")
    if(str == "aic_asd") abr_mtd <- c("chg", "asd")

    trajs <- traj_class(sets=set, str=str, abr_mtd=abr_mtd, asd_thr=asd_thr,
                        asd_chk=TRUE, type="data", showplots=showplots,
                        apriori=FALSE, run_loo=run_loo, two_bkps=two_bkps,
                        smooth_signif=smooth_signif, outplot=outplot,
                        ind_plot=ind_plot, lowwl=5, highwl="default",
                        mad_thr=3, edge_lim=5, congr_brk=5, dirname=dirname,
                        save_plot=save_plot)

    # Add linear slope, breakpoints location, magnitude:
    # trajs$best_traj <- trajs$best_traj %>%
    #   dplyr::mutate(slope = trajs$res_detail$res_lin$slope)

    traj_ts <- traj_ts %>% dplyr::bind_rows(trajs$best_traj)
    outlist[[names(df_list)[i]]] <- trajs

    if (i%%10 == 0) print(paste0(i,"/",length(df_list)))
  }

  traj_ts_full <- traj_ts %>%
    dplyr::mutate(species = simu_id %>%
                    sub("_iter01","", .))

  return(list("traj_ts_full"=traj_ts_full, "outlist"=outlist))

}



# IV) Make plots --------------------------------------------------------------


#' Make classification output plots
#'
#' @param sets List of data frame ready for analyses.
#' @param rslt Classification output data frame returned by 'fit_models'
#' function and relative to a given trajectory.
#' @param plot_class Class of trajectory corresponding the rslt.
#' @param best_traj Data frame with the output of the best fitting model,
#' trajectory, class, LOO proportion.
#' @param asd_thr Numeric threshold for as_detect method.
#' @param best_traj_loo Output from LOO process.
#' @param detection_plot Logical to plot the detection plot.
#' @param ind_plot Character to return the last plot for a given trajectory
#' (either "nch", "lin", "pol", "abt", or "best").
#' @param is_date Logical that indicates if the x axis is given as date or
#' numeric.
#'
#' @return List of plots with best fit
#'
#' @export

plot_traj_multi_abt <- function(sets, rslt, plot_class, best_traj,
                                asd_thr=NULL, best_traj_loo=NULL,
                                detection_plot=TRUE, ind_plot=NULL,
                                is_date=FALSE){

  # Initiate plots:
  if (plot_class != "abrupt") plots <-
      vector(mode = "list", length = nrow(rslt))

  if (plot_class == "abrupt") plots <- plot_bis <-
      vector(mode = "list", length = length(rslt$shifts_res))

  for (i in 1:length(sets$ts)){ # For each timeseries

    # Define y-axis name:
    if (sets$ts_type %in% c("TB", "TBbest")) ts_type <- "Biomass"
    # if (sets$ts_type %in% c("TB", "TBbest")) ts_type <-
    #     "State variable (Biomass)"
    else if (sets$ts_type %in% c("TC", "TCbest")) ts_type <- "Catch"
    else if (sets$ts_type %in% c("SProd")) ts_type <- "Surplus production"
    else if (sets$ts_type %in% c("Index", "index")) ts_type <- "Index"
    else if (sets$ts_type %in% c("R")) ts_type <- "Recruitment"
    else ts_type <- sets$ts_type

    # Define x-axis name:
    if (sets$time_type %in% c("year", "Year")) time_type <- "Time unit (year)"
    else time_type <- sets$time_type

    # Plot timeseries and model fit [smooth]:
    if (plot_class != "abrupt"){

      p <- local({
        i <- i
        ggplot(sets$ts[[i]], aes(x = X, y = Y))+
          geom_line()+
          theme_light(base_size = 7)+
          labs(x = time_type, y = ts_type)+
          expand_limits(y=0)+
          stat_function(fun=function(x){
            rslt[i,]$alpha2*as.numeric(x)^2+
              rslt[i,]$alpha1*as.numeric(x)+rslt[i,]$inter}, color="blue")+
          stat_function(fun=function(x){rslt[i,]$alpha2*as.numeric(x)^2+
              rslt[i,]$alpha1*as.numeric(x)+rslt[i,]$inter-rslt[i,]$strd},
              linetype = "dashed", color="blue")+
          stat_function(fun=function(x){rslt[i,]$alpha2*as.numeric(x)^2+
              rslt[i,]$alpha1*as.numeric(x)+rslt[i,]$inter+rslt[i,]$strd},
              linetype = "dashed", color="blue")
      })

      # Plot timeseries and model fit [abrupt]:
    } else {

      table_chg <- rslt$abt_res$chg %>% dplyr::slice(i) %>%
        dplyr::mutate(loc_brk = as.numeric(loc_brk))

      # Plot as_detect detection score:
      if(!is.null(rslt$abt_res$asd)){

        asd_ts <- rslt$shifts_res[[i]]$asd_detect %>%
          tibble::as_tibble() %>%
          tibble::rownames_to_column(var="year") %>%
          dplyr::mutate(year=sets$ts[[i]]$X[as.numeric(year)])

        # Location of shift(s):
        asd_loc <- rslt$abt_res$asd %>% dplyr::slice(i) %>%
          dplyr::pull(loc_brk) %>% strsplit(";") %>% `[[`(1) %>%
          as.numeric() %>% suppressWarnings()

        # Uncertainty around shift(s):
        asd_unc <- data.frame(xmn = rslt$shifts_res[[i]]$asd_run %>%
                                lapply(first) %>% unlist(),
                              xmx = rslt$shifts_res[[i]]$asd_run %>%
                                lapply(last) %>% unlist(),
                              ymn = rep(-Inf, length(asd_loc)),
                              ymx = rep(Inf, length(asd_loc)))

        # Add asdetect shift uncertainty:
        if (!is.na(asd_loc[1])){
          p_asd <- ggplot()+
            geom_rect(aes(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx),
                      fill="pink", alpha=0.5,
                      data=asd_unc)
        } else {
          p_asd <- ggplot()
        }

        # Detection score plot:
        p_asd <-
          p_asd+
          geom_hline(yintercept = c(-asd_thr, asd_thr),
                     col="red", linetype="dotted")+
          geom_line(data=asd_ts, aes(x=year, y=value), inherit.aes=FALSE)+
          theme_light(base_size = 7)+
          labs(y="Detection")+
          theme(plot.title = element_text(hjust = 0.5))+
          expand_limits(y=c(-1,1))+
          ggtitle("as_detect detection score")


        # Add breakpoint on plot:
        if (!is.na(asd_loc[1])) p_asd <- p_asd +
          geom_vline(xintercept = asd_loc, col="red", linetype="dashed")


      } else {
        asd_loc <- NA
      }

      # Add asdetect shift uncertainty:
      if (!is.na(asd_loc[1])){

        p <- local({
          i <- i
          ggplot()+
            geom_rect(aes(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx),
                      fill="pink", alpha=0.5,
                      data=asd_unc)
        })
      } else {
        p <- local({
          i <- i
          ggplot()
        })
      }

      # Plot timeseries with abrupt model fit:
      p <- p+
        geom_line(data=sets$ts[[i]], aes(x = X, y = Y))+
        theme_light(base_size = 7)+
        labs(x = time_type, y = ts_type)+
        expand_limits(y=0)+
        geom_vline(xintercept = table_chg$loc_brk,
                   col="blue", linetype="dashed")+
        geom_line(data = rslt$shifts_res[[i]]$chg_outlist$pred_chg,
                  aes(x=year, y=bp), col = "blue", alpha=0.7)+
        scale_colour_manual(values = rep("red", table_chg$n_brk+1))+
        theme(legend.position = "none") %>%
        suppressWarnings()

      # Plot additional breakdates [chg]:
      if ("loc_aux1_chg" %in% names(best_traj)) {

        loc_aux_chg <- c(best_traj[["loc_aux1_chg"]][i],
                         best_traj[["loc_aux2_chg"]][i])
        p <- p +
          geom_vline(xintercept = loc_aux_chg[!is.na(loc_aux_chg)],
                     col="dodgerblue4", linetype="dotted", alpha=0.5)
      }

      # Plot as_detect breakdates [asd]:
      if (!is.na(asd_loc[1])){
        p <- p +
          geom_vline(xintercept = asd_loc, col="red", linetype="dashed")
      }
    }

    # Define title according to the trajectory class:
    if(plot_class == "no_change"){
      title_part <- paste0("<br>Intercept = ",
                           format(rslt[i,]$inter, digits=2, scientific=TRUE))
    }

    if(plot_class == "linear"){
      if(rslt[i,]$slope_p_value<0.001){
        pval <- " p_val<0.001"
      } else {
        pval <- paste0(" p_val = ", round(rslt[i,]$slope_p_value, digits=3))
      }
      title_part <- paste0("<br>Slope = ",
                           format(rslt[i,]$alpha1, digits=2, scientific=TRUE),
                           pval)
    }

    if(plot_class == "quadratic"){
      if(rslt[i,]$second_order_pvalue<0.001){
        pval <- " p_val<0.001"
      } else {
        pval <- paste0(" p_val = ",
                       round(rslt[i,]$second_order_pvalue, digits=3))
      }
      title_part <- paste0("<br>2nd_order_coeff = ",
                           format(rslt[i,]$alpha2, digits=2, scientific=TRUE),
                           pval)
    }

    if(plot_class == "abrupt"){
      if(is_date==TRUE){

        title_part <- paste0(
          "<br>Breakdate(s): <span style='color:blue'>",
          lubridate::as_date(table_chg$loc_brk), "</span>",
          if("loc_aux1_chg" %in% names(best_traj)) {
            if(!(is.na(best_traj$loc_aux1_chg[i]) &
                 is.na(best_traj$loc_aux2_chg[i]))) {
              paste0(" (<span style='color:dodgerblue4'>",
                     lubridate::as_date(best_traj$loc_aux1_chg[i]),"</span>",
                     ", <span style='color:dodgerblue4'>",
                     lubridate::as_date(best_traj$loc_aux2_chg[i]),"</span>)")
            }
          },
          if(!is.na(asd_loc[1])) {
            paste0("; <span style='color:red'>",
                   paste(lubridate::as_date(asd_loc),
                         collapse=","),"</span>")},
          " Abruptness = ",
          round(table_chg$abruptness, digits=2))

      } else {

        title_part <- paste0(
          "<br>Breakdate(s): <span style='color:blue'>",
          table_chg$loc_brk, "</span>",
          if("loc_aux1_chg" %in% names(best_traj)) {
            if(!(is.na(best_traj$loc_aux1_chg[i]) &
                 is.na(best_traj$loc_aux2_chg[i]))) {
              paste0(" (<span style='color:dodgerblue4'>",
                     best_traj$loc_aux1_chg[i],"</span>",
                     ", <span style='color:dodgerblue4'>",
                     best_traj$loc_aux2_chg[i],"</span>)")
            }
          },
          if(!is.na(asd_loc[1])) {
            paste0("; <span style='color:red'>",
                   paste(asd_loc,collapse=","),"</span>")},
          " Abruptness = ",
          round(table_chg$abruptness, digits=2))

      }
    }

    # Complement title if LOO performed:
    if (!is.null(best_traj_loo)){

      # Add title and LOO points (for smooth fit) [LOO]:
      if (plot_class != "abrupt"){

        p <- p +
          ggtitle(paste0(
            "<b><span style = 'color:#000000;'>",
            stringr::str_to_title(
              sub("_"," ",rslt[i,]$best_class))," ",
            rslt[i,]$trend,"</span></b>",
            "<br>AICc = ", round(rslt[i,]$aic, digits=2),
            "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
            "  LOO = ", round(best_traj[[paste0("loo_", plot_class)]][i],
                              digits=2),
            "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i],
                                digits=2),
            title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10,
                                                      lineheight = 1.1,
                                                      colour="grey25"))+
          geom_point(data=best_traj_loo[[i]] %>%
                       dplyr::left_join(sets$ts[[i]], by="X") %>%
                       dplyr::filter(class==plot_class), aes(x=X, y=Y),
                     col="orange", alpha=0.7)

        # Add title and LOO points (for abrupt fit) [LOO]:
      } else {

        # Add title:
        p <- p +
          ggtitle(paste0(
            "<b><span style = 'color:#000000;'>Abrupt ",
            rslt$abt_res$chg[i,]$trend,"</span></b>",
            "<br>AICc = ", round(table_chg$aic, digits=2),
            "  wAICc = ",
            best_traj[[paste0("weight_aic_", plot_class)]][i],
            "  LOO = ",
            round(best_traj[[paste0("loo_", plot_class)]][i],
                  digits=2),
            "  NRMSE = ",
            round(best_traj[[paste0("nrmse_", plot_class)]][i],
                  digits=2),
            title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10,
                                                      lineheight = 1.1,
                                                      colour="grey25"))+

          # Add histograms for LOO breakdates:
          geom_histogram(data=best_traj_loo[[i]] %>%
                           dplyr::filter(class=="abrupt"),
                         aes(x=loc_brk_chg,
                             y=after_stat(ncount)*
                               max(sets$ts[[i]][[names(sets$ts[[i]])[2]]])/2),
                         fill="blue", alpha=.3, colour=NA, binwidth = 1)

        if(!is.null(rslt$abt_res$asd)){ # histograms for asd breakpoints
          p <- p +
            geom_histogram(data=best_traj_loo[[i]] %>%
                             dplyr::filter(class=="abrupt") %>%
                             dplyr::select(
                               dplyr::matches("^loc_brk_asd_[0-9]+$")) %>%
                             tidyr::pivot_longer(cols=dplyr::everything(),
                                                 names_to = "brk_nb",
                                                 values_to = "loc_brk_asd") %>%
                             tidyr::drop_na(),
                           aes(x=loc_brk_asd,
                               y=after_stat(ncount)*
                                 max(sets$ts[[i]][[names(sets$ts[[i]])[2]]]
                                 )/2),
                           fill="red", alpha=.3, colour=NA, binwidth = 1)
        }

        # Add LOO points:
        p <- p +
          geom_point(data=best_traj_loo[[i]] %>%
                       dplyr::left_join(sets$ts[[i]], by="X") %>%
                       dplyr::filter(class==plot_class), aes(x=X, y=Y),
                     col="orange", alpha=0.7)

      }


      # Complement title if no LOO performed:
    } else {

      # Add title (for smooth fit):
      if (plot_class != "abrupt"){

        p <- p +
          ggtitle(paste0(
            "<b><span style = 'color:#000000;'>",
            stringr::str_to_title(
              sub("_"," ",rslt[i,]$best_class))," ",
            rslt[i,]$trend, "</span></b>",
            "<br>AICc = ", round(rslt[i,]$aic, digits=2),
            "  wAICc = ", best_traj[[paste0("weight_aic_", plot_class)]][i],
            "  NRMSE = ", round(best_traj[[paste0("nrmse_", plot_class)]][i],
                                digits=2),
            title_part))+
          theme(plot.title = ggtext::element_markdown(size=10,
                                                      lineheight = 1.1,
                                                      colour="grey25"))

        # Add title (for abrupt fit):
      } else {

        p <- p +
          ggtitle(paste0(
            "<b><span style = 'color:#000000;'>Abrupt ",
            rslt$abt_res$chg[i,]$trend,"</span></b>",
            "<br>AICc = ", round(table_chg$aic, digits=2),
            "  wAICc = ", best_traj[[paste0("weight_aic_",
                                            plot_class)]][i],
            "  NRMSE = ", round(best_traj[[paste0("nrmse_",
                                                  plot_class)]][i],
                                digits=2),
            title_part)) +
          theme(plot.title = ggtext::element_markdown(size=10,
                                                      lineheight = 1.1,
                                                      colour="grey25"))

      }
    }


    # Color the background according to the best trajectory:
    if (best_traj$class[i] == plot_class){

      bkg_col <- "lightskyblue"

    } else { bkg_col <- "grey80" }

    p <- p + theme(plot.background = element_rect(fill = bkg_col))


    # Add radar plots:
    if(!is.null(ind_plot)){

      # Manage the size of the radar plots

      #   # Weight/LOO radar plot:
      #   wrad <- best_traj %>%
      #     dplyr::select(simu_id | dplyr::contains("weight")) %>%
      #     `colnames<-`(c("simu_id","nch","lin","qdr","abt")) %>%
      #     dplyr::slice(i)
      #
      #   # AICc weight radar plot:
      #   wrad_plot <-
      #     ggradar::ggradar(wrad, axis.label.size = 1.5,
      #                      grid.label.size = 0, group.point.size = 1,
      #                      group.line.width = .2,
      #                      group.colours = "red",
      #                      background.circle.transparency=0, centre.y=0,
      #                      gridline.mid.colour="grey20",
      #                      gridline.min.colour="grey20",
      #                      gridline.max.colour="grey20",
      #                      axis.line.colour="grey20", grid.line.width=0.25)+
      #     theme(
      #       plot.background = element_blank(),
      #       panel.background = element_blank(),
      #       plot.caption = element_text(hjust = 0.5, vjust = 0, size = 6))+
      #     labs(caption = "wAICc")
      #
      #   p <- cowplot::ggdraw(p)+
      #     # cowplot::draw_plot(wrad_plot, x = 0.75, y = .86,
      #     #                    width = .15, height = .15)
      #     cowplot::draw_plot(wrad_plot, x = 0.75, y = .8,
      #                               width = .3, height = .3)
      #
      #   # LOO radar plot:
      #   if (!is.null(best_traj_loo)){
      #
      #     loorad <- best_traj %>%
      #       dplyr::select(simu_id | dplyr::contains("loo")) %>%
      #       `colnames<-`(c("simu_id","nch","lin","qdr","abt")) %>%
      #       dplyr::slice(i)
      #
      #     loorad_plot <-
      #       ggradar::ggradar(loorad, axis.label.size = 1.5,
      #                        grid.label.size = 0, group.point.size = 1,
      #                        group.line.width = .2,
      #                        group.colours = "red",
      #                        background.circle.transparency=0, centre.y=0,
      #                        gridline.mid.colour="grey20",
      #                        gridline.min.colour="grey20",
      #                        gridline.max.colour="grey20",
      #                        axis.line.colour="grey20", grid.line.width=0.25)+
      #       theme(
      #         plot.background = element_blank(),
      #         panel.background = element_blank(),
      #         plot.caption = element_text(hjust = 0.5, vjust = 0, size = 6))+
      #       labs(caption = "LOO")
      #
      #     # Add to plot:
      #     p <- cowplot::ggdraw(p) +
      #       cowplot::draw_plot(loorad_plot, x = 0.85, y = .8,
      #                          width = .3, height = .3)
      #
      # }


      # Add as_detect plot:
      if (detection_plot==TRUE &
          plot_class == "abrupt" & !is.null(rslt$abt_res$asd)){

        p_bis <- cowplot::plot_grid(p,
                                    p_asd +
                                      theme(plot.background =
                                              element_rect(fill = bkg_col)),
                                    align = 'hv', ncol=1, rel_heights = c(3, 2))

        # Add to secondary plot:
        p_bis <- cowplot::ggdraw() +
          cowplot::draw_plot(p_bis, x = 0, y = 0, width = 1, height = 1)

        plot_bis[[i]] <- p_bis

      }
    }

    plots[[i]] <- p

  }

  # Either return the trajectory plot only, or
  # also with detection plot if abrupt:
  if (plot_class == "abrupt" & !is.null(rslt$abt_res$asd) & detection_plot){

    return(list("plots"=plots, "plot_bis"=plot_bis))

  } else {

    return(list("plots"=plots))
  }
}



#' Plot simulations in a simple way
#'
#' @param lib_ts Data frame with simulated timeseries.
#' @param var Column name containing timeseries values.
#' @param tbcolor Line color.
#' @param alpha Line transparency.
#' @param xname Character to name x axis.
#' @param yname Character to name y axis.
#' @param ylim Maximal value to constrain y axis.
#'
#' @return Minimal plot with TB curve.
#'
#' @export

plot_simu_simple <- function (lib_ts, var, tbcolor="black", alpha=1,
                              xname="", yname="", ylim=12.5){

  TB_plot <- ggplot(lib_ts, aes(x=year, y=get(var),
                                color=as.factor(iter), group=iter))+
    geom_line(col=tbcolor, alpha=alpha)+
    scale_x_continuous(name = xname)+
    scale_y_continuous(name = yname)+
    expand_limits(y = c(0, ylim))+
    theme_light()+
    theme(legend.position = "none")

  return(TB_plot)
}



# V) Post hoc analyses -------------------------------------------------------

#' Compute a confusion matrix
#'
#' @param best_traj Data frame with the classification output from simulated
#' timeseries, with best model and expected trajectory and class.
#'
#' @return List of confusion matrices considering either the nine trajectories
#'  or the four trajectory classes.
#'
#' @export

make_conf_mat <- function(best_traj){

  conf_mat_traj <- caret::confusionMatrix(data=best_traj$traj,
                                          reference=best_traj$expected_traj)
  conf_mat_class <- caret::confusionMatrix(data=best_traj$class,
                                           reference=best_traj$expected_class)

  return(list("conf_mat_traj" = conf_mat_traj,
              "conf_mat_class"=conf_mat_class))
}



#' Combine classification output from different length
#'
#' @param outlist_all List of "outlists" corresponding to different
#' timeseries lengths.
#'
#' @return Data frame combining classification output info from timeseries
#' of all lengths.
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
#'  corresponding to a given timeseries lengths or combine_length function.
#' @param save Logical, should the figure be saved automatically?
#' @param dirname Directory name where to save plots.
#' @param filename Character to name the file saved.
#'
#' @return Confusion matrix as heatmap.
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



#' Pool confusion matrices into one (by length)
#'
#' @param path Path to the .rds file output from traj function.
#' @param title Character to display as title.
#' @param show_legend Logical, whether to show the legend.
#'
#' @return Confusion matrix as pheatmap object.
#'
#' @export

conf_mat_pool <- function(path, title, show_legend=FALSE) {

  # Load matrices:
  outlist <- readRDS(path)
  pool_mat <- outlist$mat_list[[1]]$conf_mat_class$table

  # Pool all matrices (from different noise levels):
  for(i in 2:length(outlist$mat_list)){
    pool_mat <- pool_mat + outlist$mat_list[[i]]$conf_mat_class$table
  }
  colnames(pool_mat) <- rownames(pool_mat) <- sub("_"," ",colnames(pool_mat))

  if(show_legend){
    ftsz <- 15 ; ftsz_nb <- 20
  } else { ftsz <- 150 ; ftsz_nb <- 200 }

  # Make heatmap with proportions:
  byclass <- sum(pool_mat[,1])
  p <- pheatmap::pheatmap(pool_mat/byclass,
                          color = grDevices::colorRampPalette(
                            rev(RColorBrewer::brewer.pal(
                              n = 7, name ="RdYlBu")))(100),
                          breaks = seq(0, 1, .01),
                          cluster_rows = FALSE, cluster_cols = FALSE,
                          legend = show_legend, display_numbers = TRUE,
                          number_format = "%.2f", fontsize = ftsz,
                          number_color="black", fontsize_number = ftsz_nb,
                          main=title)
  return(p)
}



#' Make confusion matrices by timeseries length and noise combination
#'
#' @param outlist List of "outlist" from the traj_class function
#'  corresponding to  timeseries lengths or combine_length function.
#' @param save Logical, should the figure be saved automatically?
#' @param dirname Directory name where to save plots.
#' @param filename Character to name the file saved.
#'
#' @return Patchwork of cowplots.
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
#' @param outlist Data frame combining classification output info from
#' all timeseries of all lengths.
#'
#' @return Data frame with wAICc, LOO, and NRMSE in a long format.
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
