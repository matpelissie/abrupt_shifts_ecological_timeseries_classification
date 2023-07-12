###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###
#
# 18/01/2022 mathieu.pelissie@ens-lyon.fr
#
# Functions to explore the RAMLDB
#
# functions_RAMLBD.R
#
###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###-###


# NB: For all the following functions the RAMLBD should be loaded in the environment

# II) Individual stock plots -----------------------------------------------------

## II-1) Time series

#' List the types of time series available for a given stock
#'
#' @param id short stock ID character
#'
#' @return a list of types of time series available for this stock
#' @export

ts_available <- function(id){

  tsmetr_light <- tsmetrics %>%
    dplyr::arrange(tsshort) %>%
    tidyr::separate(tsshort, c("tsshort.1","tsshort.2"),"-") %>%
    dplyr::select(c(2,4)) %>%
    dplyr::rename(ts_type = tsshort.1) %>%
    tidyr::separate(tslong, c("tslong.1","tslong.2"),"[(]") %>%
    dplyr::select(c(1,2)) %>%
    dplyr::rename(tslong = tslong.1) %>%
    dplyr::group_by(ts_type) %>%
    dplyr::slice(1) %>%
    dplyr::group_by(tslong = stringi::stri_trans_totitle(tslong)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    suppressWarnings() %>%
    dplyr::bind_rows(c(ts_type="SProd", tslong="Surplus Production"))

  types <- timeseries_years_views %>%
    dplyr::filter(stockid == id)

  ts_type <- names(types[,!is.na(types)])[-c(1,2)]
  ts_type_df <- data.frame(ts_type) %>%
    dplyr::left_join(tsmetr_light, by="ts_type") %>%
    dplyr::rename(meaning = tslong)

  return(ts_type_df)
}


#' #' Get optimal number of breakpoints
#' #' from https://www.marinedatascience.co/blog/2019/09/28/comparison-of-change-point-detection-methods/
#' #'
#' #' @param id short stock ID character
#' #'
#' #' @return the optimal number of breakpoints in the time series based on the BIC
#' #' @export
#'
#' opt_bpts <- function(id) {
#'   #id = bpts_sum$RSS["BIC",]
#'   n <- length(id)
#'   lowest <- vector("logical", length = n-1)
#'   lowest[1] <- FALSE
#'   for (i in 2:n) {
#'     lowest[i] <- id[i] < id[i-1] & id[i] < id[i+1]
#'   }
#'   out <- as.integer(names(id)[lowest])
#'   # if (is.na(out)) out <- 5
#'
#'   return(out)
#' }


#' Extract a stock RAMLBD time series
#'
#' @param id short stock ID character
#' @param ts_type type of time series (TBbest, ERbest...) character
#'
#' @return a data frame of the selected time series
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


#' Extract a stock RAMLBD time series cutting at start and end specitied
#'
#' @param id short stock ID character
#' @param ts_type type of time series (TBbest, ERbest...) character
#' @param min_y start year
#' @param max_y end year
#'
#' @return a data frame of the selected time series
#' @export

extract_RAM_lim <- function(id, ts_type, min_y, max_y){

  ts <- timeseries_values_views %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(year, tidyselect::all_of(ts_type)) %>%
    dplyr::filter(year>=min_y & year<=max_y) %>%
    dplyr::mutate(scen = paste(ts_type, id, sep="_")) %>%
    dplyr::relocate(scen) %>%
    na.omit()

  return(ts)

}


#' Extract a stock RAMLBD time series and plot it
#'
#' @param id short stock ID character
#' @param ts_type type of time series (TBbest, ERbest...) character
#'
#' @return a list of three objects,
#' first a time series plot, second the data frame used to plot it,
#' then the output of the breakpoints analysis
#' @export

# account for incomplete time series (zoo::zoo?)

shifts_RAM <- function(id, ts_type){

  ts <- extract_RAM(id, ts_type)

  incomplete <- FALSE
  if(nrow(ts) != tail(ts$year,1) - head(ts$year,1) + 1){
    warning("Time series might be incomplete")
    incomplete <- TRUE
  }

  brk_out <- shifts(ts, abr_mtd=c("asd","chg"), asd_thr=0.15, asd_chk=TRUE)

  return(brk_out)
}


#' Extract a stock RAMLBD time series and plot it
#'
#' @param id short stock ID character
#' @param ts_type type of time series (TBbest, ERbest...) character
#'
#' @return a list of three objects,
#' first a time series plot, second the data frame used to plot it
#' @export

plot_RAM_bp <- function(id, ts_type){

  ts <- extract_RAM(id, ts_type)

  stock_ts <- ts %>%
    pull(ts_type) %>%
    ts(start=ts$year[1], end=tail(ts$year,1))

  metadt <- metadata %>%
    dplyr::filter(stockid == id) %>%
    dplyr::select(c(2,3,5,6,10))

  taxo <- taxonomy %>%
    dplyr::filter(scientificname == metadt[,3]) %>%
    dplyr::select(c(2,4,6,7,15))

  ts_info <- left_join(metadt, taxo, by = "scientificname")

  stocklong <- metadt$stocklong
  sciname <- ts_info$scientificname
  famname <- ts_info$family
  reg <- ts_info$region

  shifts_res <- shifts_RAM(id, ts_type)

  p <- plot_bp(ts, shifts_res = shifts_res, abr_mtd=c("asd","chg"))

  title <- paste0(stocklong, " [",id,"] – (",
                  ts$year[1],"–",tail(ts$year, 1),")\n",
                  sciname," – ",famname," – ",reg)

  pop <- p[[1]]+
    ggtitle(title)

  bp_fit <- p[[2]]+
    ggtitle(title)

  # bp_slop <- p[[3]]+
  #   ggtitle(title)

  return(list("pop"=pop, "bp_fit"=bp_fit, "ts"=ts))
}
