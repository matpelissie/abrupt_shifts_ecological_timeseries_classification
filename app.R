library("shiny")
library("leaflet")
library("tidyverse")
source("analyses/00_packages.R")
source("R/functions_trajclass.R")


# Input data --------------------------------------------------------------

color_gradient_bounded <- function(dt, column_name,
                                   gradient_colors =
                                     c("#6666FF", "#DDDDDD", "#FF6666")) {
  col_func <- grDevices::colorRampPalette(gradient_colors)
  vals <- sort(unique(dt$x$data[[column_name]]), decreasing = TRUE)
  dt %>%
    DT::formatStyle(column_name,
                    backgroundColor = DT::styleEqual(
                      vals,
                      col_func(100)[round(round(vals, 2)*99+1,0)]))
}

prep_data_any <- function(df){

  # For any empirical data (with 2 columns)
  time_time <- names(df)[1]
  ts_type <- names(df)[2]
  X <- df[,1, drop=TRUE]
  Y <- df[,2, drop=TRUE]

  # Format of the csv:
  # 2 columns:
  # 1st column with time (e.g. year)
  # 2nd column with the values of the variable (e.g. biomass)

  sets <- list("ts" = data.frame(X=X, Y=Y))

  return(list("ts"=sets, "ts_type"=ts_type, "time_type"=time_time))
}

quality_df <- readr::read_csv("data/04_app/quality_control_df.csv")


# user interface ----------------------------------------------------------

ui <- shiny::fluidPage(
  # title = "Timeseries classifier",
  shinyjs::useShinyjs(),
  shiny::fluidPage(  # Application title
    shiny::h1("Timeseries classifier"),
    shiny::h2("Load single timeseries")
  ),

  shiny::sidebarLayout(

    shiny::sidebarPanel(

      # Selector for input file:
      shiny::fileInput("file", shiny::HTML("Select a .csv file<br/>
                       1st column: time<br/>
                       2nd column: variable of interest")),


      shiny::checkboxInput(inputId = "ind_plot",
                           label = "Show only best model fit",
                           value = FALSE),

      shiny::checkboxInput(inputId = "loo",
                           label = "Compute Leave-one-out indices",
                           value = FALSE),

      shinyBS::bsCollapse(id = "cntrlC1", open = "Panel 2",
                          shinyBS::bsCollapsePanel(title = "Adjust parameters",

                                                   shiny::checkboxInput(
                                                     inputId = "use_asd",
                                                     label = "Use asdetect method",
                                                     value = TRUE),

                                                   shiny::sliderInput(
                                                     inputId = "asd_thr",
                                                     label = "asdetect detection threshold",
                                                     min = 0,  max = 1,
                                                     value = 0.15, sep="",
                                                     width='50%'),

                                                   shiny::sliderInput(
                                                     inputId = "mad_thr",
                                                     label = "asdetect MAD threshold",
                                                     min = 0,  max = 5,
                                                     value = 3, sep="",
                                                     width='50%'),

                                                   shiny::sliderInput(
                                                     inputId = "edge_lim",
                                                     label = "Limit to edges",
                                                     min = 0,  max = 10,
                                                     value = 5, sep="",
                                                     width='50%'),

                                                   shiny::sliderInput(
                                                     inputId = "congr_brk",
                                                     label = "Congruence limit between brekpoints",
                                                     min = 0,  max = 10,
                                                     value = 5, sep="",
                                                     width='50%'),

                                                   htmlOutput("dateSelector", inline = TRUE),
                                                   style = "info"
                          )
      ),

      # # Sources
      p(shiny::h3(tags$strong("Sources"))),

      p(h6(shiny::HTML('Pélissié M., Devictor V. & Dakos V. 2024. A systematic approach for detecting abrupt shifts in ecological timeseries.',
           '<i>Biological Conservation</i> 290: 110429.',
           '<a href="https://www.sciencedirect.com/science/article/pii/S000632072300530X">https://doi.org/10.1016/j.biocon.2023.110429</a>.'))),

      p(h6("The code to run classification on several timeries is available",
           a(href="https://github.com/matpelissie/abrupt_shifts_ecological_timeseries_classification", "here")))

    ), # end sidebarPanel

    # Plot map
    shiny::mainPanel(

      h2("Model fitting"),

      shiny::plotOutput("plot"),

      h6(shiny::HTML("wAICc (AICc weight) expresses how much the best model performs in terms of AIC compared to the others. <br/>
         LOO (Leave-one-out score) indicates how robust is the model choice to the removal of individual data points. <br/>
         NRMSE (Normalized root mean square error) gives the ratio between the variation not explained by the model and the overall variation.")),

      h2("Quality of fit"),

      shiny::plotOutput("qul"),

      h6("The plots above show the quality scores of the classification compared with the distribution of those from simulated timeseries correctly (in blue) and incorrectly (in red) classified."),

      h2("Quality of fit - table"),

      DT::dataTableOutput("tbl")

      # To add as input
      # Possibility to show asdetect detection timeseries
      # Possibility to change the moving window on which the classification is done

      # To add as output
      # When breakpoint: shape before/after
      # Relative change size

      # Warning when fit is not good


    ) # end mainPanel
  ) # end sidebarLayout
) # end fluidPage


# server function ---------------------------------------------------------


server <- function(input, output, session) {

  # Import input file:
  input_file <- shiny::reactive({
    if (is.null(input$file)) {
      return("")
    }
    readr::read_csv(file = input$file$datapath)
  })


  # Trajectory classification:
  traj_class_out <- shiny::reactive({

    # render only if there is data available
    shiny::req(input_file())

    data <- input_file()

    set <- prep_data_any(df=data)

    if (input$ind_plot==TRUE){
      display <- "best"
    } else {
      display <- NULL
    }

    if (input$use_asd==TRUE){
      abr_mtd_usd <- c("chg", "asd")
      str <- "aic_asd"
    } else {
      abr_mtd_usd <- c("chg")
      str <- "aic"
    }

    classif <- traj_class(sets=set, str=str, abr_mtd=abr_mtd_usd,
                          asd_thr=input$asd_thr, asd_chk=TRUE, type="data",
                          showplots=TRUE, apriori=FALSE, run_loo=input$loo,
                          edge_lim=input$edge_lim, congr_brk=input$congr_brk,
                          mad_thr=input$mad_thr,
                          save_plot=FALSE, two_bkps=TRUE, smooth_signif=TRUE,
                          outplot=TRUE, ind_plot=display, dirname=NULL)

    if(input$use_asd==TRUE) classif$best_traj <-
      classif$best_traj %>%
      dplyr::select(-nrmse_asd)


    return(classif)

  })

  # Output plot :
  output$plot <- shiny::renderPlot({

    # render only if there is data available
    shiny::req(input_file())

    if (input$ind_plot==TRUE){

      traj_class_out()$class_plot$plots[[1]]

    } else {

      traj_class_out()$class_plot
    }

  })


  # Output table:
  tbl <- shiny::reactive({

    shiny::req(traj_class_out)

    bestclass <- traj_class_out()$best_traj$class

    if(input$loo == TRUE){

      tb <- traj_class_out()$best_traj %>%
        dplyr::select(dplyr::contains("weight")) %>%
        t() %>% cbind(
          traj_class_out()$best_traj %>%
            dplyr::select(dplyr::contains("loo")) %>%
            t()
        ) %>% cbind(
          traj_class_out()$best_traj %>%
            dplyr::select(dplyr::contains("nrmse"))%>%
            t()
        ) %>%
        `colnames<-`(c("wAICc", "LOO", "NRMSE")) %>%
        `rownames<-`(c("no change", "linear", "quadratic", "abrupt")) %>%
        as_tibble(rownames = "class") %>%
        dplyr::mutate(best_class = ifelse(class==bestclass, 1, 0)) %>%
        dplyr::mutate(wAICc = round(wAICc, digits=3),
                      LOO = round(LOO, digits=3),
                      NRMSE = round(NRMSE, digits=3)) %>%
        DT::datatable() %>%
        DT::formatStyle(
          "best_class",
          target = 'row',
          backgroundColor = DT::styleEqual(c(1, 0), c('lightgreen', 'white'))) %>%
        color_gradient_bounded(column_name="wAICc", c("grey95","#6666FF")) %>%
        color_gradient_bounded(column_name="LOO", c("grey95","#6666FF")) %>%
        color_gradient_bounded(column_name="NRMSE", c("#6666FF","grey95"))

    } else {

      tb <- traj_class_out()$best_traj %>%
        dplyr::select(dplyr::contains("weight")) %>%
        t() %>% cbind(
          traj_class_out()$best_traj %>%
            dplyr::select(dplyr::contains("nrmse")) %>%
            t()
        ) %>%
        `colnames<-`(c("wAICc", "NRMSE")) %>%
        `rownames<-`(c("no change", "linear", "quadratic", "abrupt")) %>%
        as_tibble(rownames = "class") %>%
        dplyr::mutate(best_class = ifelse(class==bestclass, 1, 0)) %>%
        dplyr::mutate(wAICc = round(wAICc, digits=3),
                      NRMSE = round(NRMSE, digits=3)) %>%
        DT::datatable() %>%
        DT::formatStyle(
          "best_class",
          target = 'row',
          backgroundColor = DT::styleEqual(c(1, 0), c('lightgreen', 'white'))) %>%
        color_gradient_bounded(column_name="wAICc", c("grey95","#6666FF")) %>%
        color_gradient_bounded(column_name="NRMSE", c("#6666FF","grey95"))

    }

    return(tb)

  })

  output$tbl <- DT::renderDataTable({ tbl() })


  qul <- shiny::reactive({

    shiny::req(traj_class_out)

    bestclass <- traj_class_out()$best_traj$class

    if(input$loo == TRUE){

      quality_scores <- traj_class_out()$best_traj %>%
        dplyr::select(paste0(c("loo_", "nrmse_","weight_aic_"), bestclass)) %>%
        `colnames<-`(c("LOO","NRMSE", "wAICc")) %>%
        tidyr::pivot_longer(cols = dplyr::everything(),
                            names_to = "metric", values_to = "vals")

      quality_inds <- ggplot(quality_df %>%
                               dplyr::filter(class_tested==class,
                                             class==bestclass) %>%
                               dplyr::mutate(correct = ifelse(class==expected_class,
                                                              TRUE, FALSE)),
                             aes(x = value, fill = correct)) +
        geom_density(alpha = 0.4, position="identity", adjust=3) +
        theme_classic() +
        expand_limits(y=c(0,1))+
        facet_wrap(nrow=1, facets=vars(metric), scales="free_y")+
        theme(legend.position = "none")+
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
        scale_colour_manual(values = c("TRUE" = "#00bfc4", "FALSE" = "#f8766d")) +
        geom_vline(data = quality_scores, aes(xintercept = vals),
                   colour="red", linewidth=1)

    } else {

      quality_scores <- traj_class_out()$best_traj %>%
        dplyr::select(paste0(c("nrmse_","weight_aic_"), bestclass)) %>%
        `colnames<-`(c("NRMSE", "wAICc")) %>%
        tidyr::pivot_longer(cols = dplyr::everything(),
                            names_to = "metric", values_to = "vals")

      quality_inds <- ggplot(quality_df %>%
                               dplyr::filter(class_tested==class,
                                             class==bestclass,
                                             metric != "LOO") %>%
                               dplyr::mutate(correct = ifelse(class==expected_class,
                                                              TRUE, FALSE)),
                             aes(x = value, fill = correct)) +
        geom_density(alpha = 0.4, position="identity", adjust=3) +
        theme_classic() +
        expand_limits(y=c(0,1))+
        facet_wrap(nrow=1, facets=vars(metric), scales="free_y")+
        theme(legend.position = "none")+
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed")) +
        scale_colour_manual(values = c("TRUE" = "#00bfc4", "FALSE" = "#f8766d")) +
        geom_vline(data = quality_scores, aes(xintercept = vals),
                   colour="red", linewidth=1)

    }

    return(quality_inds)

  })

  output$qul <- shiny::renderPlot({ qul() })

} # end server


# Launch app --------------------------------------------------------------

shiny::shinyApp(ui=ui, server=server)
