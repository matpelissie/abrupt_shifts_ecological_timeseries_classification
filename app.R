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
  ts_type <- names(df)[2]
  X <- df[,1, drop=TRUE]
  Y <- df[,2, drop=TRUE]

  # Format of the csv:
  # 2 columns:
  # 1st column with time (e.g. year)
  # 2nd column with the values of the variable (e.g. biomass)

  sets <- list("ts" = data.frame(X=X, Y=Y))

  return(list("ts"=sets, "ts_type"=ts_type))
}


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


      shiny::sliderInput(inputId = "asd_thr",
                         label = "asdetect threshold", min = 0,
                         max = 1, value = 0.15, sep="", width='30%'),


      h5(shiny::HTML("wAICc: AICc weight <br/>
         LOO: Leave-one-out score  <br/>
         NRMSE: Normalized root mean square error")),


      # Sources
      p(shiny::h3(tags$strong("Sources"))),

      p("..."),

    ), # end sidebarPanel

    # Plot map
    shiny::mainPanel(

      shiny::plotOutput("plot"),

      DT::dataTableOutput("tbl")

      # To add as input
      # Expandable panel with asd_thr, asd or not,

      # To add as output
      # When breakpoint: shape before/after


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

    classif <- traj_class(sets=set, str="aic_asd", abr_mtd=c("chg", "asd"),
                          asd_thr=input$asd_thr, asd_chk=FALSE, type="data",
                          showplots=TRUE, apriori=FALSE, run_loo=input$loo,
                          save_plot=FALSE, two_bkps=TRUE, smooth_signif=TRUE,
                          outplot=TRUE, ind_plot=display, dirname=NULL)
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
              dplyr::select(dplyr::contains("nrmse") &
                              !"nrmse_asd") %>%
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
              dplyr::select(dplyr::contains("nrmse") &
                              !"nrmse_asd") %>%
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


}


# Launch app --------------------------------------------------------------

shiny::shinyApp(ui=ui, server=server)
