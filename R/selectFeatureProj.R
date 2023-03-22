#' @title Select feature when projection strategy is employed for the
#'  case where features are shared with multiple omics profiles
#'
#' @param HOSVD HOSVD 
#' @param Multi list of omics profiles, row: sample, column: feature
#' @param cond list of conditions for individual omics profiles
#' @param de initial value for optimization of standard deviation
#' @param p0 Threshold P-value
#' @param breaks The number of bins of histogram of P-values
#' @param input_all The number of selected feature. if null, interactive mode
#' is activated 
#'
#' @return list composed of logical vector that represent which features are selected and p-values
#' @export
#'
#' @examples
#' library(TDbasedUFE)
#' Multi <- list(matrix(runif(1000),10),matrix(runif(1000),10),
#' matrix(runif(1000),10),matrix(runif(1000),10))
#' Z <- prepareTensorfromList(Multi,10L)
#' Z <- aperm(Z,c(2,1,3))
#' Z <- PrepareSummarizedExperimentTensor(feature =as.character(1:10),
#'                                       sample=array("",1),value=Z)
#' HOSVD <- computeHosvd(Z)
#' cond <- rep(list(rep(1:2,each=5)),4)
#' index <- selectFeatureProj(HOSVD,Multi,cond,de=0.1,input_all=2)
selectFeatureProj <-
  function(HOSVD, Multi, cond, de = 1e-4, p0 = 0.01, breaks = as.integer(100),
           input_all = NULL) {
    # Augument check
    stopifnot("`HOSVD` must be a list." = is.list(HOSVD))
    stopifnot("`Multi` must be a list." = is.list(Multi))
    stopifnot("`de` must be a numeric." = is.numeric(de))
    stopifnot("`p0` must be a numeric." = is.numeric(p0))
    stopifnot("`breaks` must be a integer." = is.integer(breaks))
    stopifnot("`input_all` must be a vector." = is.vector(input_all) |
      is.null(input_all))
    #
    interact <- FALSE
    Multi_list <- lapply(
      Multi,
      function(x) {
        data.matrix(x) %*% data.matrix(HOSVD$U[[1]])
      }
    )
    if (is.null(input_all)) {
      interact <- TRUE
      j <- 1
      ui <- fluidPage(
        sidebarLayout(
          sidebarPanel(
            h3("Select one with prefarable dependence"),
            actionButton(inputId = "action", label = "Next"),
            actionButton(inputId = "prev", label = "Prev"),
            actionButton(inputId = "select", label = "Select")
          ),
          mainPanel(
            plotOutput("plot")
          )
        )
      )
      server <- function(input, output) {
        observeEvent(input$action, {
          if (j < dim(HOSVD$U[[1]])[2]) j <<- j + 1
        })
        observeEvent(input$prev, {
          if (j != 1) {
            j <<- j - 1
          }
        })
        observeEvent(input$select, {
          input_all <<- j
          stopApp()
        })
        output$plot <- renderPlot({
          input$action
          input$prev
          par(mfrow = c(length(cond), 1))
          par(mai = c(0.3, 0.2, 0.2, 0.2))
          for (i in seq_along(cond))
          {
            boxplot(Multi_list[[i]][, j] ~ cond[[i]],
              main = paste(j, i, sep = "-")
            )
            abline(0, 0, col = 2, lty = 2)
          }
          par(mfrow = c(1, 1))
        })
      }
      app <- shinyApp(ui, server)
      runApp(app)
      input_all <- j
    } else {
      par(mfrow = c(length(cond), 1))
      par(mai = c(0.3, 0.2, 0.2, 0.2))
      for (i in seq_along(cond))
      {
        boxplot(Multi_list[[i]][, input_all] ~ cond[[i]],
          main = paste(input_all, i, sep = "-")
        )
        abline(0, 0, col = 2, lty = 2)
      }
      par(mfrow = c(1, 1))
    }
    th <- function(sd, breaks, p0) {
      P2 <- pchisq((u / sd)^2, 1, lower.tail = FALSE)
      hc <- hist(1 - P2, breaks = breaks, plot = FALSE)
      return(sd(hc$count[seq_len(sum(hc$breaks
      < 1 - min(P2[p.adjust(P2, "BH") > p0])))]))
    }
    u <- HOSVD$U[[1]][, input_all]
    sd <- optim(de, function(x) {
      th(x, breaks, p0)
    },
    control = list(warn.1d.NelderMead = FALSE)
    )$par
    sd1 <- seq(0.1 * sd, 2 * sd, by = 0.1 * sd)
    th0 <- apply(matrix(sd1, ncol = 1), 1, function(x) {
      th(x, breaks, p0)
    })
    P2 <- pchisq((u / sd)^2, 1, lower.tail = FALSE)
    ui <- fluidPage(
      sidebarLayout(
        sidebarPanel(
          actionButton(inputId = "action", label = "Next")
        ),
        mainPanel(
          plotOutput("plot")
        )
      )
    )
    server <- function(input, output) {
      observeEvent(input$action, {
        stopApp()
      })
      output$plot <- renderPlot({
        input$action
        par(mfrow = c(1, 2))
        plot(sd1, th0, type = "o")
        arrows(sd, max(th0), sd, min(th0), col = 2)
        hist(1 - P2, breaks = breaks)
        par(mfrow = c(1, 1))
      })
    }
    app <- shinyApp(ui, server)
    if (interact) {
      runApp(app)
    } else {
      par(mfrow = c(1, 2))
      plot(sd1, th0, type = "o")
      arrows(sd, max(th0), sd, min(th0), col = 2)
      hist(1 - P2, breaks = breaks)
      par(mfrow = c(1, 1))
    }
    index <- p.adjust(P2, "BH") < p0
    index_all <- list(index = index, p.value = P2)
    return(index_all)
  }


