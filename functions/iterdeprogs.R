#' getIterDEResultsUI
#' Creates a panel to visualize DE results
#'
#' @param id, namespace id
#' @return panel
#' @examples
#'     x <- getIterDEResultsUI("batcheffect")
#'
#' @export
#'
getIterDEResultsUI<- function (id) {
  ns <- NS(id)
  list(
    fluidRow(
      shinydashboard::box(title = "Results",
                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                          fluidRow(
                            column(12,
                                   uiOutput(ns("IterDEResults"))
                            ),
                            actionButtonDE("goMain", "Go to Main Plots", styleclass = "primary")
                          )
      ),
      shinydashboard::box(title = "Membership Scores",
                          solidHeader = T, status = "info",  width = 12, collapsible = TRUE,
                          fluidRow(
                            column(12,
                                   uiOutput(ns("MembershipScores"))
                            )
                          )
      )
    )
  )
}