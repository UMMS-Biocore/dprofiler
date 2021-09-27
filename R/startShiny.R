#' startDprofiler
#'
#' Starts the Dprofiler to be able to run interactively.
#'
#' @note \code{startDprofiler}
#' @return the app
#'
#' @examples
#'     startDprofiler()
#'
#' @export
#'
startDprofiler <- function(){
  if (interactive()) {
    #the upload file size limit is 30MB
    options( shiny.maxRequestSize = 30 * 1024 ^ 2, warn = -1,
             shiny.sanitize.errors = TRUE)
    # addResourcePath(prefix = "demo", directoryPath =
    #                   system.file("extdata", "demo", 
    #                               package = "dprofiler"))
    # addResourcePath(prefix = "www", directoryPath =
    #                   system.file("extdata", "www", 
    #                               package = "dprofiler"))
    environment(dprofilerServer) <- environment()
    
    app <- shinyApp( ui = shinyUI(dprofilerUI),
                     server = shinyServer(dprofilerServer))
    runApp(app)
  }
}

