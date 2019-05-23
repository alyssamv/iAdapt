#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "iAdapt")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing iAdapt", call. = FALSE)
  }
   shiny::runApp(appDir, display.mode = "normal")
}



