#' getTabUpdateJS
#' 
#' prepmenu tab and discovery menu tab updates, adapted from debrowser::getTabUpdateJS()
#'
#' @return the JS for tab updates
#'
#' @examples
#'     x<- getTabUpdateJS()
#'
#' @export
getTabUpdateJS<-function(){
  tags$script(HTML( "
                      $(function() {
                      //hide buttons on entrance
                      $('.sidebar-menu > ').css('display', 'none');
                      $('.sidebar-menu > :nth-child(1)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(2)').css('display', 'inline');
                      $(document).on('click', '#load-Filter', function () {
                      $('.sidebar-menu > :nth-child(2)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(3)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(4)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#lcf-goDE', function () {
                      $('.sidebar-menu > :nth-child(4)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#batcheffect-goDE', function () {
                      $('.sidebar-menu > :nth-child(4)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(5)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#deresults-gotodeconvolute', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#lcf-gotodeconvolute', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#batcheffect-gotodeconvolute', function () {
                      $('.sidebar-menu > :nth-child(5)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(6)').css('display', 'none');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#deresults-gotoprofile', function () {
                      $('.sidebar-menu > :nth-child(6)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#lcf-gotoprofile', function () {
                      $('.sidebar-menu > :nth-child(6)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      $(document).on('click', '#batcheffect-gotoprofile', function () {
                      $('.sidebar-menu > :nth-child(6)').css('display', 'inline');
                      $('.sidebar-menu > :nth-child(7)').css('display', 'none');
                      });
                      
                      $(document).on('click', '#deresults-summarisescores', function () {
                      $('.sidebar-menu > :nth-child(7)').css('display', 'inline');
                      });
                      $(document).on('click', '#deconvolute-summarisescores', function () {
                      $('.sidebar-menu > :nth-child(7)').css('display', 'inline');
                      });
                      $(document).on('click', '#profiling-summarisescores', function () {
                      $('.sidebar-menu > :nth-child(7)').css('display', 'inline');
                      });
                      
                      })
                      "))
}

###
## Auxiliary Support ####
###

#' getIconLabel
#' 
#' creates a label with information icon
#'
#' @param label label
#' @param message message
#'     
#' @export
#' 
getIconLabel <- function(label = NULL, message = NULL){
  
  icon_label <- tags$span(label, 
                          tags$i(
                            class = "glyphicon glyphicon-info-sign", 
                            style = "color:#0072B2;",
                            title = message
                          )
  )
  return(icon_label)
}


#' hideMulti
#' 
#' hides multiple shiny elements given a vector of input.id's
#'
#' @param IDs 
#'
#' @export
#'
hideMulti <- function(IDs){
  for(id in IDs){
    shinyjs::hide(id)
  }
} 

#' showMulti
#' 
#' shows multiple shiny elements given a vector of input.id's
#'
#' @param IDs 
#'
#' @export
#'
showMulti <- function(IDs){
  for(id in IDs){
    shinyjs::show(id)
  }
} 

#' RDS_from_web
#' 
#' download url from web and read as rds
#'
#' @param url 
#'
#' @export
#'
RDS_from_web <- function(url) {
  
  tempFile_location<- tempfile()
  download.file(url, tempFile_location)
  b <- readRDS(paste0(tempFile_location,"/countdata.rds"))
  file.remove(tempFile_location)
  b
}