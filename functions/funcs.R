venn.diagram <- function (x, filename, height = 3000, width = 3000, resolution = 500, 
                          imagetype = "tiff", units = "px", compression = "lzw", 
                          na = "stop", main = NULL, sub = NULL, main.pos = c(0.5, 
                                                                             1.05), main.fontface = "plain", main.fontfamily = "serif", 
                          main.col = "black", main.cex = 1, main.just = c(0.5, 
                                                                          1), sub.pos = c(0.5, 1.05), sub.fontface = "plain", 
                          sub.fontfamily = "serif", sub.col = "black", 
                          sub.cex = 1, sub.just = c(0.5, 1), category.names = names(x), 
                          force.unique = TRUE, print.mode = "raw", sigdigs = 3, 
                          direct.area = FALSE, area.vector = 0, hyper.test = FALSE, 
                          total.population = NULL, lower.tail = TRUE, ...) 
{
  time.string = gsub(":", "-", gsub(" ", "_", as.character(Sys.time())))
  if (direct.area) {
    if (1 == length(area.vector)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = area.vector[1], 
                                                 category = list.names, ind = FALSE, ...)
    }
    if (3 == length(area.vector)) {
      grob.list <- VennDiagram::draw.pairwise.venn(area1 = area.vector[1], 
                                                   area2 = area.vector[2], cross.area = area.vector[3], 
                                                   category = category.names, ind = FALSE, print.mode = print.mode, 
                                                   sigdigs = sigdigs, ...)
    }
  }
  else {
    if (force.unique) {
      for (i in 1:length(x)) {
        x[[i]] <- unique(x[[i]])
      }
    }
    if ("none" == na) {
      x <- x
    }
    else if ("stop" == na) {
      for (i in 1:length(x)) {
        if (any(is.na(x[[i]]))) {
          flog.error("NAs in dataset", call. = FALSE, 
                     name = "VennDiagramLogger")
          stop("NAs in dataset", call. = FALSE)
        }
      }
    }
    else if ("remove" == na) {
      for (i in 1:length(x)) {
        x[[i]] <- x[[i]][!is.na(x[[i]])]
      }
    }
    else {
      flog.error("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"", 
                 name = "VennDiagramLogger")
      stop("Invalid na option: valid options are \"none\", \"stop\", and \"remove\"")
    }
    if (0 == length(x) | length(x) > 5) {
      flog.error("Incorrect number of elements.", 
                 call. = FALSE, name = "VennDiagramLogger")
      stop("Incorrect number of elements.", call. = FALSE)
    }
    if (1 == length(x)) {
      list.names <- category.names
      if (is.null(list.names)) {
        list.names <- ""
      }
      grob.list <- VennDiagram::draw.single.venn(area = length(x[[1]]), 
                                                 category = list.names, ind = FALSE, ...)
    }
    else if (2 == length(x)) {
      grob.list <- draw.pairwise.venn(area1 = length(x[[1]]), 
                                      area2 = length(x[[2]]), cross.area = length(intersect(x[[1]], 
                                                                                            x[[2]])), category = category.names, ind = FALSE, 
                                      print.mode = print.mode, sigdigs = sigdigs, ...)
    }
    else {
      flog.error("Invalid size of input object", 
                 name = "VennDiagramLogger")
      stop("Invalid size of input object")
    }
  }
  if (length(x) == 2 & !is.null(total.population) & hyper.test) {
    val.p = calculate.overlap.and.pvalue(x[[1]], x[[2]], 
                                         total.population, lower.tail = lower.tail)
    if (is.null(sub)) {
      sub = paste0("p = ", signif(val.p[3], digits = 2))
    }
    else {
      sub = paste0(sub, ", p = ", signif(val.p[3], 
                                         digits = 2))
    }
  }
  if (!is.null(sub)) {
    grob.list <- add.title(gList = grob.list, x = sub, pos = sub.pos, 
                           fontface = sub.fontface, fontfamily = sub.fontfamily, 
                           col = sub.col, cex = sub.cex)
  }
  if (!is.null(main)) {
    grob.list <- add.title(gList = grob.list, x = main, pos = main.pos, 
                           fontface = main.fontface, fontfamily = main.fontfamily, 
                           col = main.col, cex = main.cex)
  }
  if (!is.null(filename)) {
    current.type <- getOption("bitmapType")
    if (length(grep("Darwin", Sys.info()["sysname"]))) {
      options(bitmapType = "quartz")
    }
    else {
      options(bitmapType = "cairo")
    }
    if ("tiff" == imagetype) {
      tiff(filename = filename, height = height, width = width, 
           units = units, res = resolution, compression = compression)
    }
    else if ("png" == imagetype) {
      png(filename = filename, height = height, width = width, 
          units = units, res = resolution)
    }
    else if ("svg" == imagetype) {
      svg(filename = filename, height = height, width = width)
    }
    else {
      flog.error("You have misspelled your 'imagetype', please try again", 
                 name = "VennDiagramLogger")
      stop("You have misspelled your 'imagetype', please try again")
    }
    grid.draw(grob.list)
    dev.off()
    options(bitmapType = current.type)
    return(1)
  }
  return(grob.list)
}

getLoadingMsg <- function (output = NULL) 
{
  addResourcePath(prefix = "www", directoryPath = system.file("extdata", 
                                                              "www", package = "debrowser"))
  imgsrc_full <- "www/images/loading_start.gif"
  imgsrc_small <- "www/images/loading.gif"
  
  a <- list(tags$head(tags$style(type = "text/css", "\n            #loadmessage {\n            position: fixed;\n            top: 0px;\n            left: 0px;\n            width: 100%;\n            height: 100%;\n            padding: 5px 0px 5px 0px;\n            text-align: center;\n            font-weight: bold;\n            font-size: 100%;\n            color: #000000;\n            opacity: 0.8;\n            z-index: 100;\n            }\n            #loadmessage_small {\n            position: fixed;\n            left: 50%;\n            transform: translateX(-50%);\n            top: 50px;\n            text-align: center;\n            opacity: 0.8;\n            z-index: 999999;\n            }\n                             ")), 
            conditionalPanel(condition = paste0("$('html').hasClass('shiny-busy')", 
                                                "& input.MenuItems=='CondSelect'"), 
                             tags$div(id = "loadmessage", tags$img(src = imgsrc_full))), 
            conditionalPanel(condition = paste0("$('html').hasClass('shiny-busy')", 
                                                "& !(input.MenuItems=='CondSelect')"), 
                             tags$div(id = "loadmessage_small", tags$img(src = imgsrc_small))))
}

getLoadingMsg1 <- function (output = NULL) 
{
  addResourcePath(prefix = "www", directoryPath = system.file("extdata", 
                                                              "www", package = "debrowser"))
  imgsrc_full <- "www/images/loading_start.gif"
  imgsrc_small <- "www/images/loading.gif"
  
  a <- list(tags$head(tags$style(type = "text/css", "\n            #loadmessage {\n            position: fixed;\n            top: 0px;\n            left: 0px;\n            width: 100%;\n            height: 100%;\n            padding: 5px 0px 5px 0px;\n            text-align: center;\n            font-weight: bold;\n            font-size: 100%;\n            color: #000000;\n            opacity: 0.8;\n            z-index: 100;\n            }\n            #loadmessage_small {\n            position: fixed;\n            left: 50%;\n            transform: translateX(-50%);\n            top: 50px;\n            text-align: center;\n            opacity: 0.8;\n            z-index: 999999;\n            }\n                             ")), 
            conditionalPanel(condition = paste0("$('html').hasClass('shiny-busy')"),                  
            tags$div(id = "loadmessage", tags$img(src = imgsrc_full))))
}

getLoadingMsg2 <- function (output = NULL) 
{
  addResourcePath(prefix = "www", directoryPath = system.file("extdata", 
                                                              "www", package = "debrowser"))
  imgsrc_full <- "www/images/loading_start.gif"
  imgsrc_small <- "www/images/loading.gif"
  
  a <- list(tags$head(tags$style(type = "text/css", "\n            #loadmessage {\n            position: fixed;\n            top: 0px;\n            left: 0px;\n            width: 100%;\n            height: 100%;\n            padding: 5px 0px 5px 0px;\n            text-align: center;\n            font-weight: bold;\n            font-size: 100%;\n            color: #000000;\n            opacity: 0.8;\n            z-index: 100;\n            }\n            #loadmessage_small {\n            position: fixed;\n            left: 50%;\n            transform: translateX(-50%);\n            top: 50px;\n            text-align: center;\n            opacity: 0.8;\n            z-index: 999999;\n            }\n                             ")), 
            conditionalPanel(condition = paste0("$('html').hasClass('shiny-busy')"),   
            tags$div(id = "loadmessage", tags$img(src = imgsrc_small))))
}
