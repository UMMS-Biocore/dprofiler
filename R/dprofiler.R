#' \code{dprofiler} package
#'
#' Dprofiler
#'
#' @docType package
#' @name dprofiler
#' @importFrom dplyr %>%
#' @importFrom purrr %||%
NULL

utils::globalVariables(c("Conds", "ExpressionSet", "Samples","Scores", 'dc', 'demodata', 'demoprofdata',
                         'demoscdata', 'exprs', 'fData', 'grid.draw', 'grid.newpage', 'nnls', 'pData', 'pData<-',
                           'profileConds', 'silhouette', 'x', 'y',"waiting_screen","CellType","nCount_integratedRNA_norm"),
                       add = FALSE)