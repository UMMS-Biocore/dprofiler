#' \code{dprofiler} package
#'
#' Dprofiler
#'
#' @docType package
#' @name dprofiler
#' 
NULL

utils::globalVariables(c("Conds", "ExpressionSet", "Samples","Scores", 'dc', 'demodata', 'demoprofdata',
                         'demoscdata', 'exprs', 'fData', 'grid.draw', 'grid.newpage', 'nnls', 'pData', 'pData<-',
                           'profileConds', 'silhouette', 'x', 'y',"CellType","nCount_integratedRNA_norm","Condition",
                         "DGEList","calcNormFactors"),
                       package = "dprofiler", add = FALSE)
# utils::globalVariables(c("."))