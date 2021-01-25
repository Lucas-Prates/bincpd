#' @title
#' Fits a binary CPD model according to the method chosen
#'
#' @description
#' The main function of the package, wrapping up all methods that perform binary
#' change point detection. It calls different functions depending on the method
#' argument passed. In order how to pass arguments correctly for each method,
#' their documentations are provided, albeit they can not be directly called.
#' The function returns a S3 object of the type bincpd. Some methods have more
#' or less attributes associated with their fitting procedure, but all of them
#' returns at least:
#' \itemize{
#'  \item[changepoints] a list containing the set of estimated change points;
#'  \item[probabilities}] a list containing the estimated probability parameters
#'  for each block;
#'  \item[loss] the final loss evaluated on the entire data set for the
#'  returned model;
#'  \item[n_cp] number of change points estimated.
#' }
#'
#' @param data_matrix Data frame or matrix containing the data set codified
#' in 0 or 1. The data can contain several rows, the methods use all data to
#' estimate more accurately the change point locations.
#'
#' @param method The method that will be used to fit the model. The current
#' implemented models are:
#' \itemize{
#'  \item[\link[=compute_hierseg]{hierseg}] Hierarchical segmentation, also
#'  known as binary segmentation;
#'  \item[\link[=compute_dynseg]{'dynseg'}] Dynamical programming segmentation;
#'  \item[\link[=compute_cvdynseg]{'cvdynseg'}] Dynseg with cross validation
#'  on regularization constant;
#'  \item[\link[=compute_cvseg]{'cvseg'}] Dynseg with cross validation on the
#'  number of change points;
#'  \item[\link[=compute_fusedlasso]{'fusedlasso'}] Fused Lasso segmentation.
#' }
#' @param arguments list containing the arguments that will be used to fit the
#' selected method. Each method might contain unique parameter.
#' @export
fit_model <- function(data_matrix, method = "hierseg", arguments = list()) {

  IMPLEMENTED_METHODS <- c("hierseg", "dynseg", "cvdynseg",
                           "cvseg", "fusedlasso")

  if ( !(method %in% IMPLEMENTED_METHODS) ) {
    stop("The 'method' argument provided is not implemented!")
  }

  methodcall_name <- paste0("compute_", method) # package function name
  fit_arguments   <- c(list(data_matrix = data.matrix(data_matrix)), arguments)
  model <- do.call(methodcall_name, fit_arguments)

  model$metadata <- list(method = method, arguments = model$arguments,
                         n = nrow(data_matrix), m = ncol(data_matrix),
                         columns = colnames(data_matrix))
  model$arguments <- NULL
  class(model) <- "bincpd"

  return(model)

}
