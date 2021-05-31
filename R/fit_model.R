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
fit_model <- function(data_matrix, method = "hierseg", arguments = list(),
                      bootstrap = FALSE, boot_sample_size = 100) {

  IMPLEMENTED_METHODS <- c("hierseg", "dynseg", "cvdynseg",
                           "cvseg", "fusedlasso")

  if ( !(method %in% IMPLEMENTED_METHODS) ) {
    stop("Error! The 'method' argument provided is not implemented!")
  }

  methodcall_name <- paste0("compute_", method) # package function name
  fit_arguments   <- c(list(data_matrix = data.matrix(data_matrix)), arguments)
  model <- do.call(methodcall_name, fit_arguments)

  n <- nrow(data_matrix)
  m <- ncol(data_matrix)

  # Bootstrap computation
  if(bootstrap){
    cp_freq <- rep(0, m) # Probability of a column being detected as a cp
    sym_diff_samp <- rep(0, boot_sample_size)
    rand_samp <- rep(0, boot_sample_size)
    ncp_samp <- rep(0, boot_sample_size)

    for(i in 1:boot_sample_size){
      boot_samp <- sample(n, n, replace = TRUE)
      fit_arguments_boot <- c(list(data_matrix = data.matrix(data_matrix[boot_samp, ])), arguments)
      model_boot <- do.call(methodcall_name, fit_arguments_boot)
      cp_freq[model_boot$changepoints] = cp_freq[model_boot$changepoints] + 1
      sym_diff_samp[i] = compute_symdiff(model$changepoints,
                                         model_boot$changepoints)
      rand_samp[i] = compute_rand(model$changepoints,
                                  model_boot$changepoints, m)
      ncp_samp[i] = length(model_boot$changepoints)

    }
    bootstrap_info <- list(b_samples = boot_sample_size,
                           cp_freq = cp_freq,
                           sym_diff_samp = sym_diff_samp,
                           rand_samp = rand_samp,
                           ncp_samp = ncp_samp)

    model$bootstrap_info = bootstrap_info
  }

  model$metadata <- list(method = method, arguments = model$arguments,
                         n = n, m = m,
                         columns = colnames(data_matrix))
  model$arguments <- NULL
  class(model) <- "bincpd"

  return(model)

}

#' @title
#' Plot change point blocks with probabilties
#'
#' @description
#'
#' Method of the generic \link[=plot]{plot} for bincpd objects.
#'
#' @export
plot.bincpd <- function(object, ...){
  m <- object$metadata$m
  method_name <- object$metadata$method
  n_cp <- object$n_cp
  block_ends <- c(1, object$changepoints, m)
  probs <- c(object$probabilities, object$probabilities[length(object$probabilities)])
  plot(block_ends, probs, type = "s",
       xlab = "",
       ylab = "Block probability",
       main = paste0("Block plot for ", method_name)
       )
}
