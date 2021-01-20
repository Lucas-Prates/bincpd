#' @title
#' Segments the data using dynamical programming algorithm
#'
#' @description
#' Computes the exact solution of the regularized loss optimization problem,
#' providing change point locations and the probabilities of each blocks
#'
#' @param data_matrix Data matrix to perform change point analysis
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for performance issues on large data sets
#' @param lambda Penalization constant
#' @param pen_func The penalization function used for the computation
#' of the estimator. A user specified function can be provided,
#' with function signature/arguments is function(left_index, right_index, n, m),
#' which represents the loss for the block interval i:j with left_index.
#' The default loss is the BIC loss \eqn{-log(n)(|cp|+1)}.
compute_dynseg <- function(data_matrix,
                           segthr = NULL,
                           lambda = 1,
                           pen_func = bic_loss_hs){

  n <- nrow(data_matrix)
  m <- ncol(data_matrix)
  if (is.null(segthr)) {segthr <- m - 1}

  # Penalization function that will be called in .cpp extension
  ds_pen_function <- function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, n, m) )
  }

  ds_output <- compute_dynseg_cpp(data_mat = data_matrix,
                                  ncol = m,
                                  segthr = segthr,
                                  pen_func = ds_pen_function)

  arguments <- list(segthr = segthr,
                    lambda = lambda,
                    pen_func = pen_func)

  model_info <- list(changepoints = ds_output[[1]],
                     probabilities = ds_output[[2]],
                     loss = ds_output[[3]],
                     n_cp = ds_output[[4]],
                     arguments = arguments)

  return(model_info)

}
