#' @title
#' Segments the data using dynamical programming with cross validation on the
#' regularization constant \eqn{\lambda}
#'
#' @description
#' Minimizes the k-fold CV mean negative regularized loss.
#' For each fold, the dynseg algorithm is run and we compute the best estimator
#' best estimator for each penalization constant and their respective log
#' likelihoods. In the end, the algorithm is rerun in the whole data set using
#' and select the best estimator from the best number of change points.
#'
#' @param data_matrix Data matrix to perform change point analysis
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for performance issues on large data sets
#' @param lambda_set The set of penalization constants that will be used in the
#' k-fold CV.
#' @param n_folds Number of folds used in the k-fold CV
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, n, m),
#'   where the left_index:right_index is the integer interval, n the sample
#'   size and m the number of variables/columns.
compute_cvdynseg <- function(data_matrix,
                             segthr = NULL,
                             lambda_set = c(0.1, 1, 10),
                             n_folds = 5,
                             pen_func = bic_loss_hs) {
  n <- nrow(data_matrix)
  m <- ncol(data_matrix)
  if (is.null(segthr)) {segthr <- m - 1}
  else if ((segthr >= m)||(segthr < 0))  stop(paste0("Argument segthr must be",
                                                     "between 0 and m."))

  # A function builder passed to the cpp in order to build the correct
  # regularization functions for each value of lambda
  create_pen_func <- function(lambda, n, m){
    pen_func_aux <- function(left_index, right_index){
      return(lambda*pen_func(left_index, right_index, n, m))
    }
    return(pen_func_aux)
  }

  cvseg_output <- compute_cvdynseg_cpp(data_matrix,
                                       segthr,
                                       lambda_set,
                                       n_folds,
                                       m,
                                       n,
                                       create_pen_func)

  arguments <- list(segthr = segthr,
                    n_folds = n_folds,
                    lambda_set = lambda_set,
                    pen_func = pen_func)

  model_info <- list(changepoints = cvseg_output[[1]],
                     probabilities = cvseg_output[[2]],
                     loss = cvseg_output[[3]],
                     n_cp = cvseg_output[[4]],
                     best_lambda = cvseg_output[[5]],
                     loglike_cv = cvseg_output[[6]],
                     lambda_set = lambda_set,
                     arguments = arguments)

  return(model_info)
}
