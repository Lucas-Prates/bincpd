#' @title
#' Segments the data using dynamical programming with cross validation on the
#' number of segments
#'
#' @description
#' Minimizes the k-fold CV mean negative log likelihood computed, no
#' penalization term included. For each fold, the dynseg algorithm is run and we
#' compute the best estimator best estimator for each number of change points
#' and their respective log likelihoods. In the end, the algorithm is rerun in
#' the whole data set using and select the best estimator from the best number
#' of change points.
#'
#' @param data_matrix Data matrix to perform change point analysis
#' @param n_folds Number of folds used in the k-fold CV
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for performance issues on large data sets
compute_cvseg <- function(data_matrix,
                          segthr = NULL,
                          n_folds = 5) {
  n <- nrow(data_matrix)
  m <- ncol(data_matrix)
  if (is.null(segthr)) {segthr <- m - 1}
  else if ((segthr >= m)||(segthr < 0))  stop(paste0("Argument segthr must be",
                                                    "between 0 and m."))

  cvseg_output <- compute_cvseg_cpp(data_matrix,
                             segthr,
                             n_folds,
                             m,
                             n)

  arguments <- list(segthr = segthr, n_folds = n_folds)

  model_info <- list(changepoints = cvseg_output[[1]],
                     probabilities = cvseg_output[[2]],
                     loss = cvseg_output[[3]],
                     n_cp = cvseg_output[[4]],
                     loglike_cv = cvseg_output[[5]],
                     arguments = arguments)

  return(model_info)

}
