#' Segments the data using hierarchical algorithm
#'
#' @description
#' Uses binary splitting in a greedy solution approach to obtain the
#' log likelihood, change points and probabilities
#' as a function of the number of segments
#'
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, n, m),
#'   where the left_index:right_index is the integer interval, n the sample
#'   size and m the number of variables/columns.
compute_hierseg <- function(data_matrix,
                            lambda = 1,
                            pen_func = bic_loss_hs) {
  n <- nrow(data_matrix)
  m <- ncol(data_matrix)
  hs_pen_function <- function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, n, m) )
  }

  hs_output <- hierarchical_algorithm_cpp(data_matrix = data_matrix,
                                          n = n,
                                          m = m,
                                          pen_func = hs_pen_function)

  outlist <- list(changepoints = hs_output[[1]],
                 probabilities = hs_output[[2]],
                 loss = hs_output[[3]])

}
