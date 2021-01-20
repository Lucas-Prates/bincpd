#' @title
#' Segments the data using hierarchical algorithm
#'
#' @description
#' Uses binary splitting in a greedy solution of the regularized loss
#' optimization problem to obtain the change points and probabilities of each
#' block
#'
#' @param data_matrix Data matrix to perform change point analysis
#' @param lambda Penalization constant
#' @param pen_func A penalization function defined i integer intervals
#'   The function signature should be pen(left_index, right_index, n, m),
#'   where the left_index:right_index is the integer interval, n the sample
#'   size and m the number of variables/columns.
compute_hierseg <- function(data_matrix,
                            lambda = 1,
                            pen_func = bic_loss_hs) {
  n <- nrow(data_matrix)
  m <- ncol(data_matrix)

  # Penalization function that will be called in .cpp extension
  hs_pen_function <- function(left_index, right_index) {
    return( lambda * pen_func(left_index, right_index, n, m) )
  }

  hs_output <- compute_hierseg_cpp(data_matrix = data_matrix,
                                   n = n,
                                   m = m,
                                   hs_pen_function)

  arguments <- list(lambda = lambda,
                    pen_func = pen_func)

  model_info <- list(changepoints = hs_output[[1]],
                  probabilities = hs_output[[2]],
                  loss = hs_output[[3]],
                  n_cp = length(hs_output[[1]]),
                  arguments = arguments)

  return(model_info)

}
