#----------------------- Utility functions -------------------------------------
#' NA vector for rBBM
#'
#' Returns a vector containing NA or 1's that will be multiplied by the data
#' matrix to generate NA's in the BBM.
#'
#' @noRd
generate_na <- function(n_samples, prob){
  sample = stats::rbinom(n_samples, size = 1, prob = 1 - prob)
  sample[sample == 0] = NA
  return(sample)
}

#----------------------- Main function -----------------------------------------
#' @title
#' Sampler for the Binary CP Block Model
#'
#' @description
#' Creates a \eqn{n \times m} matrix with \eqn{k} change points. In between
#' change points, the random variables are i.i.d. Bernoullis with the same
#' probability parameter.
#'
#' @param n Sample size
#' @param m Array length
#' @param k Number of change points.  The number of blocks is \eqn{k + 1}. It is
#' overridden if changepoints is non-NULL.
#' @param prob_vec Vector \eqn{k + 1} dimensional containing the probability
#' parameter of each block. If NULL, the parameters are sampled independently
#' from a uniform distribution
#' @param changepoints A \eqn{k} dimensional vector containing the
#'  change point locations. The change points are between 1 and \eqn{m-1}. If
#'  NULL, the change points are sampled uniformly in \eqn{[1, m-1]}
#' @param  prob_NA Probability of each entry of being NA
#' @export
rBBM <- function(n = 100,
                 m = 50,
                 k = 1,
                 prob_vec = NULL,
                 changepoints = NULL,
                 prob_NA = NULL) {

  # Setup variables if NULL or check for input errors if the user specified the
  # arguments
  #---------
  if (is.null(changepoints)) {
    if ((k >= m)||(k < 0)){
      stop(paste0("Input error! The number of change points k must be between ",
                  "0 and m."))
    }
    changepoints <- sort(c(0, sample(1:(m-1), k, replace = FALSE), m))
  }

  # Check if the change point vector provided is valid
  # and append auxiliary change points
  else {

    if ((any(changepoints <= 0)) || (any(changepoints >= m))) {
      stop("Input error! Change point vector entries must vary from 1 to m-1.")
    }

    k <- length(changepoints)
    # Auxiliary change points for sampling
    changepoints <- c(0, sort(changepoints), m)

  }

  if (is.null(prob_vec)) {
    prob_vec <- stats::runif(k+1)
  }

  # Check if the probability vector provided is valid.
  else{
    if ((any(prob_vec < 0)) || (any(prob_vec > 1))) {
      stop("Input error! Probability vector entries must vary from 0 to 1.")
    }
  }

  # Check if probability vector and change point vector have compatible sizes
  # The -2 on the LHS is due to the append of 2 extra auxiliary change points!
  if ( (length(changepoints) - 2) != (length(prob_vec) - 1)) {
    stop(paste0("Input error! Length of the change point vector must be equal ",
                "to the length of the probability vector - 1."))
  }
  #---------

  # Generate data
  data_matrix <- matrix(0, nrow = n, ncol = m)
  for (i in 1:(k + 1)) {
    # A block starts at changepoints[i] + 1 and ends at changepoints[i+1]
    interval <- (changepoints[i] + 1):(changepoints[i + 1])
    n_samples <- length(interval) * n
    data_matrix[, interval] <- stats::rbinom(n_samples,
                                          size = 1 ,
                                          prob = prob_vec[i])
    if (!is.null(prob_NA)) {
      data_matrix[, interval] <- data_matrix[, interval] * generate_na(n_samples,
                                                                 prob_NA)
    }
  }

  # Put relevant information in a list and remove auxiliary change points
  data_list <- list(data_matrix = data.matrix(data_matrix),
                    changepoints = changepoints[-c(1, k+2)],
                    probabilities = prob_vec)
  return(data_list)
}
