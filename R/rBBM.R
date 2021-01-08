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

#' Sampler for the Binary CP Block Model
#'
#' Creates a \eqn{n \times m} matrix with \eqn{k} change points. In between
#' change points, the random variables are i.i.d. Bernoullis with the same
#' probability parameter.
#' @param n Sample size
#' @param m Array length
#' @param k Number of change points.  The number of blocks is \eqn{k + 1}.
#' @param prob_vec Vector \eqn{k + 1} dimensional containing the probability
#' parameter of each block. If NULL, the parameters are sampled independently
#' from a uniform distribution
#' @param changepoints A \eqn{k} dimensional vector containing the
#'  change point locations. The change points are between 1 and \eqn{m-1}. If
#'  NULL, the change points are sampled uniformly in \eqn{[1, m-1]}
#' @param  prob_NA Probability of each entry of being NA
#' @export
rBBM <- function(n = 100,
                 m = 1000,
                 k = 10,
                 prob_vec = NULL,
                 changepoints = NULL,
                 prob_NA = NULL) {

  if (is.null(changepoints)) {
    changepoints = sort(c(0, sample(1:(m-1), k, replace = FALSE), m))
  }

  # Check if the change point vector provided is valid
  # and append auxiliary change points
  else {

    if ((any(changepoints <= 0)) || (any(changepoints >= m))) {
      stop("Change point vector entries must vary from 1 to m-1.")
    }

    # Auxiliary change points for sampling
    changepoints = c(0, sort(changepoints), m)
  }

  if (is.null(prob_vec)) {
    prob_vec  = stats::runif(k+1)
    data_mat = matrix(0, nrow = n, ncol = m)
  }

  # Check if the probability vector provided is valid.
  else{
    if ((any(prob_vec < 0)) || (any(prob_vec > 1))) {
      stop("Probability vector entries must vary from 0 to 1.")
    }
  }

  for (i in 1:(k + 1)) {
    # A block starts at changepoints[i] + 1 and ends at changepoints[i+1]
    interval = (changepoints[i] + 1):(changepoints[i + 1])
    n_samples = length(interval) * n
    data_mat[, interval] = stats::rbinom(n_samples,
                                         size = 1 ,
                                         prob = prob_vec[i])
    if (!is.null(prob_NA)) {
      data_mat[, interval] = data_mat[, interval] * generate_na(n_samples,
                                                                prob_NA)
    }
  }

  # Put relevant information in a list and remove auxiliary change points
  data_list = list(data_mat = data.matrix(data_mat),
                   changepoints = changepoints[-c(1, k+2)],
                   probabilities = prob_vec)
  return(data_list)
}
