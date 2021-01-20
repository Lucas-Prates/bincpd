# Requires: CVXR - R package for Disciplined Convex Optimization

#----------------------- Utility functions -------------------------------------
# Get the number of 0's and 1's in each column. Entries with Na's are removed.
get_nobs <- function(data_matrix){

    n1_vec <- as.numeric(apply(data_matrix, 2, sum, na.rm = TRUE))
    n0_vec <- as.numeric(apply( (1 - data_matrix), 2, sum, na.rm = TRUE))
    nobs_vec <- list(n0_vec = n0_vec, n1_vec = n1_vec)

    return(nobs_vec)

}

#----------------------- Main function -----------------------------------------
#' @title
#' Computes the Fused Lasso estimator.
#'
#' @description
#' The estimator is parameter that minimizes the fused lasso loss defined by
#' \eqn{-l(p; X) + J(n) \times lambda \times sum_{i=1}^m|p_i - p_{i+1}|},
#'  where l is the log likelihood, and the penalization is on the difference of
#'  the probability parameters
#' @param data_matrix Data matrix to perform change point analysis
#' @param lambda Penalization constant
#' @param sample_pen_func The penalization function J(n) multiplying
#' the fused regularization term
#' @param  precision Precision to calculate the change points.
#' If \eqn{|p_{i+1} - p_{i}| > precision}, then \eqn{i} is considered
#' a change point.
#' @param MLE_prob If TRUE, the probability parameter will be estimated using
#' Maximum likehood approach, which is the mean of the block. So only
#' the change points will be determined by the L1 loss.
#' @param check_dcp a flag to check if the optimization problem satisfies the
#' DCP rules. Checking the rules is expensive, and since our problem satisfies
#' the rules, the default is FALSE.
compute_fusedlasso <- function(data_matrix,
                               lambda = 1,
                               sample_pen_func = log,
                               precision = 1e-3,
                               MLE_prob = FALSE,
                               check_dcp = FALSE
                               ){

  nobs <- get_nobs(data_matrix)
  n0_vec <- nobs$n0_vec
  n1_vec <- nobs$n1_vec
  n <- nrow(data_matrix)
  m <- ncol(data_matrix)

  ## CVXR setup
  prob_vec <- CVXR::Variable(m)

  # components of objective functions
  loglike <- -sum(n0_vec*log(1-prob_vec) + n1_vec*log(prob_vec)) #/n
  fused_loss <- sample_pen_func(n)*(lambda*sum(abs(CVXR::diff(prob_vec))))
  #

  objective <- CVXR::Minimize(loglike + fused_loss)
  problem <- CVXR::Problem(objective)
  prob_data <- CVXR::get_problem_data(problem, solver = "ECOS")
  ##

  # Solving the problem using CVXR
  # Solve checking dcp rules
  if(check_dcp) {result <- CVXR::solve(problem)}

  # Solve the problem without checking dcp rules, which is faster!
  # CVXR syntax changed a few things after 1.x, so we solve for any version
  else{
    if (utils::packageVersion("CVXR") > "0.99-7") {
      ECOS_dims <- CVXR::ECOS.dims_to_solver_dict(prob_data$data[["dims"]])
    } else {
      ECOS_dims <- prob_data$data[["dims"]]
    }

    solver_output <- ECOSolveR::ECOS_csolve(c = prob_data$data[["c"]],
                                            G = prob_data$data[["G"]],
                                            h = prob_data$data[["h"]],
                                            dims = ECOS_dims,
                                            A = prob_data$data[["A"]],
                                            b = prob_data$data[["b"]])
    if (utils::packageVersion("CVXR") > "0.99-7") {
      result <- CVXR::unpack_results(problem, solver_output,
                               prob_data$chain, prob_data$inverse_data)
    } else {
      result <- CVXR::unpack_results(problem, "ECOS", solver_output)
    }
  }

  # summarize results
  probabilities <- as.numeric(result$getValue(prob_vec))
  fused_cp <-  abs(diff(probabilities))
  fused_cp <- which(fused_cp > precision)
  fused_block_begin = 1 + c(0, fused_cp)
  fused_block_end <- c(fused_cp, m)

  probabilities <- probabilities[c(1, fused_cp + 1)]


  if(MLE_prob){
    for(i in 1:length(probabilities)){
      probabilities[i] <- mean(data_matrix[, fused_block_begin[i]:fused_block_end[i]],
                               na.rm = TRUE)
    }
  }


  arguments <- list(lambda = lambda,
                    sample_pen_func = sample_pen_func,
                    precision = precision,
                    MLE_prob = MLE_prob,
                    check_dcp = check_dcp)

  model_info <- list(changepoints = fused_cp,
                     probabilities = probabilities,
                     loss = result$value,
                     n_cp = length(fused_cp),
                     CVXR_output = list(solver = result$solver,
                                        num_iter = result$num_iter,
                                        solve_time = result$solve_time,
                                        setup_time = result$setup_time,
                                        status = result$status
                                        ),
                     arguments = arguments
                     )

  return(model_info)

}
