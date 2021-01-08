# This script requires CVXR for Disciplined Convex Optimization
# check_dcp is a flag to check for the DCP rules.
# Checking the rules is a bit expensive,
# so since we know that our problem follows the rules,
# we usually don't check them
# especially on the cv_fused_cvxr algorithm
### Auxiliar functions

# Get the number of 0's and 1's in each column
# NA's are removed
get_nobs__ <- function(data){

    n1_vec <- as.numeric(apply(data, 2, sum, na.rm = TRUE))
    n0_vec <- as.numeric(apply( (1 - data), 2, sum, na.rm = TRUE))
    nobs_vec <- list(n0_vec = n0_vec, n1_vec = n1_vec)
    return(nobs_vec)

}

# Given the probability vector of the blocks and the change points,
# returns the probability vector for all positions

build_full_prob__ <- function(ncol, changepoint, prob){

  cp_begin <- 1 + c(0, changepoint)
  cp_end <- c(changepoint, ncol)
  full_prob <- rep(0, ncol)

  for(i in 1:length(cp_begin)) full_prob[cp_begin[i]:cp_end[i]] <- prob[i]

  return(full_prob)
}


###-------------
###Main function
#' Computes the Fused Lasso estimator.
#'
#' The estimator is the argmin of -l(p; X) + lambda*sum_{i=1}^m|p_i - p_{i+1}|,
#'  where l is the log likelihood, and the penalization is on the difference of
#'  the probability parameters
#' @param lambda Penalization constant
#'
#' @param  fused_prec Precision to calculate the change points
#' if |p_{i+1} - p_{i}| > fused_prec, then {i} is considered a change point.
#'
#' @param MLE_prob If TRUE, the probability parameter will be estimated using
#' Maximum Likehood approach (in thi case, the mean of the block). So only
#' the breaks will be determined by the L1 loss.
#' @export
fit_fused_cvxr <- function(data, lambda = 1, fused_prec = 1e-4,
                           check_dcp = FALSE, MLE_prob = FALSE){

  nobs <- get_nobs__(data)
  n0_vec <- nobs$n0_vec
  n1_vec <- nobs$n1_vec
  n = nrow(data); m = ncol(data)
  prob_vec <- CVXR::Variable(m)

  loglike = -sum(n0_vec*log(1-prob_vec) + n1_vec*log(prob_vec))/n
  fused_loss = (lambda*sum(abs(CVXR::diff(prob_vec))))

  # objective_func <- function(prob_vec){
  #   return(comp_fused_loglike(nobs$n0_vec, nobs$n1_vec, prob_vec, lambda = lambda))
  # }

  objective <- CVXR::Minimize(loglike + fused_loss)

  problem <- CVXR::Problem(objective)
  #result <- solve(problem)

  prob_data <- CVXR::get_problem_data(problem, solver = "ECOS")

  if(check_dcp) result <- CVXR::solve(problem)#, ignore_dcp = TRUE)

  # solve the problem without checking dcp (way faster!)
  # CVXR twinked a few things after 1.x, so we solve for any version
  else{
    if (packageVersion("CVXR") > "0.99-7") {
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
    if (packageVersion("CVXR") > "0.99-7") {
      result <- CVXR::unpack_results(problem, solver_output,
                               prob_data$chain, prob_data$inverse_data)
    } else {
      result <- CVXR::unpack_results(problem, "ECOS", solver_output)
    }
  }

  prob <- as.numeric(result$getValue(prob_vec)) # prob_vec %>% result$getValue() %>% as.numeric()
  fused_cp <-  abs(diff(prob))  # (prob %>% diff() %>% abs())
  fused_cp <- which(fused_cp > fused_prec)
  fused_block_begin = 1 + c(0, fused_cp)
  fused_block_end = c(fused_cp, m)

  prob <- prob[c(1, fused_cp + 1)]

  fused_metadata = list(n_samples = n,
                        n_col = m,
                        colnames = colnames(data),
                        method = 'fused',
                        fused_result = result)

  if(MLE_prob){
    for(i in 1:length(prob))
      prob[i] <- mean(data[, fused_block_begin[i]:fused_block_end[i]],
                      na.rm = TRUE)
  }

  outlist <- list(metadata = fused_metadata,
                  changepoints = fused_cp,
                  block_begin = fused_block_begin,
                  block_end = fused_block_end,
                  prob = prob)

  return(outlist)

}


fit_cvfused_cvxr <- function(data, lambda_set = seq(0.5, 10, by = 1),
                             nfolds = 5, fused_prec = 1e-4, seed = 10,
                             omit_progress = FALSE, MLE_prob = FALSE){

  set.seed(seed)
  init_time <- Sys.time()
  n_row <- nrow(data)
  n_col <- ncol(data)
  loss <- rep(0, length(lambda_set))
  lambda_num_block <- rep(0, length(lambda_set))
  rand_ind <- sample(1:n_row, rep = FALSE)
  n_cv <- floor(n_row/nfolds)

  for(i in 1:nfolds){

    if(i < nfolds) test_ind <- rand_ind[( (i-1)*n_cv + 1):(i*n_cv) ]
    else test_ind <- rand_ind[((i-1)*n_cv + 1):n_row]

    for(j in 1:length(lambda_set)){

      fit_fused <- fit_fused_cvxr(data[-test_ind, ], lambda = lambda_set[j],
                                  fused_prec = fused_prec)
      full_prob <- build_full_prob__(ncol = n_col, changepoint = fit_fused$changepoint,
                                  prob = fit_fused$prob)
      loss[j] <- loss[j] +
                    + comp_loss__(data = data[test_ind, ], full_prob = full_prob)
      lambda_num_block[j] <- lambda_num_block[j] + length((fit_fused$prob))

    }

    gc()
    if(!omit_progress) cat(paste0("Progress: ", 100*i/nfolds, "%\n"))

  }

  loss <- loss/nfolds
  lambda_num_block <- lambda_num_block/nfolds
  lambda_min <- lambda_set[which.min(loss)]

  model <- fit_fused_cvxr(data = data, lambda = lambda_min,
                              fused_prec = fused_prec, MLE_prob = MLE_prob)

  fplot <- function(){
    par(mfrow = c(1, 2))
    plot(lambda_set, loss, type = "l", xlab = "lambda", ylab = "Loss")
    plot(lambda_set, lambda_num_block, type = "l",
         xlab = "lambda", ylab = "Average number of blocks")

    par(mfrow = c(1, 1))
  }

  outlist <- list(lambda_set = lambda_set, loss = loss, lambda_num_block = lambda_num_block,
                  lambda_min = lambda_min, model = model, fplot = fplot)

  rm(.Random.seed, envir=globalenv())
  end_time <- Sys.time()
  if(!omit_progress) cat(paste0("Total time taken: ",
                                round(end_time - init_time, digits = 2), " seconds\n"))
  # Frees memory used
  gc()

  return(outlist)

}


