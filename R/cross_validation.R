
# Auxiliar function that calculates the test loglikelihood for the cv test batch
compute_test_loglike__ <- function(data_mat, segments, prob_vec){
  segments = c(segments, ncol(data_mat))
  data_mat_list = apply(data_mat, 2, na.omit)
  if( typeof(data_mat_list) == "list"){
    data_mat_length = sapply(data_mat_list, length)
    data_mat_sum = sapply(data_mat_list, sum)
  }else{
    data_mat_length = apply(data_mat_list, 2, length)
    data_mat_sum = apply(data_mat_list, 2, sum)
  }
  n = length(data_mat_sum)

  cumsums1 = cumsum(data_mat_sum)
  cumlengths1 = cumsum(data_mat_length)

  seglen = length(segments)
  blocksum = vector(length = length(seglen))
  blocklength = vector(length = length(seglen))
  blocksum[1] = cumsums1[segments[1]]
  blocklength[1] = cumlengths1[segments [1]]
  if(seglen > 1){
    for (k in 2:seglen){
      blocksum[k] =  cumsums1[segments[k]] - cumsums1[segments[k-1]]
      blocklength[k] =  cumlengths1[segments[k]] - cumlengths1[segments[k-1]]
    }
  }
  loglike = 0
  for (k in 1:seglen){
    if (prob_vec[k] != 0 & prob_vec[k] != 1){
      loglike = (loglike + blocksum[k]*log(prob_vec[k])
                 + (blocklength[k] - blocksum[k])*log(1 - prob_vec[k]))
    }
  }
  return(loglike)
}

###-------------
#' Estimates the number of segments using maximum likelihood and
#' K-fold cross validation
#'
#' After splitting the data in a k-fold CV form, the estimator calculates,
#' for each number of segments k, the maximum likehood
#' estimator conditional to that number of segments on the train data.
#' Then, it calculates the likelihood on the validation data. The estimator is
#' then recalculated using the whole dataset for the number of segmetns that
#' had the highest likelihood average.
#'
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for perfomance issues on large datasets
#' @param n_folds number of k-folds in cross valition
#' @param progress_omit Omit % of progress.
#' @param time_omit Omit time taken to run the estimator.
#' @export
fit_cv_seg <- function(data_mat, dynseg = NULL, segthr = NULL, nfolds = 5,
                            seed = NULL, progress_omit = FALSE, time_omit = FALSE){
  init_time = proc.time()[3]
  if(!is.null(seed)) set.seed(seed)
  numcol = ncol(data_mat)
  if(is.null(segthr)) segthr = numcol
  else segthr = min(segthr, numcol)
  numrow = nrow(data_mat)
  unused_indexes =  1:numrow
  test_batch_size = floor(numrow/nfolds)
  loglike_avg = vector(length = segthr) #segments loglikelihood average
  for (j in 1:nfolds){
    if(j != nfolds) nfolds_ind = sample(unused_indexes, size = test_batch_size, replace = FALSE)
    else nfolds_ind = unused_indexes
    unused_indexes = setdiff(unused_indexes, nfolds_ind)
    # First obtain the segments that maximize the pml
    dm = build_matrices_cpp__(data_mat[-nfolds_ind, ])#j-th nfold data
    ds_cpp = comp_dynseg_cpp_fast__(dm[[1]], dm[[2]], numcol, segthr)# [[1]] is the loglike vec
    for (numseg in 1:segthr){
      seg_vec = sort(ds_cpp[[2]][[numseg]])
      prob_vec = ds_cpp[[3]][[numseg]]
    # We now calculate the *loglikelihood on the cv test batch*
      loglike_uns = compute_test_loglike__(data_mat[nfolds_ind, ],
                    seg_vec, prob_vec)
    #  full_prob = build_full_prob(ncol(data_mat), seg_vec, prob_vec)
    #  loglike_uns = -comp_loss(data_mat[nfolds_ind, ], full_prob = full_prob)
      loglike_avg[numseg] = loglike_avg[numseg] + loglike_uns
    }

    gc()
    if (!progress_omit) cat("\nProgress:", round(100*(j/nfolds)), "%")

  }
  loglike_avg = loglike_avg/numrow
  seg_loss_min = which.min(-loglike_avg)

  if(is.null(dynseg)) dynseg = comp_dynseg__fast__(data_mat, numcol, segthr)

  cp <- dynseg$changepoints[[seg_loss_min]]
  prob <- dynseg$prob_list[[seg_loss_min]]
  final_model = list(numseg = seg_loss_min, changepoints = cp,
                     prob = prob)


  outlist = list(segments = 1:segthr, seg_loss_avg = -loglike_avg,
                 final_model = final_model)

  end_time = proc.time()[3]

  if (!time_omit) cat("\nTime Elapsed: ",end_time-init_time, " seconds\n")
  rm(.Random.seed, envir = globalenv())
  return(outlist)
}

#' Computes the penalized maximum likelihood estimator using cross validation
#' Estimates the number of segments using maximum likelihood and
#' K-fold cross validation
#'
#' After splitting the data in a k-fold CV form, the estimator calculates
#' the penalized maximum likehood estimator, which is the
#' argmin of -l(p; X) + lambda*R(p;X), for the train data.
#' Then, it calculates the likelihood on the validation data. The estimator is
#' then recalculated using the whole dataset for the number of segmetns that
#' had the highest likelihood average.
#'
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for perfomance issues on large datasets
#' @param n_folds number of k-folds in cross valition
#' @param c_interval penalization constant set. Each value of the set will be
#' used to compute
#' @param progress_omit Omit % of progress.
#' @param time_omit Omit time taken to run the estimator.
#' @export
fit_cv_pml <- function(data_mat, dynseg = NULL, segthr = NULL, nfolds = 5,
                       c_interval = seq(.5, 10, by = .25), pen_func = NULL,
                       seed = NULL, progress_omit = FALSE, time_omit = FALSE,
                       auto_search = FALSE, search_limit = 5){

  init_time = proc.time()[3]

  numcol = ncol(data_mat)
  numrow = nrow(data_mat)
  if(!is.null(seed)) set.seed(seed)
  if(is.null(segthr)) segthr = numcol
  else segthr = min(segthr, numcol)
  auto_search_flag = TRUE
  num_searchs = 1
  c_len = length(c_interval)
  test_batch_size = floor(numrow/nfolds)
  while( (auto_search_flag == TRUE) && (num_searchs <= search_limit) ){
    unused_indexes =  1:numrow
    loglike_avg = rep(0, c_len) #segments loglikelihood average
    c_numseg_avg = rep(0, c_len)
    for (j in 1:nfolds){
      if(j != nfolds) nfolds_ind = sample(unused_indexes, size = test_batch_size, replace = FALSE)
      else nfolds_ind = unused_indexes
      unused_indexes = setdiff(unused_indexes, nfolds_ind)
      # First obtain the segments that maximize the pml
      dm = build_matrices_cpp__(data_mat[-nfolds_ind, ])#j-th nfold data
      ds_cpp = comp_dynseg_cpp_fast__(dm[[1]], dm[[2]], numcol, segthr)# [[1]] is the loglike vec
      for (i in 1:c_len){
        if(is.null(pen_func)){
          pen_ll = -ds_cpp[[1]] + c_interval[i]*((1 + (1:segthr)))*log(numrow)

        }
        else{
          pen_ll = -ds_cpp[[1]]
          for(i in 1:segthr){
            pen_ll[i] = pen_ll[i] + lambda*pen_func(n_sample, n_col,
                                              ds_cpp[[3]][[i]],
                                              ds_cpp[[2]][[i]])
          }
        }

        numseg = which.min(pen_ll)
        seg_vec = sort(ds_cpp[[2]][[numseg]])
        prob_vec = ds_cpp[[3]][[numseg]]
        # We now calculate the *loglikelihood on the cv test batch*
        loglike_uns = compute_test_loglike__(data_mat[nfolds_ind, ],
                                      seg_vec, prob_vec)
        loglike_avg[i] = loglike_avg[i] + loglike_uns
        c_numseg_avg[i] = c_numseg_avg[i] + numseg
      }

      gc()
      if (!progress_omit) cat("\nProgress:", round(100*(j/nfolds)), "%")

    }
    loglike_avg = loglike_avg/numrow
    c_numseg_avg = c_numseg_avg/nfolds
    loglike_max = which.max(loglike_avg)
    if(auto_search == FALSE) auto_search_flag = FALSE
    else{
      if((loglike_max == 1) && (loglike_avg[1] != max(loglike_avg[-1]))){
        c_interval = seq(c_interval[1]/20, c_interval[1], by = c_interval[1]/20)
        c_len = length(c_interval)
        num_searchs = num_searchs + 1
        cat(blue("\nMaximum on the lower boundary, starting new search_"))
      }
      else if((loglike_max == c_len && (loglike_avg[c_len] != max(loglike_avg[-c_len])))){
        c_interval = seq(c_interval[c_len], 5*c_interval[c_len], by = c_interval[c_len]/4)
        c_len = length(c_interval)
        num_searchs = num_searchs + 1
        cat(blue("\nMaximum on the upper boundary, starting new search_"))
      }
      else auto_search_flag = FALSE
    }
  }


  c_min <- c_interval[loglike_max] # c_min refeers to the min of the loss(-loglike)
  final_model <- fit_pml(data_mat, dynseg = dynseg, lambda = c_min,
                       pen_func = pen_func)

  outlist <- list(c_interval = c_interval, c_seg_avg = c_numseg_avg,
                  loss_avg = -loglike_avg, c_min = c_min,
                  final_model = final_model)

  end_time <- proc.time()[3]
  if (!time_omit) cat("\nTime Elapsed: ",end_time-init_time, " seconds\n")
  rm(.Random.seed, envir=globalenv())
  return(outlist)
}

