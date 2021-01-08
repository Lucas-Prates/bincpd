
build_matrices__ <- function(data_mat){
  data_mat_list = apply(data_mat, 2, na.omit)
  if( typeof(data_mat_list) == "list"){
    data_mat_length = sapply(data_mat_list, length)
    data_mat_sum = sapply(data_mat_list, sum)
  }else{
    if(nrow(data_mat) > 1){
      data_mat_length = apply(data_mat_list, 2, length)
      data_mat_sum = apply(data_mat_list, 2, sum)
    }
    else{
      data_mat_length = rep(1, length(data_mat_list))
      data_mat_sum = data_mat_list
    }
  }
  n = length(data_mat_sum)

  cumsums1 = cumsum(data_mat_sum)
  cumlengths1 = cumsum(data_mat_length)

  cumsummat  = t(outer(cumsums1, c(0, cumsums1[-n]), '-'))
  lengthsmat = t(outer(cumlengths1, c(0, cumlengths1[-n]), '-'))

  p_mat  = cumsummat/lengthsmat#p_mat[i,j] is the mle for p in the interval i:j
  loglike_mat = ifelse (p_mat == 0|p_mat == 1, 0,
                        cumsummat*log(p_mat) + (lengthsmat - cumsummat)*log(1 - p_mat))
  #loglike_mat[i,j] is the likelihood in the interval i:j
  outlist = list(loglike_mat = loglike_mat, p_mat = p_mat)
  return(outlist)
}



###---------------
# Computes the dynamical segmentation,
# obtaining the loglikelihood, cp and probabilities
# as a function of the number of segments
comp_dynseg__ <- function(data, ncol, segthr = NULL){

  if(is.null(segthr)) {segthr = ncol(data)}
  data_bm = build_matrices__(data)
  data_ds = comp_dynseg_cpp__(data_bm$loglike_mat, data_bm$p_mat, segthr = segthr)
  changepoints = list(); prob_list = list()
  n = nrow(data_ds[[2]]);
  for(i in 1:n){
    if(i != 1) changepoints[[i]] = data_ds[[2]][i, (n - i + 1):(n-1)]
    prob_list[[i]] = data_ds[[3]][i, (n - i + 1):n]
  }
  outlist = list( loglike =  data_ds[[1]], changepoints = changepoints,
                  prob_list = prob_list)
  return(outlist)

}


#' Segments the data using dynamic programming algorithm
#'
#' Uses dynamical segmentation to obtain the
#'log likelihood, changepoints and probabilities
#' as a function of the number of segments
#'
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for perfomance issues on large datasets
comp_dynseg__fast__ <- function(data, ncol = NULL, segthr = NULL){
  if(is.null(ncol)) { ncol = ncol(data)}
  if(is.null(segthr)) {segthr = ncol}
  data_bm = build_matrices_cpp__(data)
  data_ds = comp_dynseg_cpp_fast__(data_bm[[1]], data_bm[[2]],
                                   ncol, segthr = segthr)

  changepoints = list(); prob_list = list()
  outlist = list( loglike =  data_ds[[1]], changepoints = data_ds[[2]],
                  prob_list = data_ds[[3]])
  return(outlist)

}

#' Computes the penalized maximum likelihood estimator
#'
#' The estimator is the argmin of -l(p; X) + lambda*R(p;X), where l is the
#' log likelihood and R is the penalization function
#'
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for perfomance issues on large datasets
#' @param lambda Penalization constant
#' @param pen_func The penalization function used for the computation of the estimator.
#' A user specified function can be provided, and the function signature is
#' f(n, m, prob, cp)
#' The default loss is the bic_loss -log(n)log(m)(len(cp)+1).
fit_pml <- function(data, segthr = NULL, lambda = 1, pen_func = NULL,
                    dynseg = NULL){

  n_samples = nrow(data); n_col = ncol(data)
  if(is.null(dynseg)){
    dynseg = comp_dynseg__fast__(data, n_col, segthr)
  }
  segthr = length(dynseg$loglike)
  if(is.null(pen_func)){
    # equivalent to bic_loss on loss_functions.R
    pml <- -1*(dynseg$loglike - lambda*((1:segthr))*log(n_samples))
  }
  else{
    pml = -dynseg$loglike
    for(i in 1:segthr){
      pml[i] = pml[i] + lambda*pen_func(n_samples, n_col,
                                  dynseg$prob_list[[i]],
                                  dynseg$changepoints[[i]])
    }
  }

  pml_numseg = which.min(pml)
  pml_cp = dynseg$changepoints[[pml_numseg]]
  pml_block_begin = 1 + c(0, pml_cp)
  pml_block_end = c(pml_cp, n_col)
  pml_prob = dynseg$prob_list[[pml_numseg]]

  pml_metadata = list(n_samples = n_samples,
                      n_col = n_col,
                      colnames = colnames(data),
                      method = 'pml')

  outlist = list(metadata = pml_metadata,
                 numseg = pml_numseg,
                 changepoints = pml_cp,
                 block_begin = pml_block_begin,
                 block_end = pml_block_end,
                 prob = pml_prob, pml = pml)

  return(outlist)

}
