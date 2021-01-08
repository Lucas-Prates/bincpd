#' Segments the data using dynamic programming algorithm
#'
#' Uses dynamical segmentation to obtain the
#'log likelihood, changepoints and probabilities
#' as a function of the number of segments
#'
#' @param segthr Threshold on the number of block segments to fit the model.
#' This is highly recommend for perfomance issues on large datasets
compute_dynseg <- function(data,
                           segthr = NULL,
                           lambda = 1,
                           pen_func = bic_loss_hs) {

  ncol = ncol(data)
  if(is.null(segthr)) {segthr = ncol}
  data_bm = build_matrices_cpp__(data)
  data_ds = comp_dynseg_cpp_fast__(data_bm[[1]], data_bm[[2]],
                                   ncol, segthr = segthr)

  changepoints = list(); prob_list = list()
  outlist = list( loglike =  data_ds[[1]], changepoints = data_ds[[2]],
                  prob_list = data_ds[[3]])
  return(outlist)

}
