#' @title
#' Compare or evaluate model performance with respect
#' to other model or ground truth
#'
#' @description
#'
#' Compares or evaluates model estimated change point set against another model
#' or ground truth. The comparison is made using common metrics to compare
#' clusters. The metrics provided are
#'\itemize{
#'  \item["hausdorff"}] Hausdorff Distance metric (default);
#'  \item["rand"] Rand Index metric;
#'  \item["symdiff"] Symmetric difference metric.
#' }
#'
#' @param model1 The first bincpd object or list of sorted integers
#' representing the change point set.
#' @param model2 The second bincpd object or list of sorted integers
#' representing the change point set.
#' @param m The number of variables which the model was fitted on. Only needs to
#' be passed if both arguments are change point sets instead of a bincpd object.
#' @export

compare_model <- function(model1,
                          model2,
                          m = NULL){

  bincpd_flag = FALSE # flag for bincpd object

  # Model 1 setup
  if (class(model1) == "bincpd"){
    m = model1$metadata$m
    cp1 = model1$changepoints
    bincpd_flag = TRUE
  } else {cp1 = model1} # cp1 must be a change point set

  # Model 2 setup
  if (class(model2) == "bincpd"){
    m = model2$metadata$m
    cp2 = model2$changepoints
    bincpd_flag = TRUE
  } else {cp2 = model2} # cp1 must be a change point set


  if (!bincpd_flag){
    if (is.null(m)){
      stop("Error! No bincpd models were passed and m was not provided!")
    }
  }

  # Metrics computations using Rcpp
  haus = compute_hausdorff(cp1, cp2)
  rand = compute_rand(cp1, cp2, m)
  symdiff = compute_symdiff(cp1, cp2)

  metrics_list = list(haus = haus,
                      rand = rand,
                      symdiff = symdiff)

  return(metrics_list)

}
