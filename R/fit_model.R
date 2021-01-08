#'Fits a binary CPD model according to the method chosen
#'
#'
#'@export
fit_model <- function(data, method = "hierseg", arguments = list()) {

  IMPLEMENTED_METHODS <- c("hierseg", "dynseg", "cvseg",
                           "cvdynseg", "fusedlasso")

  if ( !(method %in% IMPLEMENTED_METHODS) ) {
    stop("The 'method' argument provided is not implemented!")
  }

  methodcall_name <- paste0("fit_", method) # package function name
  fit_arguments   <- c(list(data_matrix = data.matrix(data)), arguments)
  model <- do.call(methodcall_name, fit_arguments)

  class(model) <- "bincpd"
  return(model)

}
