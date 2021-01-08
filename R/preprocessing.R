### The paper "Detecting autozygosity through runs of homozygosity:
###A comparison of three autozygosity detection algorithms" uses
### GWAS (Genome Wide Association Studies) cleaning procedures,
### which consists of:
### 1) Dropping SNPs with missingness rate above 2%;
### 2) Dropping individuals with SNPs missing rate above 5%;
### 3) Dropping SNPs with MAF below 1%;
### 4) Dropping SNPs out of HWE (Hardy Weinberg Equilibrium)
###    with p-squared < 0.0001

# Auxiliar function which handles individuals missing rate and snps missing rate
handle_missing__ <- function(snps_data, thr_mis_snps = 0.01, thr_mis_ind = 0.05){
  n <- nrow(snps_data)
  aux_missing <- (apply(is.na(snps_data), 2, sum)/n) > thr_mis_snps
  new_snps_data <- snps_data
  if(any(aux_missing)){
    new_snps_data <- new_snps_data[, - which( aux_missing)]
  }
  cat("Removing", sum(aux_missing), "SNPs that had missing rates above the",
      thr_mis_snps, "threshold.\n", sep = ' ')


  # m <- ncol(new_snps_data)
  # aux_missing <- (apply(is.na(new_snps_data), 1, sum)/m) > thr_mis_ind
  # if(any(aux_missing)){
  #   new_snps_data <- new_snps_data[- which( aux_missing), ]
  # }
  return(new_snps_data)
}


# Auxiliar function that process snps before applying
# the changepoint detection algorithmn
# the cutoff is related to the missing rate for the sequencing
smooth_snps__ <- function(snps, radius = 10, cutoff = 0.95){

  len_snps <- length(snps)

  smoothed_snps <- rep(0, len_snps)
  for(i in 1:radius){
    smoothed_snps[i] <- 1*(mean(snps[i:(2*radius + i - 1)],
                                na.rm = TRUE) >= cutoff)
    k = len_snps - i +1
    smoothed_snps[k] <- 1*(mean(snps[(k - 2*radius + 1):k],
                                na.rm = TRUE) >= cutoff)
  }
  for(i in (radius+1):(len_snps - radius)){
    smoothed_snps[i] <- 1*(mean(snps[(i - radius):(i + radius)],
                             na.rm = TRUE) >= cutoff)
  }

  return(smoothed_snps)
}


#' Data prep for ROH island analysis
#' Description: Process SNPs data so it is ready for model fitting.
#' The preprocessing consists of computing the mean in a neighboorhood
#' of the specified radius (total size = 2*radius + 1)
#' and checking if it is above a cutoff (default 95%).
#' If that is the case, we output 1 to that row and SNP position in a new
#' dataframe with the same dimensions as the original, otherwise we output 0.
#' @export
preprocess <- function(snps_raw, radius = 15,
                       thr_mis_snps = 0.01, thr_mis_ind = 0.05){

  snps_mat = data.matrix(snps_raw)
  snps_mat[snps_mat == 2] = 0
  snps_mat = 1 - snps_mat  ## now 1 corresponds to homozygous SNP.

  # Not handling missing yet
  snps_mat = handle_missing__(snps_mat, thr_mis_snps, thr_mis_ind)
  processed_snps <- t(apply(X = snps_mat, MARGIN = 1, FUN = smooth_snps__,
                            radius = radius))
  colnames(processed_snps) = sub('\\_[0-9]', '', colnames(snps_mat))
  return(processed_snps)
}

# Auxiliar function that computes the cp_begin and cp_end of blocks,
# the blocks lengths and probabilities
summary_snps_seg <- function(snps_cp, snps_prob, m, min_size = 1){

  model = list()
  model$snps_prob = rep(NA, m) #-1 represents that the snps is in a block whose length is smaller than min
  model$cp_begin = c(1, snps_cp + 1)
  model$cp_end = c(snps_cp, m)
  model$block_length = model$cp_end - model$cp_begin + 1
  model$block_prob = snps_prob

  valid.intervals = (model$block_length >= min_size)

  model$cp_begin = model$cp_begin[valid.intervals]
  model$cp_end = model$cp_end[valid.intervals]
  model$block_length = model$block_length[valid.intervals]
  model$block_prob = model$block_prob[valid.intervals]
  for(i in 1:length(model$cp_begin)){
    model$snps_prob[model$cp_begin[i]:model$cp_end[i]] = model$block_prob[i]
  }

  model
}

# Receives changepoints of a model, returning the roh islands,
# number of islands and total of snps contained
# in the islands
# if the "n_isl" parameter is provided it will override the cutoff parameter
detect_roh_isl <- function(snps_cp, snps_prob, m, min_size, cutoff, n_isl = NULL){

  snps_model = summary_snps_seg(snps_cp, snps_prob, m, min_size = min_size)
  snps_nisl = length(snps_model$cp_begin)
  snps_isl = rep(0, m)
  if(!is.null(n_isl)) snps_isl_ind = which(snps_model$prob >= sort(snps_model$block_prob, decreasing = T)[n_isl])
  else snps_isl_ind = which(snps_model$block_prob >= quantile(snps_model$prob, cutoff))
  num_isl = length(snps_isl_ind)
  for(i in snps_isl_ind){
    snps_isl[snps_model$cp_begin[i]:snps_model$cp_end[i]] = 1
  }
  outlist = list()
  outlist$isl = snps_isl; outlist$num_isl = num_isl;
  outlist$tot_snps = sum(snps_isl);
  outlist
}

