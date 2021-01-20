#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
using namespace Rcpp;

//----------------------- Incomplete Matrix ------------------------------------
// Creates a data structure which represents an upper triangular matrix.
// It consists of a vector of the correct with methods to access the indexes
// using double notation as with a matrix.
// TODO: Reformulate this into a class!

inline NumericVector create_im(int ncol){
  int total_size = ncol*(ncol+1)/2;
  NumericVector x(total_size);
  return(x);
}

inline int comp_offset(const int& ncol, const int& i, const int& j){
  if(j < i){
    std::range_error("The column index j cannot be smaller \
                      than the row index i!");
  }
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  return(offset);
}

inline double getvalue_im(const NumericVector& inc_mat,
                   const int& ncol,
                   const int& i,
                   const int& j){
  if(j < i){
    std::range_error("The column index j cannot be smaller \
                      than the row index i!");
  }
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  return( inc_mat[offset + j - i] );
}

inline void setvector_row_im(NumericVector inc_mat,
                      const NumericVector& vec,
                      const int& ncol,
                      const int& i) {
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  for(int k = 0; k < ncol - i; k++){
    inc_mat[offset + k] = vec[k];
  }
}

inline NumericVector getvector_row_im(const NumericVector& inc_mat,
                               const int& ncol,
                               const int& i){
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  NumericVector row_vec(ncol-i);
  for(int k = 0; k < ncol - i; k++){
    row_vec[k] = inc_mat[offset + k];
  }
  return(row_vec);
}

//----------------------- Utility functions ------------------------------------

// Function to retrieve values from Incomplete Matrix structure
// Given the indexes from a normal Matrix, it retrieves the value for the IM

//calculates the minimum of a vector before the cutoff
inline double partial_min(const NumericVector& vec, const int& cutoff){
  double min_vec = vec[0];
  for (int i = 1; i < cutoff; i++) if(min_vec > vec[i]) min_vec = vec[i];
  return min_vec;
}

//calculates the index of the minimum of a vector before the cutoff
inline int partial_indmin(const NumericVector& vec, const int& cutoff){
  double min_vec = vec[0];
  int indmin = 0;
  for (int i = 1; i < cutoff; i++){
    if (min_vec > vec[i]){
      min_vec = vec[i];
      indmin = i;
    }
  }
  return indmin;
}

// Build DP recursion matrix
//
// Given the data matrix, the function constructs the incomplete matrices of
// the loss function and the probabilities. The element [i, j], where i <= j
// always holds, is the loss and empiric probability of the block i:j
// These matrices are essential for the Dynamical Programming algorithm and
// for the CVSEG function.
inline List build_ds_matrices(const NumericMatrix data_mat,
                              const Function& pen_func,
                              const int& ncol,
                              const bool& flag_reg = true,
                              const bool& flag_return_prob = true){

  int m = ncol;
  NumericVector loss_mat = create_im(m), prob_mat = create_im(m);
  NumericVector cs_obs(m), cs_len(m); // cumulative sum of observations and len
  cs_obs[0] = sum(na_omit(data_mat(_, 0)));
  cs_len[0] = na_omit(data_mat(_, 0)).size();
  for(int i = 1; i < m; i++){
    cs_obs[i] = cs_obs[i-1] + sum(na_omit(data_mat(_, i)));
    cs_len[i] = cs_len[i-1] + na_omit(data_mat(_, i)).size();
  }

  int offset;
  double aux_obs, aux_len, prob, loglike;
  for(int i = 0; i < m; i++){
    Rcpp::checkUserInterrupt();
    for(int j = i; j < m; j++){
      if(i == 0){
        aux_obs = cs_obs[j];
        aux_len = cs_len[j];
      }
      else{
        aux_obs = cs_obs[j] - cs_obs[i - 1];
        aux_len = cs_len[j] - cs_len[i - 1];
      }
      offset = comp_offset(m, i, j);
      prob = aux_obs/aux_len;
      prob_mat[offset + j - i] = prob;
      if((prob == 0) | (prob == 1)){ loss_mat[offset + j - i] = 0;}
      else{
        loglike = aux_obs*log(prob) + (aux_len - aux_obs)*log(1-prob);
        loss_mat[offset + j - i] = -loglike;
      }

      // regularization term
      if (flag_reg){
        loss_mat[offset + j - i] += *REAL(pen_func(i, j));
      }

    }
  }


  if (flag_return_prob){
    List outlist(2);
    outlist[0] = loss_mat; outlist[1] = prob_mat;
    return outlist;
  }else{
    List outlist(1);
    outlist[0] = loss_mat;
    return outlist;
  }

}

// Evaluates the log likelihood in the test data set for the parameters
// estimated in the cv training set.
inline float compute_loglike_cv(const IntegerVector& test_col_sum,
                                const IntegerVector& test_col_samples,
                                IntegerVector changepoints,
                                const NumericVector& probabilities,
                                const int& k,
                                const int& ncol){
  //Maximum likelihood probabilities
  int block_sum;
  int n_samples;
  //auxiliary cp
  changepoints.push_front(0);
  changepoints.push_back(ncol);
  float loglike_test = 0;
  for(int i = 0 ; i < (k + 1); i++){
    block_sum = 0;
    n_samples = 0;
    // if prob[i] = 0 && 1, then the loglike would be infinite. Should it be
    // allowed? For small samples, this could be a problem.
    if ((probabilities[i] != 0)&&(probabilities[i] != 1)){
      for(int c = (changepoints[i]); c <= changepoints[i + 1] - 1; c++){
        block_sum += test_col_sum[c];
        n_samples += test_col_samples[c];
      }
      loglike_test += block_sum*log(probabilities[i]);
      loglike_test += (n_samples - block_sum)*log(1 - probabilities[i]);
    }

  }

  return -loglike_test;
}

#endif
