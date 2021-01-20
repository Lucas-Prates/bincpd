#include <Rcpp.h>
#include "utils.hpp"
using namespace Rcpp;

//--------------------------- Main Function ------------------------------------

// Dynamical programming Algorithm
//
// Fast implementation of the dynamic segmentation using the Incomplete Matrix
// and returning list of lists instead of matrices
// [[Rcpp::export]]
List compute_dynseg_cpp(const NumericMatrix& data_mat,
                        const int& ncol,
                        int segthr,
                        const Function& pen_func){
  if (segthr >= ncol){segthr = ncol - 1;}

  NumericVector loss_mat = create_im(ncol), prob_mat = create_im(ncol);
  List ds_matrices(2);

  ds_matrices = build_ds_matrices(data_mat, pen_func, ncol);

  loss_mat = ds_matrices[0];
  prob_mat = ds_matrices[1];

  // Builds dynamic programming matrices
  // Pseudocode provided at master thesis.
  int offset;
  NumericMatrix zaux(segthr, ncol);
  IntegerMatrix raux(segthr - 1, ncol);
  NumericVector vec(ncol);
  zaux(0, _) = getvector_row_im(loss_mat, ncol, 0);

  int best_k = 0; // best number of change points
  float min_loss = zaux(0, ncol - 1); // loss for best_k

  if (segthr >= 1){
    for (int k = 1; k <= segthr; k++){
      checkUserInterrupt();
      for (int j = k; j < ncol; j++){
        //vec = vector(length = j - k + 1)
        for (int i = k; i <= j; i++){
          offset = comp_offset(ncol, i, j);
          vec[i - k] = zaux(k - 1, i - 1) + loss_mat[offset + j - i];
        }
        zaux(k, j) = partial_min(vec, j - k + 1); // actual size of vector
                                                  // is j - k + 1
        raux(k - 1, j) = partial_indmin(vec, j - k + 1) + k - 1;
      }
      if (min_loss > zaux(k, ncol - 1)){
        min_loss = zaux(k, ncol - 1);
        best_k = k;
      }
    }
  }

  // Retrieve change points and probabilities.

  IntegerVector changepoints(best_k);
  NumericVector prob_vec(best_k + 1);
  int segf, segi;

//    changepoints[k] = ncol;
  segf = ncol - 1;
  segi = 0;

  //    p(k, segthr - 1) = prob_mat(segi + 1, segf);
  if (best_k > 0) {
    NumericVector aux_seg(best_k), aux_p(best_k+1);
    for (int j = best_k; j > 0; j--){
      segi = raux(j - 1, segf);
      offset = comp_offset(ncol, segi + 1, segf);
      prob_vec[j] = prob_mat[offset + segf -(segi + 1)];
      changepoints[j-1] = segi + 1;
      segf = segi;
    }
    offset = comp_offset(ncol, 0, segf);
    prob_vec[0] = prob_mat[offset + segf];
  }
  else{
    offset = comp_offset(ncol, 0, ncol -1);
    prob_vec[0] = prob_mat[offset + ncol - 1];
  }

  List model_info(4);
  model_info[0] = changepoints; model_info[1] = prob_vec;
  model_info[2] = min_loss; model_info[3] = best_k;

  return model_info;
}
