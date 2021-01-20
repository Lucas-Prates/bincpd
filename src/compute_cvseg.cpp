#include <Rcpp.h>

#include "utils.hpp"
using namespace Rcpp;

//------------------------ Auxiliary Function ----------------------------------

// DP algorithm for the cvseg model
// This function performs the same dynamical programming algorithm of the
// compute_dynseg. However, it only returns the log likelihood for each number
// of segments, so that we can evaluate the best estimator in cvseg.
NumericVector compute_loglike_cvseg(const NumericMatrix& data_mat,
                                    const NumericMatrix& data_mat_test,
                                    const int& nrow_test,
                                    const int& ncol,
                                    const int& segthr){

  NumericVector loss_mat = create_im(ncol), prob_mat = create_im(ncol);
  List ds_matrices(2);
  Function dummy_f("I"); // just to call build_ds_matrices
  // build ds matrices without regularization and returning probs
  ds_matrices = build_ds_matrices(data_mat,
                                  dummy_f,
                                  ncol,
                                  false);  // no penalization
  loss_mat = ds_matrices[0];
  prob_mat = ds_matrices[1];

  // Builds dynamic programming matrices
  // Pseudocode provided at master thesis.
  int offset;
  NumericMatrix zaux(segthr, ncol);
  IntegerMatrix raux(segthr - 1, ncol);
  NumericVector vec(ncol);
  zaux(0, _) = getvector_row_im(loss_mat, ncol, 0);

  if (segthr >= 1){
    for (int k = 1; k <= segthr; k++){
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
    }
  }


  NumericVector loglike_vec(segthr + 1);
  // Retrieve change points and probabilities.
  int segf, segi;

  //    changepoints[k] = ncol;
  for( int k = 0; k <= segthr; k++){
    segf = ncol - 1;
    segi = 0;
    IntegerVector changepoints(k);
    NumericVector prob_vec(k + 1);
    //    p(k, segthr - 1) = prob_mat(segi + 1, segf);
    if (k > 0) {
      NumericVector aux_seg(k), aux_p(k+1);
      for (int j = k; j > 0; j--){
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

    IntegerVector test_col_sum(ncol);
    IntegerVector test_col_samples(ncol);
    for (int i = 0; i < ncol; i++){
      test_col_sum[i] = sum(na_omit(data_mat_test(_, i)));
      test_col_samples[i] = nrow_test - sum(is_na(data_mat_test(_, i)));
    }

    loglike_vec[k] = compute_loglike_cv(test_col_sum,
                                        test_col_samples,
                                        changepoints,
                                        prob_vec,
                                        k,
                                        ncol);
  }

  //---

  return loglike_vec;

}

// Used to fit the final model, forcing the value of k selected by the cv
List compute_final_model(const NumericMatrix& data_mat,
                         const int& ncol,
                         const int& segthr){

  NumericVector loss_mat = create_im(ncol), prob_mat = create_im(ncol);
  List ds_matrices(2);
  Function dummy_f("I"); // just to call build_ds_matrices
  // build ds matrices without regularization and returning probs
  ds_matrices = build_ds_matrices(data_mat,
                                  dummy_f,
                                  ncol,
                                  false, // no penalization
                                  true); // returns probabilities!
  loss_mat = ds_matrices[0];
  prob_mat = ds_matrices[1];

  // Builds dynamic programming matrices
  // Pseudocode provided at master thesis.
  int offset;


  NumericMatrix zaux(segthr, ncol);
  // if segthr = 0, we set the array size to be 0.
  IntegerMatrix raux(std::max(segthr - 1, 0), ncol);
  NumericVector vec(ncol);

  if (segthr >= 1){
    zaux(0, _) = getvector_row_im(loss_mat, ncol, 0);
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
    }
  }

  int best_k = segthr;
  float final_loss;
  IntegerVector changepoints(best_k);
  NumericVector prob_vec(best_k + 1);

  if (best_k > 0){
    final_loss = zaux(best_k, ncol - 1);
    // Retrieve change points and probabilities.
    int segf, segi;

    //    changepoints[k] = ncol;
    segf = ncol - 1;
    segi = 0;

  //    p(k, segthr - 1) = prob_mat(segi + 1, segf);
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
    offset = comp_offset(ncol, 0, ncol - 1);
    final_loss = loss_mat[offset + ncol - 1];
    prob_vec[0] = prob_mat[offset + ncol - 1];
  }

  List model_info(4);
  model_info[0] = changepoints; model_info[1] = prob_vec;
  model_info[2] = final_loss; model_info[3] = best_k;

  return model_info;
}



//------------------------ Main Function ---------------------------------------
// [[Rcpp::export]]
List compute_cvseg_cpp(const NumericMatrix& data_mat,
                       int segthr,
                       const int & n_folds,
                       const int& ncol,
                       const int& nrow){
  if (segthr >= ncol) segthr = ncol - 1;

  NumericVector loglike_cv(segthr + 1);
  for (int i = 0; i <= segthr; i++){ loglike_cv[i] = 0;}

  // shuffled sample numbers for cv folds

  IntegerVector sample_numbers = seq(0, nrow - 1);
  sample_numbers = sample(sample_numbers, nrow,  false);
  //std::random_shuffle(sample_numbers.begin(), sample_numbers.end());

  int nrow_cv, nrow_test;
  int avg_samples_cv = ceil(nrow/n_folds);

  //---
  // Perform cross validation
  for (int k = 0; k < n_folds; k++){
    checkUserInterrupt();
    //---
    // computes the number of rows on the cv data set and test set
    if (k != n_folds - 1){
      nrow_cv = nrow - avg_samples_cv;
    }else{
      nrow_cv = k*avg_samples_cv;
    }
    nrow_test = nrow - nrow_cv;
    //---

    //---
    // builds cv data set, setting each row
    NumericMatrix data_mat_cv(nrow_cv, ncol);
    NumericMatrix data_mat_test(nrow_test, ncol);
    for(int i = 0; i < k*avg_samples_cv; i++){
      data_mat_cv(i, _) = data_mat(sample_numbers[i], _);
    }
    for(int i = (k + 1)*avg_samples_cv; i < nrow; i++){
      data_mat_cv(i-avg_samples_cv, _) = data_mat(sample_numbers[i], _);
    }

    for(int i = k*avg_samples_cv;
        i < std::min((k+1)*avg_samples_cv, nrow);
        i++){

      data_mat_test(i - (k*avg_samples_cv), _) = data_mat(sample_numbers[i], _);

    }
    //---

    // Fit CV model and compute negative log likelihood in test data.
    NumericVector loglike_aux(segthr + 1);
    loglike_aux = compute_loglike_cvseg(data_mat_cv,
                                        data_mat_test,
                                        nrow_test,
                                        ncol,
                                        segthr);

    for(int i = 0; i <= segthr; i++){
      loglike_cv[i] += loglike_aux[i];
    }

  }
  // Finished CV
  //---
  int best_k = 0;
  float min_loglike_cv = loglike_cv[0];
  for(int i = 0; i <= segthr; i++){
    loglike_cv[i] /= n_folds;
    if (loglike_cv[i] < min_loglike_cv){
      min_loglike_cv = loglike_cv[i];
      best_k = i;
    }
  }

  //---
  // Recompute final model
  List model_info = compute_final_model(data_mat,
                                        ncol,
                                        best_k);
  model_info.push_back(loglike_cv);

  return model_info;
}
