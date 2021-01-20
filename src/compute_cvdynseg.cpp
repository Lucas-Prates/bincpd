#include <Rcpp.h>

#include "utils.hpp"
#include "compute_dynseg.h"
using namespace Rcpp;


//------------------------ Main Function ---------------------------------------
// The estimator for the k-fold cross validation on the penalization set. We use
// dynamical programming to estimate the models for each lambda, and then refit
// the model in the whole data set for the best lambda.
// [[Rcpp::export]]
List compute_cvdynseg_cpp(const NumericMatrix& data_mat,
                          int segthr,
                          const NumericVector& lambda_set,
                          const int & n_folds,
                          const int& ncol,
                          const int& nrow,
                          const Function& create_pen_func){
  if (segthr >= ncol) segthr = ncol - 1;

  int len_lambda_set = lambda_set.size();
  NumericVector loglike_cv(len_lambda_set);
  for (int i = 0; i < len_lambda_set; i++){ loglike_cv[i] = 0;}

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
    List cv_model(4);
    for (int i = 0; i < len_lambda_set; i++){

      // we have to build the proper penalization function for the constant!
      Function pen_func = create_pen_func(lambda_set[i], nrow_cv, ncol);

      // Fit CV model and compute negative log likelihood in test data.
      cv_model = compute_dynseg_cpp(data_mat_cv,
                                    ncol,
                                    segthr,
                                    pen_func);

      IntegerVector test_col_sum(ncol);
      IntegerVector test_col_samples(ncol);
      for (int i = 0; i < ncol; i++){
        test_col_sum[i] = sum(na_omit(data_mat_test(_, i)));
        test_col_samples[i] = nrow_test - sum(is_na(data_mat_test(_, i)));
      }

      loglike_cv[i] += compute_loglike_cv(test_col_sum,
                                          test_col_samples,
                                          cv_model[0],
                                          cv_model[1],
                                          cv_model[3],
                                          ncol);

    }

  }
  // Finished CV
  //---

  float best_lambda = lambda_set[0];
  float min_loglike_cv = loglike_cv[0];
  for(int i = 0; i < len_lambda_set; i++){
    loglike_cv[i] /= n_folds;
    if (loglike_cv[i] < min_loglike_cv){
      min_loglike_cv = loglike_cv[i];
      best_lambda = lambda_set[i];
    }
  }

  //---
  // Recompute final model
  Function pen_func = create_pen_func(best_lambda, nrow_cv, ncol);
  List model_info = compute_dynseg_cpp(data_mat,
                                       ncol,
                                       segthr,
                                       pen_func);
  model_info.push_back(best_lambda);
  model_info.push_back(loglike_cv);

  return model_info;
}
