#include <Rcpp.h>
using namespace Rcpp;

//-----------------------------Global Variables---------------------------------

// Global variables across program run used to speed up calculations for
// recursive calls, avoid passing unnecessary arguments

Function g_pen_func("I"); // penalization function defined by user
NumericVector g_row_sum; // row sum vector of data set
NumericVector g_row_samples; // number of samples per row in data set
std::vector<int> g_changepoints; // change points detected by method
float g_loss ; // total loss function for the estimated model

//------------------------------------------------------------------------------

//----------------------------Auxiliary functions-------------------------------
// Auxiliary functions using global variables

// Compute mean in the block ranging from left index to right index
float compute_prob(const int& left_index, const int& right_index) {
  int total_sum = 0, total_samples = 0;

  for(int i = left_index - 1; i <= right_index - 1; i++) {
    total_sum += g_row_sum[i];
    total_samples += g_row_samples[i];
  }

  float prob = float(total_sum)/total_samples;
  return prob;
}


float compute_negloglike(const int& left_index,
                         const int& right_index,
                         const int& n,
                         const float& prob) {

  double loglike = 0;
  int block_size = right_index - left_index + 1;
  if( (prob != 0)&(prob != 1) ) {
    loglike = n*block_size*(prob*log(prob) + (1-prob)*log(1-prob));
  }

  return(-loglike);
}

// Computes the model sufficient statistics in order to calculate the loss
// It sets the global row_sum and row_samples to the required values
void compute_suff_stats(const NumericMatrix& data_matrix,
                        const int& n,
                        const int& m) {

  NumericVector row_sum = colSums(data_matrix, true); // second argument
                                                      // removing NA
  NumericVector row_samples(m);
  for(int i = 0; i < m; i++){
    row_samples[i] = n - sum(is_na(data_matrix(_, i)));
  }

  g_row_sum = row_sum;
  g_row_samples = row_samples;
  return;

}

// Recursive implementation of the hierarchical algorithm
// Perform calculations, selection best splitting indexes and appending to
// the global change point set.
void binary_split(const int& left_index,
                  const int& right_index,
                  const int& n,
                  const float& current_loss) {

  Rcpp::checkUserInterrupt(); // check user interruption in rcpp

  if(left_index == right_index) {
    g_loss += current_loss;
    return;
  }

  float left_loss = 0;
  float right_loss = current_loss;
  int split_index = 0;
  float new_left_loss, new_right_loss;
  float left_prob, right_prob;

  for(int j = left_index; j < right_index; j++) {
    left_prob = compute_prob(left_index, j);


    new_left_loss = compute_negloglike(left_index, j, n, left_prob);
    //negloglike + regularization
    new_left_loss = new_left_loss + *REAL(g_pen_func(left_index, j));

    right_prob = compute_prob(j + 1, right_index);

    new_right_loss = compute_negloglike(j + 1, right_index, n, right_prob);
    //negloglike + regularization
    new_right_loss = new_right_loss + *REAL(g_pen_func(j + 1, right_index));

    if(new_left_loss + new_right_loss < left_loss + right_loss){
      left_loss = new_left_loss;
      right_loss = new_right_loss;
      split_index = j;
    }


  }

  if(split_index != 0) {
    g_changepoints.push_back(split_index);
    binary_split(left_index, split_index, n, left_loss);

    binary_split(split_index + 1, right_index, n, right_loss);
  } else {
    g_loss += current_loss;
  }


}

//------------------------------------------------------------------------------


//---------------------------- Main Function -----------------------------------
// Hierarchical segmentation Algorithm cpp
//
// Obtain the solution for the optimization problem using the
// hierarchical algorithm. Does the setup for the calculations and call
// binary_split for the actual computations.
// [[Rcpp::export]]
List compute_hierseg_cpp(const NumericMatrix& data_matrix,
                         const int& n,
                         const int& m,
                         const Function& pen_func) {

  //----
  // assigns global variables
  g_loss = 0;
  g_changepoints.resize(0);
  g_pen_func = pen_func;
  compute_suff_stats(data_matrix, n, m);

  //----
  // Estimation of parameters
  float block_prob = compute_prob(1, m);
  float current_loss = compute_negloglike(1, m, n, block_prob);
  current_loss = current_loss + *REAL(g_pen_func(1, m));

  binary_split(1, m, n, current_loss);

  //----
  // setup output parameters
  std::sort(g_changepoints.begin(), g_changepoints.end());
  int number_of_blocks = g_changepoints.size() + 1;
  NumericVector probs(number_of_blocks);

  // Compute the probabilities of each block estimated
  for(int i = 0; i < number_of_blocks; i++) {
    //Start of each block
    int block_begin, block_end;
    if(i == 0) {
      block_begin = 1;
    } else {
      block_begin = g_changepoints[i-1]+1;
    }

    // End of each block
    if(i == number_of_blocks-1) {
      block_end = m;
    } else {
      block_end = g_changepoints[i];
    }

    probs[i] = compute_prob(block_begin, block_end);


  }
  //----

  List model_info(3);
  model_info[0] = g_changepoints;
  model_info[1] = probs;
  model_info[2] = g_loss;

  return model_info;
}

//------------------------------------------------------------------------------
