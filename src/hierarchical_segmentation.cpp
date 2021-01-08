#include <Rcpp.h>

using namespace Rcpp;


float __compute_prob(const int& left_index, const int& right_index,
                      const NumericVector& row_sum,
                      const NumericVector& row_samples){
  int total_sum = 0, total_samples = 0;

  for(int i = left_index - 1; i <= right_index - 1; i++){
    total_sum += row_sum[i];
    total_samples += row_samples[i];
  }

  float prob = float(total_sum)/total_samples;
  return prob;
}

float __compute_loss(const int& left_index, const int& right_index,
                     const int& n,
                     const float& prob,
                     const Function& pen_func){

  double loglike = 0;
  int block_size = right_index - left_index + 1;
  if( (prob != 0)&(prob != 1) ){
    loglike = n*block_size*(prob*log(prob) + (1-prob)*log(1-prob));
  }

  NumericVector pen_value = pen_func(left_index, right_index);

  double loss = -loglike + pen_value[0];
  return(loss);
}

// Computes the model sufficient statistics in order to calculate the loss
List __compute_suff_stats(const NumericMatrix& data_matrix,
                          const int& n, const int& m){

  NumericVector row_sum = colSums(data_matrix, true); // second argument
                                                      // removing NA
  NumericVector row_samples(m);
  for(int i = 0; i<m; i++){
    row_samples[i] = n - sum(is_na(data_matrix(_, i)));
  }

  List sufficient_stats(2);
  sufficient_stats[0] = row_sum; sufficient_stats[1] = row_samples;
  return sufficient_stats;

}

void __binary_split(const int& left_index, const int& right_index,
                    const int& n,
                    const float& current_loss,
                    const NumericVector& row_sum,
                    const NumericVector& row_samples,
                    const Function& pen_func,
                    std::vector<int>& changepoints,
                    float& loss){
  if(left_index == right_index){
    loss += current_loss;
    return;
  }

  float left_loss = 0;
  float right_loss = current_loss;
  int split_index = 0;
  float new_left_loss, new_right_loss;
  float left_prob, right_prob;

  for(int j = left_index; j < right_index; j++){
    left_prob = __compute_prob(left_index, j,
                                     row_sum,
                                     row_samples);
    new_left_loss = __compute_loss(left_index, j, n,
                                   left_prob, pen_func);

    right_prob = __compute_prob(j + 1, right_index,
                               row_sum,
                               row_samples);
    new_right_loss = __compute_loss(j + 1, right_index, n,
                                    right_prob, pen_func);

    if(new_left_loss + new_right_loss < left_loss + right_loss){
      left_loss = new_left_loss;
      right_loss = new_right_loss;
      split_index = j;
    }


  }

  if(split_index != 0){
    changepoints.push_back(split_index);
    __binary_split(left_index, split_index, n, left_loss,
                   row_sum, row_samples, pen_func,
                   changepoints, loss);
    __binary_split(split_index + 1, right_index, n, right_loss,
                   row_sum, row_samples, pen_func,
                   changepoints, loss);
  }else{
    loss += current_loss;
  }


}

// Obtain the solution for the optimization problem using the
// hierarchical algorithm.
// [[Rcpp::export]]
List hierarchical_algorithm_cpp( const NumericMatrix& data_matrix,
                                             const int& n,
                                             const int& m,
                                             const Function& pen_func){

  std::vector<int> changepoints;
  float loss = 0;

  List suff_stats = __compute_suff_stats(data_matrix, n, m);
  NumericVector row_sum = suff_stats[0];
  NumericVector row_samples = suff_stats[1];
  float block_prob = __compute_prob(1, m, row_sum, row_samples);
  float current_loss = __compute_loss(1, m, n,
                                      block_prob,
                                      pen_func);

  __binary_split(1, m, n,
                 current_loss,
                 row_sum, row_samples,
                 pen_func,
                 changepoints, loss);


  std::sort(changepoints.begin(), changepoints.end());

  int number_of_blocks = changepoints.size() + 1;
  NumericVector probs(number_of_blocks);

  // Compute the probabilities of each block estimated
  for(int i = 0; i < number_of_blocks; i++){
    //Start of each block
    int block_begin, block_end;
    if(i == 0){
      block_begin = 1;
    }else{
      block_begin = changepoints[i-1]+1;
    }

    //End of each block
    if(i == number_of_blocks-1){
      block_end = m;
    }else{
      block_end = changepoints[i];
    }

    probs[i] = __compute_prob(block_begin, block_end, row_sum, row_samples);


  }

  List output_list(3);
  output_list[0] = changepoints;
  output_list[1] = probs;
  output_list[2] = loss;

  return output_list;
}

