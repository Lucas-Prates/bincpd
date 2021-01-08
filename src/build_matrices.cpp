#include <Rcpp.h>

using namespace Rcpp;

NumericVector create_im(int ncol){
  int total_size = ncol*(ncol+1)/2;
  NumericVector x(total_size);
  return(x);
}

int comp_offset__(int ncol, int i, int j){
  if(j < i){
    std::range_error("The column index j cannot be smaller than the row index i!");
  }
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  return(offset);
}

// [[Rcpp::export]]
List build_matrices_cpp__(NumericMatrix data_mat){
  int n = data_mat.nrow();
  int m = data_mat.ncol();
  NumericVector ll_mat = create_im(m), prob_mat = create_im(m);
  NumericVector cs_obs(m), cs_len(m); // cumulative sum of observations and len
  cs_obs[0] = sum(na_omit(data_mat(_, 0)));
  cs_len[0] = na_omit(data_mat(_, 0)).size();
  for(int i = 1;i < m;i++){
    cs_obs[i] = cs_obs[i-1] + sum(na_omit(data_mat(_, i)));
    cs_len[i] = cs_len[i-1] + na_omit(data_mat(_, i)).size();
  }

  int offset;
  double aux_obs, aux_len, prob, ll;
  for(int i = 0; i < m; i++){
    for(int j = i; j < m; j++){
      if(i == 0){
        aux_obs = cs_obs[j];
        aux_len = cs_len[j];
      }
      else{
        aux_obs = cs_obs[j] - cs_obs[i - 1];
        aux_len = cs_len[j] - cs_len[i - 1];
      }
      offset = comp_offset__(m, i, j);
      prob = aux_obs/aux_len;
      prob_mat[offset + j - i] = prob;
      if((prob == 0) | (prob == 1)){ ll_mat[offset + j - i] = 0;}
      else{
        ll = aux_obs*log(prob) + (aux_len - aux_obs)*log(1-prob);
        ll_mat[offset + j - i] = ll;
      }
//      prob_mat[i] = aux_prob;
//      ll_mat[i] = aux_ll;
    }
  }

  List outlist(2);
  outlist[0] = ll_mat; outlist[1] = prob_mat;
  return(outlist);
}



