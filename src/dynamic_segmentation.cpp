#include <Rcpp.h>
using namespace Rcpp;

// Function to retrieve values from Incomplete Matrix structure
// Given the indexes from a normal Matrix, it retrieves the value for the IM

int comp_offset(int ncol, int i, int j){
  if(j < i){
    std::range_error("The column index j cannot be smaller than the row index i!");
  }
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  return(offset);
}

double getvalue_im(NumericVector inc_mat, int ncol, int i, int j){
  if(j < i){
    std::range_error("The column index j cannot be smaller than the row index i!");
  }
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  return( inc_mat[offset + j - i] );
}

void setvector_row_im(NumericVector inc_mat, NumericVector vec, int ncol, int i){
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  for(int k = 0; k < ncol - i; k++){
    inc_mat[offset + k] = vec[k];
  }
}

NumericVector getvector_row_im(NumericVector inc_mat, int ncol, int i){
  int offset = ncol*(ncol+1)/2 - (ncol-i)*((ncol-i)+1)/2;
  NumericVector row_vec(ncol-i);
  for(int k = 0; k < ncol - i; k++){
    row_vec[k] = inc_mat[offset + k];
  }
  return(row_vec);
}


//calculates the maximum of a vector before the cutoff
double partial_max(NumericVector vec, int cutoff){
  double max_vec = vec[0];
  for (int i = 1; i < cutoff; i++) if(max_vec < vec[i]) max_vec = vec[i];
  return max_vec;
}
//calculates the index of the maximum of a vector before the cutoff
int partial_indmax(NumericVector vec, int cutoff){
  double max_vec = vec[0];
  int indmax = 0;
  for (int i = 1; i < cutoff; i++){
    if (max_vec < vec[i]){
      max_vec = vec[i];
      indmax = i;
    }
  }
  return indmax;
}

// Dynamical Segmentation Algorithm to maximize the likelihood function
// [[Rcpp::export]]
List comp_dynseg_cpp__(NumericMatrix loglike_mat, NumericMatrix p_mat, int segthr = 100){
  int n = loglike_mat.ncol();
  if (segthr > n) segthr = n;
  NumericMatrix zaux(segthr, n), raux(segthr - 1, n);
  NumericVector vec(n);
  zaux(0, _) = loglike_mat(0, _);
  if (segthr > 1){
    for (int k = 1; k < segthr; k++){
      for (int j = k; j < n; j++){
        //vec = vector(length = j - k + 1)
        for (int i = k; i <= j; i++) vec[i - k] = zaux(k - 1, i - 1) + loglike_mat(i, j);
        zaux(k, j) = partial_max(vec, j - k + 1); //actual size of vec is j - k + 1
          raux(k - 1, j) = partial_indmax(vec, j - k + 1) + k - 1;
      }
    }
  }

  NumericVector loglike_vec = zaux(_, n - 1);
  NumericMatrix segments(segthr, segthr);
  NumericMatrix p(segthr, segthr);
  int segf, segi;

  for (int k = 0; k < segthr; k++){
    segments(k, segthr - 1) = n;
    segf = n - 1;
    segi = 0;
//    p(k, segthr - 1) = p_mat(segi + 1, segf);
    if (k > 0) {
      for (int j = k; j > 0; j--){
        // segf = segments[[k]][length(segments[[k]])]
        segi = raux(j - 1, segf);
        p(k, segthr - 1 - (k - j) ) = p_mat(segi + 1, segf);
        segments(k, segthr - 1 - (k + 1 - j) ) = segi + 1;
        //The +1 is a correction because the indexes in cpp start with 0
        // p[[k]] = c(0,p[[k]])
        segf = segi;
      }
      p(k, segthr - 1 - k) = p_mat(0, segf);
      // p[[k]] = c(p.fit[1,segf],p[[k]])
    }
    else p(k, segthr - 1) = p_mat(0, n - 1);
  }
   List outlist(3);
    outlist[0] = loglike_vec; outlist[1] = segments; outlist[2] = p;
    return outlist;
}

// Fast implementation of the dynamic segmentation using the Incomplete Matrix
// and returning list of lists instead of matrices
// [[Rcpp::export]]
List comp_dynseg_cpp_fast__(NumericVector loglike_mat, NumericVector p_mat,
                            int ncol, int segthr = 100){
  if (segthr > ncol) segthr = ncol;
  int offset;
  NumericMatrix zaux(segthr, ncol), raux(segthr - 1, ncol);
  NumericVector vec(ncol);
  zaux(0, _) = getvector_row_im(loglike_mat, ncol, 0);
  if (segthr > 1){
    for (int k = 1; k < segthr; k++){
      for (int j = k; j < ncol; j++){
        //vec = vector(length = j - k + 1)
        for (int i = k; i <= j; i++){
          offset = comp_offset(ncol, i, j);
          vec[i - k] = zaux(k - 1, i - 1) + loglike_mat[offset + j - i];
        }
        zaux(k, j) = partial_max(vec, j - k + 1); //actual size of vec is j - k + 1
        raux(k - 1, j) = partial_indmax(vec, j - k + 1) + k - 1;
      }
    }
  }

  NumericVector loglike_vec = zaux(_, ncol - 1);
  List segments(segthr), p(segthr);
  int segf, segi;

  for (int k = 0; k < segthr; k++){
//    segments[k] = ncol;
    segf = ncol - 1;
    segi = 0;
    //    p(k, segthr - 1) = p_mat(segi + 1, segf);
    if (k > 0) {
      NumericVector aux_seg(k), aux_p(k+1);
      for (int j = k; j > 0; j--){
        // segf = segments[[k]][length(segments[[k]])]
        segi = raux(j - 1, segf);
//        p(k, segthr - 1 - (k - j) ) = getvalue_im( p_mat, ncol, segi + 1, segf);
        offset = comp_offset(ncol, segi + 1, segf);
        aux_p[j] = p_mat[offset + segf -(segi + 1)];
        aux_seg[j-1] = segi + 1;
//        segments(k, segthr - 1 - (k + 1 - j) ) = segi + 1;
        //The +1 is a correction because the indexes in cpp start with 0
        // p[[k]] = c(0,p[[k]])
        segf = segi;
      }
      offset = comp_offset(ncol, 0, segf);
      aux_p[0] = p_mat[offset + segf];
      p[k] = aux_p;
      segments[k] = aux_seg;
//      p(k, segthr - 1 - k) = getvalue_im(p_mat, ncol, 0, segf);
      // p[[k]] = c(p.fit[1,segf],p[[k]])
    }
    else{
      offset = comp_offset(ncol, 0, ncol -1);
      p[k] = p_mat[offset + ncol - 1];
    }
  }
  List outlist(3);
  outlist[0] = loglike_vec; outlist[1] = segments; outlist[2] = p;
  return outlist;
}
