# n: sample size
# m: vector size
# prob: probability vector
# cp: changepoint set
bic_loss = function(n, m, p, cp){
  return((1+length(cp))*log(n))
}

sqrt_loss_ds = function(n, m, p, cp){
  return((1+length(cp))*sqrt(n))
}

bic_loss_hs = function(left_index, right_index, n, m){
  loss = log(n)
  return(loss)
}

sqrt_loss_hs = function(left_index, right_index, n, m){
  loss = sqrt(n)
  return(loss)
}

exp_loss = function(left_index, right_index, n, m){
  loss = log(n)*(1+exp(-(right_index - left_index)))
  return(loss)
}

sqrt_n_exp_loss = function(left_index, right_index, n, m){
  loss = sqrt(n)*(1+exp(-(right_index - left_index)))
  return(loss)
}
