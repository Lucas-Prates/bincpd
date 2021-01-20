bincpd: a package for binary change point detection
================

The `bincpd` package implements statistical models based on regularized
likelihood to perform offline change point detection on binary data. In
most change point problems, we usually work with a single time series
which we wish to segment. In `bincpd`, however, we are interested in
working with several time series data sampled from the same
distribution.

This vignette provides a brief explanation of the statistical model
considered and the estimators implemented in the package.

# Installation

A github repository contains the source files for the package. To
install it using the
[devtools](https://cran.r-project.org/web/packages/devtools/index.html)
package, run the following command.

``` r
devtools::install_github("package-repository-url")
```

# Statistical Model and estimators

The likelihood used to fit the models is based on the distribution of a
vector with \(m\) with variables, c change points and \(c+1\) blocks.
The \(i\)-th block has a probability parameter \(p_i\), and the entries
on the blocks are independent Bernoulli with the parameter of the block.

# Example

## Simulation Setup

Here we provide an example of usage of the package on simulated data.
The function `rBBM` generates data according to the model proposed.
Consider \(50\) data samples from a vector of size \(20\) separated in
\(4\) blocks: \(1:5\), \(6:10\), \(11:15\) and \(16:20\). We assign the
probability vector \((0.3, 0.7, 0.9, 0.1)\) to represent the probability
of each of the blocks. Here is the code for sampling from the proposed
setup.

``` r
library(bincpd)
set.seed(42) 
sim_info <- rBBM(n = 50, m = 20, prob_vec = c(0.3, 0.7, 0.9, 0.1), changepoints = c(5, 10, 15))
```

Notice that we have \(4\) blocks, but only \(3\) change points, which
are precisely the end point of each block, with the exception of the
last block. The functions returns us the `data_matrix`, the
`changepoints` and the probability vector `prob_vec`. The last two are
useful when we do not specify these arguments, and hence the change
points and probabilities are simulated uniformly before simulating the
data vector.

## Estimators

A single function called `fit_model` wraps up all methods. The function
takes \(3\) arguments: the `data_matrix`, the data that change point
analysis will be done; the `method`, a string specifying the method to
be fitted; and `arguments`, a list of arguments that are specific to
that method. Check the documentation to see which arguments can be
passed to which methods.

Five methods are provided: `hierseg`, `dynseg`, `cvdynseg`, `cvseg` and
`fusedlasso`. All methods were implemented using
[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) to speed
up calculations. Below we show the code to fit each of the models using
their default arguments.

``` r
sim_data  <- sim_info$data_matrix
# Using default arguments for all methods.
hs_model  <- fit_model(data_matrix = sim_data, method = "hierseg") 
ds_model  <- fit_model(data_matrix = sim_data, method = "dynseg")
cvds_model  <- fit_model(data_matrix = sim_data, method = "cvdynseg")
cvseg_model <- fit_model(data_matrix = sim_data, method = "cvseg")
fused_model <- fit_model(data_matrix = sim_data, method = "fusedlasso")
```

Every model returns an S3 object of the class *bincpd*. Different models
might return more or less information about their fit
procedure/arguments, but all of them always contain at least four
attributes:

  - *changepoints*: a list containing the set of estimated change
    points;
  - *probabilities*: a list containing the estimated probability
    parameters for each block;
  - *loss*: the final loss evaluated on the entire data set for the
    returned model;
  - *n\_cp*: number of change points estimated.

Below, we show the change point set and probabilities estimated by each
of the methods.

    #> [1] "Summary of hierseg"
    #> Changepoints: 5 10 15
    #> Probabilities: 0.328 0.736 0.9 0.116
    #> Loss: 489.1267
    #> [1] "Summary of dynseg"
    #> Changepoints: 5 10 15
    #> Probabilities: 0.328 0.736 0.9 0.116
    #> Loss: 489.1267
    #> [1] "Summary of cvdynseg"
    #> Changepoints: 5 10 15 19
    #> Probabilities: 0.328 0.736 0.9 0.14 0.02
    #> Loss: 488.0979
    #> [1] "Summary of cvseg"
    #> Changepoints: 4 5 10 15 19
    #> Probabilities: 0.355 0.22 0.736 0.9 0.14 0.02
    #> Loss: 467.908
    #> [1] "Summary of fusedlasso"
    #> Changepoints: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
    #> Probabilities: 0.4209291 0.2709076 0.3637832 0.34 0.2492841 0.7290924 0.7120819 0.72 0.7724995 0.749388 0.9067716 0.8768924 0.929784 0.9 0.88 0.14 0.1160525 0.14 0.1590682 0.02165782
    #> Loss: 472.1623

Albeit all models report the loss attained by the estimated model, they
might not be comparable between different methods. For instance, the
*hierseg* method and the *fusedlasso* have different loss functions, so
loss comparison might not make sense.

Notice that the `fusedlasso` has over segmented the array. The reason is
that the penalization constant default \(\lambda = 1\) was too low. Here
we exemplify how to pass arguments to a method, changing the \(\lambda\)
arguments and also the `precision` argument. The precision is the
threshold used by the `fusedlasso` to detect change points. When two
consecutive probabilities differ by more than the precision, a change
point is added to that location.

``` r
arguments <- list(lambda = 10, precision = 1e-2)
fused_model <- fit_model(data_matrix = sim_data, method = "fusedlasso", arguments = arguments)
```

The summary of the new fitted model is shown below.

    #> [1] "Summary of fusedlasso"
    #> Changepoints:
    #> 5 8 10 15 19
    #> Probabilities:
    #> 0.364236 0.7199909 0.76 0.8629977 0.140014 0.07283751
    #> Loss:
    #> 472.1623
