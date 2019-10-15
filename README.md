# bcd
Block coordinate descent for group lasso

This package uses block coordinate descent to fit common GLMs with a group lasso penalty. The implementation is in C++ with wrapper functions using Rcpp. The package is capable of fitting models with overlapping groups (with implicit rather than explicit duplication of the design matrix columns for memory saving). The package can also be used to fit sparse group lasso type models using the overlapping group lasso framework. Run help("fit_bcd") in R for some examples of using the package.
