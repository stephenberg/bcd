# bcd
Block coordinate descent for group lasso

Package fits common GLMs with group lasso penalty by coordinate descent. Implemented in C++ via Rcpp. Capable of fitting models with overlapping groups (with implicit rather than explicit duplication of the design matrix columns for memory saving). Can also fit sparse group lasso type models using the overlapping group lasso framework. Run help("fit_bcd") in R for some examples of using the package.
