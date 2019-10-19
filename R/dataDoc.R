#' Design matrix for tests and examples
#'
#' A 1000 by 50 matrix with 1's for the intercept in the 
#' first column and numeric values in the remaining columns. 
#' Used in tests and examples
#'
#' @format A matrix with 53940 rows and 10 variables:
"X"

#' True coefficients for tests
#'
#' Length 50 coefficient vector for all tests except
#' multinomial and multiresponse Gaussian. The 11th 
#' through 30th coefficients are all 0. The remainder
#' are nonzero.
#'
#' @format Length 50 numeric vector
"beta"

#' True coefficients for multiresponse tests
#' 
#' 50 by 3 coefficient matrix in the tests for
#' multiresponse Gaussian and multinomial tests. 
#' The coefficients in rows 11-30 are all 0.
#'
#' @format 50 by 3 numeric matrix
"beta_multiresponse"

#' Grouping used for the example tests
#' 
#' Grouping used for the examples. Length 4 list of 
#' integer vectors.
#'
#' @format Length 4 list of integer vectors.
"grouping"

#' Number of categories for multiresponse examples
#' 
#' Integer: k=3.
#'
#' @format integer
"k"

#' Number of samples in examples
#' 
#' Integer: n=1000.
#'
#' @format integer
"n"

#' Number of predictors in examples
#' 
#' Integer: p=50.
#'
#' @format integer
"p"

#' penaltyFactor
#' 
#' Length 4 vector containing penalization level for groups
#' in the examples: c(0,1,1,1). Intercept group unpenalized.
#'
#' @format Length 4 numeric vector
"penaltyFactor"

#' Sample weights
#' 
#' Length 1000 numeric vector of sample weights for the examples.
#'
#' @format Length 1000 numeric vector
"sampleWeights"

#' Binary response vector
#' 
#' Length 1000 vector of 0 and 1 values for logistic 
#' regression example.
#'
#' @format Length 1000 integer vector
"y_binary"

#' Multinomial response example
#' 
#' Length 1000 vector with values 1, 2, and 3 for the multinomial
#' regression example.
#'
#' @format Length 1000 numeric vector
"y_multinomial"

#' Multiresponse Gaussian response matrix
#' 
#' 1000 by 3 numeric matrix containing responses for the multiresponse 
#' Gaussian example
#'
#' @format 1000 by 3 numeric matrix
"y_multiresponse"

#' Gaussian response vector
#' 
#' Length 1000 numeric vector containing responses for the 
#' linear regression example.
#'
#' @format Length 1000 numeric vector
"y_gaussian"

#' Poisson response vector
#' 
#' Length 1000 integer vector containing responses for the Poisson regression example.
#'
#' @format Length 1000 numeric vector
"y_count"

#' Linear regression reference coefficients
#' 
#' Linear regression reference coefficients for test without sample weights.
#' 
#' @format Length 100 list of 50 by 1 numeric matrices
"reference_Linear"

#' Linear regression reference coefficients
#' 
#' Linear regression reference coefficients for test with sample weights.
#' 
#' @format Length 100 list of 50 by 1 numeric matrices
"reference_Linear_Weighted"

#' Logistic regression reference coefficients
#' 
#' Logistic regression reference coefficients for test without sample weights.
#' 
#' @format Length 100 list of 50 by 2 numeric matrices
"reference_Logistic"

#' Logistic regression reference coefficients
#' 
#' Logistic regression reference coefficients for test with sample weights.
#' 
#' @format Length 100 list of 50 by 2 numeric matrices
"reference_Logistic_Weighted"

#' Multinomial regression reference coefficients
#' 
#' Multinomial regression reference coefficients for test without sample weights.
#' 
#' @format Length 100 list of 50 by 3 numeric matrices
"reference_Multinomial"

#' Multinomial regression reference coefficients
#' 
#' Multinomial regression reference coefficients for test with sample weights.
#' 
#' @format Length 100 list of 50 by 3 numeric matrices
"reference_Multinomial_Weighted"

#' Multiresponse Gaussian reference coefficients
#' 
#' Multiresponse Gaussian reference coefficients for test without sample weights.
#' 
#' @format Length 100 list of 50 by 3 numeric matrices
"reference_Multiresponse"

#' Multiresponse Gaussian reference coefficients
#' 
#' Multiresponse Gaussian reference coefficients for test with sample weights.
#' 
#' @format Length 100 list of 50 by 3 numeric matrices
"reference_Multiresponse_Weighted"

#' Poisson regression reference coefficients
#' 
#' Poisson regression reference coefficients for test without sample weights.
#' 
#' @format Length 100 list of 50 by 1 numeric matrices
"reference_Poisson"

#' Poisson regression reference coefficients
#' 
#' Poisson regression reference coefficients for test with sample weights.
#' 
#' @format Length 100 list of 50 by 1 numeric matrices
"reference_Poisson_Weighted"




