#include <Rcpp.h>
#include <RcppEigen.h>
#include "MultinomialModel.h"
using Rcpp::Rcout;
// 
// [[Rcpp::export]]
Rcpp::List multinomialDense(int nGroups_
                                                          ,Eigen::VectorXi groupStart_
                                                          ,Eigen::VectorXi groupEnd_
                                                          ,Eigen::VectorXi groupColumns_
                                                          ,double tol_
                                                          ,Eigen::VectorXd penaltyFactor_
                                                          ,int maxit_
                                                          ,Eigen::MatrixXd response_
                                                          ,Eigen::MatrixXd X_
                                                          ,int k_
                                                          ,int nLambda_
                                                          ,double lambdaMinRatio_
                                                          ,double eigenValueTolerance_
                                                          ,bool scale_
                                                          ,bool useLambda_
                                                          ,Eigen::VectorXd lambda_
                                                          ,bool useWeights_
                                                          ,Eigen::VectorXd sampleWeights_
                                                          ,bool useDevTol_
                                                          ,double devTol_){
  int p_=X_.cols();
  int n_=X_.rows();
  
  MultinomialModel<Eigen::MatrixXd>* multModel;
  
  if (!useWeights_){
    multModel=new MultinomialModel<Eigen::MatrixXd>(p_
                                           ,n_
                                           ,k_
                                           ,nGroups_
                                           ,groupStart_
                                           ,groupEnd_
                                           ,groupColumns_
                                           ,tol_
                                           ,penaltyFactor_
                                           ,1
                                           ,response_
                                           ,X_
                                           ,tol_
                                           ,maxit_
                                           ,eigenValueTolerance_
                                           ,scale_);
  }
  else{
    multModel=new MultinomialModel<Eigen::MatrixXd>(p_
                                                      ,n_
                                                      ,k_
                                                      ,nGroups_
                                                      ,groupStart_
                                                      ,groupEnd_
                                                      ,groupColumns_
                                                      ,tol_
                                                      ,penaltyFactor_
                                                      ,1
                                                      ,response_
                                                      ,X_
                                                      ,tol_
                                                      ,maxit_
                                                      ,eigenValueTolerance_
                                                      ,scale_
                                                      ,sampleWeights_);
  }
  return(solveModel(multModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));
}

// [[Rcpp::export]]
Rcpp::List multinomialSparse(int nGroups_
                                                ,Eigen::VectorXi groupStart_
                                                ,Eigen::VectorXi groupEnd_
                                                ,Eigen::VectorXi groupColumns_
                                                ,double tol_
                                                ,Eigen::VectorXd penaltyFactor_
                                                ,int maxit_
                                                ,Eigen::MatrixXd response_
                                                ,Eigen::SparseMatrix<double> X_
                                                ,int k_
                                                ,int nLambda_
                                                ,double lambdaMinRatio_
                                                ,double eigenValueTolerance_
                                                ,bool scale_
                                                ,bool useLambda_
                                                ,Eigen::VectorXd lambda_
                                                ,bool useWeights_
                                                ,Eigen::VectorXd sampleWeights_
                                                ,bool useDevTol_
                                                ,double devTol_){
  int p_=X_.cols();
  int n_=X_.rows();
  
  MultinomialModel<Eigen::SparseMatrix<double> >* multModel;
  
  if (!useWeights_){
    multModel=new MultinomialModel<Eigen::SparseMatrix<double> >(p_
                                                      ,n_
                                                      ,k_
                                                      ,nGroups_
                                                      ,groupStart_
                                                      ,groupEnd_
                                                      ,groupColumns_
                                                      ,tol_
                                                      ,penaltyFactor_
                                                      ,1
                                                      ,response_
                                                      ,X_
                                                      ,tol_
                                                      ,maxit_
                                                      ,eigenValueTolerance_
                                                      ,scale_);
  }
  else{
    multModel=new MultinomialModel<Eigen::SparseMatrix<double> >(p_
                                                      ,n_
                                                      ,k_
                                                      ,nGroups_
                                                      ,groupStart_
                                                      ,groupEnd_
                                                      ,groupColumns_
                                                      ,tol_
                                                      ,penaltyFactor_
                                                      ,1
                                                      ,response_
                                                      ,X_
                                                      ,tol_
                                                      ,maxit_
                                                      ,eigenValueTolerance_
                                                      ,scale_
                                                      ,sampleWeights_);
  }
  
  return(solveModel(multModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));

}
