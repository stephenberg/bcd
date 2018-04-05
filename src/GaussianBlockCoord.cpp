#include <Rcpp.h>
#include <RcppEigen.h>
#include "Model.h"
using Rcpp::Rcout;

// [[Rcpp::export]]
Rcpp::List multiResponseGaussianDense(int nGroups_
                                                          ,Eigen::VectorXi groupStart_
                                                          ,Eigen::VectorXi groupEnd_
                                                          ,Eigen::VectorXi groupColumns_
                                                          ,double tol_
                                                          ,Eigen::VectorXd penaltyFactor_
                                                          ,int maxit_
                                                          ,Eigen::MatrixXd response_
                                                          ,Eigen::MatrixXd X_
                                                          ,double boundConstant_
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
  
  Model<Eigen::MatrixXd>* quadModel;
  
  if (!useWeights_){
    quadModel=new Model<Eigen::MatrixXd>(p_
                                     ,n_
                                     ,k_
                                     ,nGroups_
                                     ,groupStart_
                                     ,groupEnd_
                                     ,groupColumns_
                                     ,tol_
                                     ,penaltyFactor_
                                     ,maxit_
                                     ,response_
                                     ,X_
                                     ,1
                                     ,tol_
                                     ,1
                                     ,eigenValueTolerance_
                                     ,scale_);
  }
  else{
    quadModel=new Model<Eigen::MatrixXd>(p_
                                           ,n_
                                           ,k_
                                           ,nGroups_
                                           ,groupStart_
                                           ,groupEnd_
                                           ,groupColumns_
                                           ,tol_
                                           ,penaltyFactor_
                                           ,maxit_
                                           ,response_
                                           ,X_
                                           ,1
                                           ,tol_
                                           ,1
                                           ,eigenValueTolerance_
                                           ,scale_
                                           ,sampleWeights_);
  }  
  return(solveModel(quadModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));
}


// [[Rcpp::export]]
Rcpp::List multiResponseGaussianSparse(int nGroups_
                                                          ,Eigen::VectorXi groupStart_
                                                          ,Eigen::VectorXi groupEnd_
                                                          ,Eigen::VectorXi groupColumns_
                                                          ,double tol_
                                                          ,Eigen::VectorXd penaltyFactor_
                                                          ,int maxit_
                                                          ,Eigen::MatrixXd response_
                                                          ,Eigen::SparseMatrix<double> X_
                                                          ,double boundConstant_
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
  
  Model<Eigen::SparseMatrix<double> >* quadModel;
  
  if (!useWeights_){
    quadModel=new Model<Eigen::SparseMatrix<double> >(p_
                                           ,n_
                                           ,k_
                                           ,nGroups_
                                           ,groupStart_
                                           ,groupEnd_
                                           ,groupColumns_
                                           ,tol_
                                           ,penaltyFactor_
                                           ,maxit_
                                           ,response_
                                           ,X_
                                           ,1
                                           ,tol_
                                           ,1
                                           ,eigenValueTolerance_
                                           ,scale_);
  }
  else{
    quadModel=new Model<Eigen::SparseMatrix<double> >(p_
                                           ,n_
                                           ,k_
                                           ,nGroups_
                                           ,groupStart_
                                           ,groupEnd_
                                           ,groupColumns_
                                           ,tol_
                                           ,penaltyFactor_
                                           ,maxit_
                                           ,response_
                                           ,X_
                                           ,1
                                           ,tol_
                                           ,1
                                           ,eigenValueTolerance_
                                           ,scale_
                                           ,sampleWeights_);
  }
  return(solveModel(quadModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));
}
