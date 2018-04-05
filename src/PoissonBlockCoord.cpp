#include <Rcpp.h>
#include <RcppEigen.h>
#include "PoissonModel.h"
using Rcpp::Rcout;
// 
// [[Rcpp::export]]
Rcpp::List poissonDense(int nGroups_
                                                ,Eigen::VectorXi groupStart_
                                                ,Eigen::VectorXi groupEnd_
                                                ,Eigen::VectorXi groupColumns_
                                                ,double tol_
                                                ,Eigen::VectorXd penaltyFactor_
                                                ,int maxit_
                                                ,Eigen::MatrixXd response_
                                                ,Eigen::MatrixXd X_
                                                ,int nLambda_
                                                ,double lambdaMinRatio_
                                                ,double eigenValueTolerance_
                                                ,bool scale_
                                                ,bool useLambda_
                                                ,Eigen::VectorXd lambda_
                                                ,bool useWeights_
                                                ,Eigen::VectorXd sampleWeights_
                                                ,Eigen::MatrixXd offset_
                                                ,bool useDevTol_
                                                ,double devTol_){
  int p_=X_.cols();
  int n_=X_.rows();
  
  PoissonModel<Eigen::MatrixXd>* poisModel;
  
  if (!useWeights_){
    poisModel=new PoissonModel<Eigen::MatrixXd>(p_
                                                      ,n_
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
                                                      ,offset_);
  }
  else{
    poisModel=new PoissonModel<Eigen::MatrixXd>(p_
                                                      ,n_
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
                                                      ,offset_
                                                      ,sampleWeights_);
  }
  return(solveModel(poisModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));
}

// [[Rcpp::export]]
Rcpp::List poissonSparse(int nGroups_
                                            ,Eigen::VectorXi groupStart_
                                            ,Eigen::VectorXi groupEnd_
                                            ,Eigen::VectorXi groupColumns_
                                            ,double tol_
                                            ,Eigen::VectorXd penaltyFactor_
                                            ,int maxit_
                                            ,Eigen::MatrixXd response_
                                            ,Eigen::SparseMatrix<double> X_
                                            ,int nLambda_
                                            ,double lambdaMinRatio_
                                            ,double eigenValueTolerance_
                                            ,bool scale_
                                            ,bool useLambda_
                                            ,Eigen::VectorXd lambda_
                                            ,bool useWeights_
                                            ,Eigen::VectorXd sampleWeights_
                                            ,Eigen::MatrixXd offset_
                                            ,bool useDevTol_
                                            ,double devTol_){
  int p_=X_.cols();
  int n_=X_.rows();
  
  PoissonModel<Eigen::SparseMatrix<double> >* poisModel;
  
  if (!useWeights_){
    poisModel=new PoissonModel<Eigen::SparseMatrix<double> >(p_
                                                  ,n_
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
                                                  ,offset_);
  }
  else{
    poisModel=new PoissonModel<Eigen::SparseMatrix<double> >(p_
                                                  ,n_
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
                                                  ,offset_
                                                  ,sampleWeights_);
  }
  return(solveModel(poisModel,nLambda_,lambdaMinRatio_,useLambda_,lambda_,useWeights_,useDevTol_,devTol_));
}