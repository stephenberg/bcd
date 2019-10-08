#include <RcppEigen.h>
#include <cstdlib>
#include <Eigen/SparseQR>
#include "QuadraticBlockCoord.h"

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi intVector;
typedef Eigen::MatrixXd matrix;

using Rcpp::Rcout;
using namespace Eigen;
using namespace Rcpp;


template<typename XMatrix>
class Model
{
protected:
  double nullDeviance; //deviance under the null model
  double saturatedDeviance; //deviance under the saturated model
  double currentDeviance; //current deviance
  double previousDeviance; //previous deviance
  bool useDeviance; //use deviance to determine convergence
  double devianceTolerance;
  
  double outerTol; //tolerance for changes in coefficients between different 
                   //quadratic approximation steps
                   
  matrix response;  //the observed response
  matrix workingResponse; //the irls working response
  matrix mu;
  
  int maxitOuter; //number of outer quadratic approximation steps that are permitted

  vec lambda; //lambda vector
  double lambdaMin; //minimum lambda value
                    
  vec beta_prev_outer;//vector before solving the quadratic coordinate descent problem
  vec beta_outer;//vector after solving the quadratic coordinate descent problem
  std::vector<matrix> betaMat;
  bool gaussian;
  
public:
  
  QuadraticBlockCoord<XMatrix> model; //object that will handle most of the solving details
  
  //constructor
  Model(int p_
          ,int n_
          ,int k_
          ,int nGroups_
          ,intVector groupStart_
          ,intVector groupEnd_
          ,intVector groupColumns_
          ,double innerTol_
          ,vec penaltyFactor_
          ,int maxitInner_
          ,matrix response_
          ,XMatrix X_
          ,double boundConstant_
          ,double outerTol_
          ,int maxitOuter_
          ,double eigenValueTolerance_
          ,bool scale_):model(p_
                                         ,n_
                                         ,k_
                                         ,nGroups_
                                         ,groupStart_
                                         ,groupEnd_
                                         ,groupColumns_
                                         ,innerTol_
                                         ,penaltyFactor_
                                         ,maxitInner_
                                         ,response_
                                         ,X_
                                         ,boundConstant_
                                         ,eigenValueTolerance_
                                         ,scale_){
    response=response_;
    workingResponse.setZero(n_,k_);
    mu.setZero(n_,k_);
    outerTol=outerTol_;
    maxitOuter=maxitOuter_;
    beta_prev_outer=model.beta;
    beta_outer=beta_prev_outer;
  }

  //constructor
  Model(int p_
          ,int n_
          ,int k_
          ,int nGroups_
          ,intVector groupStart_
          ,intVector groupEnd_
          ,intVector groupColumns_
          ,double innerTol_
          ,vec penaltyFactor_
          ,int maxitInner_
          ,matrix response_
          ,XMatrix X_
          ,double boundConstant_
          ,double outerTol_
          ,int maxitOuter_
          ,double eigenValueTolerance_
          ,bool scale_
          ,vec sampleWeights_):model(p_
                                ,n_
                                ,k_
                                ,nGroups_
                                ,groupStart_
                                ,groupEnd_
                                ,groupColumns_
                                ,innerTol_
                                ,penaltyFactor_
                                ,maxitInner_
                                ,response_
                                ,X_
                                ,boundConstant_
                                ,eigenValueTolerance_
                                ,scale_
                                ,sampleWeights_){
    response=response_;
    workingResponse.setZero(n_,k_);
    mu.setZero(n_,k_);
    outerTol=outerTol_;
    maxitOuter=maxitOuter_;
    beta_prev_outer=model.beta;
    beta_outer=beta_prev_outer;
  }
  
  
  //return solution path
  Eigen::MatrixXd getBetaMat(){
    return(betaMat);
  }
  
  //given linear predictor values in model and response, set working response
  void computeWorkingResponse(){
    workingResponse.setZero(model.n,model.k);
    computeMeans();
    
    //1/boundConstant
    workingResponse=workingResponse+(1/model.boundConstant)*(response-mu)+model.linPred;
  }
  
  //compute mean for each observation
  virtual void computeMeans(){
    mu=model.linPred;
  }
  
  virtual void computeHoldoutMeans(matrix& mu_, matrix& linPred_,intVector& fold){
    mu_=linPred_;
  }
  
  virtual double computeHoldoutDeviance(matrix& mu_,intVector& fold_, bool useWeights_, vec& sampleWeights_){
    double deviance=0;
    for (int nInd=0;nInd<fold_.size();nInd++){
      int nCurrent=fold_(nInd);
      if (useWeights_){
        deviance=deviance+sampleWeights_(nCurrent)*(response.block(nCurrent,0,1,this->model.k)-mu_.block(nInd,0,1,this->model.k)).array().pow(2).sum();
      }
      else{
        deviance=deviance+(response.block(nCurrent,0,1,this->model.k)-mu_.block(nInd,0,1,this->model.k)).array().pow(2).sum();
      }
    }
    return(deviance);
  }
  
  //compute current deviance
  virtual double computeDeviance(){
    double deviance=0;
    deviance=(response-model.linPred).array().pow(2).sum();
    return(deviance);  
  }
  
  //get means
  matrix getMeans(){
    return(mu);
  }
  
  //get beta vector
  matrix getBeta(bool compressedBeta){
    return(model.getRescaledBeta(compressedBeta));
  }
  
  //get working response
  matrix getWorkingResponse(){
    return(workingResponse);
  }
  
  //get response
  matrix getresponse(){
    return(response);
  }
  
  //set response
  // void setresponse(matrix newresponse){
  //   response=newresponse;
  // }
  // 
  //check convergence
  bool checkConverged(){
    bool converged=true;
    
    //we may only need to check the groups 
    //that are active at the end of the current quadratic fitting
    //for each lambda, this set only gets bigger
    
    //for now, check all of the groups
    for (int i=0;i<model.active.nActive;i++){
      int p_group=model.active.getGroup(i);
      converged=checkConvergedGroup(p_group);
      if (!converged){
        return(converged);
      }
    }
    return(converged);
  }
  
  bool checkConvergedGroup(int p_group){
    bool converged=true;
    Map<MatrixXd> beta_p_outer=model.getSlice(p_group,beta_outer);
    Map<MatrixXd> beta_p_outer_prev=model.getSlice(p_group,beta_prev_outer);
    model.checkConvergedGroup(beta_p_outer,beta_p_outer_prev,p_group,converged,outerTol);
    return(converged);
  }
  
  //solve unpenalized problem
  //basically irls for the unpenalized covariates
  int solveUnpenalized(){
    bool converged;
    int iterations;
    for (iterations=0;iterations<maxitOuter;iterations++){
      model.updateLinearPredictorFromResiduals();
      computeWorkingResponse();
      model.setResponse(workingResponse);
      model.updateResiduals();
      beta_prev_outer=model.beta;
      model.solveUnpenalized();
      beta_outer=model.beta;
      converged=checkConverged();
      if (converged){
        model.updateLinearPredictorFromResiduals();
        computeWorkingResponse();
        model.setResponse(workingResponse);
        model.updateResiduals();
        break;
      }
    }
    return(iterations);
  }
  
  //compute maximum lambda value
  double computeMaxLambda(){
    solveUnpenalized();
    return(model.computeMaxLambda());
  }
  
  //set lambda path
  void setLambda(vec lambda_){
    lambda=lambda_;
  }
  
  //compute lambda path
  vec computeLambdaPath(int pathLength, double lambdaMinRatio){
    return(model.computeLambdaPath(pathLength,lambdaMinRatio));
  }
  
  //solve the penalized multinomial regression problem at a certain lambda value
  int solve(double lambda_){
    bool converged;
    int iterations;
    for (iterations=0;iterations<maxitOuter;iterations++){
      beta_prev_outer=model.beta;
      
      model.updateLinearPredictorFromResiduals();
      computeWorkingResponse();
      model.setResponse(workingResponse);
      model.updateResiduals();

      model.solve(lambda_);
      
      beta_outer=model.beta;
      converged=checkConverged();
      if (converged){
        break;
        
      }
    }
    return(iterations);
  }
  
  //using the lambda path computed by computeLambdaPath(),
  //solve the coordinate descent problem for each lambda value
  void solvePath(){
    for (int i=0;i<lambda.size();i++){
      solve(lambda(i));
    }
  }
  
  vec getLambdaSequence(){
    return(lambda);
  }
};

template<typename XMatrix>
Rcpp::List solveModel(Model<XMatrix>* outerModel, int nLambda_, double lambdaMinRatio_, bool useLambda_, Eigen::VectorXd lambda_,bool useWeights_,bool useDevTol, double devTol_){
  
  Eigen::VectorXd lambda;
  double nullDev=0;
  outerModel->solveUnpenalized();
  if (useDevTol){
    nullDev=outerModel->computeDeviance();
  }
  if (useLambda_){
    lambda=lambda_;
    nLambda_=lambda.size();
  }
  else{
    lambda=outerModel->model.computeLambdaPath(nLambda_,lambdaMinRatio_);
  }
  
  Eigen::VectorXd deviance;
  deviance.setZero(nLambda_);

  std::vector<Eigen::MatrixXd> betaMat;
  
  for (int i=0;i<nLambda_;i++){
    outerModel->solve(lambda(i));
    deviance(i)=outerModel->computeDeviance();
    if (useDevTol){
      if ((nullDev-deviance(i))>devTol_*nullDev){
        break;
      }
    }
    betaMat.push_back(outerModel->getBeta(true));
  }

  if (useWeights_){
    return(Rcpp::List::create(Rcpp::Named("beta")=betaMat,
                              Rcpp::Named("lambda")=lambda,
                              Rcpp::Named("deviance")=deviance,
                              Rcpp::Named("sampleWeights")=outerModel->model.sampleWeights
    ));
  }
  else{
    Eigen::VectorXd sampleWeights;
    sampleWeights.setOnes(outerModel->model.n);
    return(Rcpp::List::create(Rcpp::Named("beta")=betaMat,
                              Rcpp::Named("lambda")=lambda,
                              Rcpp::Named("deviance")=deviance,
                              Rcpp::Named("sampleWeights")=sampleWeights
    ));
  }
}
