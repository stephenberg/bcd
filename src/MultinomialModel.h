#include <RcppEigen.h>
#include <cstdlib>
#include <Eigen/SparseQR>
#include "Model.h"

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi intVector;
typedef Eigen::MatrixXd matrix;

using Rcpp::Rcout;
using namespace Eigen;
using namespace Rcpp;

template<typename XMatrix>
class MultinomialModel: public Model<XMatrix>{

public:
  
  //constructor
  MultinomialModel(int p_
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
          ,double outerTol_
          ,int maxitOuter_
          ,double eigenValueTolerance_
          ,bool scale_):Model<XMatrix>(p_
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
                                    ,0.5
                                    ,outerTol_
                                    ,maxitOuter_
                                    ,eigenValueTolerance_
                                    ,scale_){
  }
  
  //constructor
  MultinomialModel(int p_
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
                     ,double outerTol_
                     ,int maxitOuter_
                     ,double eigenValueTolerance_
                     ,bool scale_
                     ,vec sampleWeights_):Model<XMatrix>(p_
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
                                                    ,0.5
                                                    ,outerTol_
                                                    ,maxitOuter_
                                                    ,eigenValueTolerance_
                                                    ,scale_
                                                    ,sampleWeights_
                                                    ){
                     }
  
  //compute mean for each observation
  virtual void computeMeans(){

    for (int nInd=0;nInd<this->model.n;nInd++){
      double expMax=this->model.linPred.block(nInd,0,1,this->model.k).maxCoeff();
      double expSum=(this->model.linPred.block(nInd,0,1,this->model.k).array()-expMax).exp().sum();
      this->mu.block(nInd,0,1,this->model.k)=(this->model.linPred.block(nInd,0,1,this->model.k).array()-expMax).exp()/expSum;
    }
  }
  
  //compute mean for each observation
  virtual void computeMeans(matrix& mu_, matrix& linPred_, intVector& fold){
    mu_.setZero(linPred_.rows(),linPred_.cols());
    int n_fold=linPred_.rows();
    for (int nInd=0;nInd<n_fold;nInd++){
      double expMax=linPred_.block(nInd,0,1,this->model.k).maxCoeff();
      double expSum=(linPred_.block(nInd,0,1,this->model.k).array()-expMax).exp().sum();
      mu_.block(nInd,0,1,this->model.k)=(linPred_.block(nInd,0,1,this->model.k).array()-expMax).exp()/expSum;
    }
  }
  
  //compute deviance
  virtual double computeDeviance(){
    double deviance=0;
    for (int nInd=0;nInd<this->model.n;nInd++){
      for (int kInd=0;kInd<this->model.k;kInd++){
        if (this->model.useWeights){
          deviance=deviance+(this->response)(nInd,kInd)*std::log((this->mu)(nInd,kInd))*this->model.sampleWeights(nInd);
        }
        else{
          deviance=deviance+(this->response)(nInd,kInd)*std::log((this->mu)(nInd,kInd));
        }
      }
    }
    deviance=-2*deviance;
    
    return(deviance);
  }
  
  virtual double computeHoldoutDeviance(matrix& mu_,intVector& fold_, bool useWeights_, vec& sampleWeights_){
    double deviance=0;
    for (int nInd=0;nInd<fold_.size();nInd++){
      int nCurrent=fold_(nInd);
      for (int kInd=0;kInd<this->model.k;kInd++){
        if (useWeights_){
          deviance=deviance+(this->response)(nCurrent,kInd)*std::log(mu_(nInd,kInd))*sampleWeights_(nCurrent);
        }
        else{
          deviance=deviance+(this->response)(nCurrent,kInd)*std::log(mu_(nInd,kInd));
        }
      }
    }
    deviance=-2*deviance;
    return(deviance);
  }
  
};