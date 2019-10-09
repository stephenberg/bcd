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
class PoissonModel: public Model<XMatrix>{
protected:
  matrix offset;
public:
  
  //constructor
  PoissonModel(int p_
                     ,int n_
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
                     ,matrix offset_):Model<XMatrix>(p_
                                                                ,n_
                                                                ,1
                                                                ,nGroups_
                                                                ,groupStart_
                                                                ,groupEnd_
                                                                ,groupColumns_
                                                                ,innerTol_
                                                                ,penaltyFactor_
                                                                ,maxitInner_
                                                                ,response_
                                                                ,X_
                                                                ,1
                                                                ,outerTol_
                                                                ,maxitOuter_
                                                                ,eigenValueTolerance_
                                                                ,scale_){
                       offset=offset_;
                     }
  
  //constructor
  PoissonModel(int p_
                 ,int n_
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
                 ,matrix offset_
                 ,vec sampleWeights_):Model<XMatrix>(p_
                                                        ,n_
                                                        ,1
                                                        ,nGroups_
                                                        ,groupStart_
                                                        ,groupEnd_
                                                        ,groupColumns_
                                                        ,innerTol_
                                                        ,penaltyFactor_
                                                        ,maxitInner_
                                                        ,response_
                                                        ,X_
                                                        ,1
                                                        ,outerTol_
                                                        ,maxitOuter_
                                                        ,eigenValueTolerance_
                                                        ,scale_
                                                        ,sampleWeights_){
                   offset=offset_;
                 }
  
  //compute mean for each observation
  virtual void computeMeans(){
    this->mu=(offset+this->model.linPred).array().exp();
    double maxMean=0;
    for (int nInd=0;nInd<this->model.n;nInd++){
      if (this->model.useWeights){
        double temp=(this->model.sampleWeights(nInd))*this->mu(nInd,0);
        if (temp>maxMean){
          maxMean=temp;
        }
      }
      else{
        if (this->mu(nInd,0)>maxMean){
          maxMean=this->mu(nInd,0);
        }
      }
    }
    this->model.boundConstant=maxMean;
  }
  
  virtual void computeHoldoutMeans(matrix& mu_, matrix& linPred_,intVector& fold){
    mu_.setZero(linPred_.rows(),linPred_.cols());
    for (int nInd=0;nInd<fold.size();nInd++){
      mu_(nInd,0)=std::exp(offset(fold(nInd),0)+linPred_(nInd,0));
    }
  }

  virtual double computeHoldoutDeviance(matrix& mu_,intVector& fold_, bool useWeights_, vec& sampleWeights_){
    double deviance=0;
    for (int nInd=0;nInd<fold_.size();nInd++){
      int nCurrent=fold_(nInd);
      if (useWeights_){
        double temp=0;
        if ((this->response(nCurrent,0)>0) && (sampleWeights_(nCurrent)>0)){
          temp=temp+this->response(nCurrent,0)*(std::log(mu_(nInd,0))-std::log(this->response(nCurrent,0)));
        }
        temp=temp+this->response(nCurrent,0)-mu_(nInd,0);
        deviance=deviance+sampleWeights_(nCurrent)*temp;
      }
      else{
        double temp=0;
        if ((this->response(nCurrent,0)>0) && (sampleWeights_(nCurrent)>0)){
          temp=temp+this->response(nCurrent,0)*(std::log(mu_(nInd,0))-std::log(this->response(nCurrent,0)));
        }
        temp=temp+this->response(nCurrent,0)-mu_(nInd,0);
        deviance=deviance+sampleWeights_(nCurrent)*temp;      }
    }
    return(-2*deviance);
  }
  
  virtual double computeDeviance(){
    double deviance=0;
    for (int nInd=0;nInd<this->model.n;nInd++){
      double temp=0;

      if (this->model.useWeights){
        if ((this->response(nInd,0)>0) && (this->model.sampleWeights(nInd)>0)){
          temp=temp+this->response(nInd,0)*(std::log(this->mu(nInd,0))-std::log(this->response(nInd,0)));
        }
        temp=temp+this->response(nInd,0)-this->mu(nInd,0);
        deviance=deviance+this->model.sampleWeights(nInd)*temp;
      }
      else{
        if ((this->response(nInd,0)>0)){
          temp=temp+this->response(nInd,0)*(std::log(this->mu(nInd,0))-std::log(this->response(nInd,0)));
        }
        temp=temp+this->response(nInd,0)-this->mu(nInd,0);
        deviance=deviance+temp;
      }
    }
    deviance=-2*deviance;
    return(deviance);
  }
};