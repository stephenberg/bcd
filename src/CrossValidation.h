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
class CrossValidation{
protected:
  Model<XMatrix> cvModel;
  intVector foldList;
  int nFolds;
  matrix devianceMatrix;
  matrix holdoutLinPred;
  matrix holdoutMu;
  vec foldWeights;
  bool useWeights;
  vec sampleWeights;
  
  
public:
  

  matrix crossValidate(vec lambda_){
    for (int foldInd=0;foldInd<nFolds;foldInd++){
      foldWeights=sampleWeights;
      for (int nInd=0;nInd<foldList[foldInd].size();nInd++){
        foldWeights(foldList(nInd))=0;
      }
      
      int n_fold=foldList[foldInd].size();
      holdoutLinPred.resize(n_fold,cvModel.model.k);
      holdoutLinPred.setZero(n_fold,cvModel.model.k);
      holdoutMu.resize(n_fold,cvModel.model.k);
      holdoutMu.setZero(n_fold,cvModel.model.k);
      
      for (int lamInd=0;lamInd<lambda_.size();lamInd++){
        cvModel.solve(lambda_(lamInd));
        cvModel.model.computeHoldoutLinearPredictor(foldList[foldInd],holdoutLinPred);
        cvModel.computeHoldoutMeans(holdoutMu,holdoutLinPred,foldList[foldInd]);
        devianceMatrix(lamInd,foldInd)=cvModel.computeHoldoutDeviance(holdoutMu,foldList[foldInd],useWeights,sampleWeights);
      }
    }
  }
  
};