#ifndef QuadraticBlockCoord_h
#define QuadraticBlockCoord_h

#include <RcppEigen.h>
#include "stateVec.h"
#include <vector>
#include "groupIndices.h"
#include <cstdlib>
#include <Eigen/SparseQR>
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi intVector;
typedef Eigen::MatrixXd matrix;
using Rcpp::Rcout;
using namespace Eigen;


//block coordinate descent main class
template<typename XMatrix> 
class QuadraticBlockCoord
{
protected:

  int maxit;  //maximum number of iterations
  
  double tol;           // tolerance for convergence
  
  double eigenValueTolerance; //eigenvalue tolerance to determine whether
                              //group is full rank
  bool scale;
  

  
  bool converged; //convergence information

  
  //gradient information
  vec gradient;
  
  //"standardization" when there are unpenalized groups
  matrix projectedResids; //residual vector with unpenalized groups projected out
  bool anyUnpenalized; //true if there are any unpenalized groups
  matrix R1R1t; //cross product of inverse matrices of unpenalized groups
  XMatrix XUnpenalized;
  intVector rankVector;

  XMatrix X; //design matrix
  
  bool standardize; //do standardization
  std::vector<matrix> RInverseListProjected; //for standardization purposes
  matrix xtx0; //crossproduct of X with intercept group
  Eigen::MatrixXi xtx0Flag; //has the current dot product been computed?
  
  std::vector<matrix> InterceptUpdateList; //for intercept projections
  bool interceptComputed;
                    



  
public:
  bool useWeights;
  vec sampleWeights;
  double sampleWeightsSum;
  
  double boundConstant; //Diagonal(boundConstant) bounds the within-observation hessian
  
  
  matrix response; //response vector;
  matrix linPred; //current linear predictor
  groupIndices groups; //contains grouping and state information about beta vector
  stateVec active; //state vector with active coefficients
  
  //coefficient information
  vec beta; //current coefficients
  vec beta_prev; //previous coefficients
  vec rescaledBeta; //coefficients on the original scale
  
  matrix compressedBeta; //for the overlapping group lasso:
                      //return beta as a pxk matrix
                      
  //the linear predictor/residual/response formulations are most useful for exponential family models
  //predictor, residual, and response information
  matrix resids; //residual vector
                      
  matrix betaMatrix; //for returning coefficient matrix on the original scale
  
  
  //problem dimensions
  int k;      //kronecker product dimension 
  //only !=1 for multinomial version
  int p;      // dimension of X
  int n;       // number of rows
  int p_beta; //number of predictor rows in beta for the
              //for the overlapping group lasso problem
  
  //constructor
  QuadraticBlockCoord(int p_
                        ,int n_
                        ,int k_
                        ,int nGroups_
                        ,intVector groupStart_
                        ,intVector groupEnd_
                        ,intVector groupColumns_
                        ,double tol_
                        ,vec penaltyFactor_
                        ,int maxit_
                        ,matrix response_
                        ,XMatrix X_
                        ,double boundConstant_
                        ,double eigenValueTolerance_
                        ,bool scale_)
    :groups(nGroups_
              ,groupStart_
              ,groupEnd_
              ,groupColumns_
              ,p_
              ,n_
              ,penaltyFactor_)
    ,active(nGroups_)
  {
    boundConstant=boundConstant_;
    eigenValueTolerance=eigenValueTolerance_;
    scale=scale_;
    useWeights=false;
    p_beta=groups.groupColumns.size();
    rankVector.setZero(groups.nGroup);
    p=p_;
    n=n_;
    k=k_;
    tol=tol_;
    maxit=maxit_;
    converged=false;


    
    //initialize beta, gradient, residuals, projected residuals
    rescaledBeta.setZero(groups.groupColumns.size()*k);
    beta.setZero(groups.groupColumns.size()*k);
    beta_prev.setZero(groups.groupColumns.size()*k);
    gradient.setZero(groups.groupColumns.size()*k);
    resids.setZero(n,k);
    linPred.setZero(n,k);
    projectedResids.setZero(n,k);

    X=X_;
    response=response_;
    interceptComputed=false;
    unpenalizedSetUp();
    setRInverseList();
    updateResiduals();
  }
  
  //constructor
  QuadraticBlockCoord(int p_
                        ,int n_
                        ,int k_
                        ,int nGroups_
                        ,intVector groupStart_
                        ,intVector groupEnd_
                        ,intVector groupColumns_
                        ,double tol_
                        ,vec penaltyFactor_
                        ,int maxit_
                        ,matrix response_
                        ,XMatrix X_
                        ,double boundConstant_
                        ,double eigenValueTolerance_
                        ,bool scale_
                        ,vec sampleWeights_)
    :groups(nGroups_
              ,groupStart_
              ,groupEnd_
              ,groupColumns_
              ,p_
              ,n_
              ,penaltyFactor_)
    ,active(nGroups_)
  {
    eigenValueTolerance=eigenValueTolerance_;
    scale=scale_;
    p_beta=groups.groupColumns.size();
    rankVector.setZero(groups.nGroup);
    p=p_;
    n=n_;
    k=k_;
    tol=tol_;
    maxit=maxit_;
    converged=false;

    
    //initialize beta, gradient, residuals, projected residuals
    rescaledBeta.setZero(groups.groupColumns.size()*k);
    beta.setZero(groups.groupColumns.size()*k);
    beta_prev.setZero(groups.groupColumns.size()*k);
    gradient.setZero(groups.groupColumns.size()*k);
    resids.setZero(n,k);
    linPred.setZero(n,k);
    projectedResids.setZero(n,k);
    
    X=X_;
    response=response_;
    interceptComputed=false;
    unpenalizedSetUp();
    useWeights=true;
    sampleWeights=sampleWeights_;
    
    //make weights sum to number of nonzero weights
    sampleWeightsSum=0;
    for (int nInd=0;nInd<n;nInd++){
      if (sampleWeights(nInd)>0){
        sampleWeightsSum=sampleWeightsSum+1;
      }
    }
    
    sampleWeights=(sampleWeights/sampleWeights.array().sum())*sampleWeightsSum;
    setRInverseList();
    updateResiduals();
    boundConstant=boundConstant_*sampleWeights.maxCoeff();
  }
  
  //compute t(X)%*%(y-mu)
  void xtr(int p_group,matrix& residsVec){
    int groupLength=groups.getLength(p_group);
    Map<MatrixXd> grad_p=getSlice(p_group,gradient);
    for (int p1=0;p1<groupLength;p1++){
      int XIndex=groups.getColumn(p_group,p1);
      grad_p.block(p1,0,1,k)=X.col(XIndex).adjoint()*(residsVec);
    }
    grad_p=RInverseListProjected[p_group].adjoint()*grad_p;
    
    if (groups.penaltyFactor(p_group)>0){
      Map<MatrixXd> beta_p=getSlice(p_group,beta);
      grad_p=grad_p+beta_p;
    }
    
    if (groups.anyUnpenalized & (groups.penaltyFactor(p_group)>0)){
      Map<MatrixXd> beta_Intercept=getSlice(groups.unpenalizedIndex,beta);
      grad_p=grad_p-InterceptUpdateList[p_group]*beta_Intercept;
    }
    if (useWeights){
      grad_p=grad_p/sampleWeightsSum;
    }
    else{
      grad_p=grad_p/n;
    }
  }
  
  //set new responses
  void setResponse(matrix newResponses){
    interceptComputed=false;
    response=newResponses;
  }
  
  //do block update for group p_group
  void updateGroup(int p_group, double lambda,bool isActive){
    if (groups.penaltyFactor(p_group)==0){
      if (!interceptComputed){
        interceptComputed=true;
        xtr(p_group,resids);
      }
      else{
        return;
      }
    }
    else{
      xtr(p_group,resids);
    }
    
    double gradNorm=computeGroupNorm(p_group,gradient);
    if (gradNorm<(lambda/boundConstant*pow(rankVector(p_group)*k,0.5)*groups.penaltyFactor(p_group))){
      if (active.nonZeroGroups(p_group)!=0){
        converged=false;
        Map<MatrixXd> beta_p=getSlice(p_group,beta);
        Map<MatrixXd> beta_p_prev=getSlice(p_group,beta_prev);
        beta_p_prev=beta_p;
        beta_p.setZero();
        active.nonZeroGroups(p_group)=0;
        active.groupNorms(p_group)=0;

        
        //updating residuals only
        if (groups.penaltyFactor(p_group)>0){
          updateResids(p_group,X);
        }
      }
    }
    else{
      double stepSize;
      //weirdly, in poisson and logistic regression, the norm of the gradient 
      //can be zero or possibly -0, though usually only with an unpenalized intercept
      if (gradNorm!=0){
        stepSize=1-lambda/boundConstant*pow(rankVector(p_group)*k,0.5)*groups.penaltyFactor(p_group)/gradNorm;
      }
      else{
        stepSize=1;
      }
      Map<MatrixXd> beta_p=getSlice(p_group,beta);
      Map<MatrixXd> beta_p_prev=getSlice(p_group,beta_prev);
      Map<MatrixXd> grad_p=getSlice(p_group,gradient);
      beta_p_prev=beta_p;
      if (useWeights){
        beta_p=stepSize*grad_p*sampleWeightsSum;
      }
      else{
        beta_p=stepSize*grad_p*n;
      }
      checkConvergedGroup(beta_p,beta_p_prev,p_group,converged,tol);
      
      active.groupNorms(p_group)=beta_p.norm();
      active.nonZeroGroups(p_group)=1;
      if (!isActive){
        active.addActive(p_group);
      }
      
      
      //updating residuals only
      if (groups.penaltyFactor(p_group)>0){
        updateResids(p_group,X);
      }
    }
  }
  
  void checkConvergedGroup(Map<MatrixXd>& beta_p,Map<MatrixXd>& beta_p_prev,int p_group,bool& converged_,double tol_){
    int groupLength=groups.getLength(p_group);
    for (int kInd=0;kInd<k;kInd++){
      for (int p1=0;p1<groupLength;p1++){
        if (( beta_p(p1,kInd) != 0 && beta_p_prev(p1,kInd) == 0) || (beta_p(p1,kInd) == 0 && beta_p_prev(p1,kInd) != 0) ) {
          converged_=false;
        }
        if (beta_p(p1,kInd) != 0 && beta_p_prev(p1,kInd) != 0 &&
            std::abs( (beta_p(p1,kInd) - beta_p_prev(p1,kInd)) / beta_p_prev(p1,kInd) ) > tol_) {
          converged_=false;
        }
        // if (beta_p(p1,kInd) != 0 && beta_p_prev(p1,kInd) != 0 &&
        //     std::abs( (beta_p(p1,kInd) - beta_p_prev(p1,kInd)) ) > tol_) {
        //   converged_=false;
        // }
      }
    }
  }
  
  //solve a quadratic coordinate descent problem for a given lambda
  //cycle over active groups
  //if active groups converged, cycle over inactive groups
  //if no change, we're done, otherwise repeat
  int solve(double lambda_){
    int iterations;
    double lambda=lambda_;

    for (iterations=0;iterations<maxit;iterations++){
      converged=true;

      if (active.getNumActive()>0){
       updateGroups(0,active.getNumActive(),lambda,true);
      }
      if (converged){
        if (active.getNumActive()<groups.nGroup){
          updateGroups(active.getNumActive(),groups.nGroup,lambda,false);
        }
      }
      
      if (converged){
        active.setActiveToNonZero();
        break;
      }
    }
    return(iterations);
  }
  
  //update the given groups
  void updateGroups(int first,int last,double lambda_,bool isActive){
    for (int i=first;i<last; i++){
      updateGroup(active.getGroup(i),lambda_,isActive);
    }
  }
  
  //get lambda value that gets rid of penalized coefficients
  //assumes the unpenalized coefficients have been optimized and the rest are zero
  double computeMaxLambda(){
    Rcout<<"max lambda computation"<<std::endl;
    Rcout<<"resids="<<resids<<std::endl;
    
    double maxLambda=0;
    for (int p_group=0;p_group<groups.nGroup;p_group++){
      Rcout<<"Penalty factor="<<groups.penaltyFactor(p_group)<<std::endl;
      if (groups.penaltyFactor(p_group)>0){
        xtr(p_group,resids);
        double tempLambda=computeGroupNorm(p_group,gradient)/groups.penaltyFactor(p_group);

        tempLambda=tempLambda/pow(rankVector(p_group)*k,0.5)*boundConstant;
        if (tempLambda>maxLambda){
          maxLambda=tempLambda;
        }
      }
    }
    return(maxLambda);
  }
  
  //compute group norm
  double computeGroupNorm(int p_group, vec& groupedVector){
    return(getSlice(p_group,groupedVector).norm());
  }
  
  
  //in constructor, do qr decompositions for each group
  void setRInverseList(){
    
    if (groups.anyUnpenalized){
      groupQR(groups.unpenalizedIndex,true);
    }
    for (int p_group=0;p_group<groups.nGroup; p_group++){
      if (groups.penaltyFactor(p_group)==0){
      }
      else{
        if (groups.anyUnpenalized){
          int lengthUnpenalized=groups.getLength(groups.unpenalizedIndex);
          int length_pgroup=groups.getLength(p_group);

          Eigen::MatrixXd cprod;
          cprod.setZero(length_pgroup,lengthUnpenalized);

          Eigen::MatrixXd cprodGroup;
          cprodGroup.setZero(length_pgroup,length_pgroup);

          // //Xk^tX0
          for (int p1=0;p1<lengthUnpenalized;p1++){
            for (int p2=0;p2<length_pgroup;p2++){
              int col1=groups.getColumn(groups.unpenalizedIndex,p1);
              int col2=groups.getColumn(p_group,p2);
              if (xtx0Flag(col2,col1)==0){
                cprod(p2,p1)=dot(col1,col2);
                xtx0(col2,col1)=cprod(p2,p1);
                xtx0Flag(col2,col1)=1;
              }
              else{
                cprod(p2,p1)=xtx0(col2,col1);
              }
            }
          }


          //xk^txk
          for (int p1=0;p1<length_pgroup;p1++){
            for (int p2=p1;p2<length_pgroup;p2++){
              int col1=groups.getColumn(p_group,p1);
              int col2=groups.getColumn(p_group,p2);
              cprodGroup(p2,p1)=dot(col1,col2);
              cprodGroup(p1,p2)=cprodGroup(p2,p1);
            }
          }

          Eigen::MatrixXd xtx_projected;
          xtx_projected=cprodGroup-cprod*R1R1t*cprod.adjoint();

          eigenTransform(xtx_projected,p_group,true);
          InterceptUpdateList[p_group]=RInverseListProjected[p_group].adjoint()*cprod*RInverseListProjected[groups.unpenalizedIndex];
        }
        else{
          groupQR(p_group,true);
        }
      }
    }
  }
  
  //return list of R^-1 matrices
  std::vector<Eigen::MatrixXd> getRInverseList(){
    return(RInverseListProjected);
  }
  
  //do qr decomposition of group p_group
  void groupQR(int p_group,bool pushback){
    int groupLength=groups.getLength(p_group);
    Eigen::MatrixXd cprodGroup(groupLength,groupLength);
    for (int p1=0;p1<groupLength;p1++){
      for (int p2=p1;p2<groupLength;p2++){
        cprodGroup(p2,p1)=dot(groups.getColumn(p_group,p1),groups.getColumn(p_group,p2));
        cprodGroup(p1,p2)=cprodGroup(p2,p1);
      }
    }
    eigenTransform(cprodGroup,p_group,pushback);
  }
  
  void eigenTransform(matrix& crossProd,int p_group,bool pushback){
    int length=(crossProd).cols();
    SelfAdjointEigenSolver<MatrixXd> test(crossProd);
    Eigen::VectorXd eigenValues=test.eigenvalues();
    Eigen::MatrixXd eigenVectors=test.eigenvectors();
    
    //compute projected rank
    ColPivHouseholderQR<MatrixXd> rankMat(crossProd);
    rankVector(p_group)=rankMat.rank();
    
    Eigen::MatrixXd eigenTransformMat;
    eigenTransformMat.setZero(length,length);
    for (int j=0; j<length;j++){
      if (eigenValues(j)>eigenValueTolerance){
        if (scale){
          eigenTransformMat.col(j)=pow(eigenValues(j),-0.5)*eigenVectors.col(j);
        }
        else{
          eigenTransformMat.col(j)=eigenVectors.col(j);
        }
      }
      else{
        Rcout<<"Collinearity in group "<<(p_group+1)<<" after projecting out unpenalized groups."<<std::endl;
      }
    }
    if ((p_group==groups.unpenalizedIndex)){
      R1R1t=eigenTransformMat*eigenTransformMat.adjoint();
    }
    RInverseListProjected[p_group]=eigenTransformMat;
  }
  

  void updateResids(int p_group, Eigen::MatrixXd& X_){
    int groupLength=groups.getLength(p_group);
    Map<MatrixXd> beta_p=getSlice(p_group,beta);
    Map<MatrixXd> beta_p_prev=getSlice(p_group,beta_prev);
    
    //kind of bad, using the already allocated rescaled beta space to 
    //temporarily hold onto the difference in coefficients on the original scale
    //rescaled beta will actually be totally computed in getRescaledBeta
    
    Map<MatrixXd> rescaledBeta_p=getSlice(p_group,rescaledBeta);
    rescaledBeta_p=RInverseListProjected[p_group]*(beta_p-beta_p_prev);
    
    for (int p1=0;p1<groupLength;p1++){
      int XIndex=groups.getColumn(p_group,p1);
      if (!useWeights){
        resids=resids-X.col(XIndex)*rescaledBeta_p.block(p1,0,1,k);
      }
      else{
        resids=resids-sampleWeights.asDiagonal()*X.col(XIndex)*rescaledBeta_p.block(p1,0,1,k);
      }
    }
    
    //update the intercept
    if (groups.anyUnpenalized & (groups.penaltyFactor(p_group)>0)){
      Map<MatrixXd> beta_Intercept=getSlice(groups.unpenalizedIndex,beta);
      beta_Intercept=beta_Intercept-InterceptUpdateList[p_group].adjoint()*(beta_p-beta_p_prev);
    }
  }
  void updateResids(int p_group, Eigen::SparseMatrix<double>& X_){
    int groupLength=groups.getLength(p_group);
    Map<MatrixXd> beta_p=getSlice(p_group,beta);
    Map<MatrixXd> beta_p_prev=getSlice(p_group,beta_prev);
    
    //kind of bad, using the already allocated rescaled beta space to 
    //temporarily hold onto the difference in coefficients on the original scale
    //rescaled beta will actually be totally computed in getRescaledBeta
    
    Map<MatrixXd> rescaledBeta_p=getSlice(p_group,rescaledBeta);
    rescaledBeta_p=RInverseListProjected[p_group]*(beta_p-beta_p_prev);
    
    for (int p1=0;p1<groupLength;p1++){
      int XIndex=groups.getColumn(p_group,p1);
      for (int kInd=0;kInd<k;kInd++){
        sparseDenseSum(resids,kInd,X,XIndex,-rescaledBeta_p(p1,kInd));
      }
    }
    
    //update the intercept
    if (groups.anyUnpenalized & (groups.penaltyFactor(p_group)>0)){
      Map<MatrixXd> beta_Intercept=getSlice(groups.unpenalizedIndex,beta);
      beta_Intercept=beta_Intercept-InterceptUpdateList[p_group].adjoint()*(beta_p-beta_p_prev);
    }
  }
  
  void sparseDenseSum(Eigen::MatrixXd& linPred, int kInd, Eigen::SparseMatrix<double>& X_, int XIndex,double coefficient){
    Eigen::SparseMatrix<double>::InnerIterator it(X_,XIndex);
    for (it.row(); it; ++it){
      if (!useWeights){
        linPred(it.row(),kInd)=linPred(it.row(),kInd)+it.value()*coefficient;
      }
      else{
        linPred(it.row(),kInd)=linPred(it.row(),kInd)+it.value()*coefficient*sampleWeights(it.row());
      }
    }
  }
  
  
  //update residuals
  void updateResiduals(){
    resids=response-linPred;
    if (groups.anyUnpenalized){
      Map<MatrixXd> beta_Intercept=getSlice(groups.unpenalizedIndex,beta);
      resids=resids+XUnpenalized*RInverseListProjected[groups.unpenalizedIndex]*beta_Intercept;
    }
    if (useWeights){
      resids=sampleWeights.asDiagonal()*resids;
    }
  }
  
  //if sampleweights(nInd)==0, the linear predictor for nInd computed 
  //by this function will not be correct
  void updateLinearPredictorFromResiduals(){
    if (!useWeights){
      linPred=response-resids;
    }
    else{
      linPred.setZero(n,k);
      for (int nInd=0;nInd<n;nInd++){
        if (sampleWeights(nInd)>0){
          linPred.block(nInd,0,1,k)=response.block(nInd,0,1,k)-resids.block(nInd,0,1,k)/sampleWeights(nInd);
        }
      }
    }

    //add in intercept contributions
    if (groups.anyUnpenalized){
      Map<MatrixXd> beta_Intercept=getSlice(groups.unpenalizedIndex,beta);
      linPred=linPred+XUnpenalized*RInverseListProjected[groups.unpenalizedIndex]*(beta_Intercept);
    }
  }
  
  //return linear predictor
  Eigen::MatrixXd getLinearPredictor(){
    return(linPred);
  }
  
  //return rescaled coefficients
  matrix getRescaledBeta(bool returnCompressed){
    betaMatrix.setZero(p_beta,k);
    if (returnCompressed){
      compressedBeta.setZero(X.cols(),k);
    }
    
    for (int i=0; i<active.getNumActive();i++){
      int p_group=active.getGroup(i);
      if (active.nonZeroGroups(p_group)!=0){
        Map<MatrixXd> beta_p=getSlice(p_group,beta);
        int groupLength=groups.getLength(p_group);
        betaMatrix.block(groups.getStart(p_group),0,groupLength,k)=RInverseListProjected[p_group]*beta_p;
      }
      
      if (returnCompressed){
        for (int p1=0;p1<groups.getLength(p_group);p1++){
          int XIndex=groups.getColumn(p_group,p1);
          compressedBeta.block(XIndex,0,1,k)=compressedBeta.block(XIndex,0,1,k)+betaMatrix.block(groups.getStart(p_group)+p1,0,1,k);
        }
      }
    }  
    
    if (returnCompressed){
      return(compressedBeta); 
    }
    else{
      return(betaMatrix);
    }
  }
  
  //get residuals
  matrix getResiduals(){
    return(resids);
  }
  
  //return gradient (mainly for debugging)
  vec getGradient(){
    return(gradient);
  }
  
  //update unpenalized groups only
  void updateUnpenalizedGroups(){
    if (groups.anyUnpenalized){
      bool isActive=(active.groupVector(groups.unpenalizedIndex)!=0);
      updateGroup(groups.unpenalizedIndex,0,isActive);
    }
  }
  
  //solve unpenalized
  int solveUnpenalized(){
    converged=false;
    int iterations=0;
    
    while (!converged && iterations<maxit){
      //the group update will check convergence
      //if the convergence check fails, converged will be set to false
      converged=true;
      updateUnpenalizedGroups();
      iterations++;
    }
    return(iterations);
  };
  
  //return beta
  vec getBeta(){
    return(beta);
  }
  
  //deal with unpenalized groups
  void unpenalizedSetUp(){
    groups.anyUnpenalized=false;
    groups.unpenalizedIndex=-1;
    RInverseListProjected.resize(groups.nGroup);
    for (int i=0;i<groups.nGroup;i++){
      if (groups.penaltyFactor(i)==0){
        InterceptUpdateList.resize(groups.nGroup);
        groups.anyUnpenalized=true;
        groups.unpenalizedIndex=i;
        int groupLength=groups.getLength(i);
        XUnpenalized.resize(n,groupLength);
        xtx0.setZero(X.cols(),groupLength);
        xtx0Flag.setZero(X.cols(),groupLength);
        for (int p1=0;p1<groupLength;p1++){
          XUnpenalized.col(p1)=X.col(groups.getColumn(i,p1));
        }
      }
    }
  }
  
  //sparse and dense dot product methods
  double dot(int col1,int col2){
    if (!useWeights){
      return(dotUnweighted(col1,col2,X));
    }
    else{
      return(dotWeighted(col1,col2,X,sampleWeights));
    }
  }
  
  //weighted and unweighted dot products for dense matrices
  double dotUnweighted(int col1, int col2, Eigen::MatrixXd & X_){
    return(((X_).col(col1).array()*(X_).col(col2).array()).sum());
  }
  double dotWeighted(int col1, int col2, Eigen::MatrixXd & X_, Eigen::VectorXd& weights_){
    return(((X_).col(col1).array()*(X_).col(col2).array()*weights_.array()).sum());
  }
  
  //weighted and unweighted dot products for sparse matrices
  double dotUnweighted(int col1, int col2, Eigen::SparseMatrix<double> & X_){

      double dotProd=0;
      Eigen::SparseMatrix<double>::InnerIterator it1(X_,col1);
      Eigen::SparseMatrix<double>::InnerIterator it2(X_,col2);
      for (it1.row(); it1;++it1){
        while (it2.row()<it1.row() && it2){
          ++it2;
        }

        if (it2.row()==it1.row()){
          dotProd=dotProd+it1.value()*it2.value();
        }
      }
    return(dotProd);
  }
  double dotWeighted(int col1, int col2, Eigen::SparseMatrix<double> & X_, Eigen::VectorXd& weights_){
    
    double dotProd=0;
    Eigen::SparseMatrix<double>::InnerIterator it1(X_,col1);
    Eigen::SparseMatrix<double>::InnerIterator it2(X_,col2);
    for (it1.row(); it1;++it1){
      while (it2.row()<it1.row() && it2){
        ++it2;
      }
      
      if (it2.row()==it1.row()){
        dotProd=dotProd+it1.value()*it2.value()*weights_(it2.row());
      }
    }
    return(dotProd);
  }

  
  //computes lambda path: assumes unpenalized groups have been optimized
  //and all other groups are zero
  vec computeLambdaPath(int pathLength,double lambdaMinRatio){
    Eigen::VectorXd lambdaVec;
    lambdaVec.setZero(pathLength);
    double lambdaMax=computeMaxLambda();
    double logDiff=std::log(lambdaMax)-std::log(lambdaMinRatio*lambdaMax);
    double ratio=std::exp(-logDiff/(pathLength-1));
    lambdaVec(0)=lambdaMax;
    for (int i=1; i<pathLength;i++){
      lambdaVec(i)=lambdaVec(i-1)*ratio;
    }
    return(lambdaVec);
  }
  
  
  inline Map<MatrixXd> getSlice(int p_group,vec& groupedVector){
    int start=groups.getStart(p_group)*k;
    int groupLength=groups.getLength(p_group);
    return(Map<MatrixXd>(&(groupedVector(start)),groupLength,k));
  }
  
  void setBoundConstant(double newBoundConstant_){
    boundConstant=newBoundConstant_;
  }
  
  //compute linear predictor for the observations in the given fold
  //garbage
  void computeHoldoutLinearPredictor(intVector& fold, matrix& holdoutLinPred){
    int n_fold=fold.size();
    holdoutLinPred.setZero(n_fold,k);
    
    for (int i=0;i<active.getNumActive();i++){
      
      int p_group=active.getGroup(i);
      Map<MatrixXd> beta_p=getSlice(p_group,beta);
      Map<MatrixXd> rescaledBeta_p=getSlice(p_group,rescaledBeta);
      rescaledBeta_p=RInverseListProjected[p_group]*(beta_p);
      
      for (int p1=0;p1<groups.getLength(p_group);p1++){
        
        int XIndex=groups.getColumn(p_group,p1);
        
        for (int nInd=0;nInd<n_fold;nInd++){
          
          int nCurrent=fold(nInd);
          holdoutLinPred.block(nInd,0,1,k)=holdoutLinPred.block(nInd,0,1,k)+X.block(nCurrent,XIndex,1,1)*rescaledBeta_p.block(p1,0,1,k);
        
        }
      } 
    }
  }
};

#endif
