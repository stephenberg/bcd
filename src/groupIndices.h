#include <RcppEigen.h>
typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi intVector;
typedef Eigen::MatrixXi intMat;
typedef Eigen::MatrixXd matrix;
using Rcpp::Rcout;


//information re: grouping
//to avoid painful indexing 
class groupIndices
{
public:
  int nGroup;  //number of groups

  intVector groupStart; //starting points for each group in groupColumns vector
  intVector groupEnd; //ending points for each group in groupColumns vector
  
  intVector groupColumns; //vector with the columns for each group
                          //the starts and breaks of each group are in
                          //groupStart and groupEnd
                          
  intVector groupLengths; //vector with lengths of each group
                          
  int p; //number of columns in X matrix
  int n; //number of rows in X matrix
  
  vec penaltyFactor;
  
  
  bool anyUnpenalized; //are any groups unpenalized?
  int unpenalizedIndex; //index of unpenalized group, if any
                        //-1 if no unpenalized group
  
  groupIndices(int nGroup_, intVector groupStart_, intVector groupEnd_, intVector groupColumns_, int p_,int n_,vec penaltyFactor_):
    nGroup(nGroup_), groupStart(groupStart_), groupEnd(groupEnd_), groupColumns(groupColumns_),p(p_), n(n_), penaltyFactor(penaltyFactor_)
  {
    groupLengths.setZero(nGroup);
    for (int i=0; i<nGroup;i++){
      groupLengths(i)=groupEnd[i]-groupStart[i]+1;
    }
  }

  
  //get first column (in groupColumns (X matrix)) of p_group-th group
  int getStart(int p_group){
    return(groupStart[p_group]);
  }
  
  //get last column ""
  int getEnd(int p_group){
    return(groupEnd[p_group]);
  }
  
  //get number of columns of X matrix in the p_group-th group
  int getLength(int p_group){
    return(groupLengths(p_group));
  }
  
  //get column in the X matrix of the groupPosition-th column of group 
  //p_group
  int getColumn(int p_group, int groupPosition){
    return(groupColumns(getStart(p_group)+groupPosition));
  }
};
