#include <RcppEigen.h>

typedef Eigen::VectorXd vec;
typedef Eigen::VectorXi intVector;
typedef Eigen::MatrixXd matrix;
using Rcpp::Rcout;



//Class stores active set in manner described below

//example object state:
//groupVector=1 3 2 0 8 4 5 6 7
//nActive=4
//Thus, groups 1, 3, 2, and 0 are active, while groups 8, 4, 5, 6, 7 are not active

//This storage format simulates a list and avoids resizing the vector containing the active set.
class stateVec
{
protected:
  
public:
  
  //constructor
  stateVec(int nGroups_):nGroups(nGroups_),nActive(0),groupVector(nGroups),groupIndices(nGroups),nonZeroGroups(nGroups),groupNorms(nGroups)
  {
    for (int i=0; i<nGroups; i++){
      groupVector(i)=i;
      groupIndices(i)=i;
      nonZeroGroups(i)=0;
      groupNorms(i)=0;
    }
  }
  
  int nGroups;
  int nActive;
  intVector groupVector;
  intVector groupIndices;
  intVector nonZeroGroups;
  vec groupNorms;
  
  //mostly for debugging: return the current state of the active/inactive sets
  intVector getGroups(){
    return(groupVector);
  }
  
  //the starting index of the active set is always 0
  int getActiveStart(){
    return(0);
  }
  
  //the ending index of the active set is always nActive-1
  int getNumActive(){
    return(nActive);
  }
  
  //if iterating through the active (inactive) set, return the group 
  //corresponding to the current index
  int getGroup(int index_){
    return(groupVector(index_));
  }
  

  
  //add group to active set
  void addActive(int newGroup){
    // Rcout<<"new group: "<<newGroup<<std::endl;
    // Rcout<<"start"<<std::endl;
    // Rcout<<"nActive"<<getNumActive()<<std::endl;
    // Rcout<<"active"<<std::endl;
    // Rcout<<getGroups()<<std::endl;
    //groupVector(group)=groupVector(nActive);
    //groupVector(nActive)=group;
    //nActive=nActive+1;
    
    
    //new code
    //if the next group is already in the right spot, no need to switch entries
    if (groupVector(nActive)==newGroup){
      nActive=nActive+1;
    }
    //otherwise switch entries
    else{
      int currentIndex=groupIndices(newGroup);
      int currentOccupant=groupVector(nActive);
      //(int newOccupant=newGroup;)
      groupVector(nActive)=newGroup;
      groupVector(currentIndex)=currentOccupant;
      groupIndices(currentOccupant)=currentIndex;
      groupIndices(newGroup)=nActive;
      nActive=nActive+1;
    }
    // Rcout<<"end"<<std::endl;
    // Rcout<<"nActive"<<getNumActive()<<std::endl;
    // Rcout<<"active"<<std::endl;
    // Rcout<<getGroups()<<std::endl;
  }
  
  void reset(){
    for (int i=0; i<nGroups; i++){
      groupVector(i)=i;
      groupIndices(i)=i;
      groupNorms(i)=0;
      nonZeroGroups(i)=0;
    }
    nActive=0;
  }
  
  void setActiveToNonZero(){
    nActive=0;
    for (int i=0; i<nGroups; i++){
      groupVector(i)=i;
      groupIndices(i)=i;
    }
    for (int i=0; i<nGroups;i++){
      if (groupNorms(i)>0){
        addActive(i);
      }
    }
  }

};