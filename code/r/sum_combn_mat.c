#include <R.h>
#include <Rinternals.h>
#include <float.h>

static inline int iterate_indices(int*, int, int);

SEXP sum_combn_mat(SEXP XSXP, SEXP nrowsXSXP, SEXP numCombSXP, SEXP numOutSXP){
  
  double *X = REAL(XSXP);
  int nrowsX = *INTEGER(nrowsXSXP);
  
  int numComb = *INTEGER(numCombSXP);
  int numOutMin = INTEGER(numOutSXP)[0];
  int numOutMax = INTEGER(numOutSXP)[1];
  
  SEXP indicesMinSXP = PROTECT(Rf_allocVector(INTSXP, numOutMin * numComb));
  SEXP indicesMaxSXP = PROTECT(Rf_allocVector(INTSXP, numOutMax * numComb));
  
  int *indicesWorking = (int*) R_alloc(numComb, sizeof(int));
  int *indicesMin = INTEGER(indicesMinSXP);
  int *indicesMax = INTEGER(indicesMaxSXP);
  
  for(int ind = 0; ind < numComb; ++ind){
    indicesWorking[ind] = ind;
  }
  
  SEXP minSumSXP = PROTECT(Rf_allocVector(REALSXP, numOutMin));
  SEXP maxSumSXP = PROTECT(Rf_allocVector(REALSXP, numOutMax));
  
  double *minSum = REAL(minSumSXP);
  double *maxSum = REAL(maxSumSXP);
  
  double minSumWorking = DBL_MAX;
  double maxSumWorking = 0;
  
  for(int ind = 0; ind < numOutMin; ++ind){
    minSum[ind] = DBL_MAX;
  }
  for(int ind = 0; ind < numOutMax; ++ind){
    maxSum[ind] = 0;
  }
  
  do{
    
    double sum = 0;
    for(int ind1 = 0; ind1 < numComb-1; ++ind1){
      for(int ind2 = ind1+1; ind2 < numComb; ++ind2){
        int ind3 = indicesWorking[ind2];
        sum += X[ind3*(ind3-1)/2+indicesWorking[ind1]];
      }
    }
    
    if(sum < minSumWorking){
      
      int ind1;
      for(ind1 = numOutMin-1; minSum[ind1-1] > sum && ind1 > 0; --ind1){
        for(int ind2 = 0; ind2 < numComb; ++ind2){
          indicesMin[ind1*numComb+ind2] = indicesMin[(ind1-1)*numComb+ind2];
        }
        minSum[ind1] = minSum[ind1-1];
      }
      
      for(int ind2 = 0; ind2 < numComb; ++ind2){
        indicesMin[ind1*numComb+ind2] = indicesWorking[ind2]+1;
      }
      minSum[ind1] = sum;
      minSumWorking = minSum[numOutMin-1];
    }
    
    if(sum > maxSumWorking){
      
      int ind1;
      for(ind1 = numOutMax-1; maxSum[ind1-1] < sum && ind1 > 0; --ind1){
        for(int ind2 = 0; ind2 < numComb; ++ind2){
          indicesMax[ind1*numComb+ind2] = indicesMax[(ind1-1)*numComb+ind2];
        }
        maxSum[ind1] = maxSum[ind1-1];
      }
      
      for(int ind2 = 0; ind2 < numComb; ++ind2){
        indicesMax[ind1*numComb+ind2] = indicesWorking[ind2]+1;
      }
      maxSum[ind1] = sum;
      maxSumWorking = maxSum[numOutMax-1];
    }
    
  } while(iterate_indices(indicesWorking, numComb, nrowsX));
  
  SEXP listOutSXP = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP listOutMinSXP = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP listOutMaxSXP = PROTECT(Rf_allocVector(VECSXP, 2));
  
  SET_VECTOR_ELT(listOutSXP, 0, listOutMinSXP);
  SET_VECTOR_ELT(listOutSXP, 1, listOutMaxSXP);
  
  SET_VECTOR_ELT(listOutMinSXP, 0, indicesMinSXP);
  SET_VECTOR_ELT(listOutMinSXP, 1, minSumSXP);
  
  SET_VECTOR_ELT(listOutMaxSXP, 0, indicesMaxSXP);
  SET_VECTOR_ELT(listOutMaxSXP, 1, maxSumSXP);
  
  UNPROTECT(7);
  
  return listOutSXP;
}

static inline int iterate_indices(int *indices, int length, int indexMax){
  
  int ind1;
  for(ind1 = length-1; indices[ind1] > indexMax-length+ind1-1 && ind1 > 0; --ind1);
  
  int indTemp = indices[ind1]+1;
  for(int ind2 = ind1; ind2 < length; ++ind2){
    indices[ind2] = indTemp;
    ++indTemp;
  }
  
  return ind1 > 0 || indTemp < indexMax;
}