#include <R.h>
#include <Rinternals.h>
#include <float.h>

SEXP sum_combn_mat(SEXP XSXP, SEXP kSXP){
  
  double *X = REAL(XSXP);
  
  int nrowsX = Rf_nrows(XSXP);
  
  int k = *INTEGER(kSXP);
  
  int *indices = (int*) R_alloc(k, sizeof(int));
  
  SEXP indicesMaxSXP = PROTECT(Rf_allocVector(INTSXP, k));
  int *indicesMax = INTEGER(indicesMaxSXP);
  
  double maxSum = 0;
  
  for(int temp = 0; temp < k; ++temp){
    indices[temp] = temp;
  }
  
  while(1){
    
    double sum = 0;
    for(int temp1 = 0; temp1 < k-2; ++temp1){
      for(int temp2 = temp1+1; temp2 < k-1; ++temp2){
        sum += X[indices[temp1]*nrowsX + indices[temp2]];
      }
    }
    
    if(sum > maxSum){
      maxSum = sum;
      for(int temp = 0; temp < k; ++temp){
        indicesMax[temp] = indices[temp]+1;
      }
    }
    
    int flag;
    for(int temp1 = k-1; temp1 >= 0; --temp1){
      flag = 0;
      indices[temp1] += 1;
      if(indices[temp1] < nrowsX-(k-1-temp1)){
        for(int temp2 = temp1+1; temp2 < k; ++temp2){
          indices[temp2] = indices[temp1] + (temp2 - temp1);
        }
        break;
      } else {
        flag = 1;
      }
    }
    
    if(flag == 1){
      break;
    }
  }
  
  SEXP listOutSXP = PROTECT(Rf_allocVector(VECSXP, 2));
  SEXP maxSumSXP = PROTECT(Rf_allocVector(REALSXP, 1));
  *REAL(maxSumSXP) = maxSum;
  
  SET_VECTOR_ELT(listOutSXP, 0, indicesMaxSXP);
  SET_VECTOR_ELT(listOutSXP, 1, maxSumSXP);
  
  UNPROTECT(3);
  
  return listOutSXP;
}
