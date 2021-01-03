#include <R.h>
#include <Rinternals.h>
#include <math.h>

static inline double cor_boot_vec(double*, double*, int*, int*, int, int);

SEXP get_cor(SEXP XSXP, SEXP repsSXP, SEXP indicesSXP){
  
  double *X = REAL(XSXP);
  int nrowsX = Rf_nrows(XSXP);
  
  int *reps = INTEGER(repsSXP);
  int *indices = INTEGER(indicesSXP);
  
  SEXP corSXP = PROTECT(Rf_allocVector(REALSXP, 1));
  *REAL(corSXP) = cor_boot_vec(X, X, reps, indices, Rf_xlength(repsSXP), nrowsX);
  
  UNPROTECT(1);
  return corSXP;
}

static inline double cor_boot_vec(double *X, double *Y, int *reps, int *indices, int length, int total){
  
  double sumX = 0;
  double sumY = 0;
  double sumXX = 0;
  double sumYY = 0;
  double sumXY = 0;
  
  switch(length){
    case 1 :
      {
        double tempX = X[indices[0]-1];
        double tempY = Y[indices[0]-1];
        int tempReps = reps[0];
        
        sumX = tempX*tempReps;
        sumY = tempY*tempReps;
        sumXX = tempX*tempX*tempReps;
        sumYY = tempY*tempY*tempReps;
        sumXY = tempX*tempY*tempReps;
      }
      break;
      
    case 2 :
      {
        double tempX1 = X[indices[0]-1];
        double tempX2 = X[indices[1]-1];
        double tempY1 = Y[indices[0]-1];
        double tempY2 = Y[indices[1]-1];
        int tempReps1 = reps[0];
        int tempReps2 = reps[1];
        
        sumX = tempX1*tempReps1+tempX2*tempReps2;
        sumY = tempY1*tempReps1+tempY2*tempReps2;
        sumXX = tempX1*tempX1*tempReps1+tempX2*tempX2*tempReps2;
        sumYY = tempY1*tempY1*tempReps1+tempY2*tempY2*tempReps2;
        sumXY = tempX1*tempY1*tempReps1+tempX2*tempY2*tempReps2;
      }
      break;
      
    default :
      {
        int levels = 0;
        for(int temp = length-1; temp > 0; ++levels){
          temp /= 2;
        }
        
        int indicesUpper[levels];
        
        double sumsX[levels-1];
        double sumsY[levels-1];
        double sumsXX[levels-1];
        double sumsYY[levels-1];
        double sumsXY[levels-1];
        
        indicesUpper[0] = length-1;
        
        int ind1 = 1;
        for(int indexLower = 0; indexLower < length; indexLower = indicesUpper[levels-1]+1){
          
          int ind2;
          for(ind2 = ind1; indicesUpper[ind2-1] > indexLower+2 && ind2 < levels; ++ind2){
            indicesUpper[ind2] = (indicesUpper[ind2-1]+indexLower)/2;
          }
          for(; ind2 < levels; ++ind2){
            indicesUpper[ind2] = indicesUpper[ind2-1];
          }
          
          if(indexLower < indicesUpper[levels-1]){
            
            double tempX1 = X[indices[indexLower]-1];
            double tempX2 = X[indices[indexLower+1]-1];
            double tempY1 = Y[indices[indexLower]-1];
            double tempY2 = Y[indices[indexLower+1]-1];
            int tempReps1 = reps[indexLower];
            int tempReps2 = reps[indexLower+1];
            
            sumX = tempX1*tempReps1+tempX2*tempReps2;
            sumY = tempY1*tempReps1+tempY2*tempReps2;
            sumXX = tempX1*tempX1*tempReps1+tempX2*tempX2*tempReps2;
            sumYY = tempY1*tempY1*tempReps1+tempY2*tempY2*tempReps2;
            sumXY = tempX1*tempY1*tempReps1+tempX2*tempY2*tempReps2;
          }
          else{
            
            double tempX = X[indices[indexLower]-1];
            double tempY = Y[indices[indexLower]-1];
            int tempReps = reps[indexLower];
            
            sumX = tempX*tempReps;
            sumY = tempY*tempReps;
            sumXX = tempX*tempX*tempReps;
            sumYY = tempY*tempY*tempReps;
            sumXY = tempX*tempY*tempReps;
          }
          
          int ind1;
          for(ind1 = levels-1; indicesUpper[ind1] == indicesUpper[ind1-1] && ind1 > 0; --ind1){
            
            sumX += sumsX[ind1-1];
            sumY += sumsY[ind1-1];
            sumXX += sumsXX[ind1-1];
            sumYY += sumsYY[ind1-1];
            sumXY += sumsXY[ind1-1];
            
            sumsX[ind1-1] = 0;
            sumsY[ind1-1] = 0;
            sumsXX[ind1-1] = 0;
            sumsYY[ind1-1] = 0;
            sumsXY[ind1-1] = 0;
          }
          
          if(ind1 == 0){
            break;
          }
          
          sumsX[ind1-1] = sumX;
          sumsY[ind1-1] = sumY;
          sumsXX[ind1-1] = sumXX;
          sumsYY[ind1-1] = sumYY;
          sumsXY[ind1-1] = sumXY;
      }
      break;
  }
  
  return sumXY;
}