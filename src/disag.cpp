#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::NumericVector matrixSubcol(Rcpp::NumericMatrix x, int i0, int i1, int icol){

  // Determine the length of the output
  int n = i1-i0+1;

  // Create an output vector
  Rcpp::NumericVector out(n);

  // Loop through the matrix
  for(int z = 0; z < n; ++z) {
    out(z) = x(i0+z,icol);
  }

  return out;
}

// [[Rcpp::export]]
void rcpp_rprintf(NumericVector v){
  // printing values of all the elements of Rcpp vector
  for(int i=0; i<v.length(); ++i){
    Rprintf("the value of v[%i] : %f \n", i, v[i]);
  }
}


// [[Rcpp::export]]
IntegerVector order_(NumericVector x) {
  NumericVector sorted = clone(x).sort();

  return match(sorted, x);
}

// [[Rcpp::export]]
bool any_sug(LogicalVector x){
  // Note the use of is_true to return a bool type
  return is_true(any(x == TRUE));
}


// [[Rcpp::export]]
bool all_sug(LogicalVector x) {
  // Note the use of is_true to return a bool type.
  return is_true(all(x == TRUE));
}


// [[Rcpp::export]]
int find_row(NumericMatrix m, NumericVector v) {
  // a function that returns the index of the row which is equal to the vector vec
  for (int i = 0; i < m.nrow(); i++) {
    if (all_sug(m(i,_) == v)) {
      return i;
    }
  }
  return -1;
}

// [[Rcpp::export]]
NumericVector getrmseiTemp(int Npdt,
                           IntegerVector mObs,
                           IntegerVector mSim,
                           NumericMatrix Yobs24,
                           NumericMatrix Ysim24,
                           int i){
  // getrmseiTemp returns the vector of analog scores for the simulated day i (nTobs value)
  // - Npdt: the number of time steps at the disaggregated scale (e.g. 24 from days to hours)
  // - mObs: index of the month (1..12) for the matrix of daily observations
  // - mSim: index of the month (1..12) for the matrix of daily simulations
  // - Yobs24: Matrix of observed temperatures day periods: nTobs [days] x nStat [stations]
  // - Ysim24: Matrix of simulated temperatures day periods: nTsim [days] x nStat [stations]
  // - i is a simulation day

  // Input/Output
  int nTobs = Yobs24.nrow();
  int nTsim = Ysim24.nrow();

  // locals
  bool notSameClass;
  NumericVector rmseI(nTobs);

  // if this is the first or last simulated steps, we just compute scores on the
  // current temperatures
  if(i== 0 || i==(nTsim-1)) {
    for (int j = 0; j < nTobs; j++) {
      //same month
      notSameClass = (mSim(i)/=mObs(j));
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      if(is_true(any(is_na(Yobs24(j,_)))) || notSameClass) {
        rmseI(j) = 1E30;
      } else {
        //absolute differences between temperature for this day
        rmseI(j) = sum(abs(Ysim24(i,_)-Yobs24(j,_)));
      }
    }
  }else{

    //discard first and last elements
    rmseI(0) = 2E30;
    rmseI(nTobs-1) = 2E30;

    for (int j = 1; j < (nTobs-1); j++) {
      //same month
      notSameClass = (mSim(i)/=mObs(j));
      rmseI(j) = 0;
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      for(int lag = -1; lag<2;lag++){

        //if any observed value is missing in the observed prec. field or if the months do not correspond
        //we discard this observed day
        if(is_true(any(is_na(Yobs24(j+lag,_)))) || notSameClass) {
          rmseI(j) = 1E30;
        } else {
          //absolute differences between temperature for this day
          rmseI(j) = rmseI(j) + sum(abs(Ysim24(i+lag,_)-Yobs24(j+lag,_)));
        }
      }
    }
  }
  return rmseI;
}

// [[Rcpp::export]]
List disagTempMOF(int Npdt,
                  IntegerVector mObs,
                  IntegerVector mSim,
                  NumericMatrix YobsXX,
                  NumericMatrix Yobs24,
                  NumericMatrix Ysim24){
  // disagTempMOF returns disaggregaed values
  // - Npdt: the number of time steps at the disaggregated scale (e.g. 24 from days to hours)
  // - mObs: index of the month (1..12) for the matrix of daily observations
  // - mSim: index of the month (1..12) for the matrix of daily simulations
  // - YobsXX: Matrix of observed temperatures at fine time scales: nTobs [days*Npdt] x nStat [stations]
  // - Yobs24: Matrix of observed temperatures day periods: nTobs [days] x nStat [stations]
  // - Ysim24: Matrix of simulated temperatures day periods: nTsim [days] x nStat [stations]

  // Input/Output
  int nTobs = Yobs24.nrow();
  int nStat = Yobs24.ncol();
  int nTsim = Ysim24.nrow();

  // outputs
  // Matrix of disag. simulated amounts: nTsim*Npdt [days] x nStat [stations]
  NumericMatrix Ysim(nTsim*Npdt,nStat);
  //Matrix indicating how it has been disaggregated
  NumericVector codeDisag(nTsim);

  // Locals
  NumericVector YobsXXD(Npdt);
  NumericVector rmseI(nTobs);
  IntegerVector indBestRMSEI(nTobs);
  int indBestFieldI;
  int jXX; // index of the observed day

  //For each simulated 3-day period
  for (int i = 0; i < nTsim; i++) {

    rmseI = getrmseiTemp(Npdt,mObs,mSim,Yobs24,Ysim24,i);

    // get index of minimum RMSE values
    // !!! order_ starts from 1, which means that indBestFieldI is shifted by 1
    // in terms of index of the matrices
    indBestRMSEI = order_(rmseI)-1;
    indBestFieldI = indBestRMSEI[0];

    //index just before period in the observed XXH matrix
    jXX = indBestFieldI*Npdt-1;

    //update code of disaggregation
    codeDisag(i) = indBestFieldI+1;

    //======= disaggregate =====
    for (int k = 0; k < nStat; k++) {

      // simulated temperatures at the fine scale are observed
      // fine temperatures for the analog day, shifted to match the
      // simulated daily temperature
      YobsXXD = matrixSubcol(YobsXX,jXX+1,jXX+Npdt,k);

      // loop over the Npdt time steps
      for (int ii = 0; ii < Npdt; ii++) {
        Ysim(i*Npdt+ii,k) = YobsXXD(ii) - sum(YobsXXD)/YobsXXD.length() + Ysim24(i,k);
      }
    }
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("codeDisag") = codeDisag,
                                      Rcpp::Named("Ysim") = Ysim);
  return res;
}