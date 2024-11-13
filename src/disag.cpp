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
NumericVector getrmseiPrec(int Npdt,
                           IntegerVector mObs,
                           IntegerVector mSim,
                           IntegerVector cObs,
                           IntegerVector cSim,
                           NumericMatrix Yobs24,
                           NumericMatrix Ysim24,
                           int i,
                           int nLagScore){

  // getrmseiPrec returns the vector of analog scores for the simulated day i (nTobs value)
  // - Npdt: the number of time steps at the disaggregated scale (e.g. 24 from days to hours)
  // - mObs: index of the month (1..12) for the matrix of daily observations
  // - mSim: index of the month (1..12) for the matrix of daily simulations
  // - cObs: class (=1,2,3,4) of precip for observations
  // - cSim: class (=1,2,3,4) of precip for simulations
  // - Yobs24: Matrix of observed temperatures day periods: nTobs [days] x nStat [stations]
  // - Ysim24: Matrix of simulated temperatures day periods: nTsim [days] x nStat [stations]
  // - nLagScore: lag=1
  // - i is a simulation day
  // - nLagScore is an integer

  // Input/Output
  int nTobs = Yobs24.nrow();
  int nStat = Yobs24.ncol();

  // Locals
  const double& naVal=-9999.;
  double rmseIJ;
  NumericVector adimObs(nStat);
  NumericVector adimSim(nStat);
  NumericVector rmseI(nTobs);
  bool notSameClass;


  //======= Step 1: Minimize Score for this selection
  //sum of square differences for the closest fields on two time steps

  //for the first time steps
  if(i < nLagScore) {

    for (int j = 0; j < nTobs; j++) {

      //same month and class
      notSameClass = (mSim(i)/=mObs(j)) || (cSim(i)/=cObs(j));
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      if((any_sug(Yobs24(j,_) == naVal)) || notSameClass) {
        rmseI(j) = 1E30;
      } else {
        //absolute differences between adimensioned precipitation for this day
        if(sum(Ysim24(i,_))==0) {
          adimSim = Ysim24(i,_);
        } else {
          adimSim = Ysim24(i,_)/sum(Ysim24(i,_));
        }

        if(sum(Yobs24(j,_))==0) {
          adimObs = Yobs24(j,_);
        } else {
          adimObs = Yobs24(j,_)/sum(Yobs24(j,_));
        }
        rmseI(j) = sum(abs(adimSim-adimObs));
      }

    }
  } else {

    //discard first elements
    for (int j = 0; j < nLagScore; j++) {
      rmseI(j) = 2E30;
    }


    //for the next days, compute score
    for (int j = nLagScore; j < nTobs; j++) {

      //same month and class
      notSameClass = (mSim(i)!=mObs(j)) || (cSim(i)!=cObs(j));
      //if any observed value is missing in the observed prec. field or if the months do not correspond
      //we discard this observed day
      if(any_sug(Yobs24(j,_) == naVal) || notSameClass) {
        rmseIJ = 3E30;
      } else {
        //absolute differences between adimensioned precipitation for this day
        if(sum(Ysim24(i,_))==0) {
          adimSim = Ysim24(i,_);
        } else {
          adimSim = Ysim24(i,_)/sum(Ysim24(i,_));
        }

        if(sum(Yobs24(j,_))==0) {
          adimObs = Yobs24(j,_);
        } else {
          adimObs = Yobs24(j,_)/sum(Yobs24(j,_));
        }
        rmseIJ = sum(abs(adimSim-adimObs));
      }


      // add differences for the previous days, just non na values
      // jDay is the index in the daily matrix just before the current
      // 3-day step. We  roll over the nLag previous daily steps starting
      // from jDay to jDay-nLag+1 (iLag = 0,...,(nLag-1))
      for (int iLag = 1; iLag <= nLagScore; iLag++) {
        if(any_sug(Yobs24(j-iLag,_)==naVal)) {
          rmseIJ = 4E30;
          break;
        } else {
          //absolute differences between adimensioned precipitation for this day
          if(sum(Ysim24(i-iLag,_))==0) {
            adimSim = Ysim24(i-iLag,_);
          } else {
            adimSim = Ysim24(i-iLag,_)/sum(Ysim24(i-iLag,_));
          }

          if(sum(Yobs24(j-iLag,_))==0) {
            adimObs = Yobs24(j-iLag,_);
          } else {
            adimObs = Yobs24(j-iLag,_)/sum(Yobs24(j-iLag,_));
          }

          rmseIJ = rmseIJ + sum(abs(adimSim-adimObs));
        }
      }

      rmseI(j) = rmseIJ;
    }
  }
  return rmseI;
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




// [[Rcpp::export]]
List disagPrecMOF(int Npdt,
                  IntegerVector mObs,
                  IntegerVector mSim,
                  IntegerVector cObs,
                  IntegerVector cSim,
                  NumericMatrix YobsXX,
                  NumericMatrix Yobs24,
                  NumericMatrix Ysim24,
                  int nLagScore){
  // a function that disagregate daily precipitations with subdaily precipitation
  // For each daily simulations at n gauges, we find the closest field in observations for which
  // subdaily precipitation structures are available.

  // - Npdt: the number of time steps at the disaggregated scale (e.g. 24 from days to hours)
  // - mObs: index of the month (1..12) for the matrix of daily observations
  // - mSim: index of the month (1..12) for the matrix of daily simulations
  // - cObs: class (=1,2,3,4) of precip for observations
  // - cSim: class (=1,2,3,4) of precip for simulations
  // - YobsXX: Matrix of observed temperatures at fine time scales: nTobs [days*Npdt] x nStat [stations]
  // - Yobs24: Matrix of observed temperatures day periods: nTobs [days] x nStat [stations]
  // - Ysim24: Matrix of simulated temperatures day periods: nTsim [days] x nStat [stations]
  // - nLagScore: number of preceeding days used to calcule the scores, can be 0, 1 or 2

  // Input/Output
  int nTobs = Yobs24.nrow();
  int nStat = Yobs24.ncol();
  int nTsim = Ysim24.nrow();

  // outputs
  NumericMatrix Ysim(nTsim*Npdt,nStat); // Matrix of disag. simulated amounts: nTsim*Npdt [days] x nStat [stations]
  NumericMatrix codeDisag(nTsim,nStat); //Matrix indicating how it has been disaggregated

  // Locals
  const int& nBestField=10;
  const double& naVal=-9999.;
  double rmseIJ;
  NumericVector adimObs(nStat);
  NumericVector adimSim(nStat);
  NumericVector YobsXXD(Npdt);
  NumericVector rmseI(nTobs);
  IntegerVector indBestRMSEI(nTobs);
  IntegerVector indBestFieldI(nBestField);
  bool notSameClass;
  int jXX; //index of the beginning of the period in the observed YobsXX matrix

  //For each simulated 3-day period
  for (int i = 0; i < nTsim; i++) {

    // get vector of rmse (similarity observed and simulated fields)
    rmseI = getrmseiPrec(Npdt,mObs,mSim,cObs,cSim,Yobs24,Ysim24,i,nLagScore);

    // get index of minimum RMSE values
    // !!! order_ starts from 1, which means that indBestFieldI is shifted by 1
    // in terms of index of the precipitation matrices
    indBestRMSEI = order_(rmseI)-1;
    indBestFieldI = indBestRMSEI[Range(0,nBestField-1)];

    // Look at the different case and disaggregate
    for (int k = 0; k < nStat; k++) {

      //initialise code
      codeDisag(i,k) = naVal;


      //////////case 1: no occurrence for this period and this station, nothing to disaggregate.
      if(Ysim24(i,k)==0) {
        codeDisag(i,k) = -1000;
        for (int ii = 0; ii < Npdt; ii++) {
          Ysim(i*Npdt+ii,k) = 0;
        }

        //////////case 2: look at the closest fields, if there is a positive
        // precipitation, we take the observed temporal structure
      } else {

        for (int j = 0; j < nBestField; j++) {


          //index of the beginning of the period in the observed XXH matrix
          jXX = indBestFieldI(j)*Npdt;


          //check if there is an observed precitation for the selected observed
          // field and for this station
          double Yobs24jk = Yobs24(indBestFieldI(j),k);
          if(Yobs24jk>0) {

            //update code of disaggregation
            codeDisag(i,k) = j+1;

            // simulated temperatures at the fine scale are observed
            // fine temperatures for the analog day, shifted to match the
            // simulated daily temperature
            YobsXXD = matrixSubcol(YobsXX,jXX,jXX+Npdt-1,k);

            // disaggregation
            for (int ii = 0; ii < Npdt; ii++) {
              Ysim(i*Npdt+ii,k) = YobsXXD(ii)*Ysim24(i,k)/sum(YobsXXD);
            }

            //codes to analyse how large subdaily intensities are disaggregated
            if(any_sug(matrixSubcol(Ysim,i*Npdt,i*Npdt+Npdt-1,k)>max(YobsXX(_,k)))){
              codeDisag(i,k) = codeDisag(i,k) + 10000;
            }
            //get out of the loop "for (int j = 0; j < nBestField; j++)"
            break;
          }
        }
      }

      //////////case 3: if we did not find similar structure then we find, for
      // this station, days with similar amounts
      if(codeDisag(i,k)==naVal) {
        //for the first time step
        if(i < nLagScore) {
          for (int j = 0; j < nTobs; j++) {
            //same month and class
            notSameClass = (mSim(i)/=mObs(j)) || (cSim(i)/=cObs(j));
            //if any observed value is missing in the observed prec. field or if the months do not correspond
            //we discard this observed day
            if((Yobs24(j,k) == naVal) || notSameClass) {
              rmseI(j) = 1E30;
            } else {
              rmseI(j) = abs(Ysim24(i,k)-Yobs24(j,k));
            }
          }
        } else {
          //discard first elements
          for (int j = 0; j < nLagScore; j++) {
            rmseI(j) = 1E30;
          }

          //for the next days, compute score
          for (int j = nLagScore; j < nTobs; j++) {
            //same month and class
            notSameClass = (mSim(i)!=mObs(j)) || (cSim(i)!=cObs(j));
            //if any observed value is missing in the observed prec. field or if the months do not correspond
            //or if there is no observed precip for this day
            //we discard this observed day
            if(((Yobs24(j,k) == naVal) || (Yobs24(j,k) == 0)) || notSameClass) {
              rmseI(j) = 1E30;
            } else {
              rmseIJ = abs(Ysim24(i,k)-Yobs24(j,k));
              //add differences for the previous days, but not NA values
              if(any_sug(matrixSubcol(YobsXX,j-nLagScore,j-1,k)==naVal)) {
                rmseIJ = 1E30;
              } else {
                for (int iLag = 0; iLag < nLagScore; iLag++) {
                  rmseIJ = rmseIJ + abs(Ysim24(i-iLag-1,k)-Yobs24(j-iLag-1,k));
                }
              }
              rmseI(j) = rmseIJ;
            }
          }
        }

        // get index of minimum RMSE values
        // !!! order_ starts from 1, which means that indBestFieldI is shifted by 1
        // in terms of index of the precipitation matrices
        indBestRMSEI = order_(rmseI)-1;
        indBestFieldI = indBestRMSEI[Range(0,nBestField-1)];

        codeDisag(i,k) = 2000;
        int rndDay = floor(Rcpp::runif(1)[0]*10);
        int jDay = indBestFieldI(rndDay)*Npdt;
        YobsXXD = matrixSubcol(YobsXX,jDay,jDay+Npdt-1,k);

        // disaggregation
        for (int ii = 0; ii < Npdt; ii++) {
          Ysim(i*Npdt+ii,k) = YobsXXD(ii)*Ysim24(i,k)/sum(YobsXXD);
        }

        //codes to analyse how large 3-day intensities are disaggregated
        if(any_sug(matrixSubcol(Ysim,i*Npdt,i*Npdt+Npdt-1,k)>max(YobsXX(_,k)))){
          codeDisag(i,k) = codeDisag(i,k) + 10000;
        }
      }
    }
  }

  Rcpp::List res = Rcpp::List::create(Rcpp::Named("codeDisag") = codeDisag,
                                      Rcpp::Named("Ysim") = Ysim);
  return res;
}

