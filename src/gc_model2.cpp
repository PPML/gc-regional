#include <RcppArmadillo.h>
#include "gc_model.h"
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
Rcpp::List gcSim2(
  List  x, //params
  double dt, // timestep
  List cm, //mixing
  int nYrs, //model run time
  NumericVector initPop, //initial population
  NumericVector n_sa, //size of sexually active population
  arma::mat births, //birth rate
  arma::mat births_sa, //births p_sexually active
  arma::mat births_nsa, //births p_not sexually active
  arma::mat aging, //aging rate
  arma::mat aging_nsa, //aging p.not sexually active
  bool debug,
  bool quarterly,
  bool increase_risk_all_pops
) {
  int i = 4; //number of subpopulations
  int j = 2; //number of sexes
  int k = 2; //number of sexual activity groups
  int l = 2; //number of age groups
  int dim = i*j*k*l; //total number of population groups in output matrix
  NumericVector b = as<NumericVector>(x["b"]); //transmission rate
  NumericVector pS_param = as<NumericVector>(x["symp"]); //probability susceptible
  NumericVector gamma_param = as<NumericVector>(x["gamma"]); //dur infectious, symptomatic
  NumericVector delta_param = as<NumericVector>(x["delta"]); //dur infectious, asymptomatic
  NumericMatrix rep_symp = as<NumericMatrix>(x["rep.symp"]); //prob reported, symptomatic
  NumericVector rep_asymp = as<NumericVector>(x["rep.asymp"]); //prob reported, asymptomatic
  NumericMatrix screen = as<NumericMatrix>(x["screen"]); //screening rate
  NumericVector c_msm = as<NumericVector>(x["behav"]); //rr transmit in msm (updated annually)
  NumericVector c_hetero = as<NumericVector>(x["behav.hetero"]); // increase in transmission for heterosexuals


  NumericVector psi(dim); //screening rate (updated annually)
  NumericVector rep_s(dim); //reporting prob symptomatic (updated annually)

  NumericVector cmlow1 = as<NumericVector>(cm["cm.low1"]); //contact matrix low activity, black
  NumericVector cmlow2 = as<NumericVector>(cm["cm.low2"]); //contact matrix low activity, other
  NumericVector cmlow3 = as<NumericVector>(cm["cm.low3"]); //contact matrix low activity, Hispanic
  NumericVector cmlow4 = as<NumericVector>(cm["cm.low4"]); //contact matrix low activity, MSM
  NumericVector cmhigh1 = as<NumericVector>(cm["cm.high1"]); //contact matrix high activity, black
  NumericVector cmhigh2 = as<NumericVector>(cm["cm.high2"]); //contact matrix high activity, other
  NumericVector cmhigh3 = as<NumericVector>(cm["cm.high3"]); //contact matrix high activity, Hispanic
  NumericVector cmhigh4 = as<NumericVector>(cm["cm.high4"]); //contact matrix high activity, MSM
  arma::mat birth_rate = births; //birth rate
  arma::mat aging_rate = aging; //aging rate
  arma::vec S(dim); //susceptible
  arma::vec Y(dim); //infectious, symptomatic
  arma::vec Z(dim); //infectious, asymptomatic
  arma::vec NSA(dim); //not sexually active
  arma::vec D(dim); //diagnosed and reported cases
  arma::vec R(dim); //diagnosed and reported symptomatic cases
  arma::vec CUMINC(dim); //cumulative incidence
  arma::vec gamma(dim);
  arma::vec delta(dim);
  arma::vec pS(dim);
  arma::vec lambda(dim); //force of infection
  arma::vec Pop(dim*8); //population matrix
  int out_rows;
  int out_cols;
  if (debug) {
    out_rows = nYrs*52*7;
  } else if (quarterly) {
    out_rows = nYrs*4+1;
  } else {
    out_rows = nYrs;
  }
  if (quarterly) {
    out_cols = dim*8+2;
  } else {
    out_cols = dim*8+1;
  }
  int counter_for_quarterly = 0;
  NumericMatrix outputs(out_rows, out_cols);
  NumericMatrix out(out_rows, out_cols);
  int increase_risk_all_pops_int;
  if (increase_risk_all_pops) {
    increase_risk_all_pops_int = 1;
  } else {
    increase_risk_all_pops_int = 0;
  }

  //update parameters for males
  for(int n=0; n<((i-1)*k*l); n++) {
    gamma[n] = gamma_param[0];
    delta[n] = delta_param[0];
    pS[n] = pS_param[0];
  }

  //update parameters for MSM
  for(int n=((i-1)*k*l); n<(i*k*l); n++) {
    gamma[n] = gamma_param[1];
    delta[n] = delta_param[1];
    pS[n] = pS_param[1];
  }

  //update parameters for females
  for(int n=((i*k*l)); n<(2*i*k*l); n++) {
    gamma[n] = gamma_param[2];
    delta[n] = delta_param[2];
    pS[n] = pS_param[2];
  }


  // In order to add the day loop, we make dt 7 times smaller
  dt = dt / 7;

  /// year loop
  for(int y=0; y<nYrs; y++){
    for(int w=0; w<52; w++){
      for (int d=0; d<7; d++){

        //update reporting and transmission rr
        double c_rr = c_msm[y];
        double c_rr_hetero = c_hetero[y];
        double rep_a = rep_asymp[y];

        if(y==0 && w==0 && d==0){   ///populate the model at t=0
        for(int n=0; n<dim*8; n++){
          Pop[n] = initPop[n];
        }
        }

        else{
          for(int n=0; n<dim; n++) { //calculate population size at start of each week
            S[n] = Pop[n];
            Y[n] = Pop[n+dim];
            Z[n] = Pop[n+2*dim];
            NSA[n] = Pop[n+3*dim];
            D[n] = Pop[n+4*dim];
            R[n] = Pop[n+5*dim];
            CUMINC[n] = Pop[n+6*dim];
            psi[n] = screen(y,n);  //update screening and reporting
            rep_s[n] = rep_symp(y,n);
          }


          //calculate force of infection for each sex and subpopulation
          for(int n=0; n<4; n++) { //males subpop 1
            int m = n;
            lambda[n] =
              ((b[0]) * (
                cmlow1[m]     * (Y[16] + Z[16]) / (n_sa[16]) + //F, low AC, age=1, subpop1
                cmlow1[m+4]   * (Y[18] + Z[18]) / (n_sa[18]) + //F, low AC, age=2, subpop1
                cmlow1[m+8]   * (Y[20] + Z[20]) / (n_sa[20]) +//F, low AC, age=1, subpop2
                cmlow1[m+12]  * (Y[22] + Z[22]) / (n_sa[22]) + //F, low AC, age=2, subpop2
                cmlow1[m+16]  * (Y[24] + Z[24]) / (n_sa[24]) + //F, low AC, age=1, subpop3
                cmlow1[m+20]  * (Y[26] + Z[26]) / (n_sa[26])) + //F, low AC, age=2, subpop3
               (b[0]  + (1.0-b[0])*c_rr_hetero*increase_risk_all_pops_int) * (
                cmhigh1[m]    * (Y[17] + Z[17]) / (n_sa[17]) +  //F, high AC, age=1, subpop1
                cmhigh1[m+4]  * (Y[19] + Z[19]) / (n_sa[19]) +  //F, high AC, age=2, subpop1
                cmhigh1[m+8]  * (Y[21] + Z[21]) / (n_sa[21]) + //F, high AC, age=1, subpop2
                cmhigh1[m+12] * (Y[23] + Z[23]) / (n_sa[23]) + //F, high AC, age=2, subpop2
                cmhigh1[m+16] * (Y[25] + Z[25]) / (n_sa[25]) + //F, high AC, age=1, subpop3
                cmhigh1[m+20] * (Y[27] + Z[27]) / (n_sa[27]))); //F, high AC, age=2, subpop3
          }

          for(int n=4; n<8; n++) {   //males subpop 2
            int m = n-4;
            lambda[n] =
              ((b[0]) * (
                cmlow2[m]     * (Y[16] + Z[16]) / (n_sa[16]) + //F, low AC, age=1, subpop1
                cmlow2[m+4]   * (Y[18] + Z[18]) / (n_sa[18]) + //F, low AC, age=2, subpop1
                cmlow2[m+8]   * (Y[20] + Z[20]) / (n_sa[20]) + //F, low AC, age=1, subpop2
                cmlow2[m+12]  * (Y[22] + Z[22]) / (n_sa[22]) + //F, low AC, age=2, subpop2
                cmlow2[m+16]  * (Y[24] + Z[24]) / (n_sa[24]) + //F, low AC, age=1, subpop3
                cmlow2[m+20]  * (Y[26] + Z[26]) / (n_sa[26])) + //F, low AC, age=2, subpop3
              (b[0] + (1.0-b[0])*c_rr_hetero*increase_risk_all_pops_int) * (
                cmhigh2[m]    * (Y[17] + Z[17]) / (n_sa[17]) + //F, high AC, age=1, subpop1
                cmhigh2[m+4]  * (Y[19] + Z[19]) / (n_sa[19]) + //F, high AC, age=2, subpop1
                cmhigh2[m+8]  * (Y[21] + Z[21]) / (n_sa[21]) + //F, high AC, age=1, subpop2
                cmhigh2[m+12] * (Y[23] + Z[23]) / (n_sa[23]) + //F, high AC, age=2, subpop2
                cmhigh2[m+16] * (Y[25] + Z[25]) / (n_sa[25]) + //F, high AC, age=1, subpop3
                cmhigh2[m+20] * (Y[27] + Z[27]) / (n_sa[27]))); //F, high AC, age=2, subpop3
          }

          for(int n=8; n<12; n++) {   //males subpop 3
            int m = n-8;
            lambda[n] =
              ((b[0]) * (
                cmlow3[m]     * (Y[16] + Z[16]) / (n_sa[16]) + //F, low AC, age=1, subpop1
                cmlow3[m+4]   * (Y[18] + Z[18]) / (n_sa[18]) + //F, low AC, age=2, subpop1
                cmlow3[m+8]   * (Y[20] + Z[20]) / (n_sa[20]) + //F, low AC, age=1, subpop2
                cmlow3[m+12]  * (Y[22] + Z[22]) / (n_sa[22]) + //F, low AC, age=2, subpop2
                cmlow3[m+16]  * (Y[24] + Z[24]) / (n_sa[24]) + //F, low AC, age=1, subpop3
                cmlow3[m+20]  * (Y[26] + Z[26]) / (n_sa[26])) + //F, low AC, age=2, subpop3
              (b[0] + (1.0-b[0])*c_rr_hetero*increase_risk_all_pops_int) * (
                cmhigh3[m]    * (Y[17] + Z[17]) / (n_sa[17]) + //F, high AC, age=1, subpop1
                cmhigh3[m+4]  * (Y[19] + Z[19]) / (n_sa[19]) + //F, high AC, age=2, subpop1
                cmhigh3[m+8]  * (Y[21] + Z[21]) / (n_sa[21]) + //F, high AC, age=1, subpop2
                cmhigh3[m+12] * (Y[23] + Z[23]) / (n_sa[23]) + //F, high AC, age=2, subpop2
                cmhigh3[m+16] * (Y[25] + Z[25]) / (n_sa[25]) + //F, high AC, age=1, subpop3
                cmhigh3[m+20] * (Y[27] + Z[27]) / (n_sa[27]))); //F, high AC, age=2, subpop3
          }

          for(int n=12; n<16; n++) {   //males subpop 4
            int m = n-12;
            lambda[n] =
              ((b[0]) *  (
                cmlow4[m]     * (Y[16] + Z[16]) / (n_sa[16]) + //F, low AC, age=1, subpop1
                cmlow4[m+4]   * (Y[18] + Z[18]) / (n_sa[18]) + //F, low AC, age=2, subpop1
                cmlow4[m+8]   * (Y[20] + Z[20]) / (n_sa[20]) + //F, low AC, age=1, subpop2
                cmlow4[m+12]  * (Y[22] + Z[22]) / (n_sa[22]) + //F, low AC, age=2, subpop2
                cmlow4[m+16]  * (Y[24] + Z[24]) / (n_sa[24]) + //F, low AC, age=1, subpop3
                cmlow4[m+20]  * (Y[26] + Z[26]) / (n_sa[26])) + //F, low AC, age=2, subpop3
              (b[0] + (1.0 - b[0])*c_rr_hetero*increase_risk_all_pops_int) *  (
                cmhigh4[m]    * (Y[17] + Z[17]) / (n_sa[17]) + //F, high AC, age=1, subpop1
                cmhigh4[m+4]  * (Y[19] + Z[19]) / (n_sa[19]) +  //F, high AC, age=2, subpop1
                cmhigh4[m+8]  * (Y[21] + Z[21]) / (n_sa[21]) + //F, high AC, age=1, subpop2
                cmhigh4[m+12] * (Y[23] + Z[23]) / (n_sa[23]) + //F, high AC, age=2, subpop2
                cmhigh4[m+16] * (Y[25] + Z[25]) / (n_sa[25]) + //F, high AC, age=1, subpop3
                cmhigh4[m+20] * (Y[27] + Z[27]) / (n_sa[27]))) +

                ((b[2] + (1.0-b[2])*c_rr)*
                   (cmlow4[m+24]  * (Y[12] + Z[12]) / (n_sa[12]) + //M, low AC, age=1, subpop4
                    cmlow4[m+28]  * (Y[14] + Z[14]) / (n_sa[14])) + //M, low AC, age=2, subpop4
                (b[2] + (1.0-b[2])*c_rr)*(
                    cmhigh4[m+24] * (Y[13] + Z[13]) / (n_sa[13]) +  //M, high AC, age=1, subpop4
                    cmhigh4[m+28] * (Y[15] + Z[15]) / (n_sa[15]))); //M, high AC, age=2, subpop4
          }


          for(int n=16; n<20; n++) {   //females subpop 1
            int m = n;
            lambda[n] =
              ((b[1]) * (
                cmlow1[m+16]     * (Y[0]  + Z[0])  /(n_sa[0])  + //M, low AC, age=1, subpop1
                cmlow1[m+4+16]   * (Y[2]  + Z[2])  /(n_sa[2])  + //M, low AC, age=2, subpop1
                cmlow1[m+8+16]   * (Y[4]  + Z[4])  /(n_sa[4])  + //M, low AC, age=1, subpop2
                cmlow1[m+12+16]  * (Y[6]  + Z[6])  /(n_sa[6])  + //M, low AC, age=2, subpop2
                cmlow1[m+16+16]  * (Y[8]  + Z[8])  /(n_sa[8])  + //M, low AC, age=1, subpop3
                cmlow1[m+20+16]  * (Y[10] + Z[10]) /(n_sa[10]) + //M, low AC, age=2, subpop3
                cmlow1[m+24+16]  * (Y[12] + Z[12]) /(n_sa[12]) + //M, low AC, age=1, subpop4
                cmlow1[m+28+16]  * (Y[14] + Z[14]) /(n_sa[14])) + //M, low AC, age=2, subpop4
              (b[1] + (1.0-b[1])*c_rr_hetero*increase_risk_all_pops_int) * (
                cmhigh1[m+16]    * (Y[1]  + Z[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
                cmhigh1[m+4+16]  * (Y[3]  + Z[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
                cmhigh1[m+8+16]  * (Y[5]  + Z[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
                cmhigh1[m+12+16] * (Y[7]  + Z[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
                cmhigh1[m+16+16] * (Y[9]  + Z[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
                cmhigh1[m+20+16] * (Y[11] + Z[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
                cmhigh1[m+24+16] * (Y[13] + Z[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
                cmhigh1[m+28+16] * (Y[15] + Z[15]) / (n_sa[15]))); //M, high AC, age=2, subpop4

          }

          for(int n=20; n<24; n++) {   //females subpop 2
            int m = n-4;
            lambda[n] =
              ((b[1]) * (
                  cmlow2[m+16]     * (Y[0]  + Z[0])  /(n_sa[0])  + //M, low AC, age=1, subpop1
                  cmlow2[m+4+16]   * (Y[2]  + Z[2])  /(n_sa[2])  + //M, low AC, age=2, subpop1
                  cmlow2[m+8+16]   * (Y[4]  + Z[4])  /(n_sa[4])  + //M, low AC, age=1, subpop2
                  cmlow2[m+12+16]  * (Y[6]  + Z[6])  /(n_sa[6])  + //M, low AC, age=2, subpop2
                  cmlow2[m+16+16]  * (Y[8]  + Z[8])  /(n_sa[8])  + //M, low AC, age=1, subpop3
                  cmlow2[m+20+16]  * (Y[10] + Z[10]) /(n_sa[10]) + //M, low AC, age=2, subpop3
                  cmlow2[m+24+16]  * (Y[12] + Z[12]) /(n_sa[12]) + //M, low AC, age=1, subpop4
                  cmlow2[m+28+16]  * (Y[14] + Z[14]) /(n_sa[14])) + //M, low AC, age=2, subpop4
              (b[1] + (1.0-b[1])*c_rr_hetero*increase_risk_all_pops_int) * (
                  cmhigh2[m+16]    * (Y[1]  + Z[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
                  cmhigh2[m+4+16]  * (Y[3]  + Z[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
                  cmhigh2[m+8+16]  * (Y[5]  + Z[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
                  cmhigh2[m+12+16] * (Y[7]  + Z[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
                  cmhigh2[m+16+16] * (Y[9]  + Z[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
                  cmhigh2[m+20+16] * (Y[11] + Z[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
                  cmhigh2[m+24+16] * (Y[13] + Z[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
                  cmhigh2[m+28+16] * (Y[15] + Z[15]) / (n_sa[15]))); //M, high AC, age=2, subpop4
          }

          for(int n=24; n<28; n++) {   //females subpop 3
            int m = n-8;
            lambda[n] =
              ((b[1]) *  (
                  cmlow3[m+16]     * (Y[0]  + Z[0])  /(n_sa[0])  + //M, low AC, age=1, subpop1
                  cmlow3[m+4+16]   * (Y[2]  + Z[2])  /(n_sa[2])  + //M, low AC, age=2, subpop1
                  cmlow3[m+8+16]   * (Y[4]  + Z[4])  /(n_sa[4])  + //M, low AC, age=1, subpop2
                  cmlow3[m+12+16]  * (Y[6]  + Z[6])  /(n_sa[6])  + //M, low AC, age=2, subpop2
                  cmlow3[m+16+16]  * (Y[8]  + Z[8])  /(n_sa[8])  + //M, low AC, age=1, subpop3
                  cmlow3[m+20+16]  * (Y[10] + Z[10]) /(n_sa[10]) + //M, low AC, age=2, subpop3
                  cmlow3[m+24+16]  * (Y[12] + Z[12]) /(n_sa[12]) + //M, low AC, age=1, subpop4
                  cmlow3[m+28+16]  * (Y[14] + Z[14]) /(n_sa[14])) + //M, low AC, age=2, subpop4
              (b[1] + (1.0-b[1])*c_rr_hetero*increase_risk_all_pops_int) *  (
                  cmhigh3[m+16]    * (Y[1]  + Z[1])  / (n_sa[1])  + //M, high AC, age=1, subpop1
                  cmhigh3[m+4+16]  * (Y[3]  + Z[3])  / (n_sa[3])  + //M, high AC, age=2, subpop1
                  cmhigh3[m+8+16]  * (Y[5]  + Z[5])  / (n_sa[5])  + //M, high AC, age=1, subpop2
                  cmhigh3[m+12+16] * (Y[7]  + Z[7])  / (n_sa[7])  + //M, high AC, age=2, subpop2
                  cmhigh3[m+16+16] * (Y[9]  + Z[9])  / (n_sa[9])  + //M, high AC, age=1, subpop3
                  cmhigh3[m+20+16] * (Y[11] + Z[11]) / (n_sa[11]) + //M, high AC, age=2, subpop3
                  cmhigh3[m+24+16] * (Y[13] + Z[13]) / (n_sa[13]) + //M, high AC, age=1, subpop4
                  cmhigh3[m+28+16] * (Y[15] + Z[15]) / (n_sa[15]))); //M, high AC, age=2, subpop4

          }

          // population aging and movement in/out of sexually active group
          arma::vec agingS = mv_mult(aging_rate,S);
          arma::vec agingY = mv_mult(aging_rate,Y);
          arma::vec agingZ = mv_mult(aging_rate,Z);
          arma::mat agingNSAmat = mat_mult(aging_nsa, aging_rate);
          arma::vec agingNSA = mv_mult(agingNSAmat, NSA); //aging from NSA to NSA
          arma::mat agingSNSAmat = aging_rate - mat_mult(aging_nsa, aging_rate);
          arma::vec agingSNSA = mv_mult(agingSNSAmat, NSA); //aging from NSA to SA

          arma::mat birthsSmat = mat_mult(birth_rate, births_sa);
          arma::vec Ntot = S+Y+Z+NSA;
          arma::vec birthsS = mv_mult(birthsSmat,Ntot); //births into S
          arma::mat birthsNSAmat = mat_mult(birth_rate, births_nsa);
          arma::vec birthsNSA = mv_mult(birthsNSAmat,Ntot); //births into NSA

          // Initialize the equations and integration
          for(int n=0; n<dim; n++) {
            Pop[n] += (-lambda[n]*S[n] + gamma[n]*Y[n] + delta[n]*Z[n] + psi[n]*Z[n]  + agingS[n] + agingSNSA[n] + birthsS[n]) * dt ;  // dS/dt
            Pop[n+dim] += (pS[n]*lambda[n]*S[n] - gamma[n]*Y[n] + agingY[n]) * dt; // dY/dt
            Pop[n+dim*2] += ( (1-pS[n])*lambda[n]*S[n] - delta[n]*Z[n] - psi[n]*Z[n] + agingZ[n] ) * dt; //dZ/dt
            Pop[n+dim*3] += (agingNSA[n] + birthsNSA[n]) * dt; //dNSA/dt
            Pop[n+dim*4] += (rep_s[n]*gamma[n]*Y[n] + rep_a*psi[n]*Z[n]) * dt; //dD/dt
            Pop[n+dim*5] += (rep_s[n]*gamma[n]*Y[n] ) * dt; //dR/dt
            Pop[n+dim*6] += (lambda[n]*S[n])*dt ; //dCUMINC/dt
            Pop[n+dim*7] += ( psi[n]*( S[n] + Z[n]) ) * dt; //dScr/dt
          }
        }

        // fill results table with annual outputs
        if (debug) {
          outputs(y*52*7+w*7+d, 0) = y*52*7+w*7+d;
          for (int n=0; n<dim*8; n++) {
            outputs(y*52*7+w*7+d, 1+n) = Pop[n];
          }
        } else if (quarterly && (((w == 0 || w == 12 || w == 25 || w == 38) && d == 0) || (y == nYrs - 1 && w == 51 && d == 6))) {
          outputs(counter_for_quarterly, 0) = y;
          outputs(counter_for_quarterly, 1) = w + 1;
          for (int n = 0; n<dim*8; n++) {
            outputs(counter_for_quarterly, 2+n) = Pop[n];
          }
          counter_for_quarterly += 1;
        } else if (!quarterly && w==32 && d==0) {
          outputs(y,0) = y; //year
          for(int n=0; n<dim*8; n++) {
            outputs(y,1+n) = Pop[n]; //
          }
        }

      } // end of day loop
    } /// end of week loop
  } /// end of year loop

  // for(int m=0; m<nYrs; m++) {
    //   for(int n=0; n<dim*8+1; n++) {
      //     out(m,n) = outputs(m,n);
      //   }
    // }

  return Rcpp::List::create(
    Rcpp::Named("out") = outputs
  );
}
