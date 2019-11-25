gcSimR <- function(x,
                   dt,
                   cm,
                   nYrs,
                   initPop,
                   n_sa,
                   births,
                   births_sa,
                   births_nsa,
                   aging,
                   aging_nsa,
                   debug = F) {
  # __Set Up Memory Allocation and Variables ----
  i = 4 #number of subpopulations
  j = 2 #number of sexes
  k = 2 #number of sexual activity groups
  l = 2 #number of age groups
  dim = i*j*k*l #total number of population groups in output matrix
  b = (x[["b"]]) #transmission rate
  pS_param = (x[["symp"]]) #probability susceptible
  gamma_param = (x[["gamma"]]) #dur infectious, symptomatic
  delta_param = (x[["delta"]]) #dur infectious, asymptomatic
  rep_symp = (x[["rep.symp"]]) #prob reported, symptomatic
  rep_asymp = (x[["rep.asymp"]]) #prob reported, asymptomatic
  screen = (x[["screen"]]) #screening rate
  c_msm = (x[["behav"]]) #rr transmit in msm (updated annually)

  psi <- numeric(length = dim)
  rep_s <- numeric(length = dim)
  cmlow1 = (cm[["cm.low1"]]) #contact matrix low activity, black
  cmlow2 = (cm[["cm.low2"]]) #contact matrix low activity, other
  cmlow3 = (cm[["cm.low3"]]) #contact matrix low activity, Hispanic
  cmlow4 = (cm[["cm.low4"]]) #contact matrix low activity, MSM
  cmhigh1 = (cm[["cm.high1"]]) #contact matrix high activity, black
  cmhigh2 = (cm[["cm.high2"]]) #contact matrix high activity, other
  cmhigh3 = (cm[["cm.high3"]]) #contact matrix high activity, Hispanic
  cmhigh4 = (cm[["cm.high4"]]) #contact matrix high activity, MSM

  birth_rate = births #birth rate
  aging_rate = aging #aging rate
  S = numeric(dim) #susceptible
  Y <- numeric(dim) #infectious, symptomatic
  Z <- numeric(dim) #infectious, asymptomatic
  NSA <- numeric(dim) #not sexually active
  D <- numeric(dim) #diagnosed and reported cases
  R <- numeric(dim) #diagnosed and reported symptomatic cases
  CUMINC <- numeric(dim) #cumulative incidence
  gamma <- numeric(dim)
  delta <- numeric(dim)
  pS <- numeric(dim)
  lambda <- numeric(dim) #force of infection
  Pop <- numeric(dim*8) #population matrix
  outputs <- matrix(nrow = nYrs, ncol = dim*8+1)
  if (debug)
    outputs <- matrix(nrow = nYrs*52*7, ncol = dim*8+1)

  for(n in (1:((i-1)*k*l))) {
    gamma[[n]] = gamma_param[[1]]
    delta[[n]] = delta_param[[1]]
    pS[[n]] = pS_param[[1]]
  }

  for(n in ((i-1)*k*l+1):(i*k*l)) {
    gamma[n] = gamma_param[2]
    delta[n] = delta_param[2]
    pS[n] = pS_param[2]
  }

  for(n in ((i*k*l)+1):(2*i*k*l)) {
    gamma[n] = gamma_param[3]
    delta[n] = delta_param[3]
    pS[n] = pS_param[3]
  }
  dt = dt / 7

  # __Year Loop ----
  for (y in 1:nYrs) {
    # __Week Loop ----
    for (w in 1:52) {
      # __Day Loop ----
      for (d in 1:7) {
        c_rr = c_msm[y]
        rep_a = rep_asymp[y]

        if (y == 1 && w == 1 && d == 1) {
          for (n in 1:(dim*7))
            Pop[n] = initPop[n]
        } else {
          for (n in 1:dim) {
            S[n] = Pop[n]
            Y[n] = Pop[n+dim]
            Z[n] = Pop[n+2*dim]
            NSA[n] = Pop[n+3*dim]
            D[n] = Pop[n+4*dim]
            R[n] = Pop[n+5*dim]
            CUMINC[n] = Pop[n+6*dim]
            psi[n] = screen[y,n]  #update screening and reporting
            rep_s[n] = rep_symp[y,n]
          }

          # __Calculate Force of Infection ----
          for(n in 1:4) { #males subpop 1
            m = n;
            lambda[n] =  (b[1] * (cmlow1[m]     * (Y[16+1] + Z[16+1]) / (n_sa[16+1]) + #F, low AC, age=1, subpop1
                                    cmlow1[m+4]   * (Y[18+1] + Z[18+1]) / (n_sa[18+1]) + #F, low AC, age=2, subpop1
                                    cmlow1[m+8]   * (Y[20+1] + Z[20+1]) / (n_sa[20+1]) +#F, low AC, age=1, subpop2
                                    cmlow1[m+12]  * (Y[22+1] + Z[22+1]) / (n_sa[22+1]) + #F, low AC, age=2, subpop2
                                    cmlow1[m+16]  * (Y[24+1] + Z[24+1]) / (n_sa[24+1]) + #F, low AC, age=1, subpop3
                                    cmlow1[m+20]  * (Y[26+1] + Z[26+1]) / (n_sa[26+1]) + #F, low AC, age=2, subpop3
                                    cmhigh1[m]    * (Y[17+1] + Z[17+1]) / (n_sa[17+1]) +  #F, high AC, age=1, subpop1
                                    cmhigh1[m+4]  * (Y[19+1] + Z[19+1]) / (n_sa[19+1]) +  #F, high AC, age=2, subpop1
                                    cmhigh1[m+8]  * (Y[21+1] + Z[21+1]) / (n_sa[21+1]) + #F, high AC, age=1, subpop2
                                    cmhigh1[m+12] * (Y[23+1] + Z[23+1]) / (n_sa[23+1]) + #F, high AC, age=2, subpop2
                                    cmhigh1[m+16] * (Y[25+1] + Z[25+1]) / (n_sa[25+1]) + #F, high AC, age=1, subpop3
                                    cmhigh1[m+20] * (Y[27+1] + Z[27+1]) / (n_sa[27+1]))); #F, high AC, age=2, subpop3
          }

          for(n in 5:8) {   #males subpop 2
            m = n-4;
            lambda[n] =   (b[1] * (cmlow2[m]     * (Y[16+1] + Z[16+1]) / (n_sa[16+1]) + #F, low AC, age=1, subpop1
                                     cmlow2[m+4]   * (Y[18+1] + Z[18+1]) / (n_sa[18+1]) + #F, low AC, age=2, subpop1
                                     cmlow2[m+8]   * (Y[20+1] + Z[20+1]) / (n_sa[20+1]) + #F, low AC, age=1, subpop2
                                     cmlow2[m+12]  * (Y[22+1] + Z[22+1]) / (n_sa[22+1]) + #F, low AC, age=2, subpop2
                                     cmlow2[m+16]  * (Y[24+1] + Z[24+1]) / (n_sa[24+1]) + #F, low AC, age=1, subpop3
                                     cmlow2[m+20]  * (Y[26+1] + Z[26+1]) / (n_sa[26+1]) + #F, low AC, age=2, subpop3
                                     cmhigh2[m]    * (Y[17+1] + Z[17+1]) / (n_sa[17+1]) + #F, high AC, age=1, subpop1
                                     cmhigh2[m+4]  * (Y[19+1] + Z[19+1]) / (n_sa[19+1]) + #F, high AC, age=2, subpop1
                                     cmhigh2[m+8]  * (Y[21+1] + Z[21+1]) / (n_sa[21+1]) + #F, high AC, age=1, subpop2
                                     cmhigh2[m+12] * (Y[23+1] + Z[23+1]) / (n_sa[23+1]) + #F, high AC, age=2, subpop2
                                     cmhigh2[m+16] * (Y[25+1] + Z[25+1]) / (n_sa[25+1]) + #F, high AC, age=1, subpop3
                                     cmhigh2[m+20] * (Y[27+1] + Z[27+1]) / (n_sa[27+1]))); #F, high AC, age=2, subpop3
          }

          for(n in 9:12) {   # males subpop 3
            m = n-8;
            lambda[n] =   (b[1] * ( cmlow3[m]     * (Y[16+1] + Z[16+1]) / (n_sa[16+1]) + #F, low AC, age=1, subpop1
                                      cmlow3[m+4]   * (Y[18+1] + Z[18+1]) / (n_sa[18+1]) + #F, low AC, age=2, subpop1
                                      cmlow3[m+8]   * (Y[20+1] + Z[20+1]) / (n_sa[20+1]) + #F, low AC, age=1, subpop2
                                      cmlow3[m+12]  * (Y[22+1] + Z[22+1]) / (n_sa[22+1]) + #F, low AC, age=2, subpop2
                                      cmlow3[m+16]  * (Y[24+1] + Z[24+1]) / (n_sa[24+1]) + #F, low AC, age=1, subpop3
                                      cmlow3[m+20]  * (Y[26+1] + Z[26+1]) / (n_sa[26+1]) + #F, low AC, age=2, subpop3
                                      cmhigh3[m]    * (Y[17+1] + Z[17+1]) / (n_sa[17+1]) + #F, high AC, age=1, subpop1
                                      cmhigh3[m+4]  * (Y[19+1] + Z[19+1]) / (n_sa[19+1]) + #F, high AC, age=2, subpop1
                                      cmhigh3[m+8]  * (Y[21+1] + Z[21+1]) / (n_sa[21+1]) + #F, high AC, age=1, subpop2
                                      cmhigh3[m+12] * (Y[23+1] + Z[23+1]) / (n_sa[23+1]) + #F, high AC, age=2, subpop2
                                      cmhigh3[m+16] * (Y[25+1] + Z[25+1]) / (n_sa[25+1]) + #F, high AC, age=1, subpop3
                                      cmhigh3[m+20] * (Y[27+1] + Z[27+1]) / (n_sa[27+1]))); #F, high AC, age=2, subpop3
          }

          for(n in 13:16) {   # males subpop 4
            m = n-12;
            lambda[n] =  ((b[1] *  (cmlow4[m]     * (Y[16+1] + Z[16+1]) / (n_sa[16+1]) + #F, low AC, age=1, subpop1
                                      cmlow4[m+4]   * (Y[18+1] + Z[18+1]) / (n_sa[18+1]) + #F, low AC, age=2, subpop1
                                      cmlow4[m+8]   * (Y[20+1] + Z[20+1]) / (n_sa[20+1]) + #F, low AC, age=1, subpop2
                                      cmlow4[m+12]  * (Y[22+1] + Z[22+1]) / (n_sa[22+1]) + #F, low AC, age=2, subpop2
                                      cmlow4[m+16]  * (Y[24+1] + Z[24+1]) / (n_sa[24+1]) + #F, low AC, age=1, subpop3
                                      cmlow4[m+20]  * (Y[26+1] + Z[26+1]) / (n_sa[26+1]) + #F, low AC, age=2, subpop3
                                      cmhigh4[m]    * (Y[17+1] + Z[17+1]) / (n_sa[17+1]) + #F, high AC, age=1, subpop1
                                      cmhigh4[m+4]  * (Y[19+1] + Z[19+1]) / (n_sa[19+1]) +  #F, high AC, age=2, subpop1
                                      cmhigh4[m+8]  * (Y[21+1] + Z[21+1]) / (n_sa[21+1]) + #F, high AC, age=1, subpop2
                                      cmhigh4[m+12] * (Y[23+1] + Z[23+1]) / (n_sa[23+1]) + #F, high AC, age=2, subpop2
                                      cmhigh4[m+16] * (Y[25+1] + Z[25+1]) / (n_sa[25+1]) + #F, high AC, age=1, subpop3
                                      cmhigh4[m+20] * (Y[27+1] + Z[27+1]) / (n_sa[27+1]))) +
                            ((b[3]*c_rr)*
                               (cmlow4[m+24]  * (Y[12+1] + Z[12+1]) / (n_sa[12+1]) + #M, low AC, age=1, subpop4
                                  cmlow4[m+28]  * (Y[14+1] + Z[14+1]) / (n_sa[14+1]) + #M, low AC, age=2, subpop4
                                  cmhigh4[m+24] * (Y[13+1] + Z[13+1]) / (n_sa[13+1]) +  #M, high AC, age=1, subpop4
                                  cmhigh4[m+28] * (Y[15+1] + Z[15+1]) / (n_sa[15+1])))); #M, high AC, age=2, subpop4
          }


          for(n in 17:20) {   #females subpop 1
            m = n;
            lambda[n] =  (b[2] * (cmlow1[m+16]     * (Y[0+1]  + Z[0+1])  /(n_sa[0+1])  + #M, low AC, age=1, subpop1
                                    cmlow1[m+4+16]   * (Y[2+1]  + Z[2+1])  /(n_sa[2+1])  + #M, low AC, age=2, subpop1
                                    cmlow1[m+8+16]   * (Y[4+1]  + Z[4+1])  /(n_sa[4+1])  + #M, low AC, age=1, subpop2
                                    cmlow1[m+12+16]  * (Y[6+1]  + Z[6+1])  /(n_sa[6+1])  + #M, low AC, age=2, subpop2
                                    cmlow1[m+16+16]  * (Y[8+1]  + Z[8+1])  /(n_sa[8+1])  + #M, low AC, age=1, subpop3
                                    cmlow1[m+20+16]  * (Y[10+1] + Z[10+1]) /(n_sa[10+1]) + #M, low AC, age=2, subpop3
                                    cmlow1[m+24+16]  * (Y[12+1] + Z[12+1]) /(n_sa[12+1]) + #M, low AC, age=1, subpop4
                                    cmlow1[m+28+16]  * (Y[14+1] + Z[14+1]) /(n_sa[14+1]) + #M, low AC, age=2, subpop4


                                    cmhigh1[m+16]    * (Y[1+1]  + Z[1+1])  / (n_sa[1+1])  + #M, high AC, age=1, subpop1
                                    cmhigh1[m+4+16]  * (Y[3+1]  + Z[3+1])  / (n_sa[3+1])  + #M, high AC, age=2, subpop1
                                    cmhigh1[m+8+16]  * (Y[5+1]  + Z[5+1])  / (n_sa[5+1])  + #M, high AC, age=1, subpop2
                                    cmhigh1[m+12+16] * (Y[7+1]  + Z[7+1])  / (n_sa[7+1])  + #M, high AC, age=2, subpop2
                                    cmhigh1[m+16+16] * (Y[9+1]  + Z[9+1])  / (n_sa[9+1])  + #M, high AC, age=1, subpop3
                                    cmhigh1[m+20+16] * (Y[11+1] + Z[11+1]) / (n_sa[11+1]) + #M, high AC, age=2, subpop3
                                    cmhigh1[m+24+16] * (Y[13+1] + Z[13+1]) / (n_sa[13+1]) + #M, high AC, age=1, subpop4
                                    cmhigh1[m+28+16] * (Y[15+1] + Z[15+1]) / (n_sa[15+1]))); #M, high AC, age=2, subpop4

          }

          for(n in 21:24) {   #females subpop 2
            m = n-4;
            lambda[n] =  (b[2] * (cmlow2[m+16]     * (Y[0+1]  + Z[0+1])  /(n_sa[0+1])  + #M, low AC, age=1, subpop1
                                    cmlow2[m+4+16]   * (Y[2+1]  + Z[2+1])  /(n_sa[2+1])  + #M, low AC, age=2, subpop1
                                    cmlow2[m+8+16]   * (Y[4+1]  + Z[4+1])  /(n_sa[4+1])  + #M, low AC, age=1, subpop2
                                    cmlow2[m+12+16]  * (Y[6+1]  + Z[6+1])  /(n_sa[6+1])  + #M, low AC, age=2, subpop2
                                    cmlow2[m+16+16]  * (Y[8+1]  + Z[8+1])  /(n_sa[8+1])  + #M, low AC, age=1, subpop3
                                    cmlow2[m+20+16]  * (Y[10+1] + Z[10+1]) /(n_sa[10+1]) + #M, low AC, age=2, subpop3
                                    cmlow2[m+24+16]  * (Y[12+1] + Z[12+1]) /(n_sa[12+1]) + #M, low AC, age=1, subpop4
                                    cmlow2[m+28+16]  * (Y[14+1] + Z[14+1]) /(n_sa[14+1]) + #M, low AC, age=2, subpop4

                                    cmhigh2[m+16]    * (Y[1+1]  + Z[1+1])  / (n_sa[1+1])  + #M, high AC, age=1, subpop1
                                    cmhigh2[m+4+16]  * (Y[3+1]  + Z[3+1])  / (n_sa[3+1])  + #M, high AC, age=2, subpop1
                                    cmhigh2[m+8+16]  * (Y[5+1]  + Z[5+1])  / (n_sa[5+1])  + #M, high AC, age=1, subpop2
                                    cmhigh2[m+12+16] * (Y[7+1]  + Z[7+1])  / (n_sa[7+1])  + #M, high AC, age=2, subpop2
                                    cmhigh2[m+16+16] * (Y[9+1]  + Z[9+1])  / (n_sa[9+1])  + #M, high AC, age=1, subpop3
                                    cmhigh2[m+20+16] * (Y[11+1] + Z[11+1]) / (n_sa[11+1]) + #M, high AC, age=2, subpop3
                                    cmhigh2[m+24+16] * (Y[13+1] + Z[13+1]) / (n_sa[13+1]) + #M, high AC, age=1, subpop4
                                    cmhigh2[m+28+16] * (Y[15+1] + Z[15+1]) / (n_sa[15+1]))); #M, high AC, age=2, subpop4
          }

          for(n in 25:28) {   #females subpop 3
            m = n-8;
            lambda[n] =  (b[1] *  ( cmlow3[m+16]     * (Y[0+1]  + Z[0+1])  /(n_sa[0+1])  + #M, low AC, age=1, subpop1
                                      cmlow3[m+4+16]   * (Y[2+1]  + Z[2+1])  /(n_sa[2+1])  + #M, low AC, age=2, subpop1
                                      cmlow3[m+8+16]   * (Y[4+1]  + Z[4+1])  /(n_sa[4+1])  + #M, low AC, age=1, subpop2
                                      cmlow3[m+12+16]  * (Y[6+1]  + Z[6+1])  /(n_sa[6+1])  + #M, low AC, age=2, subpop2
                                      cmlow3[m+16+16]  * (Y[8+1]  + Z[8+1])  /(n_sa[8+1])  + #M, low AC, age=1, subpop3
                                      cmlow3[m+20+16]  * (Y[10+1] + Z[10+1]) /(n_sa[10+1]) + #M, low AC, age=2, subpop3
                                      cmlow3[m+24+16]  * (Y[12+1] + Z[12+1]) /(n_sa[12+1]) + #M, low AC, age=1, subpop4
                                      cmlow3[m+28+16]  * (Y[14+1] + Z[14+1]) /(n_sa[14+1]) + #M, low AC, age=2, subpop4

                                      cmhigh3[m+16]    * (Y[1+1]  + Z[1+1])  / (n_sa[1+1])  + #M, high AC, age=1, subpop1
                                      cmhigh3[m+4+16]  * (Y[3+1]  + Z[3+1])  / (n_sa[3+1])  + #M, high AC, age=2, subpop1
                                      cmhigh3[m+8+16]  * (Y[5+1]  + Z[5+1])  / (n_sa[5+1])  + #M, high AC, age=1, subpop2
                                      cmhigh3[m+12+16] * (Y[7+1]  + Z[7+1])  / (n_sa[7+1])  + #M, high AC, age=2, subpop2
                                      cmhigh3[m+16+16] * (Y[9+1]  + Z[9+1])  / (n_sa[9+1])  + #M, high AC, age=1, subpop3
                                      cmhigh3[m+20+16] * (Y[11+1] + Z[11+1]) / (n_sa[11+1]) + #M, high AC, age=2, subpop3
                                      cmhigh3[m+24+16] * (Y[13+1] + Z[13+1]) / (n_sa[13+1]) + #M, high AC, age=1, subpop4
                                      cmhigh3[m+28+16] * (Y[15+1] + Z[15+1]) / (n_sa[15+1]))); #M, high AC, age=2, subpop4

          }

          # __Calculate Aging and Birth Matrices ----
          agingS <- as.numeric(aging_rate %*% S)
          agingY <- as.numeric(aging_rate %*% Y)
          agingZ <- as.numeric(aging_rate %*% Z)
          agingNSAmat <- aging_nsa %*% aging_rate
          agingNSA <- as.numeric(agingNSAmat %*% NSA)
          agingSNSAmat <- aging_rate - aging_nsa %*% aging_rate
          agingSNSA <- as.numeric(agingSNSAmat %*% NSA)

          birthsSmat <- birth_rate %*% births_sa
          Ntot <- S + Y + Z + NSA
          birthsS <- as.numeric(birthsSmat %*% Ntot)
          birthsNSAmat <- birth_rate %*% births_nsa
          birthsNSA <- as.numeric(birthsNSAmat %*% Ntot)

          # __Calculate Next Time Step with Differential Equations
          for (n in 1:dim) {
            Pop[n] = Pop[n] + (-lambda[n]*S[n] + gamma[n]*Y[n] + delta[n]*Z[n] + psi[n]*Z[n]  + agingS[n] + agingSNSA[n] + birthsS[n]) * dt   # dS/dt
            Pop[n+dim] = Pop[n+dim] + (pS[n]*lambda[n]*S[n] - gamma[n]*Y[n] + agingY[n]) * dt # dY/dt
            Pop[n+dim*2] = Pop[n+dim*2] + ( (1-pS[n])*lambda[n]*S[n] - delta[n]*Z[n] - psi[n]*Z[n] + agingZ[n] ) * dt #dZ/dt
            Pop[n+dim*3] = Pop[n+dim*3] + (agingNSA[n] + birthsNSA[n]) * dt #dNSA/dt
            Pop[n+dim*4] = Pop[n+dim*4] + (rep_s[n]*gamma[n]*Y[n] + rep_a*psi[n]*Z[n]) * dt #dD/dt
            Pop[n+dim*5] = Pop[n+dim*5] + (rep_s[n]*gamma[n]*Y[n] ) * dt #dR/dt
            Pop[n+dim*6] = Pop[n+dim*6] + (lambda[n]*S[n])*dt  #dCUMINC/dt
            Pop[n+dim*7] = Pop[n+dim*7] + ( psi[n]*( S[n] + Z[n]) ) * dt #dScr/dt
          }
        }
        if(debug) {
          iter <- 52*7*(y - 1) + 7*(w - 1) + d
          outputs[iter,1] = iter
          for (n in 1:(dim*8)) {
            outputs[iter,1+n] = Pop[n]
          }
        } else if (w==32 && d==1) {
          outputs[y,1] = y
          for (n in 1:(dim*8)) {
            outputs[y,1+n] = Pop[n]
          }
        }
      } # end of day loop
    } # end of week loop
  } # end of year loop
  return(outputs)
}
