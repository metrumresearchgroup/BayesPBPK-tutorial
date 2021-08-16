data{
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real rate[nt];
  vector<lower = 0>[nObs] cObs;
  
  // data for population model
  int<lower = 1> nSubject;
  int<lower = 1> start[nSubject];
  int<lower = 1> end[nSubject];
  real<lower = 0> weight[nSubject];
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 6;  // number of parameters
  int nIIV = 1;  // parameters with IIV
  int nCmt = 16;  // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];
  
  // objects to hold fixed parameters
    // Regional blood flows
    real CO[nSubject];         // Cardiac output (L/h) from White et al (1968)
    real QHT[nSubject];
    real QBR[nSubject];
    real QMU[nSubject];
    real QAD[nSubject];
    real QSK[nSubject];
    real QSP[nSubject];
    real QPA[nSubject];
    real QLI[nSubject];
    real QST[nSubject];
    real QGU[nSubject];
    real QHA[nSubject]; // Hepatic artery blood flow
    real QBO[nSubject];
    real QKI[nSubject];
    real QRB[nSubject];
    real QLU[nSubject];
    
    // Organs' volumes = organs' weights / organs' density
    real VLU[nSubject];
    real VHT[nSubject];
    real VBR[nSubject];
    real VMU[nSubject];
    real VAD[nSubject];
    real VSK[nSubject];
    real VSP[nSubject];
    real VPA[nSubject];
    real VLI[nSubject];
    real VST[nSubject];
    real VGU[nSubject];
    real VBO[nSubject];
    real VKI[nSubject];
    real VAB[nSubject];
    real VVB[nSubject];
    real VRB[nSubject];
    
    // partition coefficients
    real KbLU = exp(0.8334);
    real KbHT = exp(1.1205);
    real KbSK = exp(-0.5238);
    real KbSP = exp(0.3224);
    real KbPA = exp(0.3224);
    real KbLI = exp(1.7604);
    real KbST = exp(0.3224);
    real KbGU = exp(1.2026);
    real KbKI = exp(1.3171);
    
    // Other parameters
    real BP = 0.61;      // Blood:plasma partition coefficient
    real fup = 0.028;    // Fraction unbound in plasma
    real fub = fup/BP;   // Fraction unbound in blood
  
  //// fixed parameters ////
  for(i in 1:nSubject) {
    // Regional blood flows
    CO[i]  = (187.00*weight[i]^0.81)*60/1000;         // Cardiac output (L/h) from White et al (1968)
    QHT[i] = 4.0*CO[i]/100;
    QBR[i] = 12.0*CO[i]/100;
    QMU[i] = 17.0*CO[i]/100;
    QAD[i] = 5.0 *CO[i]/100;
    QSK[i] = 5.0 *CO[i]/100;
    QSP[i] = 3.0 *CO[i]/100;
    QPA[i] = 1.0 *CO[i]/100;
    QLI[i] = 25.5*CO[i]/100;
    QST[i] = 1.0 *CO[i]/100;
    QGU[i] = 14.0*CO[i]/100;
    QHA[i] = QLI[i] - (QSP[i] + QPA[i] + QST[i] + QGU[i]); // Hepatic artery blood flow
    QBO[i] = 5.0 *CO[i]/100;
    QKI[i] = 19.0*CO[i]/100;
    QRB[i] = CO[i] - (QHT[i] + QBR[i] + QMU[i] + QAD[i] + QSK[i] + QLI[i] + QBO[i] + QKI[i]);
    QLU[i] = QHT[i] + QBR[i] + QMU[i] + QAD[i] + QSK[i] + QLI[i] + QBO[i] + QKI[i] + QRB[i];

    // Organs' volumes = organs' weights / organs' density
    VLU[i] = (0.76 *weight[i]/100)/1.051;
    VHT[i] = (0.47 *weight[i]/100)/1.030;
    VBR[i] = (2.00 *weight[i]/100)/1.036;
    VMU[i] = (40.00*weight[i]/100)/1.041;
    VAD[i] = (21.42*weight[i]/100)/0.916;
    VSK[i] = (3.71 *weight[i]/100)/1.116;
    VSP[i] = (0.26 *weight[i]/100)/1.054;
    VPA[i] = (0.14 *weight[i]/100)/1.045;
    VLI[i] = (2.57 *weight[i]/100)/1.040;
    VST[i] = (0.21 *weight[i]/100)/1.050;
    VGU[i] = (1.44 *weight[i]/100)/1.043;
    VBO[i] = (14.29*weight[i]/100)/1.990;
    VKI[i] = (0.44 *weight[i]/100)/1.050;
    VAB[i] = (2.81 *weight[i]/100)/1.040;
    VVB[i] = (5.62 *weight[i]/100)/1.040;
    VRB[i] = (3.86 *weight[i]/100)/1.040;
  }
  
  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0> CLintHat;
  real<lower = 0> KbBR;
  real<lower = 0> KbMU;
  real<lower = 0> KbAD;
  real<lower = 0> KbBO;
  real<lower = 0> KbRB;
  
  // residual error
  real<lower = 0> sigma;
  
  // IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
  matrix[nIIV, nSubject] etaStd;
}

transformed parameters{
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nCmt, nt] x;
  matrix[nCmt, nCmt] K;
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubject, nIIV] thetaM; 
  
  // individual-level physiological paraeters
  real CLint[nSubject];

  // Matt's trick to use unit scale
  thetaHat[1] = CLintHat; 
  thetaM = (rep_matrix(thetaHat, nSubject) .* exp(diag_pre_multiply(omega, L * etaStd)))';

  for(i in 1:nSubject) {
    CLint[i] = thetaM[i, 1]; // CLint
    
    // Filling coefficient matrix
    K = rep_matrix(0, nCmt, nCmt);
    
    K[1,1] = -QLU[i]/KbLU/VLU[i];
    K[14,1] =  QLU[i]/KbLU/VLU[i];
    K[2,2] = -QHT[i]/KbHT/VHT[i];
    K[15,2] =  QHT[i]/KbHT/VHT[i];
    K[3,3] = -QBR[i]/KbBR/VBR[i];
    K[15,3] = QBR[i]/KbBR/VBR[i];
    K[4,4] = -QMU[i]/KbMU/VMU[i];
    K[15,4] = QMU[i]/KbMU/VMU[i];
    K[5,5] = -QAD[i]/KbAD/VAD[i];
    K[15,5] = QAD[i]/KbAD/VAD[i];
    K[6,6] = -QSK[i]/KbSK/VSK[i];
    K[15,6] = QSK[i]/KbSK/VSK[i];
    K[7,7] = -QSP[i]/KbSP/VSP[i];
    K[9,7] = QSP[i]/KbSP/VSP[i];
    K[8,8] = -QPA[i]/KbPA/VPA[i];
    K[9,8] = QPA[i]/KbPA/VPA[i];
    K[9,9] = -(CLint[i]*fub + QLI[i])/KbLI/VLI[i];
    K[15,9] = QLI[i]/KbLI/VLI[i];
    K[9,10] = QST[i]/KbST/VST[i];
    K[10,10] = -QST[i]/KbST/VST[i];
    K[9,11] = QGU[i]/KbGU/VGU[i];
    K[11,11] = -QGU[i]/KbGU/VGU[i];
    K[12,12] = -QBO[i]/KbBO/VBO[i];
    K[15,12] = QBO[i]/KbBO/VBO[i];
    K[13,13] = -QKI[i]/KbKI/VKI[i];
    K[15,13] = QKI[i]/KbKI/VKI[i];
    K[2,14] =  QHT[i]/VAB[i];
    K[3,14] =  QBR[i]/VAB[i];
    K[4,14] =  QMU[i]/VAB[i];
    K[5,14] =  QAD[i]/VAB[i];
    K[6,14] =  QSK[i]/VAB[i];
    K[7,14] =  QSP[i]/VAB[i];
    K[8,14] =  QPA[i]/VAB[i];
    K[9,14] =  QHA[i]/VAB[i];
    K[10,14] =  QST[i]/VAB[i];
    K[11,14] =  QGU[i]/VAB[i];
    K[12,14] =  QBO[i]/VAB[i];
    K[13,14] =  QKI[i]/VAB[i];
    K[14,14] = -QLU[i]/VAB[i];
    K[16,14] =  QRB[i]/VAB[i];
    K[1,15] =  QLU[i]/VVB[i];
    K[15,15] = -QLU[i]/VVB[i];
    K[15,16] = QRB[i]/KbRB/VRB[i];
    K[16,16] = -QRB[i]/KbRB/VRB[i];

    x[,start[i]:end[i]] = pmx_solve_linode(time[start[i]:end[i]],
                                          amt[start[i]:end[i]],
                                          rate[start[i]:end[i]],
                                          ii[start[i]:end[i]],
                                          evid[start[i]:end[i]],
                                          cmt[start[i]:end[i]],
                                          addl[start[i]:end[i]],
                                          ss[start[i]:end[i]],
                                          K, biovar, tlag);
           
    cHat[start[i]:end[i]] = x[15, start[i]:end[i]] / (VVB[i]*BP/1000);  // divide by subject's blood volume VVB
  }
  cHatObs = cHat[iObs];
}

model{
  // Priors
  CLintHat ~ lognormal(7.1, 0.25);
  KbBR ~ lognormal(1.1, 0.25);
  KbMU ~ lognormal(0.3, 0.25);
  KbAD ~ lognormal(2, 0.25);
  KbBO ~ lognormal(0.03, 0.25);
  KbRB ~ lognormal(0.3, 0.25);

  sigma ~ cauchy(0, 0.5);

  // Parameters for Matt's trick
  omega ~ cauchy(0, 0.5);
  L ~ lkj_corr_cholesky(1);
  to_vector(etaStd) ~ normal(0, 1);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
  matrix<lower = 0>[nCmt, nt] xPred;
  matrix[nCmt, nCmt] KPred;
  row_vector<lower = 0>[nt] cHatPred;
  row_vector<lower = 0>[nObs] cHatObsPred;
  vector<lower = 0>[nObs] cObsCond;
  row_vector<lower = 0>[nObs] cObsPred;

  // Variables for IIV  
  matrix[nIIV, nSubject] etaStdPred;
  matrix<lower = 0>[nSubject, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  // Individual-level model parameters
  real CLintPred[nSubject];
  
  rho = L * L';
  for(i in 1:nSubject) {
    for(j in 1:nIIV) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
  thetaPredM = (rep_matrix(thetaHat, nSubject) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(i in 1:nSubject) {
    CLintPred[i] = thetaPredM[i, 1]; // CLintPred
    
    // Filling coefficient matrix
    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1,1] = -QLU[i]/KbLU/VLU[i];
    KPred[14,1] =  QLU[i]/KbLU/VLU[i];
    KPred[2,2] = -QHT[i]/KbHT/VHT[i];
    KPred[15,2] =  QHT[i]/KbHT/VHT[i];
    KPred[3,3] = -QBR[i]/KbBR/VBR[i];
    KPred[15,3] = QBR[i]/KbBR/VBR[i];
    KPred[4,4] = -QMU[i]/KbMU/VMU[i];
    KPred[15,4] = QMU[i]/KbMU/VMU[i];
    KPred[5,5] = -QAD[i]/KbAD/VAD[i];
    KPred[15,5] = QAD[i]/KbAD/VAD[i];
    KPred[6,6] = -QSK[i]/KbSK/VSK[i];
    KPred[15,6] = QSK[i]/KbSK/VSK[i];
    KPred[7,7] = -QSP[i]/KbSP/VSP[i];
    KPred[9,7] = QSP[i]/KbSP/VSP[i];
    KPred[8,8] = -QPA[i]/KbPA/VPA[i];
    KPred[9,8] = QPA[i]/KbPA/VPA[i];
    KPred[9,9] = -(CLintPred[i]*fub + QLI[i])/KbLI/VLI[i];
    KPred[15,9] = QLI[i]/KbLI/VLI[i];
    KPred[9,10] = QST[i]/KbST/VST[i];
    KPred[10,10] = -QST[i]/KbST/VST[i];
    KPred[9,11] = QGU[i]/KbGU/VGU[i];
    KPred[11,11] = -QGU[i]/KbGU/VGU[i];
    KPred[12,12] = -QBO[i]/KbBO/VBO[i];
    KPred[15,12] = QBO[i]/KbBO/VBO[i];
    KPred[13,13] = -QKI[i]/KbKI/VKI[i];
    KPred[15,13] = QKI[i]/KbKI/VKI[i];
    KPred[2,14] =  QHT[i]/VAB[i];
    KPred[3,14] =  QBR[i]/VAB[i];
    KPred[4,14] =  QMU[i]/VAB[i];
    KPred[5,14] =  QAD[i]/VAB[i];
    KPred[6,14] =  QSK[i]/VAB[i];
    KPred[7,14] =  QSP[i]/VAB[i];
    KPred[8,14] =  QPA[i]/VAB[i];
    KPred[9,14] =  QHA[i]/VAB[i];
    KPred[10,14] =  QST[i]/VAB[i];
    KPred[11,14] =  QGU[i]/VAB[i];
    KPred[12,14] =  QBO[i]/VAB[i];
    KPred[13,14] =  QKI[i]/VAB[i];
    KPred[14,14] = -QLU[i]/VAB[i];
    KPred[16,14] =  QRB[i]/VAB[i];
    KPred[1,15] =  QLU[i]/VVB[i];
    KPred[15,15] = -QLU[i]/VVB[i];
    KPred[15,16] = QRB[i]/KbRB/VRB[i];
    KPred[16,16] = -QRB[i]/KbRB/VRB[i];

    xPred[,start[i]:end[i]] = pmx_solve_linode(time[start[i]:end[i]],
                                              amt[start[i]:end[i]],
                                              rate[start[i]:end[i]],
                                              ii[start[i]:end[i]],
                                              evid[start[i]:end[i]],
                                              cmt[start[i]:end[i]],
                                              addl[start[i]:end[i]],
                                              ss[start[i]:end[i]],
                                              KPred, biovar, tlag);

    cHatPred[start[i]:end[i]] = xPred[15, start[i]:end[i]] / (VVB[i]*BP/1000);
  }

   // predictions for observed data records
   cHatObsPred = cHatPred[iObs];

   for(i in 1:nObs) {
     cObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));  // individual predictions
     cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));  // population predictions
  }
}
