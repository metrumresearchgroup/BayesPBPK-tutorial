functions{
  real[] PBPKModelODE(real t, real[] x, real[] parms, real[] rdummy, int[] idummy){
    real dxdt[16];
    real WT = rdummy[1];
    
    //// fixed parameters ////
    // Regional blood flows
    real CO  = (187.00*WT^0.81)*60/1000;         // Cardiac output (L/h) from White et al (1968)
    real QHT = 4.0 *CO/100;
    real QBR = 12.0*CO/100;
    real QMU = 17.0*CO/100;
    real QAD = 5.0 *CO/100;
    real QSK = 5.0 *CO/100;
    real QSP = 3.0 *CO/100;
    real QPA = 1.0 *CO/100;
    real QLI = 25.5*CO/100;
    real QST = 1.0 *CO/100;
    real QGU = 14.0*CO/100;
    real QHA = QLI - (QSP + QPA + QST + QGU); // Hepatic artery blood flow
    real QBO = 5.0 *CO/100;
    real QKI = 19.0*CO/100;
    real QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
    real QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

    // Organs' volumes = organs' weights / organs' density
    real VLU = (0.76 *WT/100)/1.051;
    real VHT = (0.47 *WT/100)/1.030;
    real VBR = (2.00 *WT/100)/1.036;
    real VMU = (40.00*WT/100)/1.041;
    real VAD = (21.42*WT/100)/0.916;
    real VSK = (3.71 *WT/100)/1.116;
    real VSP = (0.26 *WT/100)/1.054;
    real VPA = (0.14 *WT/100)/1.045;
    real VLI = (2.57 *WT/100)/1.040;
    real VST = (0.21 *WT/100)/1.050;
    real VGU = (1.44 *WT/100)/1.043;
    real VBO = (14.29*WT/100)/1.990;
    real VKI = (0.44 *WT/100)/1.050;
    real VAB = (2.81 *WT/100)/1.040;
    real VVB = (5.62 *WT/100)/1.040;
    real VRB = (3.86 *WT/100)/1.040;

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
    
    //// parameters to be estimated ////
    real CLint = parms[1];
		real KbBR = parms[2];
		real KbMU = parms[3];
		real KbAD = parms[4];
		real KbBO = parms[5];
		real KbRB = parms[6];
    
    // model compartments
    real Lungs = x[1];
    real Heart = x[2];
    real Brain = x[3];
    real Muscles = x[4];
    real Adipose = x[5];
    real Skin = x[6];
    real Spleen = x[7];
    real Pancreas = x[8];
    real Liver = x[9];
    real Stomach = x[10];
    real Gut = x[11];
    real Bones = x[12];
    real Kidneys = x[13];
    real Arterial_Blood = x[14];
    real Venous_Blood = x[15];
    real Rest_of_Body = x[16];

    dxdt[1] = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);  //lungs
    dxdt[2] = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);  //heart
    dxdt[3] = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);  //brain
    dxdt[4] = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);  //muscles
    dxdt[5] = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);  //adipose
    dxdt[6] = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);  //skin
    dxdt[7] = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);  //spleen
    dxdt[8] = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);  //pancreas
    dxdt[9] = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP + QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU - CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;  //liver
    dxdt[10] = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);  //stomach
    dxdt[11] = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);  //gut
    dxdt[12] = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);  //bones
    dxdt[13] = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);  //kidneys
    dxdt[14] = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);  //arterial blood
    dxdt[15] = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR + QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;  //venous_blood
    dxdt[16] = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);  //rest of body

    return dxdt;
  }
}

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
  real biovar[nSubject, nCmt];
  real tlag[nSubject, nCmt];
  real VVB[nSubject];
  int len[nSubject];
  real<lower = 0> WT[nSubject, 1]; 
  int<lower = 0> idummy[nSubject, 1];
  real BP = 0.61;      // Blood:plasma partition coefficient
  
  for (i in 1:nSubject) {
    for (j in 1:nCmt) {
      biovar[i, j] = 1;
      tlag[i, j] = 0;
    }
    WT[i,1] = weight[i]; 
    idummy[i, 1] = 0;
    len[i] = end[i] - start[i] + 1;
    VVB[i] = (5.62 * weight[i]/100) / 1.040;
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
  real<lower = 0> parms[nSubject, nTheta]; // The [1] indicates the parameters are constant
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubject, nIIV] thetaM; 

  // Matt's trick to use unit scale
  thetaHat[1] = CLintHat; 
  thetaM = (rep_matrix(thetaHat, nSubject) .* exp(diag_pre_multiply(omega, L * etaStd)))';

  for(i in 1:nSubject) {
    parms[i, 1] = thetaM[i, 1]; // CLint
    parms[i, 2] = KbBR; 
    parms[i, 3] = KbMU; 
    parms[i, 4] = KbAD; 
    parms[i, 5] = KbBO; 
    parms[i, 6] = KbRB;
  }

  x = pmx_solve_group_rk45(PBPKModelODE, nCmt, len,
                           time, amt, rate, ii, evid, cmt, addl, ss,
                           parms, biovar, tlag, WT, idummy,
                           1e-6, 1e-6, 1e6);
                
  for(i in 1:nSubject) {             
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
  matrix[nCmt, nt] xPred;
  real<lower = 0> parmsPred[nSubject, nTheta];
  row_vector<lower = 0>[nt] cHatPred;
  row_vector<lower = 0>[nObs] cHatObsPred;
  vector<lower = 0>[nObs] cObsCond;
  row_vector<lower = 0>[nObs] cObsPred;

  // Variables for IIV  
  matrix[nIIV, nSubject] etaStdPred;
  matrix<lower = 0>[nSubject, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  rho = L * L';
  for(i in 1:nSubject) {
    for(j in 1:nIIV) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
  thetaPredM = (rep_matrix(thetaHat, nSubject) .* exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(i in 1:nSubject) {
    parmsPred[i, 1] = thetaPredM[i, 1]; // CL
    parmsPred[i, 2] = KbBR; 
    parmsPred[i, 3] = KbMU;
    parmsPred[i, 4] = KbAD; 
    parmsPred[i, 5] = KbBO; 
    parmsPred[i, 6] = KbRB;
  }

  xPred = pmx_solve_group_rk45(PBPKModelODE, nCmt, len,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parmsPred, biovar, tlag, WT, idummy,
                               1e-6, 1e-6, 1e6);

  for(i in 1:nSubject) {
    cHatPred[start[i]:end[i]] = xPred[15, start[i]:end[i]] / (VVB[i]*BP/1000);
  }

  // predictions for observed data records
  cHatObsPred = cHatPred[iObs];

  for(i in 1:nObs) {
    cObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));  // individual predictions
    cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));  // population predictions
  }
}
