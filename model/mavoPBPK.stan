// Mavoglurant population PBPK model
// Source: https://nlmixrdevelopment.github.io/nlmixr.examples/articles/mavoglurant.html

functions{
  real[] PBPKModelODE(real t,
			real[] x,
			real[] parms,
			real WT,
			real[] rdummy,
			int[] idummy){
    
    real dxdt[16];
    
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
    real KbSK = exp(-.5238);
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
    dxdt[14] = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);  //lungs
    dxdt[15] = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR +  
            QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD +
            QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO +
            QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;  //venous_blood
    dxdt[16] = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);  //arterial blood

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
  int<lower = 1> nSubjects;
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> weight[nSubjects];
  
  // data for priors
  real<lower = 0> CLintHatPrior;
  real<lower = 0> CLintHatPriorCV;
  real<lower = 0> KbBRPrior;
  real<lower = 0> KbMUPrior;
  real<lower = 0> KbADPrior;
  real<lower = 0> KbBOPrior;
  real<lower = 0> KbRBPrior;
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 6;  // number of parameters
  int nIIV = 1;  // parameters with IIV
  int nCmt = 16;  // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];
  
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
  matrix[nIIV, nSubjects] etaStd;
  
}

transformed parameters{
  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;
  matrix[nCmt, nt] x;
  real<lower = 0> parms[nTheta]; 
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; 

  // Matt's trick to use unit scale
  thetaHat[1] = CLintHat; 
  thetaM = (rep_matrix(thetaHat, nSubjects) .* 
             exp(diag_pre_multiply(omega, L * etaStd)))';

  for(i in 1:nSubjects) {

  real<lower = 0> ;
  real<lower = 0> ;
  real<lower = 0> ;
  real<lower = 0> ;
  real<lower = 0> ;
  
    parms[1] = thetaM[i, 1]; // CLint
    parms[2] = KbBR; 
    parms[3] = KbMU; 
    parms[4] = KbAD; 
    parms[5] = KbBO; 
    parms[6] = KbRB;






    x[start[i]:end[i]] = pmx_solve_rk45(PBPKModelODE, nCmt,
                                        time[start[i]:end[i]], 
                                        amt[start[i]:end[i]], 
                                        rate[start[i]:end[i]], 
                                        ii[start[i]:end[i]], 
                                        evid[start[i]:end[i]], 
                                        cmt[start[i]:end[i]], 
                                        addl[start[i]:end[i]], 
                                        ss[start[i]:end[i]],
                                        parms, biovar, tlag,
                                        1e-6, 1e-6, 1e6);
                             
    cHat[start[i]:end[i]] = x[2, start[i]:end[i]] / parms[3];  // divide by V1
    neutHat[start[i]:end[i]] = x[8, start[i]:end[i]] + parms[7];  // Add baseline
    
  }
  
  cHatObs = cHat[iObsPK];
  neutHatObs = neutHat[iObsPD];

}

////////////////////////





























data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real rate[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;
  
  // data for population model
  int<lower = 1> nSubjects;
  int<lower = 1> start[nSubjects];
  int<lower = 1> end[nSubjects];
  real<lower = 0> weight[nSubjects];
  
  // data for priors
  real<lower = 0> CLHatPrior;
  real<lower = 0> QHatPrior;
  real<lower = 0> V1HatPrior;
  real<lower = 0> V2HatPrior;
  real<lower = 0> kaHatPrior;
  real<lower = 0> CLHatPriorCV;
  real<lower = 0> QHatPriorCV;
  real<lower = 0> V1HatPriorCV;
  real<lower = 0> V2HatPriorCV;
  real<lower = 0> kaHatPriorCV;
  real<lower = 0> circ0HatPrior;
  real<lower = 0> circ0HatPriorCV;
  real<lower = 0> mttHatPrior;
  real<lower = 0> mttHatPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaHatPrior;
  real<lower = 0> alphaHatPriorCV;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nTheta = 9;  // number of ODE parameters
  int nIIV = 7;  // parameters with IIV
  int nCmt = 8;  // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];
  
  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  real<lower = 0> mttHat;
  real<lower = 0> circ0Hat;
  real<lower = 0> alphaHat;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  
  // IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
  matrix[nIIV, nSubjects] etaStd;
  
}

transformed parameters{
  row_vector[nt] cHat;
  row_vector[nObsPK] cHatObs;
  row_vector[nt] neutHat;
  row_vector[nObsPD] neutHatObs;
  matrix[nCmt, nt] x;
  real<lower = 0> parms[nTheta]; // The [1] indicates the parameters are constant
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] thetaHat;
  matrix<lower = 0>[nSubjects, nIIV] thetaM; 

  // Matt's trick to use unit scale
  thetaHat[1] = CLHat; 
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = mttHat;
  thetaHat[6] = circ0Hat;
  thetaHat[7] = alphaHat;
  thetaM = (rep_matrix(thetaHat, nSubjects) .* 
             exp(diag_pre_multiply(omega, L * etaStd)))';

  for(i in 1:nSubjects) {

    parms[1] = thetaM[i, 1] * (weight[i] / 70)^0.75; // CL
    parms[2] = thetaM[i, 2] * (weight[i] / 70)^0.75; // Q
    parms[3] = thetaM[i, 3] * (weight[i] / 70); // V1
    parms[4] = thetaM[i, 4] * (weight[i] / 70); // V2
    parms[5] = kaHat; // ka
    parms[6] = thetaM[i, 5]; // mtt
    parms[7] = thetaM[i, 6]; // circ0
    parms[8] = gamma;
    parms[9] = thetaM[i, 7]; // alpha

    x[start[i]:end[i]] = pmx_solve_rk45(twoCptNeutModelODE, nCmt,
                                        time[start[i]:end[i]], 
                                        amt[start[i]:end[i]], 
                                        rate[start[i]:end[i]], 
                                        ii[start[i]:end[i]], 
                                        evid[start[i]:end[i]], 
                                        cmt[start[i]:end[i]], 
                                        addl[start[i]:end[i]], 
                                        ss[start[i]:end[i]],
                                        parms, biovar, tlag,
                                        1e-6, 1e-6, 1e6);
                             
    cHat[start[i]:end[i]] = x[2, start[i]:end[i]] / parms[3];  // divide by V1
    neutHat[start[i]:end[i]] = x[8, start[i]:end[i]] + parms[7];  // Add baseline
    
  }
  
  cHatObs = cHat[iObsPK];
  neutHatObs = neutHat[iObsPD];

}

model{
  // Priors
  CLHat ~ lognormal(log(CLHatPrior), CLHatPriorCV);
  QHat ~ lognormal(log(QHatPrior), QHatPriorCV);
  V1Hat ~ lognormal(log(V1HatPrior), V1HatPriorCV);
  V2Hat ~ lognormal(log(V2HatPrior), V2HatPriorCV);
  kaHat ~ lognormal(log(kaHatPrior), kaHatPriorCV);
  sigma ~ cauchy(0, 1);
  
  mttHat ~ lognormal(log(mttHatPrior), mttHatPriorCV);
  circ0Hat ~ lognormal(log(circ0HatPrior), circ0HatPriorCV);
  alphaHat ~ lognormal(log(alphaHatPrior), alphaHatPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  // Parameters for Matt's trick
  L ~ lkj_corr_cholesky(1);
  to_vector(etaStd) ~ normal(0, 1);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities{
  matrix[nCmt, nt] xPred;
  real<lower = 0> parmsPred[nTheta];
  row_vector[nt] cHatPred;
  row_vector[nt] neutHatPred;
  vector<lower = 0>[nObsPK] cHatObsCond;
  row_vector<lower = 0>[nObsPK] cHatObsPred;
  vector<lower = 0>[nObsPD] neutHatObsCond;
  row_vector<lower = 0>[nObsPD] neutHatObsPred;

  // Variables for IIV  
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower = 0>[nSubjects, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  rho = L * L';
  for(i in 1:nSubjects) {
    for(j in 1:nIIV) {
      etaStdPred[j, i] = normal_rng(0, 1);
    }
  }
  thetaPredM = (rep_matrix(thetaHat, nSubjects) .* 
                exp(diag_pre_multiply(omega, L * etaStdPred)))';

  for(i in 1:nSubjects) {
    parmsPred[1] = thetaPredM[i, 1] * (weight[i] / 70)^0.75; // CL
    parmsPred[2] = thetaPredM[i, 2] * (weight[i] / 70)^0.75; // Q
    parmsPred[3] = thetaPredM[i, 3] * (weight[i] / 70); // V1
    parmsPred[4] = thetaPredM[i, 4] * (weight[i] / 70); // V2
    parmsPred[5] = kaHat; // ka
    parmsPred[6] = thetaPredM[i, 5]; // mtt
    parmsPred[7] = thetaPredM[i, 6]; // circ0
    parmsPred[8] = gamma; // gamma
    parmsPred[9] = thetaPredM[i, 7]; // alpha

    xPred[start[i]:end[i]] = pmx_solve_rk45(twoCptNeutModelODE, nCmt,
                                            time[start[i]:end[i]], 
                                            amt[start[i]:end[i]],
                                            rate[start[i]:end[i]],
                                            ii[start[i]:end[i]],
                                            evid[start[i]:end[i]],
                                            cmt[start[i]:end[i]],
                                            addl[start[i]:end[i]],
                                            ss[start[i]:end[i]],
                                            parmsPred, biovar, tlag,
                                            1e-6, 1e-6, 1e6);

    cHatPred[start[i]:end[i]] = xPred[2, start[i]:end[i]] / parmsPred[3]; // divide by V1
    neutHatPred[start[i]:end[i]] = xPred[8, start[i]:end[i]] + parmsPred[7]; // Add baseline
  }

  // predictions for observed data records
  cHatObsPred = cHatPred[iObsPK];
  neutHatObsPred = neutHatPred[iObsPD];

  for(i in 1:nObsPK) {
    cHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));
    cHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));
  }
  
  for(i in 1:nObsPD) {
    neutHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObs[i])), sigmaNeut));
    neutHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObsPred[i])), sigmaNeut));
  }
}


