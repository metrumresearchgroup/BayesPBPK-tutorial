// Bayesian linear popPBPK model for Theophylline

functions{
  vector calcKpPT(real logP,
  real pKa,
  real fu,
  real BP,
  matrix TC){

    real logD;
    real logDStar;
    real DStar;
    real P;
    real fut;
    vector[11] Kp;

    logD = 1.115 * logP - 1.35;
    logDStar = logD - log10(1 + 10^(pKa - TC[1,5]));
    DStar = 10^logDStar;
    P = 10^logP;
    fut = 1 / (1 + ((1 - fu) / fu) * 0.5);

    // adipose
    Kp[1] = ((DStar * (TC[1,3] + 0.3 * TC[1,4]) + (1 * (TC[1,2] + 0.7 * TC[1,4]))) / (DStar * (TC[12,3] + 0.3 * TC[12,4]) + (1 * (TC[12,2] + 0.7 * TC[12,4])))) * fu;

    // rest of tissues
    Kp[2:11] = ((P * (TC[2:11,3] + 0.3 * TC[2:11,4]) + (1 * (TC[2:11,2] + 0.7 * TC[2:11,4]))) ./ (P * (TC[12,3] + 0.3 * TC[12,4]) + (1 * (TC[12,2] + 0.7 * TC[12,4])))) * (fu/fut);

    return Kp;
  }
}

  
data{
  int<lower = 1> nSubject;
  int<lower = 1> nt;
  int<lower = 1> nObs;
  int<lower = 1> iObs[nObs];
  real<lower = 0> amt[nt];
  int<lower = 1> cmt[nt];
  int<lower = 0> evid[nt];
  int<lower = 1> start[nSubject];
  int<lower = 1> end[nSubject];
  real<lower = 0> time[nt];
  row_vector<lower = 0>[nObs] cObs;
  real<lower = 0> weight[nSubject];
  real<lower = 0> height[nSubject];
  real<lower = 0> rate[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  matrix<lower = 0>[12,5] TC;
}


transformed data{
  // physicochemical properties
  real pKaPrior = 8.81; //https://onlinelibrary.wiley.com/doi/pdf/10.1002/jps.10005 and DrugBank
  real fuPrior = 0.5;  //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1381359/pdf/brjclinpharm00045-0078.pdf
  real logPPrior = -0.02; //https://pubchem.ncbi.nlm.nih.gov/compound/Theophylline#section=Octanol-Water-Partition-Coefficient
  real BPPrior = 0.82;  //https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1884855/
  
  // prior physiological parameters info for 30 yo male
  // organ volumes
  real<lower = 0> VadMeanPrior = 18.2;
  real<lower = 0> VarMeanPrior = 0.295 * 5.6;
  real<lower = 0> VboMeanPrior = 10.5;
  real<lower = 0> VbrMeanPrior = 1.45;
  real<lower = 0> VguMeanPrior = 0.65;
  real<lower = 0> VheMeanPrior = 0.33;
  real<lower = 0> VkiMeanPrior = 0.31;
  real<lower = 0> VliMeanPrior = 1.8;
  real<lower = 0> VluMeanPrior = 0.5;
  real<lower = 0> VmuMeanPrior = 29;
  real<lower = 0> VskMeanPrior = 3.3;
  real<lower = 0> VspMeanPrior = 0.15;
  real<lower = 0> VveMeanPrior = 0.705 * 5.6;
  
  // organ volumes CV
  real<lower = 0> VadGSDPrior = 1.65;  //lognormal
  real<lower = 0> VarCVPrior = 5;
  real<lower = 0> VboCVPrior = 1;
  real<lower = 0> VbrCVPrior = 5;
  real<lower = 0> VguCVPrior = 20;
  real<lower = 0> VheCVPrior = 20;
  real<lower = 0> VkiCVPrior = 25;
  real<lower = 0> VliCVPrior = 23;
  real<lower = 0> VluGSDPrior = 1.3;  // lognormal
  real<lower = 0> VmuGSDPrior = 1.1;  // lognormal
  real<lower = 0> VskCVPrior = 45;
  real<lower = 0> VspGSDPrior = 1.5;  // lognormal
  real<lower = 0> VveCVPrior = 5;
  
  // individual volumes
  //real<lower = 0> Vad[nSubject];  //should be declared in transformed parameters block where it is assigned  
  real<lower = 0> Var[nSubject];
  //real<lower = 0> Vbo[nSubject];
  real<lower = 0> Vbr[nSubject];
  real<lower = 0> Vgu[nSubject];
  real<lower = 0> Vhe[nSubject];
  real<lower = 0> Vki[nSubject];
  real<lower = 0> Vli[nSubject];
  real<lower = 0> Vlu[nSubject];
  //real<lower = 0> Vmu[nSubject];
  real<lower = 0> Vsk[nSubject];
  real<lower = 0> Vsp[nSubject];
  real<lower = 0> Vve[nSubject];
  
  // blood flows
  real<lower = 0> COMeanPrior = 6.5*60;  //L/h
  real<lower = 0> COCVPrior = 5;
  
  real<lower = 0> Qad[nSubject];
  real<lower = 0> Qbo[nSubject];
  real<lower = 0> Qbr[nSubject];
  //real<lower = 0> Qgu[nSubject];
  real<lower = 0> Qha[nSubject];
  real<lower = 0> Qhe[nSubject];
  //real<lower = 0> Qki[nSubject];
  //real<lower = 0> Qli[nSubject];  //should be declared in transformed parameters block where it is assigned
  //real<lower = 0> Qlu[nSubject];  //should be declared in transformed parameters block where it is assigned
  real<lower = 0> Qmu[nSubject];
  real<lower = 0> Qsk[nSubject];
  real<lower = 0> Qsp[nSubject];
  
  // other physiological parameters
  real CLhepaticPrior = 0.009 * VliMeanPrior * 60;
  real CLrenalPrior = 0.0317 * VkiMeanPrior * 60;
  
  // calculate boundaries for parameters of interest; Mean +- 2*SD; 5% for flows
  real VmuUpper = exp(log(VmuMeanPrior) + 2*log(VmuGSDPrior));
  real VmuLower = exp(log(VmuMeanPrior) - 2*log(VmuGSDPrior));
  real VboUpper = VboMeanPrior + 2*(VboCVPrior*VboMeanPrior/100);
  real VboLower = VboMeanPrior - 2*(VboCVPrior*VboMeanPrior/100);
  real QguUpper = 0.16*COMeanPrior + 2*(5*0.16*COMeanPrior/100);
  real QguLower = 0.16*COMeanPrior - 2*(5*0.16*COMeanPrior/100);
  real QkiUpper = 0.19*COMeanPrior + 2*(5*0.19*COMeanPrior/100);
  real QkiLower = 0.19*COMeanPrior - 2*(5*0.19*COMeanPrior/100);
  //real logPUpper = logPPrior / 10;
  //real logPLower = logPPrior * 10;
  real fuUpper = 1;
  real fuLower = 0;
  real CLhepaticUpper = CLhepaticPrior * 100;  //in the paper CL is 0.009 /min
  real CLhepaticLower = CLhepaticPrior / 100;  //in the paper CL is 0.009 /min
  
  // other definitions
  //vector[nObs] logCObs = log(cObs);
  int<lower = 1> nRandom = 6;
  int<lower = 1> nParms = 7;
  int nCmt = 14;
  int nOrgans = 11;
  real F[nCmt];
  real tlag[nCmt];
  
  for (i in 1:nCmt) {
    F[i] = 1;
    tlag[i] = 0;
  }
  
// get individual params
for(j in 1:nSubject){
  //Vad[j] = VadMeanPrior * (height[j] / 1.76)^2;
  Var[j] = VarMeanPrior * (height[j] / 1.76)^0.75;
  //Vbo[j] = VboMeanPrior * (height[j] / 1.76)^2;
  Vbr[j] = VbrMeanPrior * (height[j] / 1.76)^0;
  Vgu[j] = VguMeanPrior * (height[j] / 1.76)^0.75;
  Vhe[j] = VheMeanPrior * (height[j] / 1.76)^0.75;
  Vki[j] = VkiMeanPrior * (height[j] / 1.76)^0.75;
  Vli[j] = VliMeanPrior * (height[j] / 1.76)^0.75;
  Vlu[j] = VluMeanPrior * (height[j] / 1.76)^0.75;
  //Vmu[j] = VmuMeanPrior * (height[j] / 1.76)^2;
  Vsk[j] = VskMeanPrior * (height[j] / 1.76)^1.6;
  Vsp[j] = VspMeanPrior * (height[j] / 1.76)^0.75;
  Vve[j] = VveMeanPrior * (height[j] / 1.76)^0.75;
  
  Qad[j] = 0.05 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qbo[j] = 0.05 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qbr[j] = 0.12 * COMeanPrior * (height[j] / 1.76)^0.75;
  //Qgu[j] = 0.16 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qhe[j] = 0.04 * COMeanPrior * (height[j] / 1.76)^0.75;
  //Qki[j] = 0.19 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qmu[j] = 0.05 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qsk[j] = 0.05 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qsp[j] = 0.03 * COMeanPrior * (height[j] / 1.76)^0.75;
  Qha[j] = 0.065 * COMeanPrior * (height[j] / 1.76)^0.75;
  //Qli[j] = Qsp[j] + Qgu[j] + Qha[j];
  //Qlu[j] = Qad[j] + Qbo[j] + Qbr[j] + Qhe[j] + Qki[j] + Qli[j] + Qmu[j] + Qsk[j];
 }
}


parameters{
  // population level parameters
  real<lower = VboLower, upper = VboUpper> VboHat;
  real<lower = VmuLower, upper = VmuUpper> VmuHat;
  real<lower = QguLower, upper = QguUpper> QguHat;
  real<lower = QkiLower, upper = QkiUpper> QkiHat;
  real<lower = CLhepaticLower, upper = CLhepaticUpper> CLhepaticHat;
  real<lower = 0> kaHat;
  
  // drug global parameters
  // real<lower = logPLower, upper = logPUpper> logP;
  real<lower = fuLower, upper = fuUpper> fu;
  
  // inter-Individual variability
  vector<lower = 0>[nRandom] omega;
  cholesky_factor_corr[nRandom] L;
  matrix[nRandom, nSubject] eta;
  
  // residual variability
  real<lower = 0> sigma;
}


transformed parameters{
  // Vector of parameter typical values -- only those with IIV
  vector<lower = 0>[nRandom] thetaHat;
  
  // Matrix of individual-level model parameters
  matrix<lower = 0>[nSubject, nRandom] theta;
  
  // variables to store predicted concentrations and amounts
  matrix[nCmt, nCmt] K;
  row_vector<lower = 0>[nt] cHat;
  matrix[nCmt, nt] x;  // rows are compartments and columns are datapoints
  
  // Individual-level physiological model parameters
  row_vector<lower = 0>[nSubject] Vbo;
  row_vector<lower = 0>[nSubject] Vmu;
  row_vector<lower = 0>[nSubject] Qgu;
  row_vector<lower = 0>[nSubject] Qki;
  row_vector<lower = 0>[nSubject] CLhepatic;
  row_vector<lower = 0>[nSubject] ka;  
  
  real<lower = 5> Vad[nSubject];
  real<lower = 0> Qli[nSubject];
  real<lower = 0> Qlu[nSubject];
  
  // drug-related parameters
  real CLrenal = CLrenalPrior;
  real pKa = pKaPrior;
  //real fu = fuPrior;
  real logP = logPPrior;
  real BP = BPPrior;
  vector<lower = 0>[nOrgans] Kp;
  real Kpad;
  real Kpbo;
  real Kpbr;
  real Kpgu;
  real Kphe;
  real Kpki;
  real Kpli;
  real Kplu;
  real Kpmu;
  real Kpsk;
  real Kpsp;
  
  Kp = calcKpPT(logP, pKa, fu, BP, TC);
  Kpad = Kp[1];
  Kpbo = Kp[2];
  Kpbr = Kp[3];
  Kpgu = Kp[4];
  Kphe = Kp[5];
  Kpki = Kp[6];
  Kpli = Kp[7];
  Kplu = Kp[8];
  Kpmu = Kp[9];
  Kpsk = Kp[11];  // skin comes after spleen in TC data
  Kpsp = Kp[10];
  
  thetaHat[1] = VmuHat;
  thetaHat[2] = VboHat;
  thetaHat[3] = QguHat;
  thetaHat[4] = QkiHat;
  thetaHat[5] = CLhepaticHat;
  thetaHat[6] = kaHat;

  theta = (rep_matrix(thetaHat, nSubject) .* 
          exp(diag_pre_multiply(omega, L * eta)))';

  for(j in 1:nSubject){
    // calculation of individual parameter values given theta and covariates
    Vmu[j] = theta[j, 1] * (height[j] / 1.76)^2;
    Vbo[j] = theta[j, 2] * (height[j] / 1.76)^2;
    Qgu[j] = theta[j, 3] * (height[j] / 1.76)^0.75;
    Qki[j] = theta[j, 4] * (height[j] / 1.76)^0.75;
    CLhepatic[j] = theta[j, 5];
    ka[j] = theta[j, 6];
    
    // calculation of other dependent individual parameter values
    Vad[j] = weight[j] - (Var[j] + Vbo[j] + Vbr[j] + Vgu[j] + Vhe[j] + Vki[j] + Vli[j] + Vlu[j] + Vmu[j] + Vsk[j] + Vsp[j] + Vve[j]);
    Qli[j] = Qsp[j] + Qgu[j] + Qha[j];
    Qlu[j] = Qad[j] + Qbo[j] + Qbr[j] + Qhe[j] + Qki[j] + Qli[j] + Qmu[j] + Qsk[j];

    // Filling coefficient matrix
    K = rep_matrix(0, nCmt, nCmt);
    
    K[1, 1] = -Qad[j] / Vad[j] / (Kpad/BP);
    K[1, 2] = Qad[j] / Var[j];
    K[2, 1] = -Qlu[j] / Var[j];
    K[2, 10] = Qlu[j] / Vlu[j] / (Kplu/BP);
    K[3, 2] = Qbo[j] / Var[j];
    K[3, 3] = -Qbo[j] / Vbo[j] / (Kpbo/BP);
    K[4, 2] = Qbr[j] / Var[j];
    K[4, 4] = -Qbr[j] / Vbr[j] / (Kpbr/BP);
    K[5, 5] = -ka[j];
    K[6, 2] = Qgu[j] / Var[j];
    K[6, 5] = ka[j];
    K[6, 6] = -Qgu[j] / Vgu[j] / (Kpgu/BP);
    K[7, 2] = Qhe[j] / Var[j];
    K[7, 7] = -Qhe[j] / Vhe[j] / (Kphe/BP);
    K[8, 2] = Qki[j] / Var[j];
    K[8, 8] = -(Qki[j] + fu * CLrenal) / Vki[j] / (Kpki/BP);
    K[9, 2] = Qha[j] / Var[j];
    K[9, 6] = Qgu[j] / Vgu[j] / (Kpgu/BP);
    K[9, 9] = -(Qli[j] + fu * CLhepatic[j]) / Vli[j] / (Kpli/BP);
    K[9, 13] = Qsp[j] / Vsp[j] / (Kpsp/BP);
    K[10, 10] = -Qlu[j] / Vlu[j] / (Kplu/BP);
    K[10, 14] = Qlu[j] / Vve[j];
    K[11, 2] = Qmu[j] / Var[j];
    K[11, 11] = -Qmu[j] / Vmu[j] / (Kpmu/BP);
    K[12, 2] = Qsk[j] / Var[j];
    K[12, 12] = -Qsk[j] / Vsk[j] / (Kpsk/BP);
    K[13, 2] = Qsp[j] / Var[j];
    K[13, 13] = -Qsp[j] / Vsp[j] / (Kpsp/BP);
    K[14, 1] = Qad[j] / Vad[j] / (Kpad/BP);
    K[14, 3] = Qbo[j] / Vbo[j] / (Kpbo/BP);
    K[14, 4] = Qbr[j] / Vbr[j] / (Kpbr/BP);
    K[14, 7] = Qhe[j] / Vhe[j] / (Kphe/BP);
    K[14, 8] = Qki[j] / Vki[j] / (Kpki/BP);
    K[14, 9] = Qli[j] / Vli[j] / (Kpli/BP);
    K[14, 11] = Qmu[j] / Vmu[j] / (Kpmu/BP);
    K[14, 12] = Qsk[j] / Vsk[j] / (Kpsk/BP);
    K[14, 14] = -Qlu[j] / Vve[j];
    
    x[,start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]],
                                          amt[start[j]:end[j]],
                                          rate[start[j]:end[j]],
                                          ii[start[j]:end[j]],
                                          evid[start[j]:end[j]],
                                          cmt[start[j]:end[j]],
                                          addl[start[j]:end[j]],
                                          ss[start[j]:end[j]],
                                          K, F, tlag);
                                          
    for(k in start[j]:end[j]){
      cHat[k] = fmax(machine_precision(), x[14, k]) / Vve[j] / BP;  // divide by BP to get plasma conc from blood conc
    }

    //cHat[start[j]:end[j]] = (x[14, start[j]:end[j]] ./ Vve[j]) ./ BP;  // divide by BP to get plasma conc from blood conc
  }
  //print(Kp);
}

model{
    // Priors
    // informative
    VmuHat ~ lognormal(log(VmuMeanPrior), log(VmuGSDPrior));
    VboHat ~ normal(VboMeanPrior, VboCVPrior*VboMeanPrior/100); 
    VboHat ~ normal(VboMeanPrior, VboCVPrior*VboMeanPrior/100);
    QguHat ~ normal(0.16*COMeanPrior, COCVPrior*0.16*COMeanPrior/100);
    QkiHat ~ normal(0.19*COMeanPrior, COCVPrior*0.19*COMeanPrior/100);
    
    // weakly informative
    kaHat ~ normal(0, 2); 
    CLhepaticHat ~ uniform(CLhepaticLower, CLhepaticUpper);
    //logP ~ uniform(logPLower, logPUpper);
    fu ~ uniform(fuLower, fuUpper);

    //  rho ~ lkj_corr(1); 
    L ~ lkj_corr_cholesky(1);
    omega ~ cauchy(0, 2);
    to_vector(eta) ~ normal(0, 1);
    
    // residual error
    sigma ~ cauchy(0, 2);
    
    // likelihood
    //logCObs ~ normal(log(cObs), sigma);
    cObs ~ lognormal(log(cHat[iObs]), sigma);  //lognormal dist
    // target += -log(fabs(cObs)); //log(Jacobian) to convert to likelihood wrt cObs
}


generated quantities{
  row_vector<lower = 0>[nt] cHatPred;
  vector<lower = 0>[nt] cObsCond;
  vector<lower = 0>[nt] cObsPred;
  matrix<lower = 0>[nSubject, nRandom] thetaPred;
  // vector[nObs] log_lik;
  
  // Individual-level model parameters
  vector<lower = 0>[nSubject] VmuPred;
  vector<lower = 0>[nSubject] VboPred;
  vector<lower = 0>[nSubject] QguPred;
  vector<lower = 0>[nSubject] QkiPred;
  vector<lower = 0>[nSubject] CLhepaticPred;
  vector<lower = 0>[nSubject] kaPred;
  
  vector<lower = 5>[nSubject] VadPred;
  vector<lower = 0>[nSubject] QliPred;
  vector<lower = 0>[nSubject] QluPred;
  
  // random effects
  matrix[nRandom, nSubject] etaPred;
  corr_matrix[nRandom] rho;
  
  matrix[nCmt, nt] xPred;
  matrix[nCmt, nCmt] KPred;
  
  rho = L * L';
  for(j in 1:nSubject){
    for(i in 1:nRandom)
      etaPred[i, j] = normal_rng(0, 1);
  }
  
  thetaPred = (rep_matrix(thetaHat, nSubject) .* 
          exp(diag_pre_multiply(omega, L * etaPred)))';
  
  for(j in 1:nSubject){
    // calculation of individual parameter values given theta and covariates
    VmuPred[j] = thetaPred[j, 1] * (height[j] / 1.76)^2;
    VboPred[j] = thetaPred[j, 2] * (height[j] / 1.76)^2;
    QguPred[j] = thetaPred[j, 3] * (height[j] / 1.76)^0.75;
    QkiPred[j] = thetaPred[j, 4] * (height[j] / 1.76)^0.75;
    CLhepaticPred[j] = thetaPred[j, 5];
    kaPred[j] = thetaPred[j, 6];
    
    // calculation of other dependent individual parameter values
    VadPred[j] = weight[j] - (Var[j] + VboPred[j] + Vbr[j] + Vgu[j] + Vhe[j] + Vki[j] + Vli[j] + Vlu[j] + VmuPred[j] + Vsk[j] + Vsp[j] + Vve[j]);
    QliPred[j] = Qsp[j] + QguPred[j] + Qha[j];
    QluPred[j] = Qad[j] + Qbo[j] + Qbr[j] + Qhe[j] + QkiPred[j] + QliPred[j] + Qmu[j] + Qsk[j];

    // Filling coefficient matrix
    KPred = rep_matrix(0, nCmt, nCmt);
    
    KPred[1, 1] = -Qad[j] / VadPred[j] / (Kpad/BP);
    KPred[1, 2] = Qad[j] / Var[j];
    KPred[2, 1] = -QluPred[j] / Var[j];
    KPred[2, 10] = QluPred[j] / Vlu[j] / (Kplu/BP);
    KPred[3, 2] = Qbo[j] / Var[j];
    KPred[3, 3] = -Qbo[j] / VboPred[j] / (Kpbo/BP);
    KPred[4, 2] = Qbr[j] / Var[j];
    KPred[4, 4] = -Qbr[j] / Vbr[j] / (Kpbr/BP);
    KPred[5, 5] = -kaPred[j];
    KPred[6, 2] = QguPred[j] / Var[j];
    KPred[6, 5] = kaPred[j];
    KPred[6, 6] = -QguPred[j] / Vgu[j] / (Kpgu/BP);
    KPred[7, 2] = Qhe[j] / Var[j];
    KPred[7, 7] = -Qhe[j] / Vhe[j] / (Kphe/BP);
    KPred[8, 2] = QkiPred[j] / Var[j];
    KPred[8, 8] = -(QkiPred[j] + fu * CLrenal) / Vki[j] / (Kpki/BP);
    KPred[9, 2] = Qha[j] / Var[j];
    KPred[9, 6] = QguPred[j] / Vgu[j] / (Kpgu/BP);
    KPred[9, 9] = -(QliPred[j] + fu * CLhepaticPred[j]) / Vli[j] / (Kpli/BP);
    KPred[9, 13] = Qsp[j] / Vsp[j] / (Kpsp/BP);
    KPred[10, 10] = -QluPred[j] / Vlu[j] / (Kplu/BP);
    KPred[10, 14] = QluPred[j] / Vve[j];
    KPred[11, 2] = Qmu[j] / Var[j];
    KPred[11, 11] = -Qmu[j] / VmuPred[j] / (Kpmu/BP);
    KPred[12, 2] = Qsk[j] / Var[j];
    KPred[12, 12] = -Qsk[j] / Vsk[j] / (Kpsk/BP);
    KPred[13, 2] = Qsp[j] / Var[j];
    KPred[13, 13] = -Qsp[j] / Vsp[j] / (Kpsp/BP);
    KPred[14, 1] = Qad[j] / VadPred[j] / (Kpad/BP);
    KPred[14, 3] = Qbo[j] / VboPred[j] / (Kpbo/BP);
    KPred[14, 4] = Qbr[j] / Vbr[j] / (Kpbr/BP);
    KPred[14, 7] = Qhe[j] / Vhe[j] / (Kphe/BP);
    KPred[14, 8] = QkiPred[j] / Vki[j] / (Kpki/BP);
    KPred[14, 9] = QliPred[j] / Vli[j] / (Kpli/BP);
    KPred[14, 11] = Qmu[j] / VmuPred[j] / (Kpmu/BP);
    KPred[14, 12] = Qsk[j] / Vsk[j] / (Kpsk/BP);
    KPred[14, 14] = -QluPred[j] / Vve[j];
    
    xPred[,start[j]:end[j]] = pmx_solve_linode(time[start[j]:end[j]],
                                      amt[start[j]:end[j]],
                                      rate[start[j]:end[j]],
                                      ii[start[j]:end[j]],
                                      evid[start[j]:end[j]],
                                      cmt[start[j]:end[j]],
                                      addl[start[j]:end[j]],
                                      ss[start[j]:end[j]],
                                      KPred, F, tlag);

    for(k in start[j]:end[j]){
      cHatPred[k] = fmax(machine_precision(), xPred[14, k]) / Vve[j] / BP;  // divide by BP to get plasma conc from blood conc
      //print(cHatPred);
      
    if(time[k] == 0){
      cObsCond[k] = 0;
      cObsPred[k] = 0;
    }else{
      cObsCond[k] = exp(normal_rng(log(cHat[k]), sigma)); # individual predictions
      cObsPred[k] = exp(normal_rng(log(cHatPred[k]), sigma)); # population predictions
    }
  }
    // cHatPred[start[j]:end[j]] = (xPred[14, start[j]:end[j]] ./ Vve[j]) ./ BP;  // divide by BP to get plasma conc from blood conc
  }
}