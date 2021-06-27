functions{
  real[] PBPKModelODE(real t, real[] x, real[] parms, real[] rdummy, int[] idummy){
    real dxdt[18];
    real EXPOSURE = rdummy[2];
    
    // Skin exposure parameters.
    AEXP = 20.0;         // Area of exposed skin (cm^2).
    VWELL = 1.0;        // Volume of material in exposure well (mL).
    
    // Basic anatomical and physiological parameters.
    real BW = 63.5 * 1000;   // Body mass (g).
    real QPUc = 1.2083;   // Cardiac output (mL/min/g^0.75).
    
    // Skin parameters.
    //real DSC = 1.0e-8;   // Fickian diffusion constant for stratum corneum (cm^2/min).
    real TSC = 0.003;    // Thickness of stratum corneum (cm).
    real TVE = 0.0075;   // Thickness of viable epidermis (cm).
    
    // Volumes (as fractions of body mass) for ...
    real FTLI = 0.026;       // ... liver.
    real FFAT = 0.136;       // ... fat.
    real FTPP = 0.487;      // ... poorly perfused tissue.
    real FTRP = 0.1861;       // ... richly perfused tissue.
    real FTABD = 0.022;      // ... arterial blood.
    real FTVBD = 0.045;      // ... venous blood.
    real FTREM = 0.09;       // ... remaining tissue.
    real FTPU = 0.0052      // ... pulmonary tissue
    
    // Partition coefficients for ...
    real HBA = 571.0;        // ... blood:air.
    real HTB = 3.49;         // ... general tissue:blood (incl. RP, PP, respiratory).
    real HLUB = 1.71;        // ... lung:blood.
    real HLB = 1.61;         // ... liver:blood.
    real HFB = 49.0;         // ... fat:blood.
    real HRPB = 2.12;        // ... richly perfused tissue:blood.
    real HPPB = 3.49;        // ... poorly perfused tissue:blood.
    real HVEB = 2.73;        // ... VE:blood (Eq. 23 of McCarley & Bunge, 2001)
    real HSCVE = 1.98;       // ... SC:VE (Eq. 22 of McCarley & Bunge, 2001)
    //real HSCJP8 = 1.0;       // ... SC:JP8.
    
    // Metabolism parameters.
    real VmaxLI = 100.8;      // Vmax in liver (nmol/min/mL).
    real KmLI = 25.0;         // Km in liver (nmol/mL).
    real VmaxLU = 0.27;       // Vmax in lung (nmol/min/mL).
    real KmLU = 82.0;        // Km in lung (nmol/mL).
    
    // Volumes (mL) of ...
    real VTLI = FTLI * BW;           // ... liver tissue.
    real VFAT = FFAT * BW;           // ... fat tissue.
    real VTSC = AEXP * TSC;              // ... (exposed) stratum corneum.
    real VTVE = AEXP * TVE;              // ... (exposed) viable epidermis.
    real VTRP = FTRP * BW - VTVE + VUR;    // ... richly perfused tissue.
    real VTPP = FTPP * BW;           // ... poorly perfused tissue.
    real VTAB = FTABD * BW;          // ... arterial blood.
    real VTVB = FTVBD * BW;          // ... venous blood.
    real VTREM = FTREM * BW;         // ... remaining tissue.
    real VTPU = FTPU * BW;           // ... pulmonary tissue
    
    // Discretization thickness for stratum corneum partial differential
    // equation discretization.
    real SCDX = TSC / 10;
    
    // Cardiac output (mL/min).
    real QPUa = QPUc * pow(BW, 0.75);
    
    // Blood flows (mL/min) to ...
    real QLI = FBLI * QPUa;          // ... liver.
    real QVE = FBRP * (VTVE / VTRP) * QPUa;  // ... viable epidermis.
    real QRP = FBRP * QPUa - QVE;    // ... richly perfused tissue.
    real QPP = FBPP * QPUa;          // ... poorly perfused tissue.
    real QFA = FBFA * QPUa;          // ... fat tissue.
    
    // Total blood flow. COMPARE with QPUa. - DFK 2/15/2018
    real QPU = QLI + QVE + QRP + QPP + QFA;

    //// parameters to be estimated ////
    real HSCJP8 = parms[1];
		real DSC = parms[2]*1e-8;
		
    // model compartments
    real APU = x[1];
    real ALI = x[2];
    real AFA = x[3];
    real ARP = x[4];
    real APP = x[5];
    real AWELL = x[6];
    real CSC01 = x[7];
    real CSC02 = x[8];
    real CSC03 = x[9];
    real CSC04 = x[10];
    real CSC05 = x[11];
    real CSC06 = x[12];
    real CSC07 = x[13];
    real CSC08 = x[14];
    real CSC09 = x[15];
    real AVE = x[16];
    real AAB = x[17];
    real AVB = x[18];
    
    // Naphthalene internal concentrations (nmol/mL) in ...
    real CPU = APU / VTPU;         // ... pulmonary tissue.
    real CvPU = CPU / HLUB;        // ... veins leaving pulmonary tissue.
    real CLI = ALI / VTLI;           // ... liver tissue.
    real CvLI = CLI / HLB;           // ... veins leaving liver tissue.
    real CFA = AFA / VFAT;           // ... fat tissue.
    real CvFA = CFA / HFB;           // ... veins leaving fat tissue.
    real CRP = ARP / VTRP;           // ... richly perfused tissue.
    real CvRP = CRP / HRPB;          // ... veins leaving richly perfused tissue.
    real CPP = APP / VTPP;           // ... poorly perfused tissue.
    real CvPP = CPP / HPPB;          // ... veins leaving poorly perfused tissue.
    real CAB = AAB / VTAB;           // ... arterial blood.
    real CVB = AVB / VTVB;           // ... venous blood.
    real CVE = VTVE > 0 ? AVE / VTVE : 0.0;           // ... viable epidermis.
    real CvVE = CVE / HVEB;          // ... veins leaving viable epidermis.
    real CWELL = (AWELL / VWELL) * EXPOSURE;     // ... skin exposure well.
    real CSC00 = (CWELL * HSCJP8) * EXPOSURE + CSC01*(1 - EXPOSURE);  // ... outer surface of SC.
    real CSC10 = CVE * HSCVE;        // ... interface of SC with VE.
    
    // ----- Metabolism -----
    real RMLPU = VTPU * VmaxLU * CvPU / (KmLU + CvPU);  // pulmonary tissue
    real RML = VTLI * VmaxLI * CvLI / (KmLI + CvLI);      // Liver
    
    // ----- Fluxes -----
    // Naphthalene fluxes (nmol/cm^2/min) at ...
    real JSC00 = -DSC * (CSC01 - CSC00) / SCDX;   // ... outer surface of SC.
    real JSC10 = -DSC * (CSC10 - CSC09) / SCDX;   // ... interface of SC with VE.

    // ----- Time rates of change of state variables (ODEs) in nmol/min  -----
    dxdt[1] = QPU * (CVB - CvPU) - RMLPU;                         // pulmonary tissue
    dxdt[2] = QLI * (CAB - CvLI) - RML;                           // liver tissue
    dxdt[3] = QFA * (CAB - CvFA);                                 // fat tissue
    dxdt[4] = QRP * (CAB - CvRP);                                 // richly perfused tissue
    dxdt[5] = QPP * (CAB - CvPP);                                 // poorly perfused tissue
    dxdt[6] = -JSC00 * AEXP;                                      // Dermal Exposure Well
    dxdt[7] = DSC * (CSC00 - 2 * CSC01 + CSC02) / (SCDX * SCDX);  // Stratum Corneum 1
    dxdt[8] = DSC * (CSC01 - 2 * CSC02 + CSC03) / (SCDX * SCDX);  // Stratum Corneum 2 CSC02
    dxdt[9] = DSC * (CSC02 - 2 * CSC03 + CSC04) / (SCDX * SCDX);  // Stratum Corneum 3
    dxdt[10] = DSC * (CSC03 - 2 * CSC04 + CSC05) / (SCDX * SCDX); // Stratum Corneum 4
    dxdt[11] = DSC * (CSC04 - 2 * CSC05 + CSC06) / (SCDX * SCDX); // Stratum Corneum 5
    dxdt[12] = DSC * (CSC05 - 2 * CSC06 + CSC07) / (SCDX * SCDX); // Stratum Corneum 6
    dxdt[13] = DSC * (CSC06 - 2 * CSC07 + CSC08) / (SCDX * SCDX); // Stratum Corneum 7
    dxdt[14] = DSC * (CSC07 - 2 * CSC08 + CSC09) / (SCDX * SCDX); // Stratum Corneum 8
    dxdt[15] = DSC * (CSC08 - 2 * CSC09 + CSC10) / (SCDX * SCDX); // Stratum Corneum 9
    dxdt[16] = JSC10 * AEXP + QVE * (CAB - CvVE);                 // viable epidermis
    dxdt[17] = QPU * (CvPU - CAB);                                // arterial blood
    dxdt[18] = QLI * CvLI + QFA * CvFA + QRP * CvRP + QPP * CvPP + QVE * CvVE - QPU * CVB;  // venous blood

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
}

transformed data{
  vector[nObs] logCObs = log(cObs);
  int nTheta = 2;  // number of parameters
  int nCmt = 18;  // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];
  real VTVB;
  
  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0> HSCJP8;
  real<lower = 0> NDSC;
  
  // residual error
  real<lower = 0> sigma;
}

transformed parameters{
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nCmt, nt] x;
  real<lower = 0> parms[nTheta];
  
  parms[1] = HSCJP8;
  parms[2] = NDSC;

  x = pmx_solve_rk45(PBPKModelODE, nCmt,
                     time, amt, rate, ii, evid, cmt, addl, ss,
                     parms, biovar, tlag,
                     1e-6, 1e-6, 1e6);
                
  cHat = x[18] / VTVB;  // divide by subject's blood volume VTVB
  
  cHatObs = cHat[iObs];
}

model{
  // Priors
  HSCJP8 ~ lognormal(log(1.0), 0.25);
  NDSC ~ lognormal(log(1.0), 0.25);

  sigma ~ cauchy(0.0, 0.5);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
}

generated quantities{
  matrix[nCmt, nt] xPred;
  real<lower = 0> parmsPred[nTheta];
  row_vector<lower = 0>[nt] cHatPred;
  row_vector<lower = 0>[nObs] cHatObsPred;
  vector<lower = 0>[nObs] cObsCond;
  row_vector<lower = 0>[nObs] cObsPred;
  
  parmsPred[1] = HSCJP8;
  parmsPred[2] = NDSC; 

  xPred = pmx_solve_group_rk45(PBPKModelODE, nCmt, len,
                               time, amt, rate, ii, evid, cmt, addl, ss,
                               parmsPred, biovar, tlag, WT, idummy,
                               1e-6, 1e-6, 1e6);

  cHatPred = xPred[18] / VTVB;

  // predictions for observed data records
  cHatObsPred = cHatPred[iObs];

  for(i in 1:nObs) {
    cObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));  // individual predictions
    cObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));  // population predictions
  }
}
