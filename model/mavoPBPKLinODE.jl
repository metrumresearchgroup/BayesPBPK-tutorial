function PBPK(p)
  ## parameters to be estimated ##
  #CLint = Float64(p[1]);
  #KbBR = Float64(p[2]);
  #KbMU = Float64(p[3]);
  #KbAD = Float64(p[4]);
  #KbBO = Float64(p[5]);
  #KbRB = Float64(p[6]);
  
  ## fixed parameters ##
  #WT = Float64(p[7]);

  CLint = p[1];
  KbBR = p[2];
  KbMU = p[3];
  KbAD = p[4];
  KbBO = p[5];
  KbRB = p[6];
  
  ## fixed parameters ##
  WT = p[7];

  # Regional blood flows
  CO  = (187.0*WT^0.81)*60/1000;         # Cardiac output (L/h) from White et al (1968)
  QHT = 4.0 *CO/100;
  QBR = 12.0*CO/100;
  QMU = 17.0*CO/100;
  QAD = 5.0 *CO/100;
  QSK = 5.0 *CO/100;
  QSP = 3.0 *CO/100;
  QPA = 1.0 *CO/100;
  QLI = 25.5*CO/100;
  QST = 1.0 *CO/100;
  QGU = 14.0*CO/100;
  QHA = QLI - (QSP + QPA + QST + QGU); # Hepatic artery blood flow
  QBO = 5.0 *CO/100;
  QKI = 19.0*CO/100;
  QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
  QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

  # Organs' volumes = organs' weights / organs' density
  VLU = (0.76 *WT/100)/1.051;
  VHT = (0.47 *WT/100)/1.030;
  VBR = (2.00 *WT/100)/1.036;
  VMU = (40.0*WT/100)/1.041;
  VAD = (21.42*WT/100)/0.916;
  VSK = (3.71 *WT/100)/1.116;
  VSP = (0.26 *WT/100)/1.054;
  VPA = (0.14 *WT/100)/1.045;
  VLI = (2.57 *WT/100)/1.040;
  VST = (0.21 *WT/100)/1.050;
  VGU = (1.44 *WT/100)/1.043;
  VBO = (14.29*WT/100)/1.990;
  VKI = (0.44 *WT/100)/1.050;
  VAB = (2.81 *WT/100)/1.040;
  VVB = (5.62 *WT/100)/1.040;
  VRB = (3.86 *WT/100)/1.040;

  # partition coefficients
  KbLU = exp(0.8334);
  KbHT = exp(1.1205);
  KbSK = exp(-0.5238);
  KbSP = exp(0.3224);
  KbPA = exp(0.3224);
  KbLI = exp(1.7604);
  KbST = exp(0.3224);
  KbGU = exp(1.2026);
  KbKI = exp(1.3171);
  
  # Other parameters
  BP = 0.61;      # Blood:plasma partition coefficient
  fup = 0.028;    # Fraction unbound in plasma
  fub = fup/BP;   # Fraction unbound in blood

  ncmt = 16;
  K = zeros(ncmt, ncmt);
      
  K[1,1] = -QLU/KbLU/VLU;
  K[14,1] =  QLU/KbLU/VLU;
  K[2,2] = -QHT/KbHT/VHT;
  K[15,2] =  QHT/KbHT/VHT;
  K[3,3] = -QBR/KbBR/VBR;
  K[15,3] = QBR/KbBR/VBR;
  K[4,4] = -QMU/KbMU/VMU;
  K[15,4] = QMU/KbMU/VMU;
  K[5,5] = -QAD/KbAD/VAD;
  K[15,5] = QAD/KbAD/VAD;
  K[6,6] = -QSK/KbSK/VSK;
  K[15,6] = QSK/KbSK/VSK;
  K[7,7] = -QSP/KbSP/VSP;
  K[9,7] = QSP/KbSP/VSP;
  K[8,8] = -QPA/KbPA/VPA;
  K[9,8] = QPA/KbPA/VPA;
  K[9,9] = -(CLint*fub + QLI)/KbLI/VLI;
  K[15,9] = QLI/KbLI/VLI;
  K[9,10] = QST/KbST/VST;
  K[10,10] = -QST/KbST/VST;
  K[9,11] = QGU/KbGU/VGU;
  K[11,11] = -QGU/KbGU/VGU;
  K[12,12] = -QBO/KbBO/VBO;
  K[15,12] = QBO/KbBO/VBO;
  K[13,13] = -QKI/KbKI/VKI;
  K[15,13] = QKI/KbKI/VKI;
  K[2,14] =  QHT/VAB;
  K[3,14] =  QBR/VAB;
  K[4,14] =  QMU/VAB;
  K[5,14] =  QAD/VAB;
  K[6,14] =  QSK/VAB;
  K[7,14] =  QSP/VAB;
  K[8,14] =  QPA/VAB;
  K[9,14] =  QHA/VAB;
  K[10,14] =  QST/VAB;
  K[11,14] =  QGU/VAB;
  K[12,14] =  QBO/VAB;
  K[13,14] =  QKI/VAB;
  K[14,14] = -QLU/VAB;
  K[16,14] =  QRB/VAB;
  K[1,15] =  QLU/VVB;
  K[15,15] = -QLU/VVB;
  K[15,16] = QRB/KbRB/VRB;
  K[16,16] = -QRB/KbRB/VRB;

  return(K)
end

function PBPK(p)
  ## parameters to be estimated ##
  CLint = p[1];
  KbBR = p[2];
  KbMU = p[3];
  KbAD = p[4];
  KbBO = p[5];
  KbRB = p[6];
  
  ## fixed parameters ##
  WT = p[7];

  # Regional blood flows
  CO  = (187.0*WT^0.81)*60/1000;         # Cardiac output (L/h) from White et al (1968)
  QHT = 4.0 *CO/100;
  QBR = 12.0*CO/100;
  QMU = 17.0*CO/100;
  QAD = 5.0 *CO/100;
  QSK = 5.0 *CO/100;
  QSP = 3.0 *CO/100;
  QPA = 1.0 *CO/100;
  QLI = 25.5*CO/100;
  QST = 1.0 *CO/100;
  QGU = 14.0*CO/100;
  QHA = QLI - (QSP + QPA + QST + QGU); # Hepatic artery blood flow
  QBO = 5.0 *CO/100;
  QKI = 19.0*CO/100;
  QRB = CO - (QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI);
  QLU = QHT + QBR + QMU + QAD + QSK + QLI + QBO + QKI + QRB;

  # Organs' volumes = organs' weights / organs' density
  VLU = (0.76 *WT/100)/1.051;
  VHT = (0.47 *WT/100)/1.030;
  VBR = (2.00 *WT/100)/1.036;
  VMU = (40.0*WT/100)/1.041;
  VAD = (21.42*WT/100)/0.916;
  VSK = (3.71 *WT/100)/1.116;
  VSP = (0.26 *WT/100)/1.054;
  VPA = (0.14 *WT/100)/1.045;
  VLI = (2.57 *WT/100)/1.040;
  VST = (0.21 *WT/100)/1.050;
  VGU = (1.44 *WT/100)/1.043;
  VBO = (14.29*WT/100)/1.990;
  VKI = (0.44 *WT/100)/1.050;
  VAB = (2.81 *WT/100)/1.040;
  VVB = (5.62 *WT/100)/1.040;
  VRB = (3.86 *WT/100)/1.040;

  # partition coefficients
  KbLU = exp(0.8334);
  KbHT = exp(1.1205);
  KbSK = exp(-0.5238);
  KbSP = exp(0.3224);
  KbPA = exp(0.3224);
  KbLI = exp(1.7604);
  KbST = exp(0.3224);
  KbGU = exp(1.2026);
  KbKI = exp(1.3171);
  
  # Other parameters
  BP = 0.61;      # Blood:plasma partition coefficient
  fup = 0.028;    # Fraction unbound in plasma
  fub = fup/BP;   # Fraction unbound in blood

  K = [-QLU/KbLU/VLU 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QLU/VVB 0.0;
       0.0 -QHT/KbHT/VHT 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QHT/VAB 0.0 0.0;
       0.0 0.0 -QBR/KbBR/VBR 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QBR/VAB 0.0 0.0;
       0.0 0.0 0.0 -QMU/KbMU/VMU 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QMU/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 -QAD/KbAD/VAD 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QAD/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 -QSK/KbSK/VSK 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QSK/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 -QSP/KbSP/VSP 0.0 0.0 0.0 0.0 0.0 0.0 QSP/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QPA/KbPA/VPA 0.0 0.0 0.0 0.0 0.0 QPA/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 QSP/KbSP/VSP QPA/KbPA/VPA -(CLint*fub+QLI)/KbLI/VLI QST/KbST/VST QGU/KbGU/VGU 0.0 0.0 QHA/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QST/KbST/VST 0.0 0.0 0.0 QST/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QGU/KbGU/VGU 0.0 0.0 QGU/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QBO/KbBO/VBO 0.0 QBO/VAB 0.0 0.0;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QKI/KbKI/VKI QKI/VAB 0.0 0.0;
       QLU/KbLU/VLU 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -QLU/VAB 0.0 0.0;
       0.0 QHT/KbHT/VHT QBR/KbBR/VBR QMU/KbMU/VMU QAD/KbAD/VAD QSK/KbSK/VSK 0.0 0.0 QLI/KbLI/VLI 0.0 0.0 QBO/KbBO/VBO QKI/KbKI/VKI 0.0 -QLU/VVB QRB/KbRB/VRB;
       0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 QRB/VAB 0.0 -QRB/KbRB/VRB]

  return(K)
end


ncmt = 16;
  K = zeros(ncmt, ncmt);
      
  K[1,1] = "-QLU/KbLU/VLU";
  K[14,1] =  "QLU/KbLU/VLU";
  K[2,2] = "-QHT/KbHT/VHT";
  K[15,2] =  "QHT/KbHT/VHT";
  K[3,3] = "-QBR/KbBR/VBR";
  K[15,3] = "QBR/KbBR/VBR";
  K[4,4] = "-QMU/KbMU/VMU";
  K[15,4] = "QMU/KbMU/VMU";
  K[5,5] = "-QAD/KbAD/VAD";
  K[15,5] = "QAD/KbAD/VAD";
  K[6,6] = "-QSK/KbSK/VSK";
  K[15,6] = "QSK/KbSK/VSK";
  K[7,7] = "-QSP/KbSP/VSP";
  K[9,7] = "QSP/KbSP/VSP";
  K[8,8] = "-QPA/KbPA/VPA";
  K[9,8] = "QPA/KbPA/VPA";
  K[9,9] = "-(CLint*fub + QLI)/KbLI/VLI";
  K[15,9] = "QLI/KbLI/VLI";
  K[9,10] = "QST/KbST/VST";
  K[10,10] = "-QST/KbST/VST";
  K[9,11] = "QGU/KbGU/VGU";
  K[11,11] = "-QGU/KbGU/VGU";
  K[12,12] = "-QBO/KbBO/VBO";
  K[15,12] = "QBO/KbBO/VBO";
  K[13,13] = "-QKI/KbKI/VKI";
  K[15,13] = "QKI/KbKI/VKI";
  K[2,14] =  "QHT/VAB";
  K[3,14] =  "QBR/VAB";
  K[4,14] =  "QMU/VAB";
  K[5,14] =  "QAD/VAB";
  K[6,14] =  "QSK/VAB";
  K[7,14] =  "QSP/VAB";
  K[8,14] =  "QPA/VAB";
  K[9,14] =  "QHA/VAB";
  K[10,14] =  "QST/VAB";
  K[11,14] =  "QGU/VAB";
  K[12,14] =  "QBO/VAB";
  K[13,14] =  "QKI/VAB";
  K[14,14] = "-QLU/VAB";
  K[16,14] =  "QRB/VAB";
  K[1,15] =  "QLU/VVB";
  K[15,15] = "-QLU/VVB";
  K[15,16] = "QRB/KbRB/VRB";
  K[16,16] = "-QRB/KbRB/VRB";



#A = DiffEqArrayOperator(K);

#=
    prob_lin = ODEProblem(A, u0, tspan)
    #@time sol_lin = solve(prob_lin, LinearExponential(), dt=0.1) 
    @time sol_lin = solve(prob_lin, LinearExponential(), tstops=0.0:0.1:25.0)
    res_lin = Array(sol_lin)

# calc conc
VVB = (5.62 *WT/100)/1.040;
BP = 0.61;      # Blood:plasma partition coefficient
PRED = res_lin[15,:]/(VVB*BP/1000.0)

plot(sol_lin.t, PRED, yscale=:log10)
scatter!(dat1.TIME[2:end], dat1.DV[2:end], yscale=:log10)
=#



#A = DiffEqArrayOperator(K);

#=
    prob_lin = ODEProblem(A, u0, tspan)
    #@time sol_lin = solve(prob_lin, LinearExponential(), dt=0.1) 
    @time sol_lin = solve(prob_lin, LinearExponential(), tstops=0.0:0.1:25.0)
    res_lin = Array(sol_lin)

# calc conc
VVB = (5.62 *WT/100)/1.040;
BP = 0.61;      # Blood:plasma partition coefficient
PRED = res_lin[15,:]/(VVB*BP/1000.0)

plot(sol_lin.t, PRED, yscale=:log10)
scatter!(dat1.TIME[2:end], dat1.DV[2:end], yscale=:log10)
=#
