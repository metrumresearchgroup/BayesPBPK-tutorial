function PBPKODE!(du, u, p, t)
    ## parameters to be estimated ##
    CLint = p[1];
    KbBR = p[2];
    KbMU = p[3];
    KbAD = p[4];
    KbBO = p[5];
    KbRB = p[6];
    
    ## fixed parameters ##
    WT = p[7];
    k₀ = p[8];

    # Regional blood flows
    CO  = (187.00*WT^0.81)*60/1000;         # Cardiac output (L/h) from White et al (1968)
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
    VMU = (40.00*WT/100)/1.041;
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
    
    # model compartments
    Lungs = u[1];
    Heart = u[2];
    Brain = u[3];
    Muscles = u[4];
    Adipose = u[5];
    Skin = u[6];
    Spleen = u[7];
    Pancreas = u[8];
    Liver = u[9];
    Stomach = u[10];
    Gut = u[11];
    Bones = u[12];
    Kidneys = u[13];
    Arterial_Blood = u[14];
    Venous_Blood = u[15];
    Rest_of_Body = u[16];

    du[1] = dLungs = QLU*(Venous_Blood/VVB - Lungs/KbLU/VLU);  #lungs
    du[2] = dHeart = QHT*(Arterial_Blood/VAB - Heart/KbHT/VHT);  #heart
    du[3] = dBrain = QBR*(Arterial_Blood/VAB - Brain/KbBR/VBR);  #brain
    du[4] = dMuscles = QMU*(Arterial_Blood/VAB - Muscles/KbMU/VMU);  #muscles
    du[5] = dAdipose = QAD*(Arterial_Blood/VAB - Adipose/KbAD/VAD);  #adipose
    du[6] = dSkin = QSK*(Arterial_Blood/VAB - Skin/KbSK/VSK);  #skin
    du[7] = dSpleen = QSP*(Arterial_Blood/VAB - Spleen/KbSP/VSP);  #spleen
    du[8] = dPancreas = QPA*(Arterial_Blood/VAB - Pancreas/KbPA/VPA);  #pancreas
    du[9] = dLiver = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP + QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU - CLint*fub*Liver/KbLI/VLI - QLI*Liver/KbLI/VLI;  #liver
    du[10] = dStomach = QST*(Arterial_Blood/VAB - Stomach/KbST/VST);  #stomach
    du[11] = dGut = QGU*(Arterial_Blood/VAB - Gut/KbGU/VGU);  #gut
    du[12] = dBones = QBO*(Arterial_Blood/VAB - Bones/KbBO/VBO);  #bones
    du[13] = dKidneys = QKI*(Arterial_Blood/VAB - Kidneys/KbKI/VKI);  #kidneys
    du[14] = dArterial_Blood = QLU*(Lungs/KbLU/VLU - Arterial_Blood/VAB);  #arterial blood
    du[15] = dVenous_Blood = k₀ + QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR + QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;  #venous_blood
    du[16] = dRest_of_Body = QRB*(Arterial_Blood/VAB - Rest_of_Body/KbRB/VRB);  #rest of body
end