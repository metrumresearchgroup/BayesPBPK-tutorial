## script to create mavoPBPK linear ODE model

## create coeff matrix

getMatrix <- function(odes){
  l <- unlist(strsplit(odes, ";"))
  cmts <- as.character(lapply(l, function(x) str_match(x, "dxdt_\\s*(.*?)\\s*=")[2]))
  ncmts <- length(cmts)
  K <- matrix("0", ncmts, ncmts)
  #Ks <- matrix("0", ncmts, ncmts)
  for(i in 1:ncmts){
    cmt <- cmts[i]
    rhs <- gsub(".*=", "", l[i])
    for(j in 1:ncmts){
      lrhs <- unlist(strsplit(rhs, "\\([^)]*\\)(*SKIP)(*F)|(?<=.)(?=[-|\\+])", perl = TRUE))
      for(k in 1:length(lrhs)){
        if(grepl(cmts[j], lrhs[k])) K[i,j] <- lrhs[k]
        K[i,j] <- gsub(cmts[j], 1, K[i,j])
        #K[i,j] <- str_replace_all(string=K[i,j], pattern=" ", repl="")
      }
    }
  }
  return(K)
}

m <- getMatrix(odes)
df <- as_tibble(m)

getMatrixText <- function(m){
  ms <- m
  for(i in 1:nrow(m)){
    for(j in 1:ncol(m)){
      flux <- m[i,j]
      ms[i,j] <- sprintf("K[%s,%s] = %s;", i, j, flux)
    }
  }
  v <- c(ms)
  v2 <- v[!(grepl("= 0", v))]
  return(v2)
}

# get matrix text
odes <- "dxdt_Lungs = QLU*Venous_Blood/VVB - QLU*Lungs/KbLU/VLU; 
dxdt_Heart = QHT*Arterial_Blood/VAB - QHT*Heart/KbHT/VHT; 
dxdt_Brain = QBR*Arterial_Blood/VAB - QBR*Brain/KbBR/VBR;  
dxdt_Muscles = QMU*Arterial_Blood/VAB - QMU*Muscles/KbMU/VMU;  
dxdt_Adipose = QAD*Arterial_Blood/VAB - QAD*Adipose/KbAD/VAD;  
dxdt_Skin = QSK*Arterial_Blood/VAB - QSK*Skin/KbSK/VSK;  
dxdt_Spleen = QSP*Arterial_Blood/VAB - QSP*Spleen/KbSP/VSP;  
dxdt_Pancreas = QPA*Arterial_Blood/VAB - QPA*Pancreas/KbPA/VPA;  
dxdt_Liver = QHA*Arterial_Blood/VAB + QSP*Spleen/KbSP/VSP + QPA*Pancreas/KbPA/VPA + QST*Stomach/KbST/VST + QGU*Gut/KbGU/VGU - (CLint*fub + QLI)*Liver/KbLI/VLI;
dxdt_Stomach = QST*Arterial_Blood/VAB - QST*Stomach/KbST/VST; 
dxdt_Gut = QGU*Arterial_Blood/VAB - QGU*Gut/KbGU/VGU; 
dxdt_Bones = QBO*Arterial_Blood/VAB - QBO*Bones/KbBO/VBO;  
dxdt_Kidneys = QKI*Arterial_Blood/VAB - QKI*Kidneys/KbKI/VKI;  
dxdt_Arterial_Blood = QLU*Lungs/KbLU/VLU - QLU*Arterial_Blood/VAB;  
dxdt_Venous_Blood = QHT*Heart/KbHT/VHT + QBR*Brain/KbBR/VBR + QMU*Muscles/KbMU/VMU + QAD*Adipose/KbAD/VAD + QSK*Skin/KbSK/VSK + QLI*Liver/KbLI/VLI + QBO*Bones/KbBO/VBO + QKI*Kidneys/KbKI/VKI + QRB*Rest_of_Body/KbRB/VRB - QLU*Venous_Blood/VVB;
dxdt_Rest_of_Body = QRB*Arterial_Blood/VAB - QRB*Rest_of_Body/KbRB/VRB;"

K <- getMatrix(odes)
Ks <- getMatrixText(K)

## update model file
modelDir <- "../model"
modelName <- "mavoPBPKLinODE"
m_raw <- readLines(file.path(modelDir, paste0(modelName, "Raw.stan")))
m <- append(m_raw[1:166], Ks)
m <- append(m, m_raw[167:230])
m <- append(m, Ks)
m <- append(m, m_raw[231:length(m_raw)])

# save new model file
write_lines(m, file = file.path(modelDir, paste0(modelName, ".stan")))



