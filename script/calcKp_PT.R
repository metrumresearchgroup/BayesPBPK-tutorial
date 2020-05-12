#calculate tissue:plasma partition coefficients based on: Poulin and Theil http://jpharmsci.org/article/S0022-3549(16)30889-9/fulltext
#the function returns a list of parameters that can be used directly to update the param() function in mrgsolve

library(dplyr)

calcKp_PT <- function(logP, pKa, fup, BP=1, type=1, dat){
  
  dat_all <- dat %>% filter(!tissue %in% c("Plasma","Adipose","RBCs"))
  
  n <- length(dat$tissue)
  Kp_all <- vector(mode = "numeric", length = n)
  
  Vwp <- dat$f_water[dat$tissue == "Plasma"]
  Vnlp <- dat$f_n_l[dat$tissue == "Plasma"]
  Vphp <- dat$f_pl[dat$tissue == "Plasma"]
  
  dat2 <- dat %>% filter(!tissue %in% c("Plasma","RBCs"))
  
  Vwt <- dat2$f_water[dat2$tissue != "Adipose"]
  Vwad <- dat2$f_water[dat2$tissue == "Adipose"]
  Vnlt <- dat2$f_n_l[dat2$tissue != "Adipose"]
  Vnlad <- dat2$f_n_l[dat2$tissue == "Adipose"]
  Vpht <- dat2$f_pl[dat2$tissue != "Adipose"]
  Vphad <- dat2$f_pl[dat2$tissue == "Adipose"]
  
  pH <- dat$pH[dat$tissue == "Adipose"]
  logD <- 1.115*logP-1.35 #logD is the olive oil:buffer(water) partition coefficient of nonionized species
  
  logD_star <- switch(type,
                      #1-neutral
                      logD,   
                      #2-monoprotic acid
                      logD-log10(1+10^(pH-pKa)),
                      #3-monoprotic base
                      logD-log10(1+10^(pKa-pH)), 
                      #4-diprotic acid
                      logD-log10(1+10^(2*pH-pKa[1]-pKa[2])),
                      #5-diprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]-2*pH)), 
                      #6-monoprotic acid monoprotic base (acid comes first)
                      logD-log10(1+10^(pKa[2]-pKa[1])),  
                      #7-triprotic acid
                      logD-log10(1+10^(3*pH-pKa[1]-pKa[2]-pKa[3])),  
                      #8-triprotic base
                      logD-log10(1+10^(pKa[1]+pKa[2]+pKa[3]-3*pH)),  
                      #9-diprotic acid monoprotic base (first two are acid)
                      logD-log10(1+10^(pH-pKa[1]-pKa[2]+pKa[3])), 
                      #10-diprotic base monoprotic acid (first one is acid)
                      logD-log10(1+10^(pKa[2]+pKa[3]-pKa[1]-pH)))       
  
  D_star <- 10^logD_star   
  Kpad <- ((D_star*(Vnlad+0.3*Vphad)+(1*(Vwad+0.7*Vphad)))/(D_star*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * fup
  
  P <- 10^logP
  fut <- 1/(1+((1-fup)/fup)*0.5)
  Kpt <- ((P*(Vnlt+0.3*Vpht)+(1*(Vwt+0.7*Vpht)))/(P*(Vnlp+0.3*Vphp)+(1*(Vwp+0.7*Vphp)))) * (fup/fut)
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kpad,Kpt))
  names(Kp) <- nms
  
  return(Kp)
  
  
}