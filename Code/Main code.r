pacman::p_load(deSolve,ggplot2,reshape,FME,dplyr,scales,gdata,data.table,ggpubr,gganimate,gifski,ggpmisc)
setwd("G:/My Drive/R project")

#### Start ####
TN.in <- read.csv(file = "G:/My Drive/R project/GitHub/Dynamic_Model/Data/TN_permeate data.csv",header = T, sep = ",")
TN.in <- TN.in[,-(1),drop=FALSE] 
# Assume gaps between two data points are linear #
C.TN.in <- approxfun(x = TN.in[,1], y = TN.in[,2], method = "linear")

#### Setting Parameters ####
InitStocks <- c(
  NCSb1 = 230,
  NCSb2 = 230,
  NCSb3 = 350,
  EQ3 = 50,
  NCSb1.N = 0,
  NCSb2.N = 0
)

Parms <- c(
  Sysstart = 10,
  Sysdown = 14,
  LMH = 5 / 60, 
  NCSCapa = 350, NCSmax = 360, 
  
  ZDens = 1.6, # Zeolite density (kg/L)
  Zmass = 175, # Zeolite mass in one bed (kg)
  GDens = 1, # GAC density (kg/L)
  Gmass = 80, # GAC mass in one bed (kg)
  
  ZNCapa = 16.17, # mg/g
  ZabsorbKL = 0.002 # Reference Zhang et al., 2020 (Langmuir)
  #Zabsorbn = 8.93, # Reference Zhang et al., 2020 (Freundlich)
  #ZabsorbKF = 5.35 # Reference Zhang et al., 2020 (Freundlich)
  
)

#### Setting Dynamic Model ####
DataMod <- function(time, stocks, parms){
  with(as.list(c(stocks, parms)),{ 
    # Auxiliary Variables
    SysOprD = ifelse(((time %/% 60)%/%24)%%7 < 4, 1, 0)
    SysOprH = ifelse(time %/% 60 %% 24 >= Sysstart & time %/% 60 %% 24 < Sysdown, 1, 0)
    BackwashR = ifelse(time %% 59 == 0, 1, 0)
    NCSb1Capa = NCSCapa - Zmass / ZDens # L
    NCSb2Capa = NCSCapa - Zmass / ZDens # L
    NCSb3Capa = NCSCapa - GDens / Gmass # L
    
    NCS.N.Capa = ZNCapa * Zmass # g of TN
    C.TN.in = C.TN.in(time) / 1000 # mg to g of TN
    NCS.B1.N.Avaib = 1 - NCSb1.N / NCS.N.Capa
    B1.Z.absorb.R = NCS.N.Capa * NCS.B1.N.Avaib * ZabsorbKL * C.TN.in / (1 + ZabsorbKL * C.TN.in) # Langmuir model 
    #B1.Z.absorb.R = ZabsorbKF * C.TN.in ^ (1 / Zabsorbn) # Freundlich model 
    # Flows
    QoutEQ2 = ifelse(BackwashR == 1, 0, (LMH * 15) * SysOprH * SysOprD)
    QoutNCSb1 = ifelse(NCSb1 >= NCSb1Capa, NCSb1 - NCSb1Capa , 0)
    QoutNCSb2 = ifelse(NCSb2 >= NCSb2Capa, NCSb2 - NCSb2Capa , 0)
    QoutNCSb3 = ifelse(NCSb3 >= NCSb3Capa, NCSb3 - NCSb3Capa , 0)
    # TN Mass
    QC.TN.in.NCB1 = C.TN.in * QoutEQ2
    QC.TN.in.NCB2 = ifelse(B1.Z.absorb.R >= QC.TN.in.NCB1 | QoutNCSb1 == 0, 0, QC.TN.in.NCB1 - B1.Z.absorb.R) # = B1.Z.absorb.R
    C.TN.in.NCB2 = ifelse(SysOprH * SysOprD  == 0 | QoutNCSb1 == 0, 0, QC.TN.in.NCB2 / QoutNCSb1)
    
    NCS.B2.N.Avaib = 1 - NCSb2.N / NCS.N.Capa
    B2.Z.absorb.R = NCS.N.Capa * NCS.B2.N.Avaib * ZabsorbKL * C.TN.in.NCB2 / (1 + ZabsorbKL * C.TN.in.NCB2) # Langmuir model 
    #B2.Z.absorb.R = ZabsorbKF * C.TN.in.NCB2 ^ (1 / Zabsorbn) # Freundlich model 
    
    QC.TN.in.NCB3 = ifelse(B2.Z.absorb.R >= QC.TN.in.NCB2 | QoutNCSb2 == 0, 0, QC.TN.in.NCB2 - B2.Z.absorb.R)
    C.TN.out = ifelse(SysOprH * SysOprD  == 0 | QoutEQ2 == 0 | QoutNCSb1 == 0 | QoutNCSb2 == 0, NA, QC.TN.in.NCB3 / QoutNCSb2)
    # Stocks
    dNCSb1_dt = QoutEQ2 - QoutNCSb1
    dNCSb2_dt = QoutNCSb1 - QoutNCSb2
    dNCSb3_dt = QoutNCSb2 - QoutNCSb3
    dEQ3_dt = QoutNCSb3
    dNCSb1.N_dt = QC.TN.in.NCB1 - QC.TN.in.NCB2
    dNCSb2.N_dt = QC.TN.in.NCB2 - QC.TN.in.NCB3
    dstocks = c(dNCSb1_dt,dNCSb2_dt,dNCSb3_dt,dEQ3_dt,dNCSb1.N_dt,dNCSb2.N_dt)
    # Output
    return(list(dstocks, SysOprD = SysOprD, SysOprH = SysOprH, BackwashR = BackwashR,
                QC.TN.in.NCB1 = QC.TN.in.NCB1, QC.TN.in.NCB2 = QC.TN.in.NCB2, QC.TN.in.NCB3 = QC.TN.in.NCB3, # Mass
                C.TN.in = C.TN.in, C.TN.in.NCB2 = C.TN.in.NCB2, C.TN.out = C.TN.out, # Conc.
                NCS.B1.N.Avaib = NCS.B1.N.Avaib, NCS.B2.N.Avaib = NCS.B2.N.Avaib,
                B1.Z.absorb.R = B1.Z.absorb.R, B2.Z.absorb.R = B2.Z.absorb.R,
                QoutEQ2 = QoutEQ2, QoutNCSb1 = QoutNCSb1, QoutNCSb2 = QoutNCSb2, QoutNCSb3 = QoutNCSb3 # Flow
    ))
  })
}

#### Setting Run ####
day <- 300
t0 = 1
tf = 60*24*day  #728640
tstep = 1
trange <- seq(t0,tf,tstep)

ModOut <- data.frame(ode(y = InitStocks, times = trange, func = DataMod, parms = Parms, method = 'euler'))
max(ModOut$EQ3)/ (day%/%7*5 + ifelse(day%%7 > 5, 5,day%%7))

#### Calibrate Model ####
CalibParms <- c(
  Sysstart = 10,
  Sysdown = 14,
  LMH = 5 / 60, 
  NCSCapa = 350, NCSmax = 360, 
  ZDens = 1.6, # Zeolite density (kg/L)
  Zmass = 175, # Zeolite mass in one bed (kg)
  GDens = 1, # GAC density (kg/L)
  Gmass = 80, # GAC mass in one bed (kg)
  ZNCapa = 20.5, #18.95,
  ZabsorbKL = 0.0007
)

CalibOutData <- data.frame(ode(y = InitStocks, times = trange, func = DataMod, parms = CalibParms, method = 'euler'))


#### Input Actual Data to Calibrate Parameters ####
actualData <- read.csv(file = "G:/My Drive/R project/GitHub/Dynamic_Model/Data/TN_postNCS data.csv",header = T, sep = ",")
actualData <- actualData[,-(1),drop=FALSE] 
colnames(actualData) <- c('time','C.TN.out')
longactualData <- melt(actualData, id.vars = 'time')

NCSeff <- subset(ModOut, select = c('time','C.TN.out'))
NCSeff$C.TN.out <- NCSeff$C.TN.out*1000
NCSeffCal <- subset(CalibOutData, select = c('time','C.TN.out'))
NCSeffCal$C.TN.out <- NCSeffCal$C.TN.out*1000
NCSeff <- NCSeff %>%
  merge(actualData, all.x = TRUE, by = "time")%>%
  merge(NCSeffCal, all.x = TRUE, by = "time")
colnames(NCSeff) <- c("time","Model","Actual","Model(Calibrated)")
NCSeff <- aggregate(NCSeff, list(rep(1:(nrow(NCSeff) %/% 60*24 + 1), each = 60*24, len = nrow(NCSeff))), mean, na.rm = TRUE)[-1]
NCSeffm<- melt(NCSeff,id.vars="time")

NCSB.N.eff2 <- ggplot(NCSeff,aes(`Model(Calibrated)`,Actual)) + geom_point(size=1)+
  geom_smooth(method = "lm", se=FALSE, color="red", linetype = "dashed", formula = y ~ x)+
  stat_poly_eq(use_label(c("R2")),rr.digits = 3) +
  annotate("text", x = 2, y = 230, label = "qm = 22.5 g/kg\nKL = 0.0164 L/mg\nkc = 0.0007 g/kg·min", hjust = 0) +
  labs(x = "Model-Calibrated (mg N/L)", y = "Actual (mg N/L)")+
  theme_bw()+theme(legend.position = "bottom",legend.title = element_blank())
NCSB.N.eff2
ggsave("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-TNeff2.jpg",plot=NCSB.N.eff2,width = 6, height = 4, dpi = 600)

NCSB.N.eff <- ggplot(subset(NCSeffm[NCSeffm$variable %in% c("Model(Calibrated)","Actual"),]),
                     aes(time/(60*24),value, color = variable)) + geom_point(size=1)+
  scale_x_continuous(breaks = seq(0, 400, by = 50), name = "time (d)") + labs(y = "Post-NCS (NH4-N mg/L)")+
  scale_color_manual(values=c("Model" = '#ade600',"Model(Calibrated)" = '#e69200',"Actual" = '#999999')) + 
  theme_bw()+theme(legend.key = element_rect(colour = "transparent", fill = "white"),
                   legend.background = element_rect(fill='transparent'),
                   legend.position = c(0.3,0.8),legend.title = element_blank())
NCSB.N.eff
ggsave("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-TNeff.jpg",plot=NCSB.N.eff,width = 6, height = 4, dpi = 600)

NCSB.N.eff <- NCSB.N.eff + geom_point(aes(group = seq_along(time))) + transition_reveal(time)
gNCSB.N.eff <- animate(NCSB.N.eff, renderer = gifski_renderer(), height = 4, width = 6, units = "in", res = 300)
anim_save("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-TNefft.gif", gNCSB.N.eff)

NCSTN <- subset(CalibOutData, select = c('time','NCSb1.N','NCSb2.N'))
colnames(NCSTN) <- c("time","Bed1","Bed2")
NCSTN <- aggregate(NCSTN, list(rep(1:(nrow(NCSTN) %/% 60*24 + 1), each = 60*24, len = nrow(NCSTN))), mean)[-1]
NCSTN<- melt(NCSTN,id.vars="time")
NCSB.N <- ggplot(NCSTN,aes(time/(60*24),value, color = variable)) + geom_point()+
  scale_x_continuous(breaks = seq(0, 400, by = 50), name = "time (d)") + labs(y = "NCS TN (g)")+
  scale_color_manual(values=c("Bed1" = '#E69F00',"Bed2" = '#999999')) + 
  theme_bw()+theme(legend.position = "bottom",legend.title = element_blank())
NCSB.N
ggsave("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-TN.jpg",plot=NCSB.N,width = 6, height = 4, dpi = 600)

NCSB.N <- NCSB.N + geom_point(aes(group = seq_along(time))) + transition_reveal(time)
gNCSB.N <- animate(NCSB.N, renderer = gifski_renderer(), height = 4, width = 6, units = "in", res = 150)
anim_save("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-TNt.gif", gNCSB.N)



#### weight of actual data to calibrate ####
modCalVars <- unique(longactualData$variable)
for(i in modCalVars) {
  longactualData$weight[longactualData$variable == i] <- mean(longactualData$value
                                                              [longactualData$variable == i],na.rm = T)
}
longactualData = longactualData[,c('variable','time','value','weight')]

CostFunction <- function(p,time,stocks,parms,yactual){
  whichpar <- names(parms)[names(parms) %in% names(p)]
  parms[whichpar] <- p[whichpar]
  ysim <- ode(
    y = stocks, 
    times = time, 
    func = DataMod, 
    parms = parms, 
    method = 'euler')
  MC <- modCost(
    model = ysim,
    obs = yactual,
    x = 'time',
    y = 'value',
    err = 'weight')
  return(MC)
}


#### Calibration factor upper and lower ####
p <- c(Zmass = 175,
       ZNCapa = 16.17,
       ZabsorbKL = 0.002)

pmin <- c(Zmass = 170,
          ZNCapa = 16.00,
          ZabsorbKL = 0.001)

pmax <- c(Zmass = 185,
          ZNCapa = 17.00,
          ZabsorbKL = 0.010)

#### Put upper and lower into function ####
DataModCalibOut <- modFit(
  f = CostFunction,
  p = p,
  lower = pmin,
  upper = pmax,
  time = trange,
  stocks = InitStocks,
  parms = Parms,
  yactual = longactualData,
  method = 'Port'
)


#### Get calibrated parameters ####
CalibParms <- Parms
CalibParms[names(DataModCalibOut$par)] <- DataModCalibOut$par
print(CalibParms)

CalibOutData <- data.frame(ode(y = InitStocks, times = trange, func = DataMod, parms = CalibParms, method = 'euler'))


#### Water Balance Model in Minute Spec ####
InitStocks <- c(
  EQ1 = 1000,
  OFT = 0,
  BioR = 2200,
  EQ2 = 20,
  NCSb1 = 230,
  NCSb2 = 230,
  NCSb3 = 350,
  EQ3 = 50,
  Chlor = 0,
  Prod = 0,
  CDTime = 0, # For Chlorination system cook time calculation
  PRTime = 0,
  MRTime = 0,
  MBWTime = 0
)

Parms <- c(
  Sysstart = 7,
  Sysdown = 19,
  FVolume = 10, # L
  LMH = 8 / 60, 
  PermRunT = 8,
  MemRelaxT = 1,
  MemBWT = 1,
  ChlorcookT = 17, # Time Base Chlorination system cook time (min)
  PipemaxF = 2000 / 60, # L/min
  PFR.R = 2500 / 60, # Pump flow rate_Reactor (L/min)
  PFR.Ch = 1500 / 60, # Pump flow rate_Chlorination (L/hr)
  PFR.BW = 100 / 60, # Pump flow rate_BackWash (L/hr)
  
  EQ1Capa = 2000, EQ1min = 500, # L
  BioRCapa = 2500, BioRmin = 1800, # L
  EQ2Capa = 50, EQ2min = 20, # L
  NCSCapa = 350, NCSmax = 360, 
  EQ3Capa = 150, EQ3min = 20,
  ChlorCapa = 50,
  ProdCapa = 150, Prodmin = 50,
  
  ZDens = 1.6, # Zeolite density (kg/L)
  Zmass = 175, # Zeolite mass in one bed (kg)
  GDens = 1, # GAC density (kg/L)
  Gmass = 80 # GAC mass in one bed (kg)
)


DataMod <- function(time, stocks, parms){
  with(as.list(c(stocks, parms)),{ 
    # Auxiliary Variables
    TUsage = sample(c(1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 1) # Around 125 usage per-day
    SysOprD = ifelse(((time %/% 60)%/%24)%%7 < 5, 1, 0)
    SysOprH = ifelse(time %/% 60 %% 24 >= Sysstart & time %/% 60 %% 24 < Sysdown, 1, 0)
    BackwashR = ifelse(time %% 59 == 0, 1, 0) # time %% 30 == 0
    
    PermRun = ifelse(SysOprH == 1 & SysOprD == 1 & PRTime < time, time + PermRunT, ifelse(PRTime == time, - time, 0))
    MRelax = ifelse(SysOprH == 1 & SysOprD == 1 & PRTime == 0, 1, ifelse(MRTime == 6, - 6, 0))
    BWRun = ifelse(SysOprH == 1 & SysOprD == 1 & MBWTime < time, time + (PermRunT + MemRelaxT) * 3, ifelse(MBWTime == time, - time, 0))
    
    NCSb1Capa = NCSCapa - Zmass / ZDens 
    NCSb2Capa = NCSCapa - Zmass / ZDens
    NCSb3Capa = NCSCapa - GDens / Gmass
    
    # Flows
    QinEQ1 = TUsage * FVolume
    QoutEQ1 = ifelse(EQ1 > EQ1min & BioR + PFR.R <= BioRCapa, PFR.R * SysOprH * SysOprD, 0)
    EQ1OF = ifelse(EQ1 > EQ1Capa, ifelse(EQ1 > EQ1Capa, EQ1 - EQ1Capa , 0), 0)
    
    QoutBioR = ifelse(BioR > BioRmin & EQ2 + (LMH * 15) <= EQ2Capa & BackwashR == 0, (LMH * 15) * SysOprH * SysOprD, 0)
    
    QoutEQ2 = ifelse(EQ2 > EQ2min & EQ3 <= EQ3Capa, ifelse(EQ2 > EQ2min, EQ2 - EQ2min , 0) * SysOprH * SysOprD, 0)
    EQ2OF = ifelse(EQ2 > EQ2Capa & QoutEQ2 < 0, ifelse(EQ2 > EQ2Capa, EQ2 - EQ2Capa , 0), 0)
    
    QoutNCSb1 = ifelse(NCSb1 >= NCSb1Capa, NCSb1 - NCSb1Capa , 0)
    QoutNCSb2 = ifelse(NCSb2 >= NCSb2Capa, NCSb2 - NCSb2Capa , 0)
    QoutNCSb3 = ifelse(NCSb3 >= NCSb3Capa, NCSb3 - NCSb3Capa , 0)
    
    QoutEQ3 = ifelse(EQ3 > EQ3min & Chlor < ChlorCapa & EQ3 > ChlorCapa, PFR.Ch * SysOprH * SysOprD, 0)
    EQ3OF = ifelse(EQ3 > EQ3Capa, ifelse(EQ3 > EQ3Capa, EQ3 - EQ3Capa , 0), 0)
    
    ChlorT = ifelse(Chlor == ChlorCapa & CDTime < time, time + ChlorcookT, ifelse(CDTime == time, - time, 0))
    
    QoutChlor = ifelse (ChlorT < 0 & Chlor == ChlorCapa, ChlorCapa * SysOprH * SysOprD, 0)
    
    QBackwash = ifelse(Prod > Prodmin, PFR.BW * SysOprH * SysOprD * BackwashR, 0)
    
    # Stocks
    dEQ1_dt = QinEQ1 - QoutEQ1 - EQ1OF
    dOFT_dt = EQ1OF + EQ3OF
    dBioR_dt = QoutEQ1 - QoutBioR + QBackwash
    dEQ2_dt = QoutBioR - QoutEQ2 - EQ2OF
    dNCSb1_dt = QoutEQ2 - QoutNCSb1 + EQ2OF
    dNCSb2_dt = QoutNCSb1 - QoutNCSb2
    dNCSb3_dt = QoutNCSb2 - QoutNCSb3
    dEQ3_dt = QoutNCSb3 - QoutEQ3 - EQ3OF
    dChlor_dt = QoutEQ3 - QoutChlor
    dProd_dt = QoutChlor - QBackwash
    
    dCDTime_dt = ChlorT
    dPRTime_dt = PermRun
    dMRTime_dt = MRelax
    dMBWTime_dt = BWRun
    
    dstocks = c(dEQ1_dt, dOFT_dt,dBioR_dt,dEQ2_dt,dNCSb1_dt,dNCSb2_dt,dNCSb3_dt,dEQ3_dt,dChlor_dt,dProd_dt,dCDTime_dt,dPRTime_dt,dMRTime_dt,dMBWTime_dt)
    
    # Output
    return(list(dstocks, SysOprD = SysOprD, SysOprH = SysOprH, MRelax = MRelax, BackwashR = BackwashR, TUsage = TUsage,
                QinEQ1 = QinEQ1, QoutEQ1 = QoutEQ1,QoutBioR = QoutBioR, QoutEQ2 = QoutEQ2,
                QoutNCSb1 = QoutNCSb1, QoutNCSb2 = QoutNCSb2, QoutNCSb3 = QoutNCSb3,
                QoutEQ3 = QoutEQ3, EQ3OF = EQ3OF, QoutChlor = QoutChlor, QBackwash = QBackwash
    ))
  })
}

day <- 2
t0 = 1; tf = 60*24*day; tstep = 1
trange = seq(t0,tf,tstep)

ModOut = data.frame(ode(y = InitStocks, times = trange, func = DataMod, parms = Parms, method = 'euler'))

ModOut_sumH <- aggregate(ModOut, list(rep(1:(nrow(ModOut) %/% 60 + 1), each = 60, len = nrow(ModOut))), sum)[-1]

sum(ModOut$TUsage)/day
max(ModOut$Prod)/ (day/7*5)

ggpl <- function(...){list(geom_point(size=0.75,shape=16,...), scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)"))}
gthem <-list(theme_bw(),theme(text = element_text(size=15),legend.title=element_blank()))


EQ1p <- ggplot() + ggpl(data = ModOut, aes(time/60,EQ1),color = "#ae6300")+labs(y = "Underground EQ1 (L)")+ylim(0,NA)+gthem
BioRp <- ggplot() + ggpl(data = ModOut, aes(time/60,BioR),color = "#40ae00")+labs(y = "BioReactor (L)")+gthem
EQ2p <- ggplot() + ggpl(data = ModOut, aes(time/60,EQ2),color = "#00ae85")+labs(y = "EQ2 (L)")+ylim(0,NA)+gthem
EQ3p <- ggplot() + ggpl(data = ModOut, aes(time/60,EQ3),color = "#00aeae")+labs(y = "EQ3 (L)")+ylim(0,NA)+gthem
Prodp <- ggplot() + ggpl(data = ModOut, aes(time/60,Prod),color = "#0071ae")+labs(y = "Product water (L)")+
  geom_label(aes(x=24, y=7500,label=paste("Average\n",round(max(ModOut$Prod)/ (2),0),"L/d")),color='#3882ff')+gthem
OFp <- ggplot() + ggpl(data = ModOut, aes(time/60,OFT),color = "#aea500")+labs(y = "Overflow water (L)")+gthem
Chlorp <- ggplot() + geom_line(data = ModOut, aes(time/60,Chlor),color = "#0057ae",size=0.1)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "Chlorination Tank (L)")+gthem
NCSb12p <- ggplot()+geom_point(data = ModOut, aes(time/60,NCSb1,color = "NCSb1"),size=0.5)+
  geom_point(data = ModOut, aes(time/60,NCSb2,color = "NCSb2"),size=0.5)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "NCS Bed (L)")+gthem
NCSb3p <-  ggplot()+geom_point(data = ModOut, aes(time/60,NCSb3,color = "NCSb3"),size=0.5)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "NCS Bed (L)")+gthem

QoutBioRp <- ggplot()+geom_line(data = ModOut, aes(time/60,QoutBioR),color = "Black",size=0.5)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "Q outMBR (L/min)")+gthem

backwp <- ggplot()+geom_line(data = ModOut, aes(time/60,QBackwash),color = "Black",size=0.5)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "Q Backwash (L/min)")+gthem

OprTp <- ggplot()+geom_line(data = ModOut, aes(time/60,SysOprH),color = "Black",size=0.5)+
  geom_line(data = ModOut, aes(time/60,SysOprD),color = "Red",size=0.5)+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+labs(y = "Operation")+gthem

Tuserp <- ggplot()+geom_line(data = ModOut_sumH, aes(1:48,TUsage),color = "Black",size=0.5)+
  scale_y_continuous(breaks = seq(0, 14, by = 2), name = "Toilet usage (people)")+
  scale_x_continuous(breaks = seq(0, 1680, by = 24), name = "time (hour)")+gthem


Systemp <- ggarrange(Tuserp,backwp,QoutBioRp,EQ1p,BioRp,EQ2p,EQ3p,Chlorp,Prodp, ncol=3, nrow = 3, align = "v")
Systemp
#ggsave("E:/Dropbox/R/Dynamic Model/Figure/Dynamic model-2d(8LMH).jpg",plot=Systemp,width = 12, height = 9, dpi = 600)


TUsageF <- function(t) {
  Usage = ifelse(t %% 24 == 1, sample(c(rbinom(n = 60, size = 1, prob = 6/60)), 1),ifelse(
    t %% 24 == 2, sample(c(rbinom(n = 60, size = 1, prob = 2/60)), 1),ifelse(
      t %% 24 == 3, sample(c(rbinom(n = 60, size = 1, prob = 1/60)), 1),ifelse(
        t %% 24 == 4, sample(c(rbinom(n = 60, size = 1, prob = 3/60)), 1),ifelse(
          t %% 24 == 5, sample(c(rbinom(n = 60, size = 1, prob = 6/60)), 1),ifelse(
            t %% 24 == 6, sample(c(rbinom(n = 60, size = 1, prob = 9/60)), 1),ifelse(
              t %% 24 == 7, sample(c(rbinom(n = 60, size = 1, prob = 13/60)), 1),ifelse(
                t %% 24 == 8, sample(c(rbinom(n = 60, size = 1, prob = 16/60)), 1),ifelse(
                  t %% 24 == 9, sample(c(rbinom(n = 60, size = 1, prob = 19/60)), 1),ifelse(
                    t %% 24 == 10, sample(c(rbinom(n = 60, size = 1, prob = 20/60)), 1),ifelse(
                      t %% 24 == 11, sample(c(rbinom(n = 60, size = 1, prob = 19/60)), 1),ifelse(
                        t %% 24 == 12, sample(c(rbinom(n = 60, size = 1, prob = 17/60)), 1),ifelse(
                          t %% 24 == 13, sample(c(rbinom(n = 60, size = 1, prob = 14/60)), 1),ifelse(
                            t %% 24 == 14, sample(c(rbinom(n = 60, size = 1, prob = 11/60)), 1),ifelse(
                              t %% 24 == 15, sample(c(rbinom(n = 60, size = 1, prob = 8/60)), 1),ifelse(
                                t %% 24 == 16, sample(c(rbinom(n = 60, size = 1, prob = 7/60)), 1),ifelse(
                                  t %% 24 == 17, sample(c(rbinom(n = 60, size = 1, prob = 8/60)), 1),ifelse(
                                    t %% 24 == 18, sample(c(rbinom(n = 60, size = 1, prob = 9/60)), 1),ifelse(
                                      t %% 24 == 19, sample(c(rbinom(n = 60, size = 1, prob = 11/60)), 1),ifelse(
                                        t %% 24 == 20, sample(c(rbinom(n = 60, size = 1, prob = 13/60)), 1),ifelse(
                                          t %% 24 == 21, sample(c(rbinom(n = 60, size = 1, prob = 15/60)), 1),ifelse(
                                            t %% 24 == 22, sample(c(rbinom(n = 60, size = 1, prob = 15/60)), 1),ifelse(
                                              t %% 24 == 23, sample(c(rbinom(n = 60, size = 1, prob = 14/60)), 1),ifelse(
                                                t %% 24 == 0, sample(c(rbinom(n = 60, size = 1, prob = 13/60)), 1)))))))))))))))))))))))))
  print(Usage)
}