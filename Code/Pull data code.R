pacman::p_load(plyr,dplyr,reshape2,stringr,lubridate,readxl,openxlsx,janitor,ggplot2,ggpmisc,ggpubr,ggsignif,gplots,data.table,magrittr,
               grid,tidyverse,rstatix,PtProcess,gtable,gridExtra,cowplot,png,tiff,pwr,MOTE,kableExtra,tinytex,
               knitr,quantmod,class,caret,gmodels,officedown,officer,flextable,car,gtools,ggrepel,lemon,rvg,gdata)
setwd("G:/My Drive/R project")

Site1_Master <- read.csv(file = "G:/My Drive/R project/All site data/Lab Master.csv",header = T, sep = ",")
Na.num<- c("Group","Date")
Site1_Master[!names(Site1_Master)%in%Na.num] <- sapply(Site1_Master[!names(Site1_Master)%in%Na.num], as.numeric )
Site1_Master$Date <- mdy(Site1_Master$Date)
Site1_Master[Site1_Master < 0] <- NaN
Site1_Master <- Site1_Master[order(Site1_Master[,3], Site1_Master[,2]), ]

#Site1_Master_New <- Site1_Master[Site1_Master$Date <= ymd("2021-2-15"), ]
Site1_Master_New <- Site1_Master[Site1_Master$Date <= ymd("2019-03-12"), ]

Site1_Master_New_P <-  Site1_Master_New %>% filter(Group == "S3")
Site1_Master_New_N <-  Site1_Master_New %>% filter(Group == "S4")


COD_P <- subset(Site1_Master_New_P, select = c('Date','COD_Avg','SCOD_Avg'))
COD_N <- subset(Site1_Master_New_N, select = c('Date','COD_Avg','SCOD_Avg'))
COD.data <- left_join(COD_P, COD_N, by = "Date")
colnames(COD.data) <- c('Date','COD.in','sCOD.in','COD.out','sCOD.out')
COD.data$Time <- COD.data$Date-COD.data$Date[1]+1

#write.csv(COD.data, file = "G:/My Drive/R project/GitHub/Dynamic_Model/Data/COD data.csv", na="NaN")
COD.data <- read.csv(file = "G:/My Drive/R project/GitHub/Dynamic_Model/Data/COD data.csv",header = T, sep = ",")

NCS.C.eff <- ggplot(COD.data,aes(COD.out.M,COD.out)) + geom_point(size=1)+
  geom_smooth(method = "lm", se=FALSE, color="red", linetype = "dashed", formula = y ~ x)+
  stat_poly_eq(use_label(c("R2")),rr.digits = 3) +
  annotate("text", x = -8, y = 295, label = "qm = 83 g/kg\nKL = 0.041 L/mg\nkc = 0.0038 g/kg·min", hjust = 0) +
  labs(x = "Model-Calibrated (mg COD/L)", y = "Actual (mg COD/L)")+
  theme_bw()+theme(legend.position = "bottom",legend.title = element_blank())
NCS.C.eff
ggsave("G:/My Drive/R project/GitHub/Dynamic_Model/Figure/Dynamic model-CODeff.jpg",plot=NCS.C.eff,width = 6, height = 4, dpi = 600)