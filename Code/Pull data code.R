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

Site1_Master_New_I <-  Site1_Master %>% filter(Group == "S1")

#Site1_Master_New_I$Time <- Site1_Master_New_I$Date-Site1_Master_New_I$Date[1]+1

TN.data <- subset(Site1_Master_New_I, select = c('Date','TS_Avg','VS_Avg','FS_Avg','TSS_Avg','COD_Avg','SCOD_Avg', 'pH_Avg', 'TP_Avg','TN_Avg','NH4_Avg','E.co_Avg'))
