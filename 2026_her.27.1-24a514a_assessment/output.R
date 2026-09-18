#------------------------------------------------------------------------------#
#output.R
#
#Author: 
#Sindre Vatnehol
#Institute of marine Research
#Norway 
#
#Description:
#Here we produce all the figures, tables, and presentation that is not used in
#the WGWIDE report but serve as auxilary information
#
#Chapters:
#  1. Setup
#  2. Data signal diagnostics
#  3. Parameter diagnostics
#  4. Residual diagnostics
#  5. Leave-out & retrospective diagnostics
#  6. Stock summary comparison
#------------------------------------------------------------------------------#



#================================================================================#
# 1. SETUP
#================================================================================#
rm(list=ls())




#------------------------------------------------------------------------------#
#load TAF
#------------------------------------------------------------------------------#
library(TAF)
library(ggplot2)
library(stockassessment)
library(reshape2)
source(taf.boot.path('rscript/utils.r'))




#------------------------------------------------------------------------------#
#Make the directory tree
#------------------------------------------------------------------------------#
mkdir("output")
mkdir("output/tables")
mkdir("output/figures")




#------------------------------------------------------------------------------#
#Previous assessments
#------------------------------------------------------------------------------#
load('boot/rscript/WGWIDE2025.Rda')
fit_wgwide2025 <- fit          
load('boot/rscript/fit_official_benchmark.Rda')
fit_benchmark <- fit          


#------------------------------------------------------------------------------#
#load current assessment
#------------------------------------------------------------------------------#
load('model/diag.Rda')
load('model/fit.Rda')

#Labels used throughout this script for the current vs. previous assessment's
#case (fit vs. fit_wgwide2025 above). Derived from each fit's terminal data
#year (max(fit$data$years) already reflects the WGWIDE meeting year, since
#the survey fleets - unlike the catch data - extend through it), so these
#update automatically next year without touching this script - only the
#loaded .Rda files (line 57 above) change.
current_case  <- paste0('WGWIDE', max(fit$data$years))
previous_case <- paste0('WGWIDE', max(fit_wgwide2025$data$years))





#================================================================================#
# 2. DATA SIGNAL DIAGNOSTICS
#================================================================================#

#------------------------------------------------------------------------------#
#Process the signal from the data
#First we grab the signal as from the fit object
#------------------------------------------------------------------------------#
out <-as.data.frame(fit$data$aux)
out$fleetName<-attr(fit$data,'fleetNames')[out$fleet]
out$logobs <- fit$data$logobs
out$weight <- fit$data$weight
out$yc <- out$year-out$age





#------------------------------------------------------------------------------#
#grab the input data for the residual catch, and plot it
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('caa.txt')) )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('caa_cv.txt')))
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2

out$fleetName[out$fleetName == 'Residual catch']<-'Catch'
p_catch<-gridExtra::grid.arrange(
  ggplot(out[out$fleet==1 & out$age<12,],aes(x=year,y=logobs,group=yc,colour=as.factor(yc)))+
    geom_line()+theme_bw()+
    ggtitle(paste(unique(out[out$fleet==1,]$fleetName),'signals'))+
    theme(legend.position = "none")+ylab('log catch')
  ,
  ggplot(cn, aes(x=log(cn), y=log(var))) +
    geom_point() +
    geom_point(data=cn[cn$year==(max(cn$year)),],   aes(colour='latest')) +
    geom_point(data=cn[cn$year==(max(cn$year)-1),], aes(colour='previous')) +
    theme_bw() +
    ggtitle(paste(unique(out[out$fleet==1,]$fleetName), 'variance')) +
    ylab('log variance') + xlab('log catch') +
    theme(
      legend.position = c(0.85, 0.15),        # 👈 x, y i 0-1 koordinater
      legend.background = element_rect(fill = "white", colour = "grey80")
    ),
  ncol=2
)
ggsave("output/figures/fig_catchsignal.png", p_catch, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting), one panel per observation-variance
#group for this fleet (see model.R's conf$keyVarObs row 1).
p_predsd_catch <- gridExtra::grid.arrange(
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=1,age=2),
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=1,age=3),
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=1,age=4),
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=1,age=5:12),
  ncol=2
)
ggsave("output/figures/fig_predsd_catch.png", p_predsd_catch, width = 10, height = 8, dpi = 300)




#------------------------------------------------------------------------------#
#Now process the IESNS survey, and plot it
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('survey_T2025-02-05.txt'))$IESNS )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('surveycv_T2025-02-05.txt'))$IESNS)
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2

p_IESNS<-gridExtra::grid.arrange(
  ggplot(out[out$fleet==2 & out$age<12,],aes(x=year,y=logobs,group=yc,colour=as.factor(yc)))+
    geom_line()+theme_bw()+ggtitle(paste(unique(out[out$fleet==2,]$fleetName),'signals'))+
    theme(legend.position = "none")+ylab('log index')
  ,
  ggplot(cn,aes(x=log(cn),y=log(var)))+geom_point()+
    geom_point(data=cn[cn$year==(max(cn$year)),],   aes(colour='latest')) +
    geom_point(data=cn[cn$year==(max(cn$year)-1),], aes(colour='previous')) +
    ggtitle(paste(unique(out[out$fleet==2,]$fleetName),'variance'))+
    theme_bw()+
    ylab('log variance')+xlab('log index')+
    theme(
      legend.position = c(0.85, 0.15),        # 👈 x, y i 0-1 koordinater
      legend.background = element_rect(fill = "white", colour = "grey80")
    ),
  ncol=2
)

ggsave("output/figures/fig_IESNSsignal.png", p_IESNS, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting), one panel per observation-variance
#group for this fleet (see model.R's conf$keyVarObs row 2).
p_predsd_IESNS <- gridExtra::grid.arrange(
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=2,age=3),
  plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=2,age=4:12),
  ncol=2
)
ggsave("output/figures/fig_predsd_IESNS.png", p_predsd_IESNS, width = 8, height = 4, dpi = 300)




#------------------------------------------------------------------------------#
#prepare the figure for BESS survey
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('survey_T2025-02-05.txt'))$recBar )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('surveycv_T2025-02-05.txt'))$recBar)
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2
cn <- cn[cn$age%in%c(2,3),]


out$fleetName[out$fleetName == 'recBar']<-'BESS'
p_BESS<-gridExtra::grid.arrange(
  ggplot(out[out$fleet==3 & out$age<12,],aes(x=year,y=logobs,group=age,colour=as.factor(age)))+
    geom_line()+theme_bw()+ggtitle('BESS signals')+
    labs(colour='Age')+ylab('log index')+geom_point()+
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", colour = "grey80")
    )
  ,
  ggplot(cn,aes(x=log(cn),y=log(var)))+geom_point()+
    geom_point(data=cn[cn$year==(max(cn$year)),],   aes(colour='latest')) +
    geom_point(data=cn[cn$year==(max(cn$year)-1),], aes(colour='previous')) +
    ggtitle('BESS variance')+
    theme_bw()+
    ylab('log variance')+xlab('log index')+
    theme(
      legend.position = c(0.85, 0.15),        # 👈 x, y i 0-1 koordinater
      legend.background = element_rect(fill = "white", colour = "grey80")
    )
  ,
  # ggplot(out[out$fleet==3 & !is.na(out$logobs),],aes(x=year,y=age,size=weight,colour=weight))+geom_point()+
  #   theme_bw()+ggtitle('BESS data weight')+
  #   theme(legend.position = "none")+
  #   scale_colour_gradient(low = "black", high = "red")
  # # ,
  ncol=2
)

ggsave("output/figures/fig_BESSsignal.png", p_BESS, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting). BESS has a single observation-
#variance group covering ages 2-3 (model.R's conf$keyVarObs row 3).
p_predsd_BESS <- plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=3,age=2:3)
ggsave("output/figures/fig_predsd_BESS.png", p_predsd_BESS, width = 5, height = 4, dpi = 300)



#------------------------------------------------------------------------------#
#prepare the plot for the IESNS recruitment
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('survey_T2025-02-05.txt'))$Bar )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('surveycv_T2025-02-05.txt'))$Bar)
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2


p_IESNSrec<-gridExtra::grid.arrange(
  ggplot(out[out$fleet==4 & out$age<12,],aes(x=year,y=logobs,group=age),colour='black')+
    geom_line()+theme_bw()+ggtitle('IESNS,bar signal')+
    theme(legend.position = "none")+ylab('log index')+geom_point()
  ,
  ggplot(cn,aes(x=log(cn),y=log(var)))+geom_point()+
    geom_point(data=cn[cn$year==(max(cn$year)),],   aes(colour='latest')) +
    geom_point(data=cn[cn$year==(max(cn$year)-1),], aes(colour='previous')) +
    ggtitle('IESNS,bar variance')+ylab('log variance')+xlab('log index')+
    theme_bw()+
    theme(
      legend.position = c(0.85, 0.15),
      legend.background = element_rect(fill = "white", colour = "grey80")
    )
  ,
  ncol=2
)

ggsave("output/figures/fig_IESNSrecsignal.png", p_IESNSrec, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting). IESNS_rec (Bar) has a single
#observation-variance group, age 2 only (model.R's conf$keyVarObs row 4).
p_predsd_IESNSrec <- plot_predsd_extweight(fit,cn$cn,cn$cv,fleet=4,age=2)
ggsave("output/figures/fig_predsd_IESNSrec.png", p_predsd_IESNSrec, width = 5, height = 4, dpi = 300)






#------------------------------------------------------------------------------#
#prepare the figure for the NASF old
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('survey_T2025-02-05.txt'))$NASF )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('surveycv_T2025-02-05.txt'))$NASF)
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2

p_NASFold <- gridExtra::grid.arrange(
  ggplot(out[out$fleet==5 & out$age<12,],aes(x=year,y=logobs,group=yc,colour=as.factor(yc)))+
    geom_line()+theme_bw()+ggtitle(paste(unique(out[out$fleet==5,]$fleetName),'signals'))+
    theme(legend.position = "none")+ylab('log index')+geom_point()
  ,
  ggplot(cn[cn$year<2015,],aes(x=log(cn),y=log(var)))+geom_point()+
    theme_bw()+
    ggtitle(paste(unique(out[out$fleet==5,]$fleetName),'variance'))+ylab('log variance')+xlab('log index')
  ,
  ncol=2
)
ggsave("output/figures/fig_NASFoldsignal.png", p_NASFold, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting), one panel per observation-variance
#group for this fleet (model.R's conf$keyVarObs row 5).
cn_old <- cn[cn$year<2015,]
p_predsd_NASFold <- gridExtra::grid.arrange(
  plot_predsd_extweight(fit,cn_old$cn,cn_old$cv,fleet=5,age=3),
  plot_predsd_extweight(fit,cn_old$cn,cn_old$cv,fleet=5,age=4),
  plot_predsd_extweight(fit,cn_old$cn,cn_old$cv,fleet=5,age=5:12),
  ncol=3
)
ggsave("output/figures/fig_predsd_NASFold.png", p_predsd_NASFold, width = 12, height = 4, dpi = 300)



#------------------------------------------------------------------------------#
#Prepare the plot for the NASF new
#------------------------------------------------------------------------------#
p_NASFnew <- gridExtra::grid.arrange(
  ggplot(out[out$fleet==6 & out$age<12,],aes(x=year,y=logobs,group=yc,colour=as.factor(yc)))+
    geom_line()+theme_bw()+ggtitle(paste(unique(out[out$fleet==6,]$fleetName),'signals'))+
    theme(legend.position = "none")+ylab('log index')
  ,
  ggplot(cn[cn$year>=2015,],aes(x=log(cn),y=log(var)))+geom_point()+
    geom_point(data=cn[cn$year>=max(cn$year),],aes(x=log(cn),y=log(var),colour='latest'))+
    geom_point(data=cn[cn$year==max(cn$year)-1,],aes(x=log(cn),y=log(var),colour='previous'))+
    theme_bw()+
    ggtitle(paste(unique(out[out$fleet==6,]$fleetName),'variance'))+ylab('log variance')+xlab('log index')+
    theme(
      legend.position = c(0.85, 0.15),        # 👈 x, y i 0-1 koordinater
      legend.background = element_rect(fill = "white", colour = "grey80")
    )
  ,
  ncol=2
)
ggsave("output/figures/fig_NASFnewsignal.png", p_NASFnew, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting), one panel per observation-variance
#group for this fleet (model.R's conf$keyVarObs row 6).
cn_new <- cn[cn$year>=2015,]
p_predsd_NASFnew <- gridExtra::grid.arrange(
  plot_predsd_extweight(fit,cn_new$cn,cn_new$cv,fleet=6,age=3),
  plot_predsd_extweight(fit,cn_new$cn,cn_new$cv,fleet=6,age=4),
  plot_predsd_extweight(fit,cn_new$cn,cn_new$cv,fleet=6,age=5:12),
  ncol=3
)
ggsave("output/figures/fig_predsd_NASFnew.png", p_predsd_NASFnew, width = 12, height = 4, dpi = 300)



#------------------------------------------------------------------------------#
#prepare the figure for the RFID
#------------------------------------------------------------------------------#
cn <- reshape2::melt(read.ices(taf.data.path('rfid_index_BioSampleID.txt'))$RFID )
names(cn)<-c('year','age','cn')
cn_cv <- reshape2::melt(read.ices(taf.data.path('rfid_cv_BioSampleID.txt'))$RFID)
names(cn_cv)<-c('year','age','cv')
cn <- plyr::join(cn,cn_cv)
cn$var <- (cn$cn*cn$cv)**2

p_RFID <- gridExtra::grid.arrange(
  ggplot(out[out$fleet==7 & out$age<12,],aes(x=year,y=logobs,group=yc,colour=as.factor(yc)))+
    geom_line()+theme_bw()+ggtitle(paste(unique(out[out$fleet==7,]$fleetName),'signal'))+
    theme(legend.position = "none")+ylab('log index')
  ,
  ggplot(cn,aes(x=log(cn),y=log(var)))+geom_point()+
    geom_point(data=cn[cn$year==(max(cn$year)),],aes(colour='lastest'))+
    geom_point(data=cn[cn$year==(max(cn$year)-1),],aes(colour='previous'))+
    theme_bw()+
    ggtitle(paste(unique(out[out$fleet==7,]$fleetName),'variance'))+ylab('log variance')+xlab('log index')+
    theme(
      legend.position = c(0.85, 0.15),        # 👈 x, y i 0-1 koordinater
      legend.background = element_rect(fill = "white", colour = "grey80")
    )
  ,
  ncol=2
)
ggsave("output/figures/fig_RFIDsignal.png", p_RFID, width = 12, height = 4, dpi = 300)

#Pred-SD relation (external weighting), one panel per observation-variance
#group for this fleet (model.R's conf$keyVarObs row 7). Age 2 is excluded
#here to match data.R's computeExternalWeight(mu_data=surveys$RFID,...),
#where surveys$RFID already has its first (age 2) column dropped.
cn_rfid <- cn[cn$age!=2,]
p_predsd_RFID <- gridExtra::grid.arrange(
  plot_predsd_extweight(fit,cn_rfid$cn,cn_rfid$cv,fleet=7,age=3:9),
  plot_predsd_extweight(fit,cn_rfid$cn,cn_rfid$cv,fleet=7,age=10:11),
  plot_predsd_extweight(fit,cn_rfid$cn,cn_rfid$cv,fleet=7,age=12),
  ncol=3
)
ggsave("output/figures/fig_predsd_RFID.png", p_predsd_RFID, width = 12, height = 4, dpi = 300)


#================================================================================#
# 3. PARAMETER DIAGNOSTICS
#================================================================================#

#------------------------------------------------------------------------------#
#prepare the figure for the cathability parameters
#This also shows the historical runs
#------------------------------------------------------------------------------#
dd <- as.data.frame(partable(fit))
dd$parname <- rownames(dd)
dd$case <- current_case
dd_out <- dd

dd <- as.data.frame(partable(fit_benchmark))
dd$parname <- rownames(dd)
dd$case <- 'benchmark'
dd_out <- rbind(dd_out,dd)

dd <- as.data.frame(partable(fit_wgwide2025))
dd$parname <- rownames(dd)
dd$case <- previous_case
dd_out <- rbind(dd_out,dd)



dd_out$facet <- 'rest'
dd_out$facet[grepl('logFpar',dd_out$parname)]<-'Catchabilty'


dd_out$parname[dd_out$parname=='logFpar_0']<-'IESNS age03'
dd_out$parname[dd_out$parname=='logFpar_1']<-'IESNS age04'
dd_out$parname[dd_out$parname=='logFpar_2']<-'IESNS age05'
dd_out$parname[dd_out$parname=='logFpar_3']<-'IESNS age06'
dd_out$parname[dd_out$parname=='logFpar_4']<-'IESNS age07'
dd_out$parname[dd_out$parname=='logFpar_5']<-'IESNS age08'
dd_out$parname[dd_out$parname=='logFpar_6']<-'IESNS age09'
dd_out$parname[dd_out$parname=='logFpar_7']<-'IESNS age10'
dd_out$parname[dd_out$parname=='logFpar_8']<-'IESNS age11-12'
dd_out$parname[dd_out$parname=='logFpar_9']<-'BESS age02'
dd_out$parname[dd_out$parname=='logFpar_10']<-'BESS age03'
dd_out$parname[dd_out$parname=='logFpar_11']<-'IESNS_rec age02'
dd_out$parname[dd_out$parname=='logFpar_12']<-'NASF_old age03'
dd_out$parname[dd_out$parname=='logFpar_13']<-'NASF_old age04'
dd_out$parname[dd_out$parname=='logFpar_14']<-'NASF_old age05'
dd_out$parname[dd_out$parname=='logFpar_15']<-'NASF_old age06'
dd_out$parname[dd_out$parname=='logFpar_16']<-'NASF_old age07'
dd_out$parname[dd_out$parname=='logFpar_17']<-'NASF_old age08-12'
dd_out$parname[dd_out$parname=='logFpar_18']<-'NASF_new age03'
dd_out$parname[dd_out$parname=='logFpar_19']<-'NASF_new age04'
dd_out$parname[dd_out$parname=='logFpar_20']<-'NASF_new age05'
dd_out$parname[dd_out$parname=='logFpar_21']<-'NASF_new age06'
dd_out$parname[dd_out$parname=='logFpar_22']<-'NASF_new age07'
dd_out$parname[dd_out$parname=='logFpar_23']<-'NASF_new age08'
dd_out$parname[dd_out$parname=='logFpar_24']<-'NASF_new age09'
dd_out$parname[dd_out$parname=='logFpar_25']<-'NASF_new age10'
dd_out$parname[dd_out$parname=='logFpar_26']<-'NASF_new age11-12'
dd_out$parname[dd_out$parname=='logFpar_27']<-'RFID age3-12'
dd_out$fleet <- sub(" .*", "", dd_out$parname)

plot_catchability <- ggplot(dd_out[dd_out$facet == 'Catchabilty' & dd_out$case != 'benchmark', ],
                            aes(x = parname, y = par, fill = case)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = par - `sd(par)`, ymax = par + `sd(par)`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Parameter", y = "Estimate", fill = "Case") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_catchability

ggsave("output/figures/fig_catchability.png", plot_catchability, width = 12, height = 8, dpi = 600)




#------------------------------------------------------------------------------#
#observation variance
#------------------------------------------------------------------------------#
dd_out$facet[grepl('logSdLogObs',dd_out$parname)]<-'ObservationVariance'

dd_out$parname[dd_out$parname=='logSdLogObs_0']<-'Catch age02'
dd_out$parname[dd_out$parname=='logSdLogObs_1']<-'Catch age03'
dd_out$parname[dd_out$parname=='logSdLogObs_2']<-'Catch age04'
dd_out$parname[dd_out$parname=='logSdLogObs_3']<-'Catch age05-12'

dd_out$parname[dd_out$parname=='logSdLogObs_4']<-'IESNS age03'
dd_out$parname[dd_out$parname=='logSdLogObs_5']<-'IESNS age04-12'

dd_out$parname[dd_out$parname=='logSdLogObs_6']<-'BESS age02-3'
dd_out$parname[dd_out$parname=='logSdLogObs_7']<-'IESNS_rec age02'

dd_out$parname[dd_out$parname=='logSdLogObs_8']<-'NASF_old age03'
dd_out$parname[dd_out$parname=='logSdLogObs_9']<-'NASF_old age04'
dd_out$parname[dd_out$parname=='logSdLogObs_10']<-'NASF_old age05-12'

dd_out$parname[dd_out$parname=='logSdLogObs_11']<-'NASF_new age03'
dd_out$parname[dd_out$parname=='logSdLogObs_12']<-'NASF_new age04'
dd_out$parname[dd_out$parname=='logSdLogObs_13']<-'NASF_new age05-12'
dd_out$parname[dd_out$parname=='logSdLogObs_14']<-'RFID age03-9'
dd_out$parname[dd_out$parname=='logSdLogObs_15']<-'RFID age10-11'
dd_out$parname[dd_out$parname=='logSdLogObs_16']<-'RFID age12'
dd_out$fleet <- sub(" .*", "", dd_out$parname)

#par/sd(par) are on the log scale (SAM's logSdLogObs); plotted here as the
#natural-scale observation SD (exp(par)) instead, with the error bar taken
#by transforming the log-scale Wald interval (exp(par +/- sd(par))), so the
#bar height reflects the actual magnitude of observation uncertainty rather
#than the log-parameter, which can look disproportionately large/small when
#its value is close to zero (i.e. SD close to 1).
dd_obsvar <- dd_out[dd_out$facet == 'ObservationVariance' & dd_out$case != 'benchmark', ]
dd_obsvar$sd_natural <- exp(dd_obsvar$par)
dd_obsvar$sd_low  <- exp(dd_obsvar$par - dd_obsvar$`sd(par)`)
dd_obsvar$sd_high <- exp(dd_obsvar$par + dd_obsvar$`sd(par)`)

plot_ObservationVariance <- ggplot(dd_obsvar,
                                   aes(x = parname, y = sd_natural, fill = case)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = sd_low, ymax = sd_high),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Parameter", y = "Observation SD (natural scale)", fill = "Case") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_ObservationVariance
ggsave("output/figures/fig_ObservationVariance.png", plot_ObservationVariance, 
       width = 12, height = 10, dpi = 600)



#------------------------------------------------------------------------------#
#and then the rest of the parameters
#------------------------------------------------------------------------------#
plot_rest <- ggplot(dd_out[dd_out$facet == 'rest', ], 
                    aes(x = parname, y = par, fill = case)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = par - `sd(par)`, ymax = par + `sd(par)`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Parameter", y = "Estimate", fill = "Case") +
  theme_bw() 
plot_rest

ggsave("output/figures/fig_rest_config.png", plot_rest, width = 8, height = 6, dpi = 600)

#Table version of dd_out (all parameter estimates + SD, WGWIDE2026 vs.
#benchmark vs. WGWIDE2025), underlying the three plots above. test.R diffs
#.txt files with identical() (no tolerance), so unrounded optimiser output
#could differ in the last decimal between runs (e.g. multi-threaded BLAS)
#and falsely fail the reproducibility check - rounded here to guard against
#that. par/sd(par) (log scale, ~O(1-10)) use round() to 3 digits, matching
#report.R's convention; exp(par)/Low/High (natural scale) can be <<1 (e.g.
#catchability), where round(...,3) would collapse small values to 0, so
#those use signif() (4 significant digits) instead.
dd_out_out <- dd_out
dd_out_out[c('par','sd(par)')] <- round(dd_out_out[c('par','sd(par)')], digits = 3)
dd_out_out[c('exp(par)','Low','High')] <- signif(dd_out_out[c('exp(par)','Low','High')], digits = 4)
write.csv(dd_out_out, file = 'output/tables/parameter_estimates.txt', row.names = FALSE)



#stockassessment::fitplot() only plots one fit at a time (base-R panels),
#which makes it hard to compare this year's fit against last year's run.
#plot.fitcompare() (boot/rscript/utils.r) reproduces the same observed vs.
#fitted extraction, but overlays both fits (coloured by case) per age panel.
#One figure per fleet (fleet numbers/names as set up in model.R); fleet 5
#(NASF_old) is left out here, same as in the case list above.
#------------------------------------------------------------------------------#
fitcompare_fleets <- c('1'='Catch','2'='IESNS','3'='BESS','4'='IESNSrec','6'='NASFnew','7'='RFID')
for(ff in names(fitcompare_fleets)){
  fleet_no <- as.numeric(ff)

  #Number of age panels facet_wrap(~age) will draw for this fleet, used to
  #scale the saved figure size so few-age fleets (e.g. IESNSrec, one panel)
  #don't get an oversized figure and many-age fleets (e.g. RFID) don't get
  #a cramped one. n2mfrow gives the same row/col layout ggplot's facet_wrap
  #falls back to.
  nage <- length(unique(fit$data$aux[fit$data$aux[,"fleet"]==fleet_no,"age"]))
  rc <- grDevices::n2mfrow(nage)

  #Presentation only wants this year's fit (no WGWIDE2025 overlay), plotted
  #in black rather than plot.fitcompare()'s default case-coloured palette.
  p_fitcompare <- plot.fitcompare(list('2026'=fit), fleet = fleet_no) +
    scale_colour_manual(values = c('2026'='black')) +
    theme(legend.position = 'none')
  #BESS (2 age panels) reads cramped at the default width - widen it beyond
  #the usual 3*rc[2]+2 formula. Floored at 8x6 so few-panel fleets (e.g.
  #IESNSrec, one panel) don't get a tiny canvas - theme_bw()'s fixed text/
  #point size then looks oversized relative to such a small plot area.
  fig_width  <- if(fitcompare_fleets[ff]=='BESS') 14 else max(3*rc[2]+2, 8)
  fig_height <- max(2.5*rc[1]+1, 6)
  ggsave(paste0("output/figures/fig_fitcompare_",fitcompare_fleets[ff],".png"), p_fitcompare,
         width = fig_width, height = fig_height, dpi = 600)
}



#================================================================================#
# 4. RESIDUAL DIAGNOSTICS
#================================================================================#

#------------------------------------------------------------------------------#
##RESIDUAL
#------------------------------------------------------------------------------#
res <- res_fit
class(res)<-'data.frame'
res<-res[!is.na(res$observation),]
res<-res[res$year!=1988,]
res$fleetNames <- attr(fit$data,'fleetNames')[res$fleet]

res$fleetNames[res$fleetNames=='Bar']<-'IESNS bar'
res$fleetNames[res$fleetNames=='recBar']<-'BESS'
res$fleetNames[res$fleetNames=='Residual catch']<-'Catch'

#Table version of res (one-step-ahead observation residuals), underlying
#fig_residual.png/fig_residualQQ.png below. All numeric columns rounded to
#3 digits, same reasoning as parameter_estimates.txt above - test.R diffs
#.txt files with identical() (no tolerance).
res_out <- res
num_cols <- sapply(res_out, is.numeric)
res_out[num_cols] <- lapply(res_out[num_cols], round, digits = 3)
write.csv(res_out, file = 'output/tables/residuals.txt', row.names = FALSE)

plot_residual <- ggplot(res, aes(x = year, y = age)) +
  geom_point(aes(size = abs(residual), color = residual > 0), alpha = 0.5) +
  scale_color_manual(
    values = c("TRUE" = "blue", "FALSE" = "red"),
    labels = c("FALSE" = "Negative", "TRUE" = "Positive"),
    name = "Residual sign"
  ) +
  scale_size_continuous(
    range = c(1, 6),   # adjust size range if you want
    breaks = c(0, 4, 6),
    name = "Residual magnitude"
  ) +
  facet_wrap(~fleetNames) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_residual
ggsave("output/figures/fig_residual.png", plot_residual, width = 12, 
       height = 8, dpi = 600)


#------------------------------------------------------------------------------#
##RESIDUAL QQplot
#------------------------------------------------------------------------------#
p<-ggplot(res,aes(sample=residual))+
  geom_qq(distribution = stats::qnorm)+
  geom_qq_line()+
  ggtitle('QQplot by fleet')+ylim(c(-5,5))+ylab('')+xlab('')+facet_wrap(~fleetNames)+
  theme_bw() 
p1<-ggplot(res,aes(sample=residual))+
  geom_qq(distribution = stats::qnorm)+
  geom_qq_line()+
  ggtitle('QQplot all fleet')+ylim(c(-5,5))+ylab('')+xlab('')+
  theme_bw() 
plot_qqres <-gridExtra::grid.arrange(p1,p,ncol=2)
ggsave("output/figures/fig_residualQQ.png", plot_qqres, width = 12, height = 8, dpi = 600)





#------------------------------------------------------------------------------#
#Process residuals
#------------------------------------------------------------------------------#
res <-fit_pro
class(res)<-'data.frame'
res$fleetName <- 'log(N)'
res[res$fleet==1,]$fleetName<-'ResidualCatch'

#Table version of res (process residuals), underlying
#fig_processresidual.png/fig_processresidualQQ.png below.
res_out <- res
num_cols <- sapply(res_out, is.numeric)
res_out[num_cols] <- lapply(res_out[num_cols], round, digits = 3)
write.csv(res_out, file = 'output/tables/process_residuals.txt', row.names = FALSE)

plot_processresiduals <- ggplot(res, aes(x = year, y = age)) +
  geom_point(aes(size = abs(residual), color = residual > 0), alpha = 0.5) +
  scale_color_manual(
    values = c("TRUE" = "blue", "FALSE" = "red"),
    labels = c("FALSE" = "Negative", "TRUE" = "Positive"),
    name = "Residual sign"
  ) +
  scale_size_continuous(
    range = c(1, 6),   # adjust size range if you want
    breaks = c(0, 4, 6),
    name = "Residual magnitude"
  ) +
  facet_wrap(~fleetName) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
plot_processresiduals
ggsave("output/figures/fig_processresidual.png", plot_processresiduals, width = 12, height = 8, dpi = 600)


#------------------------------------------------------------------------------#
#Process residuals
#------------------------------------------------------------------------------#
plot_residualqq<-ggplot(res,aes(sample=residual))+
  geom_qq(distribution = stats::qnorm)+
  geom_qq_line()+
  ggtitle('QQplot of process residual')+ylim(c(-5,5))+ylab('')+xlab('')+facet_wrap(~fleetName)+
  theme_bw() 
plot_residualqq
ggsave("output/figures/fig_processresidualQQ.png", plot_residualqq, width = 12, height = 8, dpi = 600)






#------------------------------------------------------------------------------#
#Compute the process error for the current fit, then for last year's official
#fit and the benchmark fit refitted under this year's configuration (conf),
#so all three are comparable on the same parameterisation.
#------------------------------------------------------------------------------#
pe <- reshape2::melt(get.procres.sdrep(fit))
colnames(pe)<-c('year','age','pe')
pe <- plyr::join(pe , data.frame(age = fit$conf$minAge:fit$conf$maxAge,
                                 par = exp(fit$pl$logSdLogN[fit$conf$keyVarLogN+1]))
)
pe$case <- current_case
N_tab <- reshape2::melt(ntable(fit))
colnames(N_tab)<-c('year','age','N')
pe_N <- plyr::join(pe,N_tab)


#Refit last year's data under this year's SAM configuration, so its process
#error is computed with the same conf as the current fit.
fit_test <- sam.fit(fit_wgwide2025$dat,fit$conf,par=defpar(fit_wgwide2025$dat,fit$conf))
pe_old <- reshape2::melt(get.procres.sdrep(fit_test))
colnames(pe_old)<-c('year','age','pe')
pe_old$case <- previous_case
pe_old <- plyr::join(pe_old , data.frame(age = fit_test$conf$minAge:fit_test$conf$maxAge,
                                         par = exp(fit_test$pl$logSdLogN[fit_test$conf$keyVarLogN+1]))
)
N_tab_old <- reshape2::melt(ntable(fit_test))
colnames(N_tab_old)<-c('year','age','N')
pe_N_old <- plyr::join(pe_old,N_tab_old)
pe_N <- rbind(pe_N,pe_N_old)

#Same idea, for the benchmark fit.
fit_test <- sam.fit(fit_benchmark$dat,fit$conf,par=defpar(fit_benchmark$dat,fit$conf))
pe_old <- reshape2::melt(get.procres.sdrep(fit_test))
colnames(pe_old)<-c('year','age','pe')
pe_old$case <- 'benchmark'
pe_old <- plyr::join(pe_old , data.frame(age = fit_test$conf$minAge:fit_test$conf$maxAge,
                                         par = exp(fit_test$pl$logSdLogN[fit_test$conf$keyVarLogN+1]))
)
N_tab_old <- reshape2::melt(ntable(fit_test))
colnames(N_tab_old)<-c('year','age','N')
pe_N_old <- plyr::join(pe_old,N_tab_old)
pe_N <- rbind(pe_N,pe_N_old)



pe_N$YC <- pe_N$year-pe_N$age

#Join key must include `case`, otherwise the same year class (YC) shared
#across WGWIDE2026/WGWIDE2025/benchmark causes a many-to-many join below
#(each pe_N row matched against every case's age-2 N for that YC), silently
#duplicating rows and attaching the wrong case's R.
dd<-pe_N[c('year','age','case','YC','N')]
dd<-unique(dd[dd$age==2,])
names(dd)[5]<-'R'
dd$age<-NULL
dd$year <- NULL
pe_N<-plyr::join(pe_N,dd)

#Table version of pe_N (process error by year/age), for this year's
#assessment only (WGWIDE2026) - WGWIDE2025/benchmark are kept in pe_N for
#the comparison plots below, but aren't part of the reproducibility check.
#Rounded to 3 digits, same reasoning as parameter_estimates.txt above -
#test.R diffs .txt files with identical() (no tolerance).
pe_N_out <- pe_N[pe_N$case==current_case, ]
pe_N_out[c('pe','par','N','R')] <- round(pe_N_out[c('pe','par','N','R')], digits = 3)
write.csv(pe_N_out, file = 'output/tables/process_error.txt', row.names = FALSE)

plot_processerror <- ggplot(pe_N[pe_N$age !=2, ], aes(x = year,y=pe,colour=case,fill=case)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.4, alpha = 0.4) +
  facet_wrap(~age) +
  xlab('year') +
  ylab('Process Error') +
  theme_bw() 
plot_processerror
ggsave("output/figures/fig_processError.png", plot_processerror, width = 12, height = 8, dpi = 600)

#Presentation-only variant: drop the benchmark case, keep only this year's
#fit vs. WGWIDE2025 (report.R still uses fig_processError.png above, which
#keeps the benchmark comparison).
plot_processerror_pres <- ggplot(pe_N[pe_N$age !=2 & pe_N$case!='benchmark', ], aes(x = year,y=pe,colour=case,fill=case)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.4, alpha = 0.4) +
  facet_wrap(~age) +
  xlab('year') +
  ylab('Process Error') +
  theme_bw()
ggsave("output/figures/fig_processError_presentation.png", plot_processerror_pres, width = 12, height = 8, dpi = 600)



#================================================================================#
# 5. LEAVE-OUT & RETROSPECTIVE DIAGNOSTICS
#================================================================================#







#------------------------------------------------------------------------------#
#Leave-one-fleet-out plot: LO_fit (model.R's leaveout()) drops one fleet's
#observations entirely, across all years, per run.
#------------------------------------------------------------------------------#
LO_fit$base<-attr(LO_fit,'fit')
names(LO_fit)[names(LO_fit)=='w.o. recBar']<-'BESS'
names(LO_fit)[names(LO_fit)=='w.o. Bar']<-'IESNS, bar'

#p_LO_ssb/p_retro_ssb are built once here (with the full, un-filtered
#SSB+FbarW+R table attached by plot.LO()/plot.retro() as attr(.,'data_out'))
#and reused below in fig_SSB_LR.png, so the tables aren't computed twice.
p_LO_ssb    <- plot.LO(LO_fit,what = 'SSB')
p_retro_ssb <- plot.retro(ret,what = 'SSB')

#Table versions of the leave-out and retrospective series (SSB/FbarW/R per
#case/peel/year), underlying fig_LO*.png/fig_SSB_LR.png/fig_Fbar_LR.png/
#fig_R_LR.png/fig_retro*.png below. Rounded to 3 digits, same reasoning as
#parameter_estimates.txt above - test.R diffs .txt files with identical()
#(no tolerance).
leave_out_table <- attr(p_LO_ssb,'data_out')
leave_out_table[c('Estimate','Low','High')] <- round(leave_out_table[c('Estimate','Low','High')], digits = 3)
write.csv(leave_out_table, file = 'output/tables/leave_out.txt', row.names = FALSE)

retro_table <- attr(p_retro_ssb,'data_out')
retro_table[c('Estimate','Low','High')] <- round(retro_table[c('Estimate','Low','High')], digits = 3)
write.csv(retro_table, file = 'output/tables/retro.txt', row.names = FALSE)

plot_SSB_LR<-gridExtra::grid.arrange(p_LO_ssb+ggtitle('SSB'),
                                     p_retro_ssb,ncol=2)
ggsave("output/figures/fig_SSB_LR.png", plot_SSB_LR, width = 12, height = 8, dpi = 600)


plot_SSB_LR<-gridExtra::grid.arrange(plot.LO(LO_fit,what = 'FbarW')+ggtitle('FbarW'),
                                     plot.retro(ret,what = 'FbarW'),ncol=2)
ggsave("output/figures/fig_Fbar_LR.png", plot_SSB_LR, width = 12, height = 8, dpi = 600)


plot_SSB_LR<-gridExtra::grid.arrange(plot.LO(LO_fit,what = 'R')+ggtitle('Recruitment'),
                                     plot.retro(ret,what = 'R'),ncol=2)
ggsave("output/figures/fig_R_LR.png", plot_SSB_LR, width = 12, height = 8, dpi = 600)


#Presentation-only combined figures: SSB, FbarW and Recruitment side by side
#in one panel, one for leave-out and one for the 5-year retrospective, instead
#of the three separate per-variable slides above.
plot_LO_all<-gridExtra::grid.arrange(plot.LO(LO_fit,what = 'SSB')+ggtitle('SSB')+theme(legend.position = "none"),
                                     plot.LO(LO_fit,what = 'FbarW')+ggtitle('FbarW')+theme(legend.position = "none"),
                                     plot.LO(LO_fit,what = 'R')+ggtitle('Recruitment'),ncol=3)
ggsave("output/figures/fig_LO_all.png", plot_LO_all, width = 16, height = 6, dpi = 600)

plot_retro_all<-gridExtra::grid.arrange(plot.retro(ret,what = 'SSB'),
                                        plot.retro(ret,what = 'FbarW'),
                                        plot.retro(ret,what = 'R'),ncol=3)
ggsave("output/figures/fig_retro_all.png", plot_retro_all, width = 16, height = 6, dpi = 600)



plot_LO<-gridExtra::grid.arrange(plot.LO(LO_fit,what = 'SSB')+theme(legend.position = "none"),
                                 plot.LO(LO_fit,what = 'FbarW')+theme(legend.position = "none"),
                                 plot.LO(LO_fit,what = 'R'),ncol=1)


ggsave("output/figures/fig_LO.png", plot_LO, width = 12, height = 8, dpi = 600)



#------------------------------------------------------------------------------#
#Leave-out-last-year: same idea as fig_LO.png above, but each run only drops
#the most recent year of ONE fleet (LO_lastyear, from model.R) instead of
#the whole fleet - shows how much a single new data point per fleet moves
#the terminal-year estimates.
#
#LO_lastyear's names are model.R's raw fleetNames ('IESNS','recBar','Bar',...),
#unlike LO_fit above (which leaveout() already prefixes with 'w.o. ', then
#gets BESS/IESNS-bar renamed). plot.LO() assigns colours by alphabetical
#case order, so mismatched naming here would give the same fleet a different
#colour than in fig_LO.png/fig_SSB_LR.png/etc. Apply the same convention so
#case labels and colours line up across every leave-out figure.
#------------------------------------------------------------------------------#
names(LO_lastyear)[-1] <- paste0('w.o. ',names(LO_lastyear)[-1])
names(LO_lastyear)[names(LO_lastyear)=='w.o. recBar']<-'BESS'
names(LO_lastyear)[names(LO_lastyear)=='w.o. Bar']<-'IESNS, bar'
plot_LO_lastyear<-gridExtra::grid.arrange(plot.LO(LO_lastyear,what = 'SSB')+theme(legend.position = "none"),
                                          plot.LO(LO_lastyear,what = 'R'),ncol=1)

ggsave("output/figures/fig_LO_lastyear.png", plot_LO_lastyear, width = 12, height = 8, dpi = 600)





#------------------------------------------------------------------------------#
#Make the retro plot
#------------------------------------------------------------------------------#
plot_retro <- gridExtra::grid.arrange(plot.retro(ret,what = 'SSB')+theme(legend.position = "none"),
                                      plot.retro(ret,what = 'FbarW')+theme(legend.position = "none"),
                                      plot.retro(ret,what = 'R'),ncol=1)

ggsave("output/figures/fig_retro.png", plot_retro, width = 12, height = 8, dpi = 600)

#Mohn's rho for the retrospective (ret) - currently only shown baked into
#fig_retro*.png's facet titles (see mohn() call in plot.retro(),
#boot/rscript/utils.r); written out here as its own table so the values are
#machine-readable too.
mohns_rho <- mohn(ret)
write.csv(data.frame(type = c('R','SSB','FbarW'), mohns_rho = round(as.numeric(mohns_rho)[1:3], digits = 3)),
          file = 'output/tables/mohns_rho.txt', row.names = FALSE)



#================================================================================#
# 6. STOCK SUMMARY COMPARISON
#================================================================================#

#------------------------------------------------------------------------------#
#Stock in numbers (millions) and fishing mortality at age, same precision as
#report.R's tables 4.5.1.2/4.5.1.3 (round() with no decimals / 3 decimals).
#------------------------------------------------------------------------------#
write.csv(round(ntable(fit)), file = 'output/tables/ntable.txt')
write.csv(round(faytable(fit), digits = 3), file = 'output/tables/faytable.txt')



#------------------------------------------------------------------------------#
#Make the stock summary plot
#------------------------------------------------------------------------------#
fbar <- compute_FbarW(fit)
ssb <- as.data.frame(ssbtable(fit))
ssb$year <- as.numeric(rownames(ssb))

rec <- as.data.frame(rectable(fit))
rec$year <- as.numeric(rownames(rec))

plot_status <- gridExtra::grid.arrange(
  ggplot(ssb,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_line()+geom_ribbon(alpha=0.4)+
    theme_bw() +ggtitle('SSB') ,
  
  ggplot(fbar,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_line()+geom_ribbon(alpha=0.4)+
    theme_bw() +ggtitle('FbarW'),
  
  ggplot(rec,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_line()+geom_ribbon(alpha=0.4)+
    theme_bw() +ggtitle('R'),
  
  ncol=3
)

ggsave("output/figures/fig_summary.png", plot_status, width = 16, height = 6, dpi = 600)



#------------------------------------------------------------------------------#
#Make the stock summary comparison plot - same 3-panel layout as fig_summary.png
#above, but this year's fit (current_case) overlaid against last year's
#official assessment (fit_wgwide2025).
#------------------------------------------------------------------------------#
fbar_compare <- compute_FbarW(fit_wgwide2025)
ssb_compare <- as.data.frame(ssbtable(fit_wgwide2025))
ssb_compare$year <- as.numeric(rownames(ssb_compare))
rec_compare <- as.data.frame(rectable(fit_wgwide2025))
rec_compare$year <- as.numeric(rownames(rec_compare))

colour_scale_compare <- scale_colour_manual(values = setNames(c('black','red'), c(current_case,previous_case)))
fill_scale_compare   <- scale_fill_manual(values = setNames(c('black','red'), c(current_case,previous_case)))

plot_status_compare <- gridExtra::grid.arrange(
  ggplot(ssb,aes(x=year,y=Estimate,ymin=Low,ymax=High))+
    geom_line(aes(colour=current_case))+geom_ribbon(aes(fill=current_case),alpha=0.4)+
    geom_line(data=ssb_compare,aes(colour=previous_case))+
    # geom_ribbon(data=ssb_compare,aes(fill=previous_case),alpha=0.4)+
    theme_bw()+theme(legend.position = "none")+ggtitle('SSB')+
    colour_scale_compare+fill_scale_compare,

  ggplot(fbar,aes(x=year,y=Estimate,ymin=Low,ymax=High))+
    geom_line(aes(colour=current_case))+geom_ribbon(aes(fill=current_case),alpha=0.4)+
    geom_line(data=fbar_compare,aes(colour=previous_case))+
    # geom_ribbon(data=fbar_compare,aes(fill=previous_case),alpha=0.4)+
    theme_bw()+theme(legend.position = "none")+ggtitle('FbarW')+
    colour_scale_compare+fill_scale_compare,

  ggplot(rec,aes(x=year,y=Estimate,ymin=Low,ymax=High))+
    geom_line(aes(colour=current_case))+geom_ribbon(aes(fill=current_case),alpha=0.4)+
    geom_line(data=rec_compare,aes(colour=previous_case))+
    # geom_ribbon(data=rec_compare,aes(fill=previous_case),alpha=0.4)+
    theme_bw()+ggtitle('R')+theme(legend.position = "none")+
    colour_scale_compare+fill_scale_compare,

  ncol=3
)

ggsave("output/figures/fig_summary_compare.png", plot_status_compare, width = 16, height = 6, dpi = 600)





fbar2 <- compute_FbarW(fit_wgwide2025)

ssb <- as.data.frame(ssbtable(fit))
ssb$year <- as.numeric(rownames(ssb))
ssb2 <- as.data.frame(ssbtable(fit_wgwide2025))
ssb2$year <- as.numeric(rownames(ssb2))

rec2 <- as.data.frame(rectable(fit_wgwide2025))
rec2$year <- as.numeric(rownames(rec2))
rec <- as.data.frame(rectable(fit))
rec$year <- as.numeric(rownames(rec))


ggplot(rec,aes(x=year,y=Estimate,ymin=Low,ymax=High))+
  geom_line(aes(colour='TMP'))+geom_ribbon(alpha=0.4)+
  geom_line(data=rec2,aes(colour=previous_case))+
  theme_bw() +ggtitle('R')


ggplot(ssb,aes(x=year,y=Estimate,ymin=Low,ymax=High))+
  geom_line(aes(colour='TMP'))+geom_ribbon(alpha=0.4)+
  geom_line(data=ssb2,aes(colour=previous_case))+
  theme_bw() +ggtitle('R')

