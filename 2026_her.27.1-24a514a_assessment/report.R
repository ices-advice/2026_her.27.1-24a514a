#------------------------------------------------------------------------------#
#report.R
#
#Author: 
#Sindre Vatnehol
#Institute of marine Research
#Norway
#
#
#Description:
#Here we produce all the figures and tables that is used for the report section
#in WGWIDE
#
#Input:  model/fit.Rda, model/diag.Rda (from model.R), boot/rscript/WGWIDE2025.Rda
#        (last year's official fit, used throughout for "previous vs. updated"
#        comparison figures/tables), data/tab.4.4.8.*.txt (from data.R)
#Output: report/tables/*.csv|.txt, report/figures/Fig.4.5.1.*.png
#        These are the tables/figures that go into the WGWIDE report.
#
#Chapters:
#  1. Setup
#  2. Tables
#  3. Figures
#------------------------------------------------------------------------------#



#================================================================================#
# 1. SETUP
#================================================================================#

#------------------------------------------------------------------------------#
#remove everything from memory
#------------------------------------------------------------------------------#
rm(list=ls())



#------------------------------------------------------------------------------#
#Load utils
#------------------------------------------------------------------------------#
source('boot/rscript/utils.r')
library(TAF)
library('stockassessment')
library('ggplot2')
library(gridExtra)
library(reshape2)



#------------------------------------------------------------------------------#
#load fit objects
#fit_old = last year's official assessment (fixed reference point for
#          "previous vs. updated" comparisons throughout this script)
#fit     = this year's fit, loaded via model/fit.Rda (overwrites the `fit`
#          used above for WGWIDE2025.Rda, so order of these two loads matters)
#diag.Rda supplies fit_pro, res_fit, ret, LO_fit used further down
#(ret2, the alternative-RFID-peel retro from model.R, is loaded here too
#but is not actually plotted/used anywhere in this script - see Fig.4.5.1.9)
#------------------------------------------------------------------------------#
load('boot/rscript/WGWIDE2025.Rda')
fit_old <- fit
load('model/fit.Rda')
load('model/diag.Rda')




#------------------------------------------------------------------------------#
#prepare folder structure
#------------------------------------------------------------------------------#
mkdir("report")
mkdir("report/tables")
mkdir("report/figures")



#================================================================================#
# 2. TABLES
#================================================================================#

#------------------------------------------------------------------------------#
#These 5 files (smoothed/external CV per fleet, from computeExternalWeight())
#were computed and written in data.R - this section only relocates them into
#report/tables, it does not recompute anything.
#------------------------------------------------------------------------------#
# Table 4.4.8.1. Norwegian spring-spawning herring. Coefficient of variance of catch-at-age used by SAM.
file.rename("data/tab.4.4.8.1.txt", "report/tables/tab.4.4.8.1.txt")
#Table 4.4.8.2. Norwegian spring-spawning herring. Coefficient of variance of NASF used by SAM.
file.rename("data/tab.4.4.8.2a.txt", "report/tables/tab.4.4.8.2a.txt")
file.rename("data/tab.4.4.8.2b.txt", "report/tables/tab.4.4.8.2b.txt")
#Table 4.4.8.3. Norwegian spring-spawning herring. Coefficient of variance of IESNS bar used by SAM.
file.rename("data/tab.4.4.8.3.txt", "report/tables/tab.4.4.8.3.txt")
#Table 4.4.8.4. Norwegian spring-spawning herring. Coefficient of variance of IESNS used by SAM.
file.rename("data/tab.4.4.8.4.txt", "report/tables/tab.4.4.8.4.txt")
#Table 4.4.8.5. Norwegian spring-spawning herring. Coefficient of variance of RFID used by SAM.
file.rename("data/tab.4.4.8.5.txt", "report/tables/tab.4.4.8.5.txt")




#------------------------------------------------------------------------------#
#Table 4.5.1.1. Norwegian spring-spawning herring. Parameter estimates of the
#final SAM model fit. The estimates from the benchmark 2025 is shown.
#------------------------------------------------------------------------------#
#Only built when `fit` is exactly one year ahead of fit_old (i.e. a normal
#annual update); this avoids a crash when testing the repo against an
#older/different reference fit. If the condition fails for any other reason
#(e.g. fit_old is stale, or this year's data spans more than one extra
#year), the table is skipped without a warning.
if(max(fit$data$years)==max(fit_old$data$years)+1){
  new <- as.data.frame(partable(fit)[,c(1,2)])
  names(new)<-paste(c('Estimate','sd'),max(fit$data$years))
  old <- as.data.frame(partable(fit_old)[,c(1,2)])
  names(old)<-paste(c('Estimate','sd'),max(fit_old$data$years))
  new$Parameter <- rownames(new)
  rownames(new)<-NULL
  old$Parameter <- rownames(old)
  rownames(old)<-NULL
  tab<-plyr::join(new,old)
  tab[,c(1,2,4,5)]<-round(tab[,c(1,2,4,5)],digits = 3)
  write.csv(tab[,c(3,1,2,4,5)],file = 'report/tables/4.5.1.1.csv')
}



#------------------------------------------------------------------------------#
#Table 4.5.1.2. Norwegian spring-spawning herring. Point estimates of Stock in
#numbers (millions). ntable() is a stockassessment built-in: stock size at
#age by year, straight from the fit - no further processing here.
#------------------------------------------------------------------------------#
write.csv(round(ntable(fit)),file = 'report/tables/4.5.1.2.csv')



#------------------------------------------------------------------------------#
# Table 4.5.1.3. Norwegian spring-spawning herring. Point estimates of
# fishing mortality. faytable() is a stockassessment built-in: F-at-age by
# year, straight from the fit - no further processing here.
#------------------------------------------------------------------------------#
write.csv(round(faytable(fit),digits=3),file = 'report/tables/4.5.1.3.csv')



#------------------------------------------------------------------------------#
# Table 4.5.1.4. Norwegian spring-spawning herring. Final stock summary table.
#High and low represent approximate 95% confidence limits.
#Builds one combined table from recruitment, SSB, landings and FbarW, with
#an extra "Average" row appended at the bottom.
#------------------------------------------------------------------------------#
rectab <- as.data.frame(rectable(fit))
rectab<-round(rectab[c('Estimate','High','Low')])

ssbtab <- as.data.frame(ssbtable(fit))
ssbtab<-round(ssbtab[c('Estimate','High','Low')])

fbartab <- compute_FbarW(fit)
fbartab<-fbartab[c('Estimate','High','Low')]

#Landings by year = sum over ages of (catch mean weight x catch-at-age).
#Re-reads caa.txt here (rather than pulling catch-at-age back out of `fit`)
#since fit$data$catchMeanWeight is per-age weight, not catch numbers.
#Negative caa values (SAM/ICES convention for "no catch recorded", as
#opposed to a true zero) are treated as NA so they don't pull the sum down.
caa <- stockassessment::read.ices(taf.data.path('caa.txt'))
caa[caa<0]<-NA


catch<-data.frame(catch=round(rowSums(fit$data$catchMeanWeight[,,1]*caa,na.rm = T),3))

tmp <- c(
  round(mean(rectab$Estimate)),round(mean(rectab$High)),round(mean(rectab$Low)),
  round(mean(ssbtab$Estimate)),round(mean(ssbtab$High)),round(mean(ssbtab$Low)),
  round(mean(catch$catch)),
  round(mean(fbartab$Estimate),3),round(mean(fbartab$High),3),round(mean(fbartab$Low),3)
  
)

fbartab<-rbind(round(fbartab,3),c('','',''))
catch <- rbind(catch,'')
tab<-cbind(rectab,ssbtab,catch,fbartab)
out <- rbind(tab,tmp)
rownames(out)[length(rownames(out))] = 'Average'
colnames(out)<-c('Recruitment','High','Low','SSB','High','Low','Catch','Fbar','High','Low')
write.table(out,file="report/tables/tab4.5.1.4.txt",sep=";",quote=F)



#Table 4.8.1.1 (input to short-term prediction) and Table 4.8.2.1 (catch
#options) are not produced here - they are written by forecast.R, which runs
#after this script.



#================================================================================#
# 3. FIGURES
#================================================================================#

#List of figures
# 
library(reshape2)
one <- reshape2::melt(fit_old$data$propMat[rownames(fit_old$data$propMat)%in%c(2020,2021),])
names(one) <- c('year','age','mo')
one$case <- 'assumed'
two <- reshape2::melt(fit$data$propMat[rownames(fit$data$propMat)%in%c(2020,2021),])
names(two) <- c('year','age','mo')
two$case <- 'observed'
data = rbind(one,two)
p<-ggplot(data[data$age%in%c(2:8),],aes(x=age,y=mo,colour=case))+geom_line()+
  facet_wrap(~year)+theme_bw()+ylab('maturity at age')


ggsave(filename = 'report/figures/Fig.4.4.5.1.png',p, height=10, width=20,units="cm", dpi = 600,  bg = "transparent")



#------------------------------------------------------------------------------#
# Figure 4.5.1.1
#Three-panel comparison of "this year" vs "previous year" (fit vs. fit_old)
#SAM parameter estimates: catchability (logFpar), observation variance
#(logSdLogObs), and everything else. The logFpar_N / logSdLogObs_N -> fleet
#name lookups below rely on the SAME fixed fleet order documented in
#model.R (row 1 = catch, 2 = IESNS, 3 = recBar/BESS, 4 = Bar/IESNS_rec,
#5 = NASF_old, 6 = NASF_new, 7 = RFID) - if that order ever changes, these
#labels will silently point at the wrong fleet/age.
#------------------------------------------------------------------------------#
dd <- as.data.frame(partable(fit))
dd$parname <- rownames(dd)
dd$case <- 'updated'
dd_out <- dd

dd <- as.data.frame(partable(fit_old))
dd$parname <- rownames(dd)
dd$case <- 'previous'
dd_out <- rbind(dd_out,dd)



dd_out$facet <- 'rest'
dd_out$facet[grepl('logFpar',dd_out$parname)]<-'Catchability'


dd_out$parname[dd_out$parname=='logFpar_0']<-'IESNS age 03'
dd_out$parname[dd_out$parname=='logFpar_1']<-'IESNS age 04'
dd_out$parname[dd_out$parname=='logFpar_2']<-'IESNS age 05'
dd_out$parname[dd_out$parname=='logFpar_3']<-'IESNS age 06'
dd_out$parname[dd_out$parname=='logFpar_4']<-'IESNS age 07'
dd_out$parname[dd_out$parname=='logFpar_5']<-'IESNS age 08'
dd_out$parname[dd_out$parname=='logFpar_6']<-'IESNS age 09'
dd_out$parname[dd_out$parname=='logFpar_7']<-'IESNS age 10'
dd_out$parname[dd_out$parname=='logFpar_8']<-'IESNS age 11-12'
dd_out$parname[dd_out$parname=='logFpar_9']<-'BESS age 02'
dd_out$parname[dd_out$parname=='logFpar_10']<-'BESS age 03'
dd_out$parname[dd_out$parname=='logFpar_11']<-'IESNS_rec age 02'
dd_out$parname[dd_out$parname=='logFpar_12']<-'NASF_old age 03'
dd_out$parname[dd_out$parname=='logFpar_13']<-'NASF_old age 04'
dd_out$parname[dd_out$parname=='logFpar_14']<-'NASF_old age 05'
dd_out$parname[dd_out$parname=='logFpar_15']<-'NASF_old age 06'
dd_out$parname[dd_out$parname=='logFpar_16']<-'NASF_old age 07'
dd_out$parname[dd_out$parname=='logFpar_17']<-'NASF_old age 08-12'
dd_out$parname[dd_out$parname=='logFpar_18']<-'NASF_new age 03'
dd_out$parname[dd_out$parname=='logFpar_19']<-'NASF_new age 04'
dd_out$parname[dd_out$parname=='logFpar_20']<-'NASF_new age 05'
dd_out$parname[dd_out$parname=='logFpar_21']<-'NASF_new age 06'
dd_out$parname[dd_out$parname=='logFpar_22']<-'NASF_new age 07'
dd_out$parname[dd_out$parname=='logFpar_23']<-'NASF_new age 08'
dd_out$parname[dd_out$parname=='logFpar_24']<-'NASF_new age 09'
dd_out$parname[dd_out$parname=='logFpar_25']<-'NASF_new age 10'
dd_out$parname[dd_out$parname=='logFpar_26']<-'NASF_new age 11-12'
dd_out$parname[dd_out$parname=='logFpar_27']<-'RFID age 03-12'

#fleet is derived from the full parname (e.g. "IESNS age 03" -> "IESNS")
#before the catchability parname is shortened down to just the age below,
#so facet_wrap(~fleet) further down still gets the right fleet per panel.
dd_out$fleet <- sub(" .*", "", dd_out$parname)
#Shorten the catchability parname down to just the age (e.g.
#"IESNS age 03" -> "03") for the x-axis, since each panel is already
#faceted by fleet and doesn't need the fleet name repeated on every tick.
#NB: subset on the row index vector (`idx`) rather than dd_out[dd_out$facet
#== 'Catchability'] directly - the latter is missing the comma for row
#indexing on a data.frame, so R treats it as (invalid) column selection.
idx <- dd_out$facet=='Catchability'
dd_out$parname[idx] <- sub(".*?(\\d{2}(?:-\\d{2})?)$", "\\1", dd_out$parname[idx])


plot_catchability <- ggplot(dd_out[dd_out$facet == 'Catchability', ],
                            aes(x = parname, y = par, fill = case)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = par - `sd(par)`, ymax = par + `sd(par)`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Parameter", y = "Estimate", fill = "Case") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~fleet, scales = "free_x")+ggtitle('Catchability parameter')+xlab('ages')





dd_out$facet[grepl('logSdLogObs',dd_out$parname)]<-'ObservationVariance'


dd_out$parname[dd_out$parname=='logSdLogObs_0']<-'Catch age 02'
dd_out$parname[dd_out$parname=='logSdLogObs_1']<-'Catch age 03'
dd_out$parname[dd_out$parname=='logSdLogObs_2']<-'Catch age 04'
dd_out$parname[dd_out$parname=='logSdLogObs_3']<-'Catch age 05-12'

dd_out$parname[dd_out$parname=='logSdLogObs_4']<-'IESNS age 03'
dd_out$parname[dd_out$parname=='logSdLogObs_5']<-'IESNS age 04-12'

dd_out$parname[dd_out$parname=='logSdLogObs_6']<-'BESS age 02-03'
dd_out$parname[dd_out$parname=='logSdLogObs_7']<-'IESNS_rec age 02'

dd_out$parname[dd_out$parname=='logSdLogObs_8']<-'NASF_old age 03'
dd_out$parname[dd_out$parname=='logSdLogObs_9']<-'NASF_old age 04'
dd_out$parname[dd_out$parname=='logSdLogObs_10']<-'NASF_old age 05-12'

dd_out$parname[dd_out$parname=='logSdLogObs_11']<-'NASF_new age 03'
dd_out$parname[dd_out$parname=='logSdLogObs_12']<-'NASF_new age 04'
dd_out$parname[dd_out$parname=='logSdLogObs_13']<-'NASF_new age 05-12'
dd_out$parname[dd_out$parname=='logSdLogObs_14']<-'RFID age 03-09'
dd_out$parname[dd_out$parname=='logSdLogObs_15']<-'RFID age 10-11'
dd_out$parname[dd_out$parname=='logSdLogObs_16']<-'RFID age 12'
dd_out$fleet <- sub(" .*", "", dd_out$parname)
# dd_out$parname <- sub(".*?(\\d{2}(?:-\\d{2})?)$", "\\1", dd_out$parname)   #left commented out - same shortening as the catchability panel, not applied here

#par/sd(par) are on the log scale (SAM's logSdLogObs); plotted here as the
#natural-scale observation SD (exp(par)) instead, with the error bar taken
#by transforming the log-scale Wald interval (exp(par +/- sd(par))), so the
#bar height reflects the actual magnitude of observation uncertainty rather
#than the log-parameter, which can look disproportionately large/small when
#its value is close to zero (i.e. SD close to 1).
dd_obsvar <- dd_out[dd_out$facet == 'ObservationVariance', ]
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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('observation variance')+xlab('')
plot_ObservationVariance

plot_rest <- ggplot(dd_out[dd_out$facet == 'rest', ], 
                    aes(x = parname, y = par, fill = case)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = par - `sd(par)`, ymax = par + `sd(par)`),
                width = 0.2, position = position_dodge(width = 0.9)) +
  labs(x = "Parameter", y = "Estimate", fill = "Case") +
  theme_bw() +xlab('')+ggtitle('Rest')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_rest

p <- gridExtra::grid.arrange(
  plot_catchability,
  plot_ObservationVariance,
  plot_rest,
  ncol = 1,
  heights = c(1.5, 1, 1)  # relative row heights
)

ggsave(filename = 'report/figures/Fig.4.5.1.1.png',p, height=30, width=20,units="cm", dpi = 600,  bg = "transparent")




#------------------------------------------------------------------------------#
# Figure 4.5.1.2.
#F-at-age normalised to the age with the highest F within each year (so
#each year's surface peaks at 1), shown as four 3D perspective views of the
#same F-pattern-over-time surface from different rotation angles.
#------------------------------------------------------------------------------#
FAY <- faytable(fit)
Fpattern<- t(faytable(fit))
for(i in 1:ncol(Fpattern))Fpattern[,i]<-Fpattern[,i]/max(Fpattern[,i])
Fpattern<- t(Fpattern)
png(file="report/figures/Fig.4.5.1.2.png")
par(mfrow=c(2,2),mar=c(2,2,1,1))
persp(y=2:12,x=as.numeric(rownames(Fpattern)),z=(Fpattern),shade=0.5,theta=-50,phi=20,
      ticktype="detailed",xlab="Year",ylab="Age",zlab="F-pattern",col="lightblue",zlim=c(0,1.2))
persp(y=2:12,x=as.numeric(rownames(Fpattern)),z=(Fpattern),shade=0.5,theta=-20,phi=20,
      ticktype="detailed",xlab="Year",ylab="Age",zlab="F-pattern",col="lightblue",zlim=c(0,1.2))
persp(y=2:12,x=as.numeric(rownames(Fpattern)),z=(Fpattern),shade=0.5,theta=20,phi=20,
      ticktype="detailed",xlab="Year",ylab="Age",zlab="F-pattern",col="lightblue",zlim=c(0,1.2))
persp(y=2:12,x=as.numeric(rownames(Fpattern)),z=(Fpattern),shade=0.5,theta=50,phi=20,
      ticktype="detailed",xlab="Year",ylab="Age",zlab="F-pattern",col="lightblue",zlim=c(0,1.2))
dev.off()






#------------------------------------------------------------------------------#
# Figure 4.5.1.3.
#Per-fleet plot of the SAM observation weight (from computeExternalWeight()
#in data.R) by year/age, i.e. a visual check of how much influence each
#observation actually has in the fit. Fleet numbers below (1=catch,
#2=IESNS, 3=BESS, 4=IESNS_rec, 5=NASF_old, 6=NASF_new, 7=RFID) follow the
#same fixed fleet order documented in model.R. Only plot1/2/5/6/7 (catch,
#IESNS, NASF_old, NASF_new, RFID) are actually combined into the saved
#figure below - plot3 (BESS) and plot4 (IESNS_rec) are built but not
#included in the grid.arrange() call further down.
#------------------------------------------------------------------------------#
dat <- fit$data
tmp <- as.data.frame(dat$aux)
tmp$weight <- dat$weight
tmp$log_obs <- dat$logobs


#prepare each panel
p<-plot.data.weight(tmp[tmp$fleet==1 & !is.na(tmp$log_obs),])
plot1<-p+ggtitle('Residual catch weight')+xlim(range(tmp$year))+ylim(range(tmp$age))

p<-plot.data.weight(tmp[tmp$fleet==2& !is.na(tmp$log_obs),])
plot2<-p+ggtitle('IESNS weight')+xlim(range(tmp$year))+ylim(range(tmp$age))
p<-plot.data.weight(tmp[tmp$fleet==3& !is.na(tmp$log_obs),])
plot3<-p+ggtitle('BESS weight')+xlim(range(tmp$year))+ylim(range(tmp$age))+
  scale_x_continuous(breaks=seq(1990,2050,2))
p<-plot.data.weight(tmp[tmp$fleet==4& !is.na(tmp$log_obs),])
plot4<-p+ggtitle('IESNS rec weight')+xlim(range(tmp$year))+ylim(range(tmp$age))+
  scale_x_continuous(breaks=seq(1990,2050,2))
p<-plot.data.weight(tmp[tmp$fleet==5& !is.na(tmp$log_obs),])
plot5<-p+ggtitle('NASF old weight')+
  scale_x_continuous(breaks=seq(1990,2050,2))
p<-plot.data.weight(tmp[tmp$fleet==6& !is.na(tmp$log_obs),])
plot6<-p+ggtitle('NASF new weight')+
  scale_x_continuous(breaks=seq(1990,2050,2))
p<-plot.data.weight(tmp[tmp$fleet==7& !is.na(tmp$log_obs),])
plot7<-p+ggtitle('RFID weight')+
  scale_x_continuous(breaks=seq(1990,2050,1))
#merge all panels and plot
p<-gridExtra::grid.arrange(plot1, plot2,plot5,plot6,plot7, ncol = 2)
ggsave('report/figures/Fig.4.5.1.3.png',p,width=15, height=15, units="cm", 
       device = "tiff", dpi = 300,  bg = "transparent")



#------------------------------------------------------------------------------#
# Figure 4.5.1.4.
#Standard SAM residual bubble plot (res_fit from model/diag.Rda), with the
#internal fleet names recBar/Bar relabelled to the report's fleet names
#(BESS / IESNS,recruitment) before plotting.
#------------------------------------------------------------------------------#
attributes(res_fit)$fleetNames[attributes(res_fit)$fleetNames=='recBar']<-'BESS'
attributes(res_fit)$fleetNames[attributes(res_fit)$fleetNames=='Bar']<-'IESNS, recruitment'
png(file="report/figures/Fig.4.5.1.4.png",
    width = 1200,   # width in pixels
    height = 1200,  # height in pixels
    res = 150)
plot(res_fit)
dev.off()






#-------------------------------------------------------------------------------
# Figure 4.5.1.5.
#Observed-vs-predicted scatter and QQ-plot of the (already-relabelled)
#residuals, faceted by fleet - a closer look at the same res_fit used in
#Fig.4.5.1.4, this time as observed/predicted pairs rather than a bubble plot.
#-------------------------------------------------------------------------------
class(res_fit)<-'data.frame'
res_fit<-res_fit[!is.na(res_fit$observation),]
res_fit$fleet<-attributes(res_fit)$fleetNames[res_fit$fleet]

pp1<-ggplot(res_fit,aes(x=observation,y=observation+residual))+geom_point()+
  facet_wrap(~fleet,ncol=1,scales = 'free')+
  ylab('Predicted')+xlab('Observed')+
  theme(
    panel.background = element_rect(fill = "white"),  # Set background to white
    panel.grid.major = element_blank(),  # Remove grid
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )+geom_smooth(method='lm',colour='green',se = F,size=0.5)+
  ggtitle('Observed vs Predicted')

pp2<-ggplot(res_fit,aes(sample=residual))+geom_qq(distribution = stats::qnorm)+geom_qq_line()+facet_wrap(~fleet,ncol = 1,scales = 'free_y')+
  ylab('Sample Quantiles')+xlab('Theoretical Quantiles')+
  theme(
    panel.background = element_rect(fill = "white"),  # Set background to white
    panel.grid.major = element_blank(),  # Remove grid
    panel.grid.minor = element_blank(),
    panel.border = element_blank()
  )+ggtitle('Norm Q-Q Plot')

p<-grid.arrange(pp1, pp2,ncol = 2)

ggsave('report/figures/Fig.4.5.1.5.png',p,width=25, height=30, units="cm", 
       device = "tiff", dpi = 300,  bg = "transparent")






#-------------------------------------------------------------------------------
#Figure 4.5.1.6, Process residuals
#Standard SAM plot of fit_pro (procres(), from model/diag.Rda): residuals
#of the state process (log-N random walk), as opposed to the observation
#residuals shown in Fig.4.5.1.4/4.5.1.5.
#-------------------------------------------------------------------------------
png(file="report/figures/Fig.4.5.1.6.png",
    width = 1200,   # width in pixels
    height = 1000,  # height in pixels
    res = 150)
plot(fit_pro)
dev.off()




#-------------------------------------------------------------------------------
#Figure 4.5.1.7, QQplot of process residuals
#fit_pro has 2 "fleets" here in the process-residual sense: 1 = the joint
#sample residual on log(N) (the state process itself), 2 = residual catch
#(process residual specific to the catch equation).
#-------------------------------------------------------------------------------
class(fit_pro)<-'data.frame'
fit_pro$fleet[fit_pro$fleet==1]<-'joint sample residual log(N)'
fit_pro$fleet[fit_pro$fleet==2]<-'Residual catch'
plot_residualqq<-ggplot(fit_pro,aes(sample=residual))+
  geom_qq(distribution = stats::qnorm)+
  geom_qq_line()+
  ggtitle('QQplot of process residual')+ylim(c(-5,5))+ylab('')+xlab('')+facet_wrap(~fleet)+
  theme_bw()
plot_rest
ggsave("report/figures/Fig.4.5.1.7.png", plot_residualqq, width = 12, height = 8, dpi = 600)





#------------------------------------------------------------------------------#
#Fig.4.5.1.8: process error
#pe = estimated process error (deviation from the log-N random walk) by
#year/age, joined with the fitted process SD per age (exp(logSdLogN)) so
#the bars can be interpreted relative to the model's assumed process noise.
#------------------------------------------------------------------------------#
pe <- reshape2::melt(get.procres.sdrep(fit))
colnames(pe)<-c('year','age','pe')
pe <- plyr::join(pe , data.frame(age = fit$conf$minAge:fit$conf$maxAge,
                                 par = exp(fit$pl$logSdLogN[fit$conf$keyVarLogN+1]))
)
pe$case <- 'updated'


#Refit last year's data under this year's SAM configuration (fit_test), so
#its process error can be recomputed under the same conf as the current fit
#for a like-for-like comparison. fit_test is only used for this and has no
#effect on the official fit.
fit_test <- sam.fit(fit_old$dat,fit$conf,par=defpar(fit_old$dat,fit$conf))
pe_old <- reshape2::melt(get.procres.sdrep(fit_test))
colnames(pe_old)<-c('year','age','pe')
pe_old$case <- 'previous'
pe_old <- plyr::join(pe_old , data.frame(age = fit_test$conf$minAge:fit_test$conf$maxAge,
                                         par = exp(fit_test$pl$logSdLogN[fit_test$conf$keyVarLogN+1]))
)


N_tab <- reshape2::melt(ntable(fit))
colnames(N_tab)<-c('year','age','N')
pe_N <- plyr::join(pe,N_tab)

N_tab_old <- reshape2::melt(ntable(fit_test))
colnames(N_tab_old)<-c('year','age','N')
pe_N_old <- plyr::join(pe_old,N_tab_old)

pe_N <- rbind(pe_N,pe_N_old)

#Joins in recruitment (N at age 2) per year-class (YC = year - age), purely
#so it's available to facet/label by if needed - age 2 itself is then
#excluded from the process-error plot below (recruitment is not a process
#error in the same sense as older ages' survival deviations).
pe_N$YC <- pe_N$year-pe_N$age

dd<-pe_N[c('year','age','YC','N')]
dd<-unique(dd[dd$age==2,])
names(dd)[4]<-'R'
dd$age<-NULL
dd$year <- NULL
pe_N<-plyr::join(pe_N,dd)

plot_processerror <- ggplot(pe_N[pe_N$age !=2, ], aes(x = year,y=pe,colour=case,fill=case)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.4, alpha = 0.4) +
  facet_wrap(~age) +
  xlab('year') +
  ylab('Process Error') +
  scale_fill_manual(values = c("black", "red")) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw() 


ggsave("report/figures/Fig.4.5.1.8.png", plot_processerror, width = 12, height = 8, dpi = 600)
#Fig.4.5.1.8_2 is simply the equivalent process-error figure produced in
#output.R (a 3-way comparison incl. the pre-benchmark fit), copied in here
#as an alternative/supplementary version - requires output.R to have been
#run first, otherwise this copy silently does nothing useful (file.copy()
#returns FALSE but does not error if the source is missing).
file.copy('output/figures/fig_processError.png','report/figures/Fig.4.5.1.8_2.png',overwrite = T)



#------------------------------------------------------------------------------#
#Fig.4.5.1.9.
#Retrospective (Mohn's rho-style) plot of SSB/FbarW/R over the 5 peels in
#`ret`. Uses `ret` (RFID peeled back 1 extra year) only - `ret2` (RFID
#peeled back 2 extra years, computed in model.R) is not used anywhere in
#the report/output scripts.
#------------------------------------------------------------------------------#
plot_retro <- gridExtra::grid.arrange(plot.retro(ret,what = 'SSB')+theme(legend.position = "none"),
                                      plot.retro(ret,what = 'FbarW')+theme(legend.position = "none"),
                                      plot.retro(ret,what = 'R'),ncol=1)

ggsave("report/figures/Fig.4.5.1.9.png", plot_retro, width = 8, height = 8, dpi = 600)



#------------------------------------------------------------------------------#
# Figure 4.5.1.10.
#Leave-one-fleet-out plot: SSB/FbarW/R from LO_fit (model.R's leaveout()),
#relabelled here to human-readable fleet names for the legend. Only the 2nd
#and 3rd elements are renamed (BESS, IESNS bar) - element 1 is the base fit
#(handled via LO_fit$base just above) and any other left-out fleets keep
#leaveout()'s default "w.o. <name>" label.
#------------------------------------------------------------------------------#
LO_fit$base<-attr(LO_fit,'fit')
names(LO_fit)[2]<-'w.o. BESS'
names(LO_fit)[3] <- 'w.o. IESNS bar'
plot_LO<-gridExtra::grid.arrange(plot.LO(LO_fit,what = 'SSB')+theme(legend.position = "none"),
                                 plot.LO(LO_fit,what = 'FbarW')+theme(legend.position = "none"),
                                 plot.LO(LO_fit,what = 'R'),ncol=1)
ggsave("report/figures/Fig.4.5.1.10.png", plot_LO, width = 8, height = 8, dpi = 600)



#------------------------------------------------------------------------------#
# Figure 4.5.1.11.
#Landings + SSB + FbarW + recruitment, the 4-panel stock-summary figure.
#------------------------------------------------------------------------------#
ssbtab <- as.data.frame(stockassessment::ssbtable(fit))
ssbtab$year<-as.numeric(rownames(ssbtab))
fbarw<-compute_FbarW(fit)
rectab <- as.data.frame(stockassessment::rectable(fit))
rectab$year <- as.numeric(rownames(rectab))
caa <- stockassessment::read.ices(taf.data.path('caa.txt'))
caa[caa<0]<-NA
catch<-data.frame(catch=rowSums(fit$data$catchMeanWeight[,,1]*caa,na.rm = T),
                  year = rownames(fit$data$catchMeanWeight))

p0<-ggplot(catch, aes(x=as.numeric(year), y=catch)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge())+ggtitle('Landings')+ylab('Thousands tonnes')+theme_bw()
p1<-ggplot(ssbtab,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_ribbon(colour='grey',fill='grey')+geom_line()+
  ggtitle('SSB')+ylab('Thousands tonnes')+theme_bw()
p2<-ggplot(fbarw,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_ribbon(colour='grey',fill='grey')+geom_line()+
  ggtitle('Fishing mortality')+ylab(expression(bar(F)[w5-12]))+theme_bw()
p3<-ggplot(rectab,aes(x=year,y=Estimate,ymin=Low,ymax=High))+geom_ribbon(colour='grey',fill='grey')+geom_line()+
  ggtitle('Recruitment (age2)')+ylab('Millions')+theme_bw()


p<-grid.arrange(p0,p3,p2,p1,ncol = 2,nrow=2)
ggsave(filename = 'report/figures/Fig.4.5.1.11.png',p, height=20, width=20,units="cm", dpi = 600,  bg = "transparent")




