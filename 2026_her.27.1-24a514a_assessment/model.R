#------------------------------------------------------------------------------#
#model.R
#
#Author: 
#Sindre Vatnehol
#Institute of marine Research
#Norway
#
#Here we define the paramters to be used.
#Also we run the main model, and several dignostic of this
#
#Input:  data/dat.Rda - the sam.data object produced by data.R
#Output: model/fit.Rda  - the fitted SAM model (used by output.R, report.R, forecast.R)
#        model/diag.Rda - diagnostics: process/obs residuals, simulation,
#                          jitter, leave-out and two retrospective runs
#
#Chapters:
#  1. Setup
#  2. Model configuration & fit
#  3. Diagnostics
#------------------------------------------------------------------------------#



#================================================================================#
# 1. SETUP
#================================================================================#
rm(list=ls())

library(TAF)
library(stockassessment)

mkdir("model")



#================================================================================#
# 2. MODEL CONFIGURATION & FIT
#================================================================================#

#------------------------------------------------------------------------------#
#Prepare the parameter settings as defined from WKBHERMAC (2025)
#
#defcon(dat) creates a default SAM configuration with one row per fleet in
#`dat` (row 1 = residual catch, rows 2-7 = the survey fleets, in the order
#they were added to `surveys` in data.R). For this stock that order is:
#  row 1: Residual catch (no catchability parameter - it's the reference fleet)
#  row 2: IESNS
#  row 3: recBar (BESS)
#  row 4: Bar (IESNS recruitment index)
#  row 5: NASF_old
#  row 6: NASF_new
#  row 7: RFID
#IMPORTANT: this row order is NOT re-derived from `dat` here, it is assumed,
#and there is no practical way around that (SAM's conf$key... matrices are
#plain fleet-indexed matrices with no fleet-name metadata to check against).
#If a fleet is ever added/removed/reordered in data.R, every conf$key...[n,]
#line below silently starts referring to the wrong fleet - there is no
#error, just a wrong model. The same fixed fleet order is also hardcoded a
#second time in output.R/report.R (the logFpar_0, logFpar_9, ... -> fleet
#name lookups). In practice this means: keep the fleet order in data.R's
#`surveys` list stable across years, and if it ever must change, update
#this configuration and the output.R/report.R name-lookups together.
#
#keyLogFpar: which catchability parameter (by age, per fleet-row) is shared
#across ages - e.g. row 2 (IESNS) uses a separate parameter per age up to
#age 10, then one shared parameter for the age-11/12 plus group (the
#repeated "8" at the end); -1 means "no observation/parameter for this age".
#
#keyVarObs: same idea but for observation variance groupings - e.g. row 1
#(catch) groups ages 5-12 under a single variance parameter (repeated "3").
#
#stockRecruitmentModelCode = 3: piecewise-constant stock-recruitment
#relationship (SAM's "no relationship, mean/segment-constant recruitment"
#code), as decided at WKBHERMAC (2025).
#
#fbarRange = c(5,12): the age range (5-12, i.e. including the plus group)
#used to compute the reference fishing mortality FbarW quoted in the advice.
#------------------------------------------------------------------------------#
load('data/dat.Rda')
conf<-defcon(dat)
conf$keyLogFpar[2,]<-c(-1,0,1,2,3,4,5,6,7,8,8)
conf$keyLogFpar[3,]<-c(9,10,-1,-1,-1,-1,-1,-1,-1,-1,-1)
conf$keyLogFpar[4,]<-c(11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)
conf$keyLogFpar[5,]<-c(-1,12,13,14,15,16,17,17,17,17,17)
conf$keyLogFpar[6,]<-c(-1,18,19,20,21,22,23,24,25,26,26)
conf$keyLogFpar[7,-1]<-27
conf$stockRecruitmentModelCode<-3
conf$fbarRange <- c(5,12)
conf$keyVarObs[1,] <- c(0,1,2,3,3,3,3,3,3,3,3)
conf$keyVarObs[2,-1] <- c(3,4,4,4,4,4,4,4,4,4)+1
conf$keyVarObs[3,1:2] <- c(6)
conf$keyVarObs[4,1] <- 7
conf$keyVarObs[5,-1] <- c(8,9,10,10,10,10,10,10,10,10)
conf$keyVarObs[6,-1] <- c(10,11,12,12,12,12,12,12,12,12)+1
conf$keyVarObs[7,-1] <- c(15,15,15,15,15,15,15,16,16,17)-1
#Fleet 5 (NASF_old) observation errors are modelled as AR(1)-correlated
#across ages, rather than independent (the default). This applies to
#NASF_old only, not to any other fleet.
conf$obsCorStruct[5]<-'AR'
conf$keyCorObs[5,-1]<-0
par<-defpar(dat,conf)

fit <- sam.fit(dat,conf,par)
save(fit,file='model/fit.Rda')


#0 = converged; any other value means the optimiser did not reach a proper
#optimum and every downstream result (output/report/forecast) would be
#built on an invalid fit, so we stop the pipeline here rather than continue.
print(paste0('Convergence:  ',fit$opt$convergence))
if(fit$opt$convergence != 0) stop('SAM did not converge (fit$opt$convergence != 0)')

#================================================================================#
# 3. DIAGNOSTICS
#================================================================================#
#fit_pro  = one-observation-ahead process residuals (procres)
#res_fit  = standard one-observation-ahead residuals, per fleet/age/year
#sim_fit  = self-simulation study (50 sims) - refits SAM on simulated data
#           to check the model can recover its own parameters
#jit_fit  = jitter run - refits from perturbed starting values, to check
#           the optimiser converges to the same optimum regardless of start
#LO_fit   = leave-one-fleet-out refits, to see how much each survey fleet
#           influences SSB/F/recruitment
#LO_lastyear = leave-out-last-year-of-one-fleet refits (see below)
#ret/ret2 = two retrospective runs (see below), used e.g. to compute Mohn's
#           rho in report.R/output.R
#------------------------------------------------------------------------------#
print('run procres')
fit_pro <- procres(fit)
print('run residuals')
res_fit <- residuals(fit)
print('run simfit')
sim_fit <- simstudy(fit,nsim = 50)
print('run jit')
jit_fit <- jit(fit)
print('run leavout')
LO_fit <- leaveout(fit)

print('run leave out last year')
#Unlike LO_fit (which drops an entire fleet), this drops only the most
#recent year of ONE fleet at a time and refits, keeping everything else
#untouched - tests how much each fleet's single newest data point moves the
#terminal-year estimates. 'base' = the unmodified fit, so plot.LO()
#(boot/rscript/utils.r) can compare the other runs against it directly.
fleet_names <- attr(dat,'fleetNames')
LO_lastyear <- list(base = fit)
for(f in 2:dat$noFleets){
  dat_tmp <- dat
  fleet_years <- dat_tmp$aux[dat_tmp$aux[,"fleet"]==f,"year"]
  idx <- dat_tmp$aux[,"fleet"]==f & dat_tmp$aux[,"year"]==max(fleet_years)
  dat_tmp$logobs[idx] <- NA
  par_tmp <- defpar(dat_tmp,conf)
  LO_lastyear[[fleet_names[f]]] <- sam.fit(dat_tmp,conf,par_tmp)
}

print('run retro')
#R/R2 are matrices of "final year to retain" per retro peel (rows) and per
#fleet (the 7 columns, same fleet order as conf$keyLogFpar above). Peeling
#back YY-0:4 removes 0 to 4 years of data one at a time (5 peels total).
#R and R2 differ only in the last column (RFID): R peels RFID back by the
#same 1 extra year as most fleets (YY-1)-0:4, R2 peels it back by 2 extra
#years (YY-2)-0:4 - i.e. two alternative assumptions for how much of the
#RFID series to withhold in the retrospective, run and stored side by side
#for comparison (ret vs. ret2).
YY <- max(fit$data$years)
R<-cbind((YY-1)-0:4, YY-0:4, (YY-1)-0:4, YY-0:4, YY-0:4, YY-0:4, (YY-1)-0:4)
ret <- retro(fit, year=R)
R2<-cbind((YY-1)-0:4, YY-0:4, (YY-1)-0:4, YY-0:4, YY-0:4, YY-0:4, (YY-2)-0:4)
ret2 <- retro(fit, year=R2)
save(fit_pro,res_fit,sim_fit,LO_fit,LO_lastyear,ret,ret2,file='model/diag.Rda')




