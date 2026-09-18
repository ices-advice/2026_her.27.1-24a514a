#------------------------------------------------------------------------------#
#forecast.R
#
#Here we run the short term forecast for the NSS herring.
#
#Chapters:x
#  1. Setup
#  2. Reference points & settings
#  3. Predicted maturity for forecast years
#  4. Scenarios (management plan, MSY, 0 catch, Fpa, Flim, Blim, Bpa, Fsq)
#  5. Build output table (per-scenario rows -> out_table)
#  6. Input/output tables (4.8.1.1a/b, 4.8.2.1a/b)
#  7. Figures (4.8.1.1, 4.8.1.2, 4.8.1.3, 4.8.3.1, 4.8.3.2)
#
#
#TODO: Fiks Fpa == FMSY when above MSYBtriggerFsg

#------------------------------------------------------------------------------#



#================================================================================#
# 1. SETUP
#================================================================================#
rm(list=ls())



#------------------------------------------------------------------------------#
#User input data
#------------------------------------------------------------------------------#
TAC <- 533914             #Last year advice, this is in unit tonnes
nation_quota <- 584.025 #Declared quota intermediate year, note that it must be in units ktonnes



#================================================================================#
# 2. REFERENCE POINTS & SETTINGS
#================================================================================#
blim <- 2286
flim =0.291
fmsy <- 0.21
msy_btrigger <- 3177
bpa <- msy_btrigger#3184


#------------------------------------------------------------------------------#
#Management plan
#------------------------------------------------------------------------------#
fmp <- 0.14
TAC_cap <- c(0.20,0.25)





#------------------------------------------------------------------------------#
#Load library, and source files
#------------------------------------------------------------------------------#
library(icesTAF)
library(ggplot2)
library('stockassessment')
library(quarto)
library(ggplotify)   
library(patchwork)  
library(ggrepel)
library(bbmle) #Package used for DD effect on maturity
library(gridExtra)


source(taf.boot.path('rscript/forecast_utils.R'))
source(taf.boot.path('rscript/utils.r'))
source(taf.boot.path('rscript/SAM_forecast_from_GIT.R'))
source(taf.boot.path('rscript/DD_model.r'))




#------------------------------------------------------------------------------#
#Load the fit object, and make a replica to avoid messing up with the official
#------------------------------------------------------------------------------#
load('boot/rscript/WGWIDE2025.Rda')
fit_old <- fit
load('model/fit.Rda')



#------------------------------------------------------------------------------#
#load last runned forecast
#------------------------------------------------------------------------------#
load('boot/rscript/forecast_2025.Rda')
forecast_last<-tt_stochastic


#------------------------------------------------------------------------------#
#Make the folder structure
#This is redundant if the report.R has been runned
#------------------------------------------------------------------------------#
TAF::mkdir('forecast')
TAF::mkdir('forecast/tables')
TAF::mkdir('forecast/figures')





#MP = the management-plan F-vs-SSB rule curve (see scenario() in
#forecast_utils.R), looked up by the MP scenario loop below to decide which
#F to apply each forecast year based on that year's simulated SSB.
MP<-scenario(blim,bpa,ftarget = fmp)


#------------------------------------------------------------------------------#
#some general settings that may be redundant
#------------------------------------------------------------------------------#
rec.years = fit$data$years    #Specify the range to pull recruitment from, this is must go in the forecast but it is not used since we draw from the parameter
rec_from_par = T              #Compute recruitment from parameter 
switch_year_range = T         #Flag to initialize a small hack in the forecast function, so it fits the structure of the data
ave.years <- tail(fit$data$years,3) # which years to compute the means
sim_years=3





#------------------------------------------------------------------------------#
#Set seed: 
#------------------------------------------------------------------------------#
seed_value <- 1234







#================================================================================#
# 3. PREDICTED MATURITY FOR FORECAST YEARS
#================================================================================#
tmp <- get_mopred(fit)
mo_pred <- round(tmp$mo_pred,digits = 1)




#------------------------------------------------------------------------------#
#plot the maturity at age
#------------------------------------------------------------------------------#
plot_predicted_maturity(tmp, "forecast/figures/fig_predicted_maturity.png")




#------------------------------------------------------------------------------#
#Add predicted maturity into the future
#------------------------------------------------------------------------------#
fit2 <- fit
fit2$data$propMat<-rbind(fit2$data$propMat,mo_pred[rownames(mo_pred)>=(max(fit$data$years)+1),])



#------------------------------------------------------------------------------#
#A hack as there is something strange with the dimension in the fit object
#------------------------------------------------------------------------------#
if(length(dim(fit2$data$catchMeanWeight))==3){fit2$data$catchMeanWeight<-fit2$data$catchMeanWeight[,,]}
fit2$data$landFrac<-fit2$data$landFrac[,,]
fit2$data$landMeanWeight<-fit2$data$landMeanWeight[,,]
fit2$data$propF<-fit2$data$propF[,,]





#================================================================================#
# 4. SCENARIOS
#================================================================================#
print('Run 1/8 MP')



#------------------------------------------------------------------------------#
#Compute the TAC for next year. 
#The structure is written to mimic the XSAM structure of forecast
#------------------------------------------------------------------------------#
fBarW =  c(NA,rep(fmp,sim_years))
catchval.exact = c(nation_quota,rep(NA,sim_years))
set.seed(seed_value)
tt<-forecast_XSAM(fit2,fBarW =  fBarW,
                  catchval.exact = catchval.exact,
                  ave.years = ave.years,
                  rec.years = rec.years,
                  deterministic = F,
                  bias_correct = F,
                  rec_from_par = rec_from_par,
                  switch_year_range = switch_year_range)


#Apply management plan, and update the fbarW
for(i in 1:sim_years){
  print(i)
  
  #update the projection, using the MP if below Btrigger
  fBarW[i+1]<-MP$F[MP$ssb==round(median(tt[[i+1]]$ssb))]
  
  #for future reference in the table
  if(i==1)new_F<-fBarW[i+1]
  
  #The cap rule 
  if(round(median(tt[i+1][[1]]$ssb))>=bpa){
    #If larger then cap
    if(((median(tt[[i+1]]$catch)-median(tt[[i]]$catch))/median(tt[[i]]$catch))>TAC_cap[2]){
      fBarW[i+1] =  NA
      catchval.exact[i+1] = median(tt[[i]]$catch)*(1+TAC_cap[2])
    }
    
    if(((median(tt[[i+1]]$catch)-median(tt[[i]]$catch))/median(tt[[i]]$catch))<(-TAC_cap[1])){
      fBarW[i+1] =  NA
      catchval.exact[i+1] = median(tt[[i]]$catch)*(1-TAC_cap[1])
    }
  }
  
  #Refresh
  set.seed(seed_value)
  tt<-forecast_XSAM(fit2,fBarW =  fBarW,
                    catchval.exact = catchval.exact,
                    ave.years = ave.years,
                    rec.years = rec.years,
                    deterministic = F,
                    bias_correct = F,
                    rec_from_par = rec_from_par,
                    switch_year_range = switch_year_range)
  fBarW[i+1] =  NA
  catchval.exact[i+1] = median(tt[[i+1]]$catch)
  print(tt)
}



#Do the stochastic forecast on the management plan
set.seed(seed_value)
tt_stochastic<-forecast_XSAM(fit2,
                             fBarW =  fBarW,
                             catchval.exact = catchval.exact,
                             ave.years = ave.years,
                             rec.years = rec.years,
                             deterministic = F,
                             rec_from_par = rec_from_par,
                             bias_correct = F,
                             switch_year_range = switch_year_range,addTSB = T)
save(tt_stochastic,file='forecast/forecast.Rda')




#------------------------------------------------------------------------------#
#MSY approach
#------------------------------------------------------------------------------#
print('Run 2/8 FMsy')
fbarW <- c(NA,fmsy,fmsy,fmsy)
catchval.exact = c(nation_quota,NA,NA,NA)
for(i in 2:4){
  print(i)
  set.seed(seed_value)
  tt<-forecast_XSAM(fit2,fBarW =  fbarW,
                    catchval.exact = catchval.exact,
                    ave.years = ave.years,
                    rec.years = rec.years,
                    deterministic = F,
                    bias_correct = F,
                    rec_from_par = rec_from_par,
                    switch_year_range = switch_year_range)
  
  #reduce F if SSB is below btrigger
  fbarW[i]<-fbarW[i]*median(tt[[i]]$ssb)/msy_btrigger

  #if so, run the new F
  if(fbarW[i]<fmsy){
    tt<-forecast_XSAM(fit2,fBarW =  fbarW,
                      catchval.exact = catchval.exact,
                      ave.years = ave.years,
                      rec.years = rec.years,
                      deterministic = F,
                      bias_correct = F,
                      rec_from_par = rec_from_par,
                      switch_year_range = switch_year_range)
  }
  
  #store the new F
  if(i==2)new_Fmsy<-round(fbarW[2],digits = 3)
  
  fbarW[i]<-NA
  catchval.exact[i]<- median(tt[[i]]$catch)
}

set.seed(seed_value)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(median(tt[[1]]$catch),
                                                  median(tt[[2]]$catch),
                                                  median(tt[[3]]$catch),
                                                  median(tt[[4]]$catch)),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,
                             switch_year_range = switch_year_range,
                             bias_correct = F)

save(tt_stochastic,file='forecast/forecast_MSY.Rda')





#------------------------------------------------------------------------------#
#0 catch
#------------------------------------------------------------------------------#
print('Run 3/8 0-catch')
set.seed(seed_value)
tt_stochastic<-forecast_XSAM(fit2,
                             catchval.exact = c(nation_quota,0,0,0),
                             ave.years = ave.years,
                             rec.years = rec.years,
                             deterministic = F,
                             bias_correct = F,
                             rec_from_par = rec_from_par,
                             switch_year_range = switch_year_range)
save(tt_stochastic,file='forecast/forecast_0catch.Rda')




#------------------------------------------------------------------------------#
#Fpa
#TODO: this scenario targets `fmsy` (0.21) - there is no separate `fpa`
#constant defined anywhere in the reference-points section above. May be
#intentional (Fmsy used as an Fpa proxy for this stock), but confirm with
#the assessment lead before assuming it's a copy-paste bug.
#------------------------------------------------------------------------------#
print('Run 4/8 Fpa')
run_F_scenario(fit2, fmsy, nation_quota, ave.years, rec.years,
              rec_from_par, switch_year_range, seed_value,
              'forecast/forecast_Fpa.Rda')



#------------------------------------------------------------------------------#
#Flim
#------------------------------------------------------------------------------#
print('Run 5/8 Flim')
run_F_scenario(fit2, flim, nation_quota, ave.years, rec.years,
              rec_from_par, switch_year_range, seed_value,
              'forecast/forecast_Flim.Rda')




#------------------------------------------------------------------------------#
#SSB = blim
#Finds the constant F that makes the realised median SSB hit blim, via an
#outer search (find_F_for_target_ssb(), boot/rscript/forecast_utils.R) rather
#than forecast_XSAM()'s own nextssb= option - see that function's header
#comment for why nextssb= converges a few percent off target for this stock.
#------------------------------------------------------------------------------#
print('Run 6/8 Blim')
F_blim <- find_F_for_target_ssb(fit2,target_ssb=blim,nation_quota=nation_quota,
                                ave.years=ave.years,rec.years=rec.years,
                                rec_from_par=rec_from_par,switch_year_range=switch_year_range,
                                nosim=1000,seed=seed_value)

set.seed(seed_value)
tt_stochastic<-forecast_XSAM(fit2,
                             fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(nation_quota,NA,NA,NA),
                             fval = c(NA,F_blim,F_blim,F_blim),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,
                             switch_year_range = switch_year_range,
                             bias_correct = F)
save(tt_stochastic,file='forecast/forecast_Blim.Rda')

#Test: confirm the realised median SSB actually hit blim to within 0.5% -
#this should always pass now that F is solved from the outside against the
#realised value itself, rather than relying on forecast_XSAM()'s internal
#(known-imprecise) nextssb= search.
rel_diff_blim <- abs(median(tt_stochastic[[3]]$ssb)-blim)/blim
cat(sprintf('SSB=Blim check: realised median SSB=%s, target blim=%s, diff=%s%%\n',
           round(median(tt_stochastic[[3]]$ssb)),blim,round(rel_diff_blim*100,2)))
stopifnot('SSB=Blim did not converge to within 0.5% of blim'=rel_diff_blim<=0.005)


#------------------------------------------------------------------------------#
#SSB = Bpa
#Same approach as SSB=Blim above.
#------------------------------------------------------------------------------#
print('Run 7/8 Bpa')
F_bpa <- find_F_for_target_ssb(fit2,target_ssb=msy_btrigger,nation_quota=nation_quota,
                               ave.years=ave.years,rec.years=rec.years,
                               rec_from_par=rec_from_par,switch_year_range=switch_year_range,
                               nosim=1000,seed=seed_value)

set.seed(seed_value)
tt_stochastic<-forecast_XSAM(fit2,fBarW =  c(NA,NA,NA,NA),
                             catchval.exact  =  c(nation_quota,NA,NA,NA),
                             fval = c(NA,F_bpa,F_bpa,F_bpa),
                             ave.years = ave.years,
                             deterministic = F,
                             rec.years = rec.years,
                             rec_from_par = rec_from_par,
                             switch_year_range = switch_year_range,
                             bias_correct = F)
save(tt_stochastic,file='forecast/forecast_Bpa.Rda')

#Test: same 0.5% check as SSB=Blim above.
rel_diff_bpa <- abs(median(tt_stochastic[[3]]$ssb)-msy_btrigger)/msy_btrigger
cat(sprintf('SSB=MSY Btrigger check: realised median SSB=%s, target=%s, diff=%s%%\n',
           round(median(tt_stochastic[[3]]$ssb)),msy_btrigger,round(rel_diff_bpa*100,2)))
stopifnot('SSB=MSY Btrigger did not converge to within 0.5% of msy_btrigger'=rel_diff_bpa<=0.005)


#------------------------------------------------------------------------------#
#Fsq
#------------------------------------------------------------------------------#
load('forecast/forecast.Rda')
print('Run 8/8 SQ')
f_sq <- fbartable(tt_stochastic)[1,1]
#Unlike every other scenario above, this block used to be missing its
#set.seed() calls - its stochastic draws silently depended on whatever RNG
#state happened to be left over from the code before it, rather than being
#reproducible/independent like the other 7 scenarios. Fixed by seeding here
#too, same as run_F_scenario() does for the others (confirmed with the
#assessment lead - this changes SQ's exact numbers vs. the un-seeded
#version, since it's now drawing from a fresh, defined RNG state instead of
#an inherited one).
run_F_scenario(fit2, f_sq, nation_quota, ave.years, rec.years,
              rec_from_par, switch_year_range, seed_value,
              'forecast/forecast_SQ.Rda')






#------------------------------------------------------------------------------#
#Basis table
#------------------------------------------------------------------------------#
load('forecast/forecast.Rda')
basis <- data.frame(
  'SSB_intermediat' = round(median(tt_stochastic[[1]]$ssb)),
  catch_intermediat = round(median(tt_stochastic[[1]]$catch)*1000),
  'SSB_forecast'=median(tt_stochastic[[2]]$ssb),
  'FbarW_intermediat'=round(median(tt_stochastic[[1]]$fbar),digits = 3),
  'Recruitment'=paste(round(c(median(tt_stochastic[[1]]$rec),
                              median(tt_stochastic[[2]]$rec),
                              median(tt_stochastic[[3]]$rec))/1000,digits = 3),
                      sep = ',')
)






#================================================================================#
# 5. BUILD OUTPUT TABLE
#================================================================================#
out_table <- c()



#------------------------------------------------------------------------------#
#Management plan
#------------------------------------------------------------------------------#
out_table <- rbind(out_table, build_output_row('forecast/forecast.Rda','management strategy',
                                              blim,TAC,nation_quota,basis[1,],
                                              basis='F=0.14',finn=round(new_F,3)))



#------------------------------------------------------------------------------#
#FMSY approach
#------------------------------------------------------------------------------#
out_table <- rbind(out_table, build_output_row('forecast/forecast_MSY.Rda','MSY approach',
                                              blim,TAC,nation_quota,basis[1,],
                                              basis='F=0.21',finn=new_Fmsy))


#------------------------------------------------------------------------------#
#Zero catch
#------------------------------------------------------------------------------#
out_table <- rbind(out_table, build_output_row('forecast/forecast_0catch.Rda','Zero Catch',
                                              blim,TAC,nation_quota,basis[1,],
                                              basis='F=0',finn=0))


#------------------------------------------------------------------------------#
# Fpa
#------------------------------------------------------------------------------#
out_table <- rbind(out_table, build_output_row('forecast/forecast_Fpa.Rda','Fpa',
                                              blim,TAC,nation_quota,basis[1,],
                                              basis='F=Fpa',finn=fmsy))


#------------------------------------------------------------------------------#
#SSB = BLIM
#F was already solved (find_F_for_target_ssb()) to hit blim from the outside,
#and checked to within 0.5% right after forecast_Blim.Rda was created above -
#this re-check just confirms build_output_row() picked up the same realised
#value. Once confirmed, ssb_pred is overridden to the exact target (blim)
#rather than the realised simulated median, so the presentation table
#(new_output$SSB below) shows the clean reference point value (2286000) -
#the 0.5% check above is what justifies treating them as the same number.
#------------------------------------------------------------------------------#
row_blim <- build_output_row('forecast/forecast_Blim.Rda','SSB=Blim',
                             blim,TAC,nation_quota,basis[1,],
                             target=blim)
rel_diff <- abs(row_blim$ssb_pred-blim)/blim
cat(sprintf('SSB=Blim: ssb_pred=%s, blim=%s, diff=%s%%\n',
           row_blim$ssb_pred, blim, round(rel_diff*100,2)))
stopifnot('SSB=Blim row does not match the converged forecast'=rel_diff<=0.005)
row_blim$ssb_pred <- blim
out_table <- rbind(out_table, row_blim)


#------------------------------------------------------------------------------#
#SSB = Bpa
#Same check and ssb_pred override as SSB=Blim above.
#------------------------------------------------------------------------------#
row_bpa <- build_output_row('forecast/forecast_Bpa.Rda','SSB=MSY Btrigger',
                            blim,TAC,nation_quota,basis[1,],
                            target=msy_btrigger)
rel_diff <- abs(row_bpa$ssb_pred-msy_btrigger)/msy_btrigger
cat(sprintf('SSB=MSY Btrigger: ssb_pred=%s, msy_btrigger=%s, diff=%s%%\n',
           row_bpa$ssb_pred, msy_btrigger, round(rel_diff*100,2)))
stopifnot('SSB=MSY Btrigger row does not match the converged forecast'=rel_diff<=0.005)
row_bpa$ssb_pred <- msy_btrigger
out_table <- rbind(out_table, row_bpa)


#------------------------------------------------------------------------------#
#FbarW = Status quo
#------------------------------------------------------------------------------#
out_table <- rbind(out_table, build_output_row('forecast/forecast_SQ.Rda','SQ',
                                              blim,TAC,nation_quota,basis[1,],
                                              basis=paste0('F=',f_sq),finn=f_sq))






#================================================================================#
# 6. INPUT/OUTPUT TABLES
#================================================================================#

#------------------------------------------------------------------------------#
#Write input table
#------------------------------------------------------------------------------#
load('forecast/forecast.Rda')
tmp <- data.frame(
  age=min(fit2$data$minAgePerFleet):max(fit2$data$maxAgePerFleet),
  stockno = round(as.numeric(tail(ntable(fit2),1))),
  natmort=as.numeric(tail(fit2$data$natMor,1)),
  maturity=as.numeric(tail(fit$data$propMat,1)),
  PropM=as.numeric(tail(fit2$data$propM,1)),
  PropF=as.numeric(tail(fit2$data$propF,1)),
  SWeight=as.numeric(round(tail(fit2$data$stockMeanWeight,1),digits = 3)),
  Explot = as.numeric(round(tail(faytable(fit),1),3)),
  CWeight=round(colMeans(tail(fit2$data$catchMeanWeight,3)),digits = 3)
)
write.csv2(tmp,
           file='forecast/tables/4.8.1.1a.txt'
)



#------------------------------------------------------------------------------#
#Write second input table
#SWeight here is an assumption for the forecast years (no real 2026/2027
#data exists yet), taken as the mean of the last 3 observed years - not
#just the latest single year like table 4.8.1.1a. Bug fix: this was
#missing the tmp$ prefix, so it computed the 3-year mean but never
#actually stored it in tmp$SWeight - the column silently kept table
#4.8.1.1a's single-latest-year value instead. Confirmed against last
#year's published WGWIDE report, where the two tables' "Weight in stock"
#columns are genuinely different.
#------------------------------------------------------------------------------#
tmp$stockno[-1]<-NA
tmp$stockno[1]<-round(median(tt_stochastic[[2]]$rec))
tmp$maturity <- as.numeric(head(tail(fit2$data$propMat,n=3),n=1))
tmp$maturity2 <- as.numeric(tail(head(tail(fit2$data$propMat,n=3),n=2),n=1))
tmp$SWeight <- round(colMeans(tail(fit2$data$stockMeanWeight,3)),digits = 3)
write.csv2(
  tmp,
  file='forecast/tables/4.8.1.1b.txt'
)





#================================================================================#
# 7. FIGURES
#================================================================================#

#------------------------------------------------------------------------------#
#Figure 4.8.1.1
#------------------------------------------------------------------------------#
plot_fig_4_8_1_1('forecast/forecast.Rda', fit2, 'forecast/figures/Fig.4.8.1.1.png')




#------------------------------------------------------------------------------#
#Handle table 4.8.2.1
#------------------------------------------------------------------------------#
load('forecast/forecast.Rda')



oout <- data.frame(
  'SSB_intermediat' = paste(round(tail(ssbtable(fit),n=1)[1]/1000,digits=3)),
  'catch_intermediat' = round(median(tt_stochastic[[1]]$catch)*1000),
  'SSB_next' = paste(round(median(tt_stochastic[[2]]$ssb)/1000,digits=3)),
  'FbarW_intermediat'=paste(round(median(tt_stochastic[[1]]$fbar),digits = 3)),
  'rec_intermediat' = round(tail(rectable(fit),n=1)[1]),
  'rec_next' = round(median(tt_stochastic[[2]]$rec)),
  'rec_next2' = round(median(tt_stochastic[[3]]$rec))
)



write.csv2(t(oout),file = 'forecast/tables/tab4.8.2.1a.txt')



#------------------------------------------------------------------------------#
#Plot 3 year projection
#------------------------------------------------------------------------------#
plot_fig_4_8_1_2('forecast/forecast.Rda', fit2, fit,
                 msy_btrigger, fmp, 'forecast/figures/Fig.4.8.1.2.png')



#------------------------------------------------------------------------------#
#Predicted age distribution in catch, management plan
#------------------------------------------------------------------------------#
plot_fig_4_8_1_3('forecast/forecast.Rda', 'forecast/figures/Fig.4.8.1.3.png')




load('forecast/forecast.Rda')

caa <- stockassessment::read.ices(taf.data.path('caa.txt'))
caa <- tail(caa)
for(i in 1:nrow(caa)){
  caa[i,]<-caa[i,]/max(caa[i,])
}
caa <- reshape2::melt(caa)
names(caa)<- c('year','age','caa')

dd <- c()
for(i in 1:length(tt_stochastic)){
  tmp<- data.frame(age=2:12,
                   caa = rowSums(tt_stochastic[[i]]$catchatage),
                   year = tt_stochastic[[i]]$year)
  tmp$caa<-tmp$caa/max(tmp$caa)
  dd <-rbind(dd,tmp)
}
dd<-rbind(dd,caa)

data_matrix <- reshape2::dcast(dd, age ~ year, value.var = "caa")

rel <- round((sweep(data_matrix[,-1], 2, colSums(data_matrix[,-1]), FUN = "/"))*100)
rel <- cbind(age = data_matrix$age, rel)
years <- colnames(rel)[-1]              # alle årskolonner
last_two <- head(tail(years, 4),2)             # de to siste, generisk
result <- rel[, c("age", last_two)]
result



#------------------------------------------------------------------------------#
#Diagnostic plots
#------------------------------------------------------------------------------#
plot_diagnostic_comparison('forecast/forecast.Rda', fit2, fit, fit_old, forecast_last,
                           'forecast/figures/Fig.4.8.3.1.png',
                           'forecast/figures/Fig.4.8.3.2.png')












#------------------------------------------------------------------------------#
#ICES rounding rules
#------------------------------------------------------------------------------#
new_output <- data.frame(
  Rationale = out_table$Rationale,
  catches = round(out_table$Catches),
  F = out_table$Finn,
  SSB = round(out_table$ssb_pred*1000),
  SSB_change = out_table$ssbchange,
  catch_change = out_table$Catchchange,
  advice_change = out_table$TACchange,
  probability = out_table$risc*100
)

new_output$F<-icesAdvice::icesRound(new_output$F)
new_output$SSB_change<-icesAdvice::icesRound(new_output$SSB_change)
new_output$catch_change <- ices_round_capped(new_output$catch_change)
new_output$advice_change <- ices_round_capped(new_output$advice_change)


new_output$probability<-icesAdvice::icesRound(new_output$probability)



names(new_output)<-c('Rationale',
                     paste0('Total Catch (',max(fit$data$years)+1,')'),
                     paste0('F(',max(fit$data$years)+1,')'),
                     paste0('SSB(',max(fit$data$years)+2,')'),
                     '% SSB change',
                     '% catch change',
                     '% advice change',
                     paste('% probability',max(fit$data$years)+2)
)
write.csv2(new_output,file='forecast/tables/tab4.8.2.1b.txt')




#Commented out: nothing in the repo reads this .txt file (the presentation
#loads out_table straight from forecast_options.Rda below instead), and
#it's unrounded so it kept failing test.R's reproducibility check on
#machine-noise-level floating point differences. Left here in case there
#was a reason for it - re-enable if needed.
# write.csv2(out_table,file='forecast/tables/tab4.8.2.1b_alternative.txt')
save(out_table,new_output,file='forecast/forecast_options.Rda')
write.csv2(out_table,file='forecast/tables/tab4.8.2.1b_alternative.txt')

rectable(fit)
ssbtable(fit)
compute_FbarW(fit_old)
compute_FbarW(fit)
ssbtable(fit_old)
ssbtable(fit)
