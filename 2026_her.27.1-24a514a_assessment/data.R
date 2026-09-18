#------------------------------------------------------------------------------#
#data.R
#
#Author: Sindre Vatnehol, Institute of Marine Research, Norway
#
#Reads and prepares all input data for the SAM assessment.
#
#Input:  raw text files fetched by boot.R, via taf.data.path()
#Output: data/dat.Rda         - sam.data object consumed by model.R
#        data/tab.4.4.8.*.txt - CV-weighting tables used in report.R
#        data/<parameter>.txt - every processed data source feeding into
#                                dat.Rda, written for test.R's reproducibility
#                                check (chapter 6)
#
#Chapters:
#  1. Setup
#  2. Catch data
#  3. Biological data (fractions, natural mortality, maturity, stock weight)
#  4. Survey/observation data
#  5. Weighting of data sources
#  6. Write out all data sources
#  7. Assemble SAM data object
#------------------------------------------------------------------------------#



#================================================================================#
# 1. SETUP
#================================================================================#
rm(list=ls())

library(TAF)
library('stockassessment')
source(taf.boot.path('rscript/utils.r'))
source(taf.boot.path('rscript/DD_model.r'))
source(taf.boot.path('rscript/forecast_utils.R'))
library(bbmle)

mkdir("data")

#Copy the input files into data/ so the exact files used in this run stay
#traceable/inspectable, separate from the raw bootstrap copies in boot/initial/data/.
file.copy(list.files(taf.data.path(), full.names = TRUE), "data", overwrite = TRUE)



#================================================================================#
# 2. CATCH DATA
#================================================================================#
#cn = catch-at-age (numbers), cn_cv = its CV (used below for weighting),
#cw = catch mean weight-at-age. Data responsible: Sigurvin Bjarnason (cn, cw),
#Are Salthaug (cn_cv).
cn <- read.ices(taf.data.path('caa.txt'))
cn_cv <- read.ices(taf.data.path('caa_cv.txt'))
cw <- read.ices(taf.data.path('cw.txt'))



#================================================================================#
# 3. BIOLOGICAL DATA
#================================================================================#

#lf = landing fraction; pf/pm = timing of F/M relative to spawning (used by
#SAM to time-correct mortality when computing SSB).
lf<-read.ices(taf.data.path('lf.txt'))
pf<-read.ices(taf.data.path('pf.txt'))
pm<-read.ices(taf.data.path('pm.txt'))

#NSS herring has no reported discards - dw/lw are placeholders only, needed
#because setup.sam.data() requires them structurally.
dw<-read.ices(taf.data.path('cw.txt'))
lw<-read.ices(taf.data.path('cw.txt'))

#Natural mortality, assumed fixed by age.
nm<-read.ices(taf.data.path('nm.txt'))


#--- Maturity -------------------------------------------------------------
#matprop.txt is trimmed to SAM's year/age range (1988+, ages 2-12 plus
#group). Three alternative maturity series are built below - pick which one
#is used by swapping the `prop.mature=` argument in setup.sam.data() further
#down. All three are kept at 1 decimal precision, matching matprop.txt's own
#resolution:
#
#Note: It was decided in WGWIDE 2026 to continiue the old timeseries and approach
#until the next evaluation of reference points. 
#The script underneath are kept for future reference. 
#
#  mo_original - matprop.txt as-is, i.e. Case 1/"status quo" (WKBHERMAC's
#                subjective curve for the most recent, back-calculation-
#                lagged years, taken unmodified). Used by default.
#  mo_case1    - Case 2: mo_original, with the last 4 years replaced
#                outright by the density-dependent model prediction
#                (mo_pred) - those years in mo_original are manually
#                assigned, not estimated, so they get replaced regardless
#                of whether prop_mat_new.txt happens to have a value there.
#  mo_case3    - Case 3: fully recalculated from the new maturity dataset
#                (prop_mat_new.txt) instead of matprop.txt, with any still-
#                missing years filled from a prediction. Because this case
#                swaps in a different maturity dataset entirely, the
#                density-dependent model is also refit on THAT data (not on
#                WGWIDE2025's original propMat) before predicting the
#                missing years - otherwise the prediction would still be
#                based on the very data this case is meant to replace.
#Mirrors Case 1/2/3 in WD_update_maturity.R - see that script for the full
#background/comparison. Data responsible: Erling Kåre Stenevik.
mo <- read.ices(taf.data.path('matprop.txt'))
mo <- mo[as.numeric(rownames(mo))>=1988,]
mo <- mo[,as.numeric(colnames(mo))%in%c(2:12)]
mo_original <- round(mo,digits = 1)

# 
# 
# load('boot/rscript/WGWIDE2025.Rda')
# fit_old <- fit
# rm(fit)
# 
# #sim_years must cover the gap between WGWIDE2025's terminal year and the
# #newest year in mo_original, otherwise get_mopred() has no predicted row for
# #the missing recent year(s) and the NA-fill below is a no-op for those years.
# sim_years <- max(as.numeric(rownames(mo_original))) - max(fit_old$data$years)
# 
# #--- Case 2 (mo_case1) ------------------------------------------------------
# mo_pred <- as.matrix(get_mopred(fit_old,digits = 1)$mo_pred)   #get_mopred() returns a data.frame; assigning a data.frame slice into a matrix below would silently turn it into a list
# replace_years <- tail(rownames(mo_original),4)
# stopifnot(all(replace_years%in%rownames(mo_pred)))
# 
# mo_case1 <- mo_original
# mo_case1[replace_years,] <- mo_pred[replace_years,colnames(mo_original)]
# 
# #--- Case 3 (mo_case3) -------------------------------------------------------
# mo_new <- round(read.ices(taf.data.path('prop_mat_new.txt')),digits = 1)
# mo_new <- mo_new[rownames(mo_original),colnames(mo_original)]   #align to mo_original; errors out if a year/age is missing
# 
# #Refit the density-dependent model on the new maturity data: a copy of
# #fit_old with $data$propMat swapped to mo_new (restricted to the
# #years fit_old actually has, so dimensions still match) before
# #calling get_mopred() on it.
# fit_new_data <- fit_old
# fit_new_data$data$propMat <- mo_new[rownames(mo_new)%in%rownames(fit_old$data$propMat),colnames(fit_old$data$propMat)]
# mo_pred_case3 <- as.matrix(get_mopred(fit_new_data,digits = 1)$mo_pred)
# 
# pred_years <- intersect(rownames(mo_new),rownames(mo_pred_case3))
# na_idx <- is.na(mo_new[pred_years,,drop=FALSE])   #years/ages not covered by prop_mat_new.txt
# 
# mo_case3 <- mo_new
# mo_case3[pred_years,][na_idx] <- mo_pred_case3[pred_years,colnames(mo_original),drop=FALSE][na_idx]


#--- Stock weight -----------------------------------------------------------
#west.txt trimmed the same way as matprop.txt; ages 12-15 are averaged into
#the age-12 plus group before older ages are dropped. Data responsible:
#Erling Kåre Stenevik.
sw<-read.ices(taf.data.path('west.txt'))
sw <- sw[as.numeric(rownames(sw))>=1988,]
sw<-sw[,as.numeric(colnames(sw))>=2]
sw[,as.numeric(colnames(sw))==12]<-rowMeans(sw[,as.numeric(colnames(sw))>=12])
sw <- sw[,as.numeric(colnames(sw))<=12]



#================================================================================#
# 4. SURVEY DATA
#================================================================================#
#One multi-fleet ICES-format file holds IESNS, recBar/Bar (BESS-derived
#recruitment index) and NASF (split into NASF_old/NASF_new below). RFID is
#read separately further down. Data responsible: Erling Kåre Stenevik,
#Georg Skaret, Are Salthaug.
surveys<-read.ices(taf.data.path('survey_T2025-02-05.txt'))
surveys.cv<-read.ices(taf.data.path('surveycv_T2025-02-05.txt'))


#Recalculated BESS (recBar) index, overriding the values read from
#survey_T2025-02-05.txt above, since the official index for this year was
#not yet available when this was written.
# surveys$recBar[1:21,2]<-c(36.495,16.167,5.535,2.595,1.626,0.433,0.315,1.504,1.078,
#                       5.029,1.822,9.023,1.573,2.138,6.035,0.209,0.231,0.12,0.882,
#                       32.92,43.995)
# surveys$recBar[1:21,3]<-c(0.901,6.973,1.62,6.378,3.987,1.807,0.234,0.006,1.285,0.092,
#                       6.825,3.214,3.089,3.465,1.299,6,1.816,0.36,NA,4.443,22.109)

#BESS (recBar) is set to NA in 2016, 2018 and 2022 due to known survey
#design/coverage issues in those years, and excluded from the fit.
surveys$recBar[rownames(surveys$recBar)%in% c(2016,2018,2022),]<-NA
surveys.cv$recBar[rownames(surveys.cv$recBar)%in% c(2016,2018,2022),]<-NA

#BESS only informs ages 2-3.
att_ <- attributes(surveys$recBar)
surveys$recBar<-surveys$recBar[,colnames(surveys$recBar)%in%c(2,3)]
surveys.cv$recBar<-surveys.cv$recBar[,colnames(surveys.cv$recBar)%in%c(2,3)]
attributes(surveys$recBar)$time <- att_$time
attributes(surveys$recBar)$twofirst <- att_$twofirst


#--- RFID -------------------------------------------------------------------
#Tag mark-recapture index (BioSampleID), not part of the multi-fleet file
#above. Data responsible: Aril Slotte.
index_tag1 <- read.ices(taf.data.path('rfid_index_BioSampleID.txt'))$RFID
index_tag1.cv <- read.ices(taf.data.path('rfid_cv_BioSampleID.txt'))$RFID


#--- NASF split ---------------------------------------------------------------
#NASF_old (<=2008) and NASF_new (>=2015) get separate catchabilities, since
#survey methodology changed between the two periods (no NASF data
#2009-2014). time=0.16 is the survey timing (~mid-Feb). Decision from
#WKBHERMAC (2025).
surveys$NASF_old <- surveys$NASF[as.numeric(rownames(surveys$NASF))<=2008,]
surveys$NASF_new <- surveys$NASF[as.numeric(rownames(surveys$NASF))>=2015,]
surveys.cv$NASF_old <- surveys.cv$NASF[as.numeric(rownames(surveys.cv$NASF))<=2008,]
surveys.cv$NASF_new <- surveys.cv$NASF[as.numeric(rownames(surveys.cv$NASF))>=2015,]
attributes(surveys$NASF_old)$time <- c(0.16,0.16)
attributes(surveys$NASF_new)$time <- c(0.16,0.16)
surveys$NASF<-NULL


#Add RFID as its own fleet. The first age column is dropped, since the RFID
#tag mark-recapture method does not produce a usable index for that age.
#time=0.9 corresponds to the November recapture window.
surveys$RFID <- index_tag1[,-1]
attributes(surveys$RFID)$time <- c(0.9,0.9)


#================================================================================#
# 5. WEIGHTING OF DATA SOURCES
#================================================================================#
#computeExternalWeight() (boot/rscript/utils.r) fits a Taylor mean-variance
#relationship per fleet and returns smoothed inverse-variance weights for
#SAM, rather than using the raw per-year CV directly. The tab.4.4.8.*.txt
#files written below hold the resulting *smoothed* CV (1/sqrt(weight)) and
#feed report.R's tables 4.4.8.1-5: 4.4.8.1=catch, 4.4.8.2a/b=NASF_old/new,
#4.4.8.3=IESNS bar, 4.4.8.4=IESNS, 4.4.8.5=RFID.
ttt<- computeExternalWeight(mu_data = cn,cv_data = cn_cv)
attributes(cn)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.1.txt')

ttt<- computeExternalWeight(mu_data = surveys$NASF_old,cv_data = surveys.cv$NASF_old)
attributes(surveys$NASF_old)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.2a.txt')

ttt<- computeExternalWeight(mu_data = surveys$NASF_new,cv_data = surveys.cv$NASF_new)
attributes(surveys$NASF_new)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.2b.txt')

ttt<- computeExternalWeight(mu_data = surveys$IESNS,cv_data = surveys.cv$IESNS)
attributes(surveys$IESNS)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.4.txt')

ttt<- computeExternalWeight(mu_data = surveys$Bar,cv_data = surveys.cv$Bar)
attributes(surveys$Bar)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.3.txt')

#recBar (BESS) weight is computed for use in setup.sam.data() below, but -
#unlike the other fleets - has no corresponding report table written here.
ttt<- computeExternalWeight(mu_data = surveys$recBar,cv_data = surveys.cv$recBar)
attributes(surveys$recBar)$weight = ttt$weight

ttt<- computeExternalWeight(mu_data = surveys$RFID,cv_data = index_tag1.cv[,-1])
attributes(surveys$RFID)$weight = ttt$weight
write.csv(round(1/sqrt(ttt$weight),digits = 3),file = 'data/tab.4.4.8.5.txt')



#================================================================================#
# 6. WRITE OUT ALL DATA SOURCES (for test.R reproducibility check)
#================================================================================#
#The raw ICES-format files are already copied verbatim into data/ (chapter 1),
#so those are covered by test.R as-is. This writes out every object that
#actually feeds into setup.sam.data() below - i.e. after all trimming/
#splitting/NA-masking above - so processing bugs are caught too, not just
#changes to the raw inputs. Named after the parameter, except where a
#report table (tab.4.4.8.*) already exists above for that source.
#------------------------------------------------------------------------------#
data_sources <- list(
  cn       = cn,
  cn_cv    = cn_cv,
  cw       = cw,
  lf       = lf,
  pf       = pf,
  pm       = pm,
  nm       = nm,
  dw       = dw,
  lw       = lw,
  mo       = mo,
  sw       = sw,
  IESNS       = surveys$IESNS,
  IESNS_cv    = surveys.cv$IESNS,
  recBar      = surveys$recBar,
  recBar_cv   = surveys.cv$recBar,
  Bar         = surveys$Bar,
  Bar_cv      = surveys.cv$Bar,
  NASF_old    = surveys$NASF_old,
  NASF_old_cv = surveys.cv$NASF_old,
  NASF_new    = surveys$NASF_new,
  NASF_new_cv = surveys.cv$NASF_new,
  RFID        = surveys$RFID,
  RFID_cv     = index_tag1.cv[,-1]
)
for(source_name in names(data_sources)){
  write.csv(data_sources[[source_name]], file = paste0('data/', source_name, '.txt'))
}



#================================================================================#
# 7. ASSEMBLE SAM DATA OBJECT
#================================================================================#

#Sanity check: taf.data.path() (boot/data/) should be an exact copy of
#boot/initial/data/ after taf.boot(). Fails loudly here if it isn't, rather
#than silently building on stale/wrong input.
stopifnot(identical(read.ices(taf.data.path('caa.txt')),read.ices('boot/initial/data/caa.txt')))
stopifnot(identical(read.ices(taf.data.path('caa_cv.txt')),read.ices('boot/initial/data/caa_cv.txt')))
stopifnot(identical(read.ices(taf.data.path('cw.txt')),read.ices('boot/initial/data/cw.txt')))
stopifnot(identical(read.ices(taf.data.path('lf.txt')),read.ices('boot/initial/data/lf.txt')))
stopifnot(identical(read.ices(taf.data.path('pf.txt')),read.ices('boot/initial/data/pf.txt')))
stopifnot(identical(read.ices(taf.data.path('pm.txt')),read.ices('boot/initial/data/pm.txt')))
stopifnot(identical(read.ices(taf.data.path('dw.txt')),read.ices('boot/initial/data/dw.txt')))
stopifnot(identical(read.ices(taf.data.path('lw.txt')),read.ices('boot/initial/data/lw.txt')))
stopifnot(identical(read.ices(taf.data.path('nm.txt')),read.ices('boot/initial/data/nm.txt')))
stopifnot(identical(read.ices(taf.data.path('matprop.txt')),read.ices('boot/initial/data/matprop.txt')))
stopifnot(identical(read.ices(taf.data.path('west.txt')),read.ices('boot/initial/data/west.txt')))
stopifnot(identical(read.ices(taf.data.path('survey_T2025-02-05.txt')),read.ices('boot/initial/data/survey_T2025-02-05.txt')))
stopifnot(identical(read.ices(taf.data.path('surveycv_T2025-02-05.txt')), read.ices('boot/initial/data/surveycv_T2025-02-05.txt')))
stopifnot(identical(read.ices(taf.data.path('rfid_index_BioSampleID.txt')),read.ices('boot/initial/data/rfid_index_BioSampleID.txt')))
stopifnot(identical(read.ices(taf.data.path('rfid_cv_BioSampleID.txt')),  read.ices('boot/initial/data/rfid_cv_BioSampleID.txt')))


#Final assembly: combines all series read/prepared above into the
#"sam.data" object SAM expects. This is the only object model.R depends on.
#prop.mature: swap in mo_case1 or mo_case3 here to run with an alternative
#maturity series instead of the official mo_original (see chapter 3 above).
dat<-setup.sam.data(surveys=surveys,
                    residual.fleet = cn,
                    prop.mature = mo,
                    stock.mean.weight = sw,
                    catch.mean.weight = cw,
                    dis.mean.weight = dw,
                    land.mean.weight = lw,
                    prop.f=pf,
                    prop.m = pm,
                    natural.mortality = nm,
                    land.frac = lf)
save(dat,file='data/dat.Rda')



