library(dplyr)
library(ggplot2)
library(cogxwalkr)
library(table1)
library(haven)
#library(patchwork)

setwd("//ars-data-01.sde.net.ucsf.edu/myresearchshared/mglymour_uk_biobank_shared/Jingxuan/PsyMCA2")


##################### 
#####readRDS function failed to directly read rds file from dropbox
#####you need to move the rds data to a local folder and use the path of this dataset from the local folder: change the path here. 
# pathWholeData<-"/Users/jwang30/Downloads/HRS-HCAPcrosswalkdata.rds"
df <- read_dta('HRS-HCAPcrosswalkdata.dta')

#####################
# due to the same issue of the readRDS function described above, we need to set the paths for the use in the crosswalk process, change them to your results paths accordingly
# results20path<-"/Users/wenjiecai/Documents/research/Conference/2025 psyMCA/workgroup/R codes/noboot_r1mmse_score_cog_20.RDS"
# 
# results27path<-"/Users/wenjiecai/Documents/research/Conference/2025 psyMCA/workgroup/R codes/noboot_r1mmse_score_cog_27.RDS"
# 
# results30path<-"/Users/wenjiecai/Documents/research/Conference/2025 psyMCA/workgroup/R codes/noboot_r1mmse_score_cog_30.RDS"
# 
# results35path<-"/Users/wenjiecai/Documents/research/Conference/2025 psyMCA/workgroup/R codes/noboot_r1mmse_score_cog_35.RDS"
# 
# results40path<-"/Users/wenjiecai/Documents/research/Conference/2025 psyMCA/workgroup/R codes/noboot_r1mmse_score_cog_40.RDS"

# outcomes: r1mmse_score, cog_20, cog_27, cog_30, cog_35, cog_40
#missingness in each outcome and conditional variable
sum(is.na(df$r1mmse_score))
sum(is.na(df$cog_20))
sum(is.na(df$cog_27))
sum(is.na(df$cog_30))
sum(is.na(df$cog_35))
sum(is.na(df$cog_40))
sum(is.na(df$dementiaBPV))

sum(complete.cases(df[, c("r1mmse_score", "cog_20", "cog_27", "cog_30", "cog_35", "cog_40")])) #use the 3236 for the whole sample analysis

df1<- df[complete.cases(df[, c("r1mmse_score", "cog_20", "cog_27", "cog_30", "cog_35", "cog_40")]),]


va <- c("PAGEGroup", "female", "RaceEthnicity", 
        "Educational_Attainment", "Marital", 
        "WorkForce", "QHouseholdIncome", 
        "APOEGroup", "rural")
df1[va] <- lapply(df1[va], as.factor)
df1$female<-as.factor(df1$female)
table1(~ PAGEGroup + female + RaceEthnicity + Educational_Attainment + Marital + WorkForce + QHouseholdIncome + APOEGroup + rural, 
       data=df1,
       render = function(missing=TRUE, ...) {
         table1::render.default(missing=TRUE, ...)
       } )



outcomes<- c("r1mmse_score", "cog_20", "cog_27", "cog_30", "cog_35", "cog_40")
# Get all pairs
pairs = t(combn(outcomes, 2))  # each row is (A, B)
pair_names  =  apply(pairs, 1, function(x) paste(x[1], x[2], sep = "_VS_"))
df$dementia<-df$dementiaBPV
# Loop and run crosswalk
for (i in seq_len(nrow(pairs))) {
  start_time = Sys.time()
  
  message("Running crosswalk pair", i, ": ", pairs[i, 1], " => ", pairs[i, 2])
  
  var1 = pairs[i, 1]
  var2 = pairs[i, 2]
  dementia_var = 'dementia'
  niter = 5000
  
  # Keep only the variables of interest
  data.temp  =  df %>% dplyr::select(all_of(c(var1, var2, dementia_var)))
  
  # Remove rows with missing values
  data.temp  =  data.temp %>% dplyr::filter(complete.cases(.))
  
  boot_settings  =  list(nboot = 1000, seed = 999, ncores = 3)
  cw  =  crosswalk(
    cog1 = var1,
    cog2 = var2,
    data = data.temp,
    control = boot_settings,
    condition_by = dementia_var, 
    niter = niter
  )
  
  res = list(cog1 = var1,
             cog2 = var2,
             conditional_var = dementia_var,
             cw.diffs = cw$diffs,
             cw.boot = cw$boot)
  # cw.noboot = cw$boot 
  
  end_time = Sys.time()
  elapsed = round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  
  message("Finished: ", pairs[i, 1], " => ", pairs[i, 2], 
          " | Time used: ", elapsed, " seconds")
  
  # print(paste0('boot_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  saveRDS(res,file=paste0('boot_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
}




df1$dementia<-df1$dementiaBPV

# stratified
for (i in seq_len(nrow(pairs))) {
  start_time = Sys.time()
  
  message("Running crosswalk pair", i, ": ", pairs[i, 1], " => ", pairs[i, 2])
  
  var1 = pairs[i, 1]
  var2 = pairs[i, 2]
  dementia_var = 'dementia'
  niter = 5000
  
  # Keep only the variables of interest
  data.temp  =  df1 %>% dplyr::select(all_of(c(var1, var2, dementia_var,va)))
  
  # Remove rows with missing values
  # data.temp  =  data.temp %>% dplyr::filter(complete.cases(.))
  
  data.temp.age.1 = data.temp %>% filter(PAGEGroup==1)
  data.temp.age.2 = data.temp %>% filter(PAGEGroup==2)
  data.temp.age.3 = data.temp %>% filter(PAGEGroup==3)
  data.temp.age.4 = data.temp %>% filter(PAGEGroup==4)
  
  data.temp.sex.1 = data.temp %>% filter(female==0)
  data.temp.sex.2 = data.temp %>% filter(female==1)
  
  data.temp.race.1 = data.temp %>% filter(RaceEthnicity==1)
  data.temp.race.2 = data.temp %>% filter(RaceEthnicity==2)
  data.temp.race.3 = data.temp %>% filter(RaceEthnicity==3)
  
  data.temp.educ.1 = data.temp %>% filter(Educational_Attainment==1)
  data.temp.educ.2 = data.temp %>% filter(Educational_Attainment==2)
  data.temp.educ.3 = data.temp %>% filter(Educational_Attainment==3)
  data.temp.educ.4 = data.temp %>% filter(Educational_Attainment==4)
  
  data.temp.marital.1 = data.temp %>% filter(Marital==1)
  data.temp.marital.2 = data.temp %>% filter(Marital==2)
  data.temp.marital.3 = data.temp %>% filter(Marital==3)
  data.temp.marital.4 = data.temp %>% filter(Marital==4)
  
  data.temp.work.1 = data.temp %>% filter(WorkForce==1)
  data.temp.work.2 = data.temp %>% filter(WorkForce==2)
  
  data.temp.income.1 = data.temp %>% filter(QHouseholdIncome %in% c(1,2))
  data.temp.income.2 = data.temp %>% filter(QHouseholdIncome %in% c(3,4))
  
  data.temp.apoe.1 = data.temp %>% filter(APOEGroup==1)
  data.temp.apoe.2 = data.temp %>% filter(APOEGroup==2)
  
  data.temp.rural.1 = data.temp %>% filter(rural==1)
  data.temp.rural.2 = data.temp %>% filter(rural==2)
  
  
  boot_settings  =  list(nboot = 1000, seed = 999, ncores = 3)
  cw.age.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.age.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.age.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.age.1$diffs,cw.boot = cw.age.1$boot)
  saveRDS(res.age.1,file=paste0('boot_','age1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.age.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.age.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.age.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.age.2$diffs,cw.boot = cw.age.2$boot)
  saveRDS(res.age.2,file=paste0('boot_','age2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.age.3  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.age.3,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.age.3 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.age.3$diffs,cw.boot = cw.age.3$boot)
  saveRDS(res.age.3,file=paste0('boot_','age3','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.age.4  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.age.4,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.age.4 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.age.4$diffs,cw.boot = cw.age.4$boot)
  saveRDS(res.age.4,file=paste0('boot_','age4','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.sex.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.sex.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.sex.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.sex.1$diffs,cw.boot = cw.sex.1$boot)
  saveRDS(res.sex.1,file=paste0('boot_','sex1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.sex.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.sex.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.sex.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.sex.2$diffs,cw.boot = cw.sex.2$boot)
  saveRDS(res.sex.2,file=paste0('boot_','sex2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.race.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.race.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.race.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.race.1$diffs,cw.boot = cw.race.1$boot)
  saveRDS(res.race.1,file=paste0('boot_','race1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.race.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.race.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.race.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.race.2$diffs,cw.boot = cw.race.2$boot)
  saveRDS(res.race.2,file=paste0('boot_','race2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.race.3  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.race.3,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.race.3 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.race.3$diffs,cw.boot = cw.race.3$boot)
  saveRDS(res.race.3,file=paste0('boot_','race3','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.educ.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.educ.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.educ.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.educ.1$diffs,cw.boot = cw.educ.1$boot)
  saveRDS(res.educ.1,file=paste0('boot_','educ1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.educ.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.educ.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.educ.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.educ.2$diffs,cw.boot = cw.educ.2$boot)
  saveRDS(res.educ.2,file=paste0('boot_','educ2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.educ.3  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.educ.3,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.educ.3 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.educ.3$diffs,cw.boot = cw.educ.3$boot)
  saveRDS(res.educ.3,file=paste0('boot_','educ3','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  
  cw.educ.4  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.educ.4,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.educ.4 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.educ.4$diffs,cw.boot = cw.educ.4$boot)
  saveRDS(res.educ.4,file=paste0('boot_','educ4','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  
  
  cw.marital.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.marital.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.marital.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.marital.1$diffs,cw.boot = cw.marital.1$boot)
  saveRDS(res.marital.1,file=paste0('boot_','marital1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.marital.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.marital.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.marital.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.marital.2$diffs,cw.boot = cw.marital.2$boot)
  saveRDS(res.marital.2,file=paste0('boot_','marital2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.marital.3  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.marital.3,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.marital.3 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.marital.3$diffs,cw.boot = cw.marital.3$boot)
  saveRDS(res.marital.3,file=paste0('boot_','marital3','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.marital.4  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.marital.4,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.marital.4 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.marital.4$diffs,cw.boot = cw.marital.4$boot)
  saveRDS(res.marital.4,file=paste0('boot_','marital4','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.work.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.work.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.work.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.work.1$diffs,cw.boot = cw.work.1$boot)
  saveRDS(cw.work.1,file=paste0('boot_','work1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.work.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.work.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.work.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.work.2$diffs,cw.boot = cw.work.2$boot)
  saveRDS(cw.work.2,file=paste0('boot_','work2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.income.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.income.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.income.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.income.1$diffs,cw.boot = cw.income.1$boot)
  saveRDS(res.income.1,file=paste0('boot_','income1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.income.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.income.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.income.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.income.2$diffs,cw.boot = cw.income.2$boot)
  saveRDS(res.income.2,file=paste0('boot_','income2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.apoe.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.apoe.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.apoe.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.apoe.1$diffs,cw.boot = cw.apoe.1$boot)
  saveRDS(res.apoe.1,file=paste0('boot_','apoe1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.apoe.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.apoe.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.apoe.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.apoe.2$diffs,cw.boot = cw.apoe.2$boot)
  saveRDS(res.apoe.2,file=paste0('boot_','apoe2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.rural.1  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.rural.1,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.rural.1 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.rural.1$diffs,cw.boot = cw.rural.1$boot)
  saveRDS(res.rural.1,file=paste0('boot_','rural1','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  cw.rural.2  =  crosswalk(cog1 = var1,cog2 = var2,data = data.temp.rural.2,control = boot_settings,condition_by = dementia_var, niter = niter)
  res.rural.2 = list(cog1 = var1,cog2 = var2,conditional_var = dementia_var,cw.diffs = cw.rural.2$diffs,cw.boot = cw.rural.2$boot)
  saveRDS(res.rural.2,file=paste0('boot_','rural2','_',pairs[i, 1], "_", pairs[i, 2],'.RDS'))
  
  
  end_time = Sys.time()
  elapsed = round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  
  message("Finished: ", pairs[i, 1], " => ", pairs[i, 2], 
          " | Time used: ", elapsed, " seconds")
  
}





