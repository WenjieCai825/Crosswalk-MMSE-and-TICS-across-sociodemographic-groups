Data cleaning_DropBox
================
Wenjie Cai
2025-08-14

``` r
library(haven)
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(readxl)
```

``` r
#change it to your address of the dataset folder
data_folder <- "/Users/wenjiecai/Library/CloudStorage/Dropbox/Crosswalk_paper_2_group_heterogeneity/datasets/"
```

``` r
data <- read_dta((file.path(data_folder, "HCAPHRS.dta")))
```

``` r
df<-data[!is.na(data$HCAP16WGTR),]
```

``` r
# race
table(df$RaceAndEthnicity) #1 is white
df <- df %>% mutate (White=ifelse(RaceAndEthnicity==1, 1,
                           ifelse(is.na(RaceAndEthnicity), NA,
                                 0)))

#dementia classification
table(df$vs1hcapdx) #BPV version
df <- df %>% mutate (dementiaBPV=ifelse(vs1hcapdx==3, 1,
                           ifelse(is.na(vs1hcapdx), NA,
                                 0)))
table(df$vs1hcapdxeap) #EAP version
df <- df %>% mutate (dementiaEAP=ifelse(vs1hcapdxeap==3, 1,
                           ifelse(is.na(vs1hcapdxeap), NA,
                                 0)))
table(df$Hudomiet_classification) #derived by Dr. Jones, a normal/cind/dementia three-level variable where a respondent is classified according to their highest class probability.
df <- df %>% mutate (dementiaProb=ifelse(Hudomiet_classification==3, 1,
                           ifelse(is.na(Hudomiet_classification), NA,
                                 0)))
table(df$cogfunction2016) #Langa
df <- df %>% mutate (dementiaLanga=ifelse(cogfunction2016==3, 1,
                           ifelse(is.na(cogfunction2016), NA,
                                 0)))

#education (0 is < high school)
table(df$Educational_Attainment)
df <- df %>% mutate (education=ifelse(Educational_Attainment==1, 0,
                           ifelse(is.na(Educational_Attainment), NA,
                                 1)))
#age (>75 is 1)
table(df$rage_cat)
df <- df %>% mutate (ageGroup=ifelse(rage_cat<3, 0,
                           ifelse(is.na(rage_cat), NA,
                                 1)))
```

``` r
# dataHurd<- read_dta("/Users/wenjiecai/Downloads/DementiaPredictedProbabilities/pdem_withvarnames.dta")

dataPower<-read_sas((file.path(data_folder, "hrsdementia_2024_1217.sas7bdat"))) %>% select(HHID, PN, expert_dem,HRS_year, expert_p) %>% filter(HRS_year==2016)

dataLanga<-read_dta((file.path(data_folder, "cogfinalimp_9520wide.dta")))  %>% select(hhid, pn, cogfunction2016, cogtot27_imp2016, fimrc_imp2016,fdlrc_imp2016,fser7_imp2016,fbwc20_imp2016,imrc_imp2016,dlrc_imp2016,ser7_imp2016,bwc20_imp2016 )

dataHumomiet<-read_dta((file.path(data_folder, "Dementia_HRS_2000-2016_Basic_Release1_2m.dta")))  %>% select (hhidpn, Cog)
```

``` r
df$HHIDPN<-paste(df$HHID, df$PN, sep="")
dataPower$HHIDPN<-paste(dataPower$HHID,dataPower$PN, sep="")

dataLanga$HHIDPN<-paste(dataLanga$hhid,dataLanga$pn,sep="")
df<- df %>% left_join(dataLanga, by="HHIDPN")

df2<-df %>% left_join(dataPower, by="HHIDPN") 
```

``` r
#ageGroup
df2$PAGEGroup<-case_when(df2$PAGE<70 ~ 1,
                         df2$PAGE<75 & df2$PAGE>69 ~ 2,
                         df2$PAGE<80 & df2$PAGE>74 ~ 3,
                         df2$PAGE>79 ~ 4)   #Age (years) (65-69, 70-74, 75-79, 80> years) 
df2$PAGEGroup <- factor(df2$PAGEGroup,
                        levels = 1:4,
                        labels = c("65–69", "70–74", "75–79", "80+"))
table(df2$PAGEGroup, useNA = "ifany")

#sex
table(df2$female, useNA = "ifany") 

#RaceEthnicity
df2$RaceEthnicity<-case_when(df2$black==0 & df2$hisp==0 ~ 1,  # Non-hispanic white
                             df2$black==1 & df2$hisp==0 ~ 2,  # Non-hispanic black
                             df2$hisp==1 ~ 3) # Hispanic
df2$RaceEthnicity <- factor (df2$RaceEthnicity, 
                             levels = 1:3, 
                             labels = c("Non-hispanic white", "Non-hispanic black", "Hispanic"))
table(df2$RaceEthnicity, useNA = "ifany")

#Education level
df2$Educational_Attainment <- factor (df2$Educational_Attainment, 
                                       levels = 1:4, 
                                      labels = c("<high school", "High school", "Some college", "Education beyond college"))
table(df2$Educational_Attainment, useNA = "ifany") # 1, 2, 3, 4: <high school, high school, some college, education beyond college

#Marital status
df2$Marital<-df2$PMARST
df2$Marital <- factor (df2$Marital, 
                       levels = 1:5,
                       labels = c("Married", "Separated/divorced", "Widowed", "Never married", "Unknown"))
table(df2$Marital, useNA = "ifany") #marital status: 1, 2, 3, 4, 5: married, separated/divorced, widowed, never married, unknown marital status
```

``` r
df2<- df2 %>% select(HHIDPN, White, dementiaBPV, dementiaEAP, dementiaProb, dementiaLanga, education, ageGroup, expert_dem, expert_p, Cog, cogtot27_imp2016, fimrc_imp2016,fdlrc_imp2016,fser7_imp2016,fbwc20_imp2016,imrc_imp2016,dlrc_imp2016,ser7_imp2016,bwc20_imp2016, PAGEGroup, female, RaceEthnicity, Educational_Attainment, Marital )  # Cog is Humomiet, expert_p and expert_dem are Power
```

``` r
dfC<-read_dta((file.path(data_folder, "HRS core and HCAP data merged-2.dta")))
dfC$HHIDPN<-paste(dfC$hhid,dfC$pn,sep="")
df2<- df2 %>% left_join(dfC, by="HHIDPN")

df2<-df2%>% rename("cog_27"="rand_cogtot2016", #TICS 27
                     "cog_35"="r_cogtot2016") #TICS 35
```

``` r
dfWu<-read_sas((file.path(data_folder, "cogvarsred_gdr_20230508.sas7bdat"))) %>% select(HHIDPN, memimp13)
df2$HHIDPN <- as.numeric(df2$HHIDPN)
df2<-df2 %>% left_join(dfWu, by="HHIDPN")
df2<-df2%>% rename("cog_20"="r_tr202016") 
```

``` r
dfHCAP<-read_dta((file.path(data_folder, "PsyMCA_2025_RAND_HRS_HCAP_raw.dta")))
dfHCAP$HHIDPN<-dfHCAP$hhidpn
df2<-df2 %>% left_join(dfHCAP, by="HHIDPN")

df22<-df2
names(df22) <- gsub("\\.x$", "", names(df22))
names(df22) <- gsub("\\.y$", "", names(df22))
df22 <- df22[, !duplicated(names(df22))]
```

``` r
# codebook for HCAP: https://hrs.isr.umich.edu/sites/default/files/meta/hcap/2016/codebook/hc16hp_ri.htm 
dfHCAP2<-read_dta((file.path(data_folder, "hc16hp_r.dta")))
dfHCAP2$HHIDPN<-paste(dfHCAP2$hhid,dfHCAP2$pn,sep="")
names(dfHCAP2) <- gsub("_", "", names(dfHCAP2))
dfHCAP2<- dfHCAP2 |> select(R1ORIENTTIME, RIORIENTP1, RIORIENTP2, RIORIENTP3, RIORIENTP3, RIORIENTP4, RIORIENTP5, R1ORIENTPLACE, R1REPEAT, HHIDPN)
dfHCAP2$HHIDPN <- as.numeric(dfHCAP2$HHIDPN)
df22<-df22 %>% left_join(dfHCAP2, by="HHIDPN")

names(df22) <- gsub("\\.x$", "", names(df22))
names(df22) <- gsub("\\.y$", "", names(df22))
df22 <- df22[, !duplicated(names(df22))]
```

``` r
dfapoe<-read_dta((file.path(data_folder, "apoe_serotonin_release.dta"))) |> select(apoe, hhid, pn)
dfapoe$HHIDPN<-paste(dfapoe$hhid,dfapoe$pn,sep="")
dfapoe$HHIDPN <- as.numeric(dfapoe$HHIDPN)
df22<-df22 %>% left_join(dfapoe, by="HHIDPN")
```

``` r
#workforce
table(df22$r_lbrf2016) # 1: works full time , 2: works parttime, 3: unemployed, 4: partly retired 5: retired 6: disabled 7: not in labor force
df22$WorkForce <- case_when(
  df22$r_lbrf2016 %in% c(1, 2, 4) ~ 1,
  df22$r_lbrf2016 %in% c(3, 5, 6, 7) ~ 0
)
df22$WorkForce <- factor (df22$WorkForce, 
                          levels = 0:1, 
                           labels = c("Not working", "Working"))

table(df22$WorkForce, useNA = "ifany") # 0 is not working, 1 is working




#household income
summary(df22$h_itot2016) #total household income (respondent + spouse only)
df22$QHouseholdIncome<- cut(df22$h_itot2016,
                          breaks = quantile(df22$h_itot2016, probs = seq(0, 1, 0.25), na.rm = TRUE),
                          include.lowest = TRUE,
                          labels = c("Q1", "Q2", "Q3", "Q4"))
table(df22$QHouseholdIncome, useNA = "ifany")


# summary(df22$h_itot2016) #total household income (respondent + spouse only)
# df22$QHouseholdIncome <- cut(
#   df22$h_itot2016,
#   breaks = quantile(df22$h_itot2016, probs = c(0, 0.5, 1), na.rm = TRUE),
#   include.lowest = TRUE,
#   labels = c("Low", "High")
# )
# table(df22$QHouseholdIncome, useNA = "ifany")


#APOE category
df22$APOEGroup<- case_when(df22$apoe %in% c(22, 23, 33) ~ 0, 
                           df22$apoe %in% c(24, 34, 44) ~ 1, ) 

df22$APOEGroup <- factor (df22$APOEGroup, 
                          levels=0:1,
                          labels = c("APOEe4 non-carrier", "APOEe4 carrier"))
table(df22$APOEGroup, useNA = "ifany")

#Urbanicity

df22$rural <- factor (df22$rural, 
                      levels=0:1, 
                      labels = c("Urban/suburban","Rural"))
table(df22$rural, useNA = "ifany") 
```

``` r
df22<-df22|> select(HHIDPN, White, dementiaBPV, dementiaEAP,dementiaProb,dementiaLanga,education,ageGroup,expert_dem,expert_p,Cog,black,hisp, cog_20,memimp13,cog_27,cog_35,vs5gcp_eap, vs2memsc, vs2exfsc, vs2lflsc, vs2vissc, vs2vdori1, r1mmse_score, r1tics_score, r1word_total, r1verbal_score, r1lc_score, r1bc_score, r1csid_score, r1word_dscore, r1bm_immscore, r1wlrec_totscore, r1cp_score, r1dig_score, r1cpdel_score, r1bm_delscore, r1lmb_recoscore, r1ns_score, r1rv_score, r1tma_score, r1tmb_score, r1smell_score, pmarst,workforce2016, educ_d, race,age2016, female, mar2016, rural, income2016, assets2016, cogtot27_imp2016, fimrc_imp2016,fdlrc_imp2016,fser7_imp2016,fbwc20_imp2016,imrc_imp2016,dlrc_imp2016,ser7_imp2016,bwc20_imp2016,R1ORIENTTIME, RIORIENTP1, RIORIENTP2, RIORIENTP3, RIORIENTP3, RIORIENTP4, RIORIENTP5, R1ORIENTPLACE, R1REPEAT,r_bwc202016, s_imrc2016,s_ser72016,r_scis2016, r_cact2016,s_pres2016, s_vp2016,r_dlrc2016,apoe, PAGEGroup, female, RaceEthnicity, Educational_Attainment, Marital, WorkForce, QHouseholdIncome, h_itot2016,APOEGroup,rural)

write_dta(df22, (file.path(data_folder, "/Derived/HRS-HCAPcrosswalkdata.dta")))
```

``` r
#coding for the items 
df22$date <- df22$R1ORIENTTIME # range from 0-5
df22 <- df22 %>%
  mutate(across(RIORIENTP1:RIORIENTP5, ~ case_when(. > 1 ~ NA, TRUE ~ .)))
# df22$address1<- df22$RIORIENTP5*3 # the fifth question * 3points "what is this address?" 
df22$address2<- rowSums(df22[, c("RIORIENTP1", "RIORIENTP2","RIORIENTP3")]==1, na.rm = TRUE)# the first, second and third questions  - According to Rich
df22$backCounting<- case_when(df22$r_bwc202016==2 ~ 2, # correct at the 1st try, 2 points
                             df22$r_bwc202016==1 ~ 1, # correct at the 2nd try, 1 point
                             df22$r_bwc202016==0 ~ 0)# incorrect, 0 point
df22$wordlist<-df22$imrc_imp2016 # range 0-10, same as the immediate word recall

df22$subtraction<-df22$ser7_imp2016

df22$responsive_naming<- rowSums(df22[, c("r_scis2016", "r_cact2016")]==1, na.rm = TRUE) # scissors, cactus naming, range from 0-2
df22$president<- rowSums(df22[, c("s_pres2016", "s_vp2016")]==1, na.rm = TRUE) #president, vice-president, range from 0-2
df22$delayedRecall<-df22$r_dlrc2016 # range 0-10
df22$Repetition<-case_when(df22$R1REPEAT==1 ~ 1, # correct
                           df22$R1REPEAT==5 ~ 0)  #incorrect

# cog_20 from RAND
# cog_27 from RAND
# cog_35 from RAND
# cog_30
df22$cog_30 <- rowSums(df22[, c("R1ORIENTTIME", "address2", "backCounting", "wordlist",
                                "subtraction", "responsive_naming", "president",
                                "Repetition")], na.rm = FALSE)
# cog_40 
df22$cog_40<- rowSums(df22[, c("R1ORIENTTIME", "address2", "backCounting", "wordlist",
                                "subtraction", "responsive_naming", "president",
                                "delayedRecall", "Repetition")], na.rm = FALSE)
```

``` r
write_dta(df22, (file.path(data_folder, "/Derived/HRS-HCAPcrosswalkdata.dta")))
saveRDS(df22, file.path(data_folder, "Derived", "HRS-HCAPcrosswalkdata.rds"))
```
