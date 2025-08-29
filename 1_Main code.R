########## Data analysis and plotting code for Smooth hammerhead vertebrae SIA manuscript

# Input: 
# - either "VertebraeSIA_RAW" and "Vertebrae_age_reading" data for reproducing all calculations
# - or "VertebraeSIA_Aged" data and already compiled "SEVb" and "Overlap_summary" data to only reproduce manuscript plots with already processed data

# Output:
# - Processed and filtered SIA data with estimated corresponding ages
# - SEVb
# - SEVb overlap
# - Figures and summary tables used in the manuscript and supplementary material

# Set-up #####
## Workspace and libraries
setwd("~/SCIENCE/Vertebrae SIA/R/Publication")

library(dplyr)
library(tidyverse)
library(lubridate)
library(readxl)
library(ggplot2)
library(GGally)
library(ggpubr)
library(maps)
library(mgcv)
library(gratia)
library(plotly)
library(rstatix)
library(vegan)
library(clipr)
library(Hmisc) 
library(heplots)
library(mvnormtest)
library(sf)


## load functions ####
source("Skinner et al. 2019 sppl._ellipsoidcode_functions_final_EDIT.r") #loads functions from Skinner et al. 2019 for SEV, SEVc and SEVb estimation

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Function to calculate age using total length based on VBGF with parameters from Rosa et al. 2017
calculate_age_TL <- function(Lt) {
  Linf <- 288.2  # Asymptotic length
  L0 <- 52.4     # Length at age 0
  k <- 0.09      # Growth coefficient
  
  # Calculate age t
  t <- (-1 / k) * log((Linf - Lt*0.78) / (Linf - L0))
  
  return(t)
}

# Function to calculate total length using age based on VBGF with parameters from Rosa et al. 2017
calculate_TL <- function(t) {
  Linf <- 288.2  # Asymptotic length
  L0 <- 52.4     # Length at age 0
  k <- 0.09      # Growth coefficient
  
  # Calculate age t
  Lt <- (Linf - exp(-k * t) * (Linf - L0)) / 0.78
  
  return(Lt)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# Function to calculate percentiles
percentile <- function(y, level = 0.95) {
  quantile(y, probs = c((1 - level) / 2, 1 - (1 - level) / 2))
}

# Data import #####

## import RAW data
#RAW <- read.csv("VertebraeSIA_210907.csv")
SIA <- read.csv("VertebraeSIA_Aged.csv")
sevb <- read.csv("SEVb.csv")

#over_list <- read.csv("Overlaps.csv", header=F)
over_summary <- read.csv("Overlap_summary.csv", header=T)
rownames(over_summary) <- over_summary[,1]
over_summary[,1] <- NULL


# Data transformation ####

## Age estimation of sampled sections ######
###### !!! not necessary if data "xxxx_Aged.csv" is loaded !!!
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## calculate ages for VC sections in new combined result df
SIA <- RAW
SIA$Agemin <- NA
SIA$Agemax <- NA

# calculate age at minimum drill radius
for(i in 1:nrow(SIA)){
  Fish <- filter(Age, FishCode==SIA[i, 2]) #create subset of Age data for corresponding FishID
  if (SIA$VR[i]==0 || (SIA$VR[i]-SIA$Drill.size..mm.[i]/2)<Fish$DA_mm[1]){  # VR=0 has negative drill minimum if subtracted normally so stays 0, other sections of age 0 also filtered
    SIA$Agemin[i]=0
  } else {
    row=1
    while ((SIA$VR[i]-SIA$Drill.size..mm.[i]/2)>Fish$DA_mm[row]){  #
      print(row)
      SIA$Agemin[i]=Fish$Age[row]
      row= row+1
    }
  }
}
# calculate age at maximum drill radius
for(i in 1:nrow(SIA)){
  Fish <- filter(Age, FishCode==SIA[i, 2]) #create subset of Age data for corresponding FishID
  if ((SIA$VR[i]-SIA$Drill.size..mm.[i]/2)<Fish$DA_mm[1]){  # VR=0 has negative drill minimum if subtracted normally so stays 0
    SIA$Agemax[i]=0
  } else {
    row=1
    while ((SIA$VR[i]+SIA$Drill.size..mm.[i]/2)>Fish$DA_mm[row] || Fish$DA_mm[row] == Fish$DA_mm[nrow(Fish)] ) {
      print(Fish$Age[row])
      SIA$Agemax[i]=Fish$Age[row]
      row= row+1
      if (row>nrow(Fish)){  #make sure loop breaks when reaching end of Fish data frame - not create NAs
        break
      }
    }
  }  
}
## calculate mean age and age span
SIA$Agespan <- SIA$Agemax-SIA$Agemin
SIA$Agemean <- (SIA$Agemin+SIA$Agemax)/2

## Export aged SIA data
#write.csv(SIA, "VertebraeSIA_Aged.csv")
#SIA <- read.csv("VertebraeSIA_Aged.csv")
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# CN ratio above 4 excluded
SIA <- filter(SIA, CNratio < 4)

#SIA$Sex <- ordered(SIA$Sex, levels=c("M", "F"))

# one-hot encoding of Sex to achieve un-ordered variable for GAMM
# Use model.matrix to create the one-hot encoding
one_hot <- model.matrix(~Sex-1, data=SIA)
# Convert the matrix to a data frame
one_hot_df <- as.data.frame(one_hot)
# To merge this back into original data frame
SIA <- cbind(SIA, one_hot_df)
SIA$SexM <- as.integer(SIA$SexM)
SIA$SexF <- as.integer(SIA$SexF)

SIA$FishCode <- as.factor(SIA$FishCode)
SIA$FishID <- as.factor(SIA$FishID)

# SIA data modification
SIA$Section <- as.character(SIA$Section)
SIA$Date <- as.Date(SIA$Date, format="%d/%m/%Y")
SIA$jDate <- julian(SIA$Date)
SIA$Agemean[(SIA$Section=="0" ) & (SIA$Agemean==0)] <- -1.5  # Define pre-birth samples as -1.5 and second pre-birth layer as -1
SIA$Agemean[(SIA$Section=="1") & (SIA$Agemean==0)] <- -1
SIA$Agemean[(SIA$Agemean==0.5)] <- 0
SIA$Agemin[(SIA$Agemean==-1.5)] <- -2
SIA$Agemax[(SIA$Agemean==-1.5)] <- -1
SIA$Agemin[(SIA$Agemean==-0.5)] <- -1
SIA$Agemax[(SIA$Agemean==-0.5)] <- 0

# Define lifestages 
ls <- c(0, # age at birth
        6.02, # age at lowest maturity estimate (193.7 cm TL - Nava Nava et al. 2014)
        11.8) # age at highest maturity estimate (265 cm TL - Last and Stevens 1994)

# Update Lifestage based on new estimates above
SIA <- SIA %>%
  mutate(Lifestage = case_when(
    Agemean <= ls[1] ~ "Pre-birth",
    Agemean < ls[2] ~ "Juvenile",
    Agemean < ls[3] ~ "Subadult",
    Agemean > ls[3] ~ "Adult",
    TRUE ~ Lifestage  # Default case to handle any unexpected values
  ))

SIA$Lifestage <- ordered(SIA$Lifestage, levels=c("Pre-birth", "Juvenile", "Subadult", "Adult"))

#SIA$Age <- SIA$Agemean+1.5 # centering time on first point

# Calculating mean dates of tissue deposition
SIA$tissue_date <- SIA$Date - (months(SIA$FinalAge*12)-months(SIA$Agemean*12)) # calculating age of deposited tissue and subtracting it from capture date

# Classifying capture latitude for plotting
SIA <- SIA %>%
  mutate(lat_class = case_when(
    Lat < -5  ~ "S",
    Lat >= -5 & Lat <= 5  ~ "Eq",
    Lat > 5 & Long > -30  ~ "NE",
    Lat > 5 & Long < -30  ~ "NW"
    ),
    lat_class = ordered(lat_class, levels=c("S", "Eq", "NW", "NE"))
  )

## Model subset #####
SIAM <- SIA %>%
  dplyr::select(ID, d13C, d15N, d34S, VR, Agemean, CalcTL, Sex, SexM, SexF, jDate, Lat, Long)
SIAM$Sex <- ordered(SIAM$Sex, levels=c("M", "F"))

### remove outliers for GAMMs
### C 
# Calculate z-scores for specified columns
z_scores <- scale(SIAM[, "d13C"])
# Define a threshold for outliers (e.g., 3 standard deviations)
threshold <- 3
# Create a boolean mask for outliers
outliers_mask <- apply(z_scores, 2, function(x) abs(x) > threshold)
# Remove rows with outliers
C_SIAM <- SIAM[!rowSums(outliers_mask), ]    # -> 1 Outlier removed

### N
# Calculate z-scores for specified columns
z_scores <- scale(SIAM[, "d15N"])
# Create a boolean mask for outliers
outliers_mask <- apply(z_scores, 2, function(x) abs(x) > threshold)
# Remove rows with outliers
N_SIAM <- SIAM[!rowSums(outliers_mask), ]    # -> 1 Outlier removed

### S
# Calculate z-scores for specified columns
z_scores <- scale(SIAM[, "d34S"])
# Create a boolean mask for outliers
outliers_mask <- apply(z_scores, 2, function(x) abs(x) > threshold)
# Remove rows with outliers
S_SIAM <- SIAM[!rowSums(outliers_mask), ]     # -> 0 Outliers removed



## Group Subsets ####
Edge <- SIA %>%  # Subset of oldest edge layer for each individual
  group_by(FishID)%>%
  top_n(1, Section)

Center <- SIA %>% # Subset of center layer
  filter(Section == 0)

Top <- SIA %>% # Subset of large individuals with more than 4 measurements
  filter(ID < 7)

Rest <- SIA %>% # Subset of all individuals with 4 measurements
  filter(ID > 6)

sp_ids <- c(5:7, 21)
Special <- SIA %>% # subset of particular IDs that showed different trajectories
  filter(ID %in% sp_ids)

not_special <- SIA %>% # subset of particular IDs that showed different trajectories
  filter(ID %!in% sp_ids)

SIAonly <- SIA%>% # Subset containing only life stage and isotopes
  select(Lifestage, Sex, d13C, d15N, d34S)


#### Lifestage * Sex subsets

# per sex
f <- SIA %>%
  filter(Sex=="F")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

m <- SIA %>%
  filter(Sex=="M")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

#per life stage
pb <- SIA %>%
  filter(Lifestage=="Pre-birth")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

juv <- SIA %>%
  filter(Lifestage=="Juvenile")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

sub <- SIA %>%
  filter(Lifestage=="Subadult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

ad <- SIA %>%
  filter(Lifestage=="Adult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

# per sex and life stage
fpb <- SIA %>%
  filter(Sex=="F", Lifestage=="Pre-birth")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

fjuv <- SIA %>%
  filter(Sex=="F", Lifestage=="Juvenile")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

fsub <- SIA %>%
  filter(Sex=="F", Lifestage=="Subadult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

fad <- SIA %>%
  filter(Sex=="F", Lifestage=="Adult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

mpb <- SIA %>%
  filter(Sex=="M", Lifestage=="Pre-birth")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

mjuv <- SIA %>%
  filter(Sex=="M", Lifestage=="Juvenile")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

msub <- SIA %>%
  filter(Sex=="M", Lifestage=="Subadult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

mad <- SIA %>%
  filter(Sex=="M", Lifestage=="Adult")%>%
  select(d13C, d15N, d34S)%>%
  na.omit()

# Standard ellipse volume calculations #####
## calculate SEVc ####

# using modified functions from Skinner et al. 2019

SEVc <- c(
  EstimateSEVc(f),
  EstimateSEVc(m),
  EstimateSEVc(pb),
  EstimateSEVc(juv),
  EstimateSEVc(sub),
  EstimateSEVc(ad),
  EstimateSEVc(fpb),
  EstimateSEVc(fjuv),
  EstimateSEVc(fsub),
  EstimateSEVc(fad),
  EstimateSEVc(mpb),
  EstimateSEVc(mjuv),
  EstimateSEVc(msub),
  EstimateSEVc(mad)
)

SEVc

## SEVb ####
### SEVb Bayesian determination of covariance and means of n-dimensional dataset and SEVb estimation ########
set.seed(12345)


## currently 100,000

## Fit multivariate normal distribution for groups
f_fit <- fitMVNdirect(f)
m_fit <- fitMVNdirect(m)
pb_fit <- fitMVNdirect(pb)
juv_fit <- fitMVNdirect(juv)
sub_fit <- fitMVNdirect(sub)
ad_fit <- fitMVNdirect(ad)
fpb_fit <- fitMVNdirect(fpb)
fjuv_fit <- fitMVNdirect(fjuv)
fsub_fit <- fitMVNdirect(fsub)
fad_fit <- fitMVNdirect(fad)
mpb_fit <- fitMVNdirect(mpb)
mjuv_fit <- fitMVNdirect(mjuv)
msub_fit <- fitMVNdirect(msub)
mad_fit <- fitMVNdirect(mad)

## calculate SEVb (40% coverage)
f_el <- fit3dEllipsoid(f_fit, vol.only=T, p.interval=0.4)
m_el <- fit3dEllipsoid(m_fit, vol.only=T, p.interval=0.4)
pb_el <- fit3dEllipsoid(pb_fit, vol.only=T, p.interval=0.4)
juv_el <- fit3dEllipsoid(juv_fit, vol.only=T, p.interval=0.4)
sub_el <- fit3dEllipsoid(sub_fit, vol.only=T, p.interval=0.4)
ad_el <- fit3dEllipsoid(ad_fit, vol.only=T, p.interval=0.4)
fpb_el <- fit3dEllipsoid(fpb_fit, vol.only=T, p.interval=0.4)
fjuv_el <- fit3dEllipsoid(fjuv_fit, vol.only=T, p.interval=0.4)
fsub_el <- fit3dEllipsoid(fsub_fit, vol.only=T, p.interval=0.4)
fad_el <- fit3dEllipsoid(fad_fit, vol.only=T, p.interval=0.4)
mpb_el <- fit3dEllipsoid(mpb_fit, vol.only=T, p.interval=0.4)
mjuv_el <- fit3dEllipsoid(mjuv_fit, vol.only=T, p.interval=0.4)
msub_el <- fit3dEllipsoid(msub_fit, vol.only=T, p.interval=0.4)
mad_el <- fit3dEllipsoid(mad_fit, vol.only=T, p.interval=0.4)

# code to compile sevb value iterations and save as csv
# NOT NECESSARY IF sevb loaded 

sevb <- data.frame(c(rep("F_ALL", 12000),
                     rep("M_ALL", 12000),
                     rep("ALL_PB", 12000),
                     rep("ALL_JUV", 12000),
                     rep("ALL_SUB", 12000),
                     rep("ALL_AD", 12000),
                     rep("F_PB", 12000),
                     rep("F_JUV", 12000),
                     rep("F_SUB", 12000),
                     rep("F_AD", 12000),
                     rep("M_PB", 12000),
                     rep("M_JUV", 12000),
                     rep("M_SUB", 12000),
                     rep("M_AD", 12000)),
                   c(f_el,
                     m_el,
                     pb_el,
                     juv_el,
                     sub_el,
                     ad_el,
                     fpb_el,
                     fjuv_el,
                     fsub_el,
                     fad_el,
                     mpb_el,
                     mjuv_el,
                     msub_el,
                     mad_el))
names(sevb) <- c("Group", "Volume")
write_csv(sevb, "SEVb.csv")


## SEVb Overlap #####

### automate compilation ####
## NOT NECESSARY IF "Overlap_summary" or "Overlap_all_long2.csv" loaded

overlap_all = data.frame(matrix(nrow = 0, ncol = 6)) 
colnames(overlap_all) <- c("Overlap_vol", "SEVb_vol", "Overlap_prob", "Group", "With", "Between")
#
row <- 3
col <- 2
#
for (row in 2:15){
  a <- over_list[row,1]
  for (col in 2:14){
    b <- over_list[1,col]
    ab <- over_list[row, col]
    if (ab!="") {
      system.time(overlap <- ellipsoid.Overlap(eval(parse(text=paste0(a, "_fit"))), eval(parse(text=paste0(b, "_fit"))), subdivide = 3)) # run overlap estimation - increase subdivide for more robust results
      temp <- overlap[,c(1,2,4)] # extract group A results
      temp2 <- overlap[,c(1,3,5)] # extract group B results
      over_summary[row, col] <- mean(overlap[,4]) # extract mean overlap prob of A and save in summary df
      over_summary[col, row] <- mean(overlap[,5]) # extract mean overlap prob or B and save in summary df
      colnames(temp2)<- colnames(temp)
      tempall <-  rbind(temp, temp2) # merge group A with B results (long format)
      tempall$Group <- c(rep(a, 12000), rep(b, 12000)) # Add group names
      tempall$With <- c(rep(b, 12000), rep(a, 12000)) # Add opposite names for with column
      tempall$Between <- c(rep(ab, 24000))  
      colnames(tempall) <- c("Overlap_vol", "SEVb_vol", "Overlap_prob", "Group", "Between")
      overlap_all <- rbind(overlap_all, tempall) # add to previous results
      print(ab)
    }
  }
}
write.csv(overlap_all, "Overlap_all_long.csv")
write.csv(over_summary, "Overlap_summary.csv")

### Sexes ######
m_f <- ellipsoid.Overlap(m_fit, f_fit, subdivide = 3) ##
write.csv(m_f, "m_f-overlap_new.csv")

### Life stages #####
juv_pb <- ellipsoid.Overlap(cj_fit, pb_fit, subdivide = 3)	##
write.csv(juv_pb, "juv_pb-overlap.csv")

sub_pb <- ellipsoid.Overlap(pj_fit, pb_fit, subdivide = 3)	##
write.csv(sub_pb, "sub_pb-overlap.csv")

sub_juv <- ellipsoid.Overlap(pj_fit, cj_fit, subdivide = 3)	##
write.csv(sub_juv, "sub_juv-overlap.csv")

ad_pb <- ellipsoid.Overlap(ad_fit, pb_fit, subdivide = 3) ##
write.csv(ad_pb, "ad_pb-overlap.csv")		

ad_juv <- ellipsoid.Overlap(ad_fit, cj_fit, subdivide = 3)	
write.csv(ad_juv, "ad_juv-overlap.csv")

ad_sub <- ellipsoid.Overlap(ad_fit, pj_fit, subdivide = 3)
write.csv(ad_sub, "ad_sub-overlap.csv")

### Maternal signatures ####
pb_f <- ellipsoid.Overlap(pb_fit, fjuv_fit, subdivide = 3)
write.csv(pb_f, "pb_fjuv-overlap.csv")

pb_fsub <- ellipsoid.Overlap(pb_fit, fsub_fit, subdivide = 3)
write.csv(pb_fsub, "pb_fsub-overlap.csv")

pb_fad <- ellipsoid.Overlap(pb_fit, fad_fit, subdivide = 3)
write.csv(pb_fad, "pb_fad-overlap.csv")

# add to Overlap summary table
over_summary["pb", "fjuv"] <- mean(pb_f[,4])
over_summary["fjuv", "pb"] <- mean(pb_f[,5])

over_summary["pb", "fsub"] <- mean(pb_fsub[,4])
over_summary["fsub", "pb"] <- mean(pb_fsub[,5])

over_summary["pb", "fad"] <- mean(pb_fad[,4])
over_summary["fad", "pb"] <- mean(pb_fad[,5])

# save Overlap summary again
write.csv(over_summary, "Overlap_summary.csv")


## Visualisation parameters ####

## literature values for ontogeny
onto <- c(0, # age at birth
          6.01, # age of leaving coastal nurseries
          11.8) # age of maturity (based on 200cm TL in Rosa et al. 2017 growth curve)
          
onto_c <- onto+1.5

### Colors ####

sexcols <- c("#EDDA2C", "#5F55A0")
lifecols <- c("#EEE7C3", "#86CAC1", "#48748D", "#27305E")
allcols <- c(sexcols, lifecols)
IDcolors <- c25 <- c(
  "dodgerblue2", "#E31A1C", # red
               "green4",
               "#6A3D9A", # purple
               "#FF7F00", # orange
               "black", "gold1",
               "skyblue2", "#FB9A99", # lt pink
               "palegreen2",
               "#CAB2D6", # lt purple
               "#FDBF6F", # lt orange
               "gray70", "khaki2",
               "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
               "darkturquoise", "green1", "yellow4", "yellow3",
               "darkorange4", "brown"
)

### Plot themes #####
pub_theme = theme_classic()+
  theme(#panel.grid.major = element_line(size = 0.5, color = "grey"),
    # axis.line = element_line(size = 0.7, color = "black"), 
    text = element_text(size = 14),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"))


# ????????????????????????????????????????????????????????????
# ????????????????????????????????????????????????????????????
# ????????????????????????????????????????????????????????????

# Manuscript Figures and Tables ######

## Table 1 - Lifestage summary ####
LS_sum <- SIA %>% 
  group_by(Lifestage) %>%
  summarise(n = n(), n_ID = n_distinct(FishID), 
            Cmean=mean(d13C), Csd=sd(d13C),
            Nmean=mean(d15N), Nsd=sd(d15N),
            Smean=mean(d34S), Ssd=sd(d34S))
overall_sum <- SIA %>%
  summarise(
    n = n(), n_ID = n_distinct(FishID), 
    Cmean = mean(d13C), Csd = sd(d13C),
    Nmean = mean(d15N), Nsd = sd(d15N),
    Smean = mean(d34S), Ssd = sd(d34S)
  ) %>%
  mutate(Lifestage = "Overall")
# Combine lifestage summary with the overall summary
LS_sum <- bind_rows(LS_sum, overall_sum)

write.csv(LS_sum, "Lifestage_summary.csv", row.names = F)




## Figure 2 - GAMM Age smooths####
# first best models are run for all three isotopes, then estimates and standard deviation extracted and corrected for intercepts and random effects
# individual plots were saved at 1000 x 300 and compiled manually
### d13C ####
cbest1 <- gam(formula = d13C ~ 
                # Fixed effects
                s(Agemean)+
                # Random effects 
                s(ID, by=Sex, bs="re")+
                s(Long, Lat, bs="re"),
              # data
              data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest1)

# Extract the smooth for Agemean
C_smooth_data <- smooth_estimates(cbest1, "s(Agemean)")

# Calculate the correction for the intercept:
# Base intercept
intercept <- coef(cbest1)["(Intercept)"]

# Calculate the mean random effects contributions
# Predict random effects contributions for the levels in the dataset
C_SIAM <- C_SIAM %>%
  mutate(
    re_ID = predict(cbest1, type = "terms", terms = "s(ID):SexF"),
    re_spatial = predict(cbest1, type = "terms", terms = "s(Long,Lat)")
  )

# Compute mean of random effects contributions
mean_re_ID <- mean(C_SIAM$re_ID, na.rm = TRUE)
mean_re_spatial <- mean(C_SIAM$re_spatial, na.rm = TRUE)

# Total correction for the intercept
intercept_correction <- intercept + mean_re_ID + mean_re_spatial

# Add correction to the smooth data
C_smooth_data <- C_smooth_data %>%
  mutate(est_corrected = .estimate + intercept_correction)


# Plot smooth over raw data
ggplot() +
  # Plot raw data points
  geom_vline(xintercept = ls[1], linewidth=1, alpha=1)+
  geom_vline(xintercept = ls[2], linewidth=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], linewidth=1, alpha=0.5)+
  geom_point(data = C_SIAM, aes(x = Agemean, y = d13C), alpha = 0.5, color = "black", size = 3) +
  # Add smooth term with confidence interval
  geom_line(data = C_smooth_data, aes(x = Agemean, y = est_corrected), color = "blue", linewidth = 0.8) +
  # Add confidence intervals
  geom_ribbon(data = C_smooth_data, aes(x = Agemean, ymin = est_corrected - .se, ymax = est_corrected + .se),
              fill = "blue", alpha = 0.2) +
  labs(
    x = "Age (Years)",
    y = expression(paste(delta^{13}, "C (\u2030)"))
  ) +
  pub_theme+
  theme(
    axis.title = element_text(size = 24), 
    axis.text = element_text(size = 18),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  )


### d15N ####
#best overall selected model without Long,Lat
nbest1 <- gam(formula = d15N ~ 
                # Fixed effects
                s(Agemean)+
                # Random effects 
                s(ID, bs="re")+
                s(jDate, bs="re"),
              # data
              data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nbest1)

# Extract the smooth for Agemean
n_smooth_data <- smooth_estimates(nbest1, "s(Agemean)")

# Calculate the correction for the intercept:
# Base intercept
n_intercept <- coef(nbest1)["(Intercept)"]

# Calculate the mean random effects contributions
# Predict random effects contributions for the levels in the dataset
N_SIAM <- N_SIAM %>%
  mutate(
    re_ID = predict(nbest1, type = "terms", terms = "s(ID)"),
    re_date = predict(nbest1, type = "terms", terms = "s(jDate)")
  )

# Compute mean of random effects contributions
n_mean_re_ID <- mean(N_SIAM$re_ID, na.rm = TRUE)
n_mean_re_date <- mean(N_SIAM$re_date, na.rm = TRUE)

# Total correction for the intercept
n_intercept_correction <- n_intercept + n_mean_re_ID + n_mean_re_date

# Add correction to the smooth data
n_smooth_data <- n_smooth_data %>%
  mutate(est_corrected = .estimate + n_intercept_correction)

# Plot smooth over raw data
ggplot() +
  # Plot raw data points
  geom_vline(xintercept = ls[1], linewidth=1, alpha=1)+
  geom_vline(xintercept = ls[2], linewidth=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], linewidth=1, alpha=0.5)+
  geom_point(data = N_SIAM, aes(x = Agemean, y = d15N), alpha = 0.5, color = "black", size = 3) +
  # Add smooth term with confidence interval
  geom_line(data = n_smooth_data, aes(x = Agemean, y = est_corrected), color = "blue", linewidth = 0.8) +
  # Add confidence intervals
  geom_ribbon(data = n_smooth_data, aes(x = Agemean, ymin = est_corrected - .se, ymax = est_corrected + .se),
              fill = "blue", alpha = 0.2) +
  labs(
    x = "Age (Years)",
    y = expression(paste(delta^{15}, "N (\u2030)"))
  ) +
  pub_theme+
  theme(
    axis.title = element_text(size = 24), 
    axis.text = element_text(size = 18),
    axis.text.x = element_blank(), 
    axis.title.x = element_blank()
  )

### d34S ####
sbest <- gam(formula = d34S ~ 
               # Fixed effects
               s(Agemean),
             # data
             data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sbest)

# Extract the partial effect of the spline s(Agemean)
s_partial_effect <- smooth_estimates(sbest, "s(Agemean)")

# Add the model's overall intercept
s_partial_effect$est_with_intercept <- s_partial_effect$.estimate + coef(sbest)["(Intercept)"]

# Plot smooth over raw data
ggplot() +
  # Plot raw data points
  geom_vline(xintercept = ls[1], linewidth=1, alpha=1)+
  geom_vline(xintercept = ls[2], linewidth=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], linewidth=1, alpha=0.5)+
  geom_point(data = S_SIAM, aes(x = Agemean, y = d34S), alpha = 0.5, color = "black", size = 3) +
  # Add smooth term with confidence interval
  geom_line(data = s_partial_effect, aes(x = Agemean, y = est_with_intercept), color = "blue", linewidth = 0.8) +
  geom_ribbon(data = s_partial_effect, aes(x = Agemean, ymin = est_with_intercept - .se, ymax = est_with_intercept + .se), 
              alpha = 0.2, fill = "blue") +
  labs(
    x = "Age (Years)",
    y = expression(paste(delta^{34}, "S (\u2030)"))
  ) +
  pub_theme+
  theme(
    axis.title = element_text(size = 24), 
    axis.text = element_text(size = 18)
  )

## Figure 3 - SEVb and overlap ####

#SEVb summary calculation
sevb_sub <- sevb %>%
  mutate(Group1 = Group)%>%
  separate(Group1, into = c("Sex", "Lifestage"), sep = "_")
sevb_sub <- sevb_sub[1:72000,]

# Calculate the percentiles for each group
SEVb_summary <- sevb_sub %>%
  dplyr::group_by(Group) %>%
  mutate(Group = ordered(Group, levels=rev(c("F_ALL", "M_ALL", "ALL_PB", "ALL_JUV", "ALL_SUB", "ALL_AD")))) %>%
  dplyr::summarise(
    lower_95 = percentile(Volume, 0.95)[1],
    upper_95 = percentile(Volume, 0.95)[2],
    lower_75 = percentile(Volume, 0.75)[1],
    upper_75 = percentile(Volume, 0.75)[2],
    lower_50 = percentile(Volume, 0.50)[1],
    upper_50 = percentile(Volume, 0.50)[2],
    mean = mean(Volume),
    short_mean = round(mean(Volume), 1),
    median = median(Volume)
  ) 
SEVb_summary <- SEVb_summary %>%
  mutate(group_num = 1:nrow(SEVb_summary))
levels(SEVb_summary$Group) <- c("Ad", "Sub", "Juv", "Pre", "M", "F")


sevb_bar <- ggplot(SEVb_summary, aes(x=Group, y=mean)) + 
  geom_bar(aes(fill = Group),stat = "identity", color="black", size=1) +
  geom_errorbar(aes(ymin=lower_95, ymax=upper_95),width=.3, size=1, alpha=1)+
  geom_point(aes(y = median), shape = 3, size = 2, stroke=1.5, color = "red", alpha=1) +
  geom_text(aes(label=short_mean, y=0.5), size=8)+
  pub_theme+
  labs(y = expression(SEV[b]))+
  theme(axis.title.y = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 24)
        )+
  scale_fill_manual(values = rev(allcols))+
  coord_flip()
sevb_bar
# saved at 600 x 1000

# Overlap part
# convert overlap data to long format
over_summary_short <- over_summary[1:6, 1:6]
overstats_long <- over_summary_short %>% 
  as.data.frame() %>%
  rownames_to_column("group") %>%
  pivot_longer(-c(group), names_to = "with", values_to = "overlap")
overstats_long$overlap <- round(overstats_long$overlap, 2)
overstats_long$group <- ordered(overstats_long$group, levels=rev(c("f","m", "pb", "juv", "sub","ad")))
overstats_long$with <- ordered(overstats_long$with, levels=c("f","m", "pb", "juv", "sub","ad"))

Overlap_table <- ggplot(overstats_long)+
  geom_tile(aes(x=with, y=group), color="white", fill="white")+
  geom_text(aes(x=with, y=group,label = overlap, size=overlap))+
  #scale_fill_viridis_b(option= "D", direction = -1, na.value="white")+
  scale_alpha_binned(range = c(0.5, 0.9))+
  scale_size_continuous(range=c(4, 20))+
  geom_hline(yintercept = seq(.5, 8.5, 1), size = .2) +
  #scale_fill_continuous(na.value = 'red')+
  scale_x_discrete(position = "top", labels = c("F", "M", "Pb", "Juv", "Sub", "Ad")) +
  pub_theme+
  theme(legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 26)
  )
Overlap_table

## Figure 4 - individual trajectories ####
# Trajectories plotted for each isotope, saved at 1000 x 300 and then compiled manually

Age_C <- ggplot(Special, aes(Agemean, d13C))+  
  geom_vline(xintercept = ls[1], size=1, alpha=1)+
  geom_vline(xintercept = 6.25, size=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], size=1, alpha=0.5)+
  #geom_errorbar(aes(xmin=Agemin, xmax=Agemax, color = Sex), size=1,width=.1, alpha=0.2)+
  geom_line(data = SIA, aes(group = FishID, linetype = Sex, color=lat_class), size=1.5, alpha=0.7)+
  geom_point(data = SIA, (aes(shape = Sex, color = lat_class)), size=4, alpha=0.9)+
  #geom_line(data = not_special, aes(group = FishID, linetype = Sex), color = "grey50", size=1, alpha=0.5)+
  #geom_point(data = not_special, (aes(shape = Sex)), color = "grey50", size=2, alpha=0.5)+
  #geom_point(aes(color=Sex, shape = Sex), size=3, alpha=1)+
  #geom_line(aes(color = Sex, linetype = Sex, group = FishID), size=1.5, alpha=1)+
  #scale_color_viridis_d(option='turbo')+
  #scale_color_manual(values=sexcols)+
  scale_color_viridis_d()+
  xlab("Mean Age")+
  ylab(expression(paste(delta^{13}, "C (\u2030)")))+
  pub_theme+
  theme(legend.position = "none",
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()
        )
Age_C

Age_N <- ggplot(Special, aes(Agemean, d15N))+  
  geom_vline(xintercept = ls[1], size=1, alpha=1)+
  geom_vline(xintercept = 6.25, size=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], size=1, alpha=0.5)+
  #geom_errorbar(aes(xmin=Agemin, xmax=Agemax, color = Sex), size=1,width=.1, alpha=0.2)+
  geom_line(data = SIA, aes(group = FishID, linetype = Sex, color=lat_class), size=1.5, alpha=0.8)+
  geom_point(data = SIA, (aes(shape = Sex, color = lat_class)), size=4, alpha=0.8)+
  #geom_line(data = not_special, aes(group = FishID, linetype = Sex), color = "grey50", size=1, alpha=0.5)+
  #geom_point(data = not_special, (aes(shape = Sex)), color = "grey50", size=2, alpha=0.5)+
  #geom_point(aes(color=Sex, shape = Sex), size=3, alpha=1)+
  #geom_line(aes(color = Sex, linetype = Sex, group = FishID), size=1.5, alpha=1)+
  #scale_color_viridis_d(option='turbo')+
  #scale_color_manual(values=sexcols)+
  scale_color_viridis_d()+
  xlab("Mean Age")+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  pub_theme+
  theme(legend.position = "none",
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()
        )
Age_N

Age_S <- ggplot(Special, aes(Agemean, d34S))+  
  geom_vline(xintercept = ls[1], size=1, alpha=1)+
  geom_vline(xintercept = 6.25, size=1, alpha=0.5)+
  geom_vline(xintercept = ls[3], size=1, alpha=0.5)+
  #geom_errorbar(aes(xmin=Agemin, xmax=Agemax, color = Sex), size=1,width=.1, alpha=0.2)+
  geom_line(data = SIA, aes(group = FishID, linetype = Sex, color=lat_class), size=1.5, alpha=0.8)+
  geom_point(data = SIA, (aes(shape = Sex, color = lat_class)), size=4, alpha=0.8)+
  #geom_line(data = not_special, aes(group = FishID, linetype = Sex), color = "grey50", size=1, alpha=0.5)+
  #geom_point(data = not_special, (aes(shape = Sex)), color = "grey50", size=2, alpha=0.5)+
  #geom_point(aes(color=Sex, shape = Sex), size=3, alpha=1)+
  #geom_line(aes(color = Sex, linetype = Sex, group = FishID), size=1.5, alpha=1)+
  #scale_color_viridis_d(option='turbo')+
  #scale_color_manual(values=sexcols)+
  scale_color_viridis_d()+
  xlab("Mean Age")+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  pub_theme+
  theme(legend.position = "none",
        axis.title = element_text(size = 24), 
        axis.text = element_text(size = 18),
        axis.text.y = element_text(size = 24)
        )
Age_S


# Supplementary Material - Tables and Figures #####

## Figure S 1 - Biplot of shark vertebrae SIA studies ####
# comparing of overall mean and sd with other Shark vertebra studies 

# import Vertebrae SIA review table
litSIA <- read_xlsx("C:/Users/Lenovo/Documents/SCIENCE/Vertebrae SIA/Shark vertebra SIA review.xlsx")
#litSIA <- litSIA %>%
#  filter(`RN` !=25)
litSIA <- litSIA[,c(1, 6, 13, 14, 17, 18)]
colnames(litSIA) <- c("RN", "Species", "Cmean", "Csd", "Nmean", "Nsd")
litSIA <- litSIA %>%
  mutate(Species_abbrev = str_replace(Species,
                                      "^(\\w)\\w*\\s(\\w+)$",   # capture genus initial + species
                                      "\\1. \\2"))


SIA_lit <- ggplot(litSIA)+
  geom_point(aes(Cmean, Nmean), alpha=0.6)+
  geom_errorbar(aes(x=Cmean, ymin=Nmean-Nsd, ymax=Nmean+Nsd), linewidth=0.5, width=0.05, alpha=0.6)+
  geom_errorbarh(aes(y=Nmean, xmin=Cmean-Csd, xmax=Cmean+Csd), linewidth=0.5, height=0.1, alpha=0.6)+
  geom_text(aes(Cmean , Nmean, label = RN), 
                vjust = -0.3, hjust = 1.2, size = 4, parse = TRUE) +
  geom_text(data = litSIA[25,], aes(Cmean , Nmean, label = RN), 
            vjust = -0.3, hjust = 1.2, size = 4, parse = TRUE, color = sexcols[1]) +
  #geom_text(aes(Cmean, Nmean, label = paste0("italic('", Species_abbrev, "')")), 
  #          vjust = -0.3, hjust = - 0.1, size = 2.5, parse = TRUE) +
  geom_point(data=LS_sum[5,], aes(Cmean, Nmean), color=sexcols[1])+
  geom_errorbar(data=LS_sum[5,], aes(x=Cmean, ymin=Nmean-Nsd, ymax=Nmean+Nsd), linewidth=0.8, width=0.05, alpha=1, , color=sexcols[1])+
  geom_errorbarh(data=LS_sum[5,], aes(y=Nmean, xmin=Cmean-Csd, xmax=Cmean+Csd), linewidth=0.8, height=0.1, alpha=1, , color=sexcols[1])+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)")))+
  pub_theme+
  theme(legend.position = "none")
# Save as SVG
ggsave("SIAlit.svg", plot = SIA_lit, width = 6, height = 4, device = "svg")


## Tables S 2 - S 4 - Main GAMM summaries #####
# All values were inserted manually from model summaries
### C model ####
summary(cbest1)
### N model ####
summary(nbest1)
### S model ####
summary(sbest)


## Table S 5 - Location GAMM summaries #####

# All values were inserted manually from model summaries and variance tables
### C model with Location ######
cbest2 <- gam(formula = d13C ~ 
                # Fixed effects
                s(Agemean)+
                SexM + SexF+
                s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
                # Random effects 
                s(ID, bs="re"),
              # data
              data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest2)
appraise(cbest2, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
variance_comp(cbest2)

### N model with Location ######
nbest2 <- gam(formula = d15N ~ 
                # Fixed effects
                s(Agemean)+
                s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
                # Random effects 
                s(ID, bs="re"),
              # data
              data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nbest2)
appraise(nbest2, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
variance_comp(nbest2)

### S model with Location ######
sbest <- gam(formula = d34S ~ 
               # Fixed effects
               s(Agemean),
             # data
             data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sbest)
variance_comp(sbest)


