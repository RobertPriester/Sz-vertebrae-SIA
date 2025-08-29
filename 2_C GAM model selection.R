########## d13C GAMM code for Smooth hammerhead vertebrae SIA manuscript

# Input: 
# - either "VertebraeSIA_RAW" and "Vertebrae_age_reading" data for reproducing all calculations
# - or "VertebraeSIA_Aged" data and already compiled "SEVb" and "Overlap_summary" data to only reproduce manuscript plots with already processed data

# Output:
# - Processed and filtered SIA data with estimated corresponding ages
# - d13C model subset data
# - GAMM model selection and model assessment
# - best model

## Workpace and libraries
setwd("~/SCIENCE/Vertebrae SIA/R")

library(dplyr)
library(lubridate)
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
library(fitdistrplus)
library(MuMIn)
library(broom)
source("Skinner et al. 2019 sppl._ellipsoidcode_functions_final_EDIT.r") #loads functions from Skinner et al. 2019 for SEV, SEVc and SEVb estimation

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Define a function for min-max scaling
scale_min_max <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

# 1. Data import #####

## import RAW data
#RAW <- read.csv("VertebraeSIA_210907.csv")
Age <- read.csv("Vertebrae_length-onto.csv")
#SIA <- read.csv("VertebraeSIA211103_Aged.csv")  ## if _Aged is loaded do not run ageing code
SIA <- read.csv("VertebraeSIA230214_Aged.csv")
rel_sum <- read.csv("Summary_reliso.csv")
rel_sum_age <- read.csv("relchange_age.csv")

over_list <- read.csv("Overlaps.csv", header=F)
#overlap <- read.csv("Overlap_all_long.csv", header=T)

# CN ratio above 4 excluded
SIA <- filter(SIA, CNratio < 4)

SIA$Sex <- ordered(SIA$Sex, levels=c("M", "F"))

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

## data frame modifications
Age <- filter(Age, Age!="VR") # filter VR layer to enable better looping
Age$FishCode <- as.factor(Age$FishCode)
Age$Age <- as.numeric(Age$Age)

## growth calculations
# Add growth rate of section (yearly average over integrated year rings)
SIA$growth <- NA
for(i in 1:nrow(SIA)){
  Fish <- SIA$FishCode[i]
  Ages <- Age%>%
    filter(FishCode==FishCode, Age %in% SIA$Agemin[i]:SIA$Agemax[i])
  SIA$growth[i] <- mean(Ages$Diff)
}
# condition calculation
SIA$condition <- SIA$growth - (-0.049*SIA$VR+1.0732)

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

SIA$Lifestage <- ordered(SIA$Lifestage, levels=c("Pre-birth", "Coastal Juv", "Pelagic Juv", "Adult"))

SIA$Age <- SIA$Agemean+1.5 # centering time on first point

# Calculating tissue date estimates
SIA$tissue_date <- SIA$Date - (months(SIA$FinalAge*12)-months(SIA$Agemean*12)) # calculating age of deposited tissue and subtracting it from capture date
SIA$tissue_j <- julian(SIA$tissue_date)

SIA$Lat_pos <- SIA$Lat
SIA$Lat_pos[(SIA$Lat_pos<=0)] <- SIA$Lat[(SIA$Lat<=0)]*-1

rel_sum$Lifestage <- ordered(rel_sum$Lifestage, levels=c("Pre-birth", "Coastal Juv", "Pelagic Juv", "Adult"))



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

pub_theme = theme_classic()+
  theme(#panel.grid.major = element_line(size = 0.5, color = "grey"),
    # axis.line = element_line(size = 0.7, color = "black"), 
    text = element_text(size = 14),
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14),
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"))

# 2 Data transformation ####

## 2.1 Model subset #####
SIAM <- SIA %>%
  dplyr::select(ID, d13C, d15N, d34S, VR, Agemean, CalcTL, growth, condition, Sex, SexM, SexF, jDate, Lat, Long)
SIAM$Cmm <- scale_min_max(SIAM$d13C)
SIAM$Nmm <- scale_min_max(SIAM$d15N)
SIAM$Smm <- scale_min_max(SIAM$d34S)

## 2.2 Outlier removal ####

### plot data to show outliers
C_Sex <- ggplot(SIA, aes(Sex, d13C))+
  geom_boxplot()+
  ylab(expression(paste(delta^{13}, "C (\u2030)")))+
  pub_theme
C_Sex
# remove one outlier

N_Sex <- ggplot(SIA, aes(Sex, d15N))+
  geom_boxplot()+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  pub_theme
N_Sex
# remove one outlier

S_Sex <- ggplot(SIA, aes(Sex, d34S))+
  geom_boxplot()+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  pub_theme
S_Sex

### remove outliers
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

## 2.3 Scaling ####
# creating scaled subset
C_SIAM_scale <- C_SIAM %>%
  mutate(across(where(~!is.integer(.) & is.numeric(.)), scale))

C_SIAM_mm <- C_SIAM %>%
  mutate(across(where(~!is.integer(.) & is.numeric(.)), scale_min_max))



# 3 Model construction ####
## 3.1 Assessing distribution of response variable ####

shapiro.test(C_SIAM$d13C)
# p < 0.05 - not normally distributed

hist(C_SIAM$d13C)
lines(density(SIA$d13C))

descdist(C_SIAM$Cmm, discrete = FALSE)
descdist(as.numeric(scale(C_SIAM$d13C)), discrete = FALSE)
descdist(as.vector(SIA$Cnorm), discrete = FALSE)

C_dist <- fitdist(C_SIAM$Cmm, "weibull")
plot(C_dist)
# Weibull seems best


### 3.1.1 a-priori distribution family tests in gam  ####

c00_scat <- gam(formula = d13C ~ 
             # Fixed effects
             1,
           # data
           data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c00_scat)
appraise(c00_scat)

c00_tw <- gam(formula = Cmm ~ 
                  # Fixed effects
                  1,
                # data
                data=C_SIAM, method="REML", family = tw(), select=F) # select includes double penalty approach - penalizes Null space
summary(c00_tw)
appraise(c00_tw)

c00_gamma <- gam(formula = Cmm ~ 
                # Fixed effects
                1,
              # data
              data=C_SIAM, method="REML", family = nb, select=F) # select includes double penalty approach - penalizes Null space
summary(c00_gamma)
appraise(c00_tw)


## 3.2 Basis Models #####

###  Zero model ####
c00 <- gam(formula = d13C ~ 
                  # Fixed effects
                  1,
                # data
                data=C_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Null space
summary(c00)

### Complete model####
cxx <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexF+
             SexM+
             #s(Long, Lat, k=5, bs="ds", m=c(1, 0.5))+
             # Fixed effects
             s(ID, bs="re")+
             s(jDate, bs="re")+
             s(Long, Lat, bs="re"),
           # data
           data=C_SIAM, method="REML", family=scat, select=F)
summary(cxx)
gam.check(cxx)
appraise(cxx, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.99,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(cxx)


## 3.3 Forward model selection #####
#### Fixed Effects #####
##### ** C ~ Age #####
c01 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c01)

##### *** C ~ VR #####
c02 <- gam(formula = d13C ~ 
             # Fixed effects
             s(VR, k=5),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c02)
gam.check(c02)

##### C ~ growth #####
c03 <- gam(formula = d13C ~ 
             # Fixed effects
             s(growth),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c03)

##### C ~ TL #####
c04 <- gam(formula = d13C ~ 
             # Fixed effects
             s(CalcTL),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c04)

##### C ~ Sex #####
c05 <- gam(formula = d13C ~ 
             # Fixed effects
             SexM + SexF,
           # data
           data=C_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Null space
summary(c05)

##### C ~ Date #####
c06 <- gam(formula = d13C ~ 
             # Fixed effects
             s(jDate),
           # data
           data=C_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c06)

##### C ~ Lat,Lon #####
c07 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c07)

##### C ~ ID #####
c08 <- gam(formula = d13C ~ 
             # Fixed effects
             s(ID, bs="fs"),
           # data
           data=C_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c08)

##### C ~ condition #####
c009 <- gam(formula = d13C ~ 
             # Fixed effects
             s(condition),
           # data
           data=C_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c009)

MuMIn::model.sel(c00, c01, c02, c03, c04, c05, c06, c07, c08, c009)

# -> Agemean and VR similar but Agemean ecologically more relevant - continue with Agemean


##### <<<<< #####

##### ** C ~ Agemean + Sex #####
c14 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF,
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c14)
draw(c09)

##### ** C ~ Agemean + s(Agemean by Sex) #####
c15 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF,
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c15)

##### C ~ te(Agemean by Sex) #####
c16 <- gam(formula = d13C ~ 
             # Fixed effects
             te(Agemean, by=Sex),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c16)

##### C ~ ti(Agemean by Sex) #####
c17 <- gam(formula = d13C ~ 
             # Fixed effects
             ti(Agemean, by = Sex),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c17)

##### ** C ~ ti(Agemean, ID) #####
c18 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF,
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c18)

MuMIn::model.sel(c01, c14, c15, c16, c17, c18)

##### <<<<<< #######
##### Potentially remove - Capture location should be random factor only #####

##### C ~ Agemean + Lat #####
c19 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Lat),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c19)
draw(c19)

##### C ~ Agemean + Lon #####
c20 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, Lat),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c20)
draw(c20)


##### C ~ Agemean + LatLon #####
c23 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c23)

MuMIn::model.sel(c19, c20, c23)

##### C ~ Agemean + Sex + LatLon #####
c24 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c24)

##### *** C ~ Agemean + s(Agemean by Sex) + LatLon #####
c25 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c25)
gam.check(c25)
appraise(c25, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(c25, residuals=T)

##### C ~ ti(Agemean, ID) #####
c26 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c26)

##### C ~ te(Agemean, Long, Lat) #####
c27 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             ti(Agemean, Long, Lat)+
             SexM+SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c27)


MuMIn::model.sel(c01, c14, c23, c24, c25, c26, c27)

##### <<<<<< ######
##### C ~ Agemean + s(Agemean by Sex) + LatLon  + condition#####
c28 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
             s(condition),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c28)

MuMIn::model.sel(c01, c25, c27)

# best fixed effect models
# without lat: C ~ s(Agemean)
cbest1 <- c01

# with lat long: C ~ s(Agemean) + s(Agemean by Sex) + SexM + SexF + s(LatLon)
cbest2 <- c25

#### Random effects ######

## Testing best previous models (cbest1, cbest2) with different remaining random effects
##### C01 + re(ID) #####
cr11 <- gam(formula = d13C ~ 
             # Fixed effects
             s(Agemean)+
             # Random effects 
             s(ID, bs="re"),
           # data
           data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr11)

##### * C01 + re(ID, by Sex) #####
cr12 <- gam(formula = d13C ~ 
               # Fixed effects
               s(Agemean)+
               # Random effects 
               s(ID, by=Sex, bs="re"),
             # data
             data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr12)

##### C01 + re(jDate) #####
cr13 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(jDate, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr13)

##### C01 + re(Long Lat) #####
cr131 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(Long, Lat, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr131)

MuMIn::model.sel(cbest1, cr11, cr12, cr13, cr131)
# best: cr131
# best most simple: cr131

##### C01 + re(ID, by Sex) + re(jDate) #####
cr14 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, by=Sex, bs="re")+
              s(jDate, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr14)

##### C01 + re(Long, Lat) + re(jDate) #####
cr141 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(Long, Lat, bs="re")+
              s(jDate, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr141)

##### ** C01 + re(ID, by Sex) + re(Lat, Lon) #####
cr15 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              SexF + SexM+
              # Random effects 
              s(ID, by=Sex, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr15)

#####  C01 + re(ID, by Sex) re(jDate) + re(Lat, Lon) #####
cr16 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, by=Sex, bs="re")+
              s(jDate, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr16)

MuMIn::model.sel(c01, cr131, cr141, cr12, cr13, cr14, cr15, cr16)
# best: cr16
# best most simple: cr15

##### <<<<<<<< #######

##### ** Cbest2 + re(ID) #####
cr21 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr21)

##### Cbest2 + re(ID, by Sex) #####
cr22 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, by=Sex, bs="re"),
            # data
            data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr22)

##### Cbest2 + re(jDate) #####
cr23 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(jDate, bs="re"),
            # data
            data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr23)

MuMIn::model.sel(cbest2, cr21, cr22, cr23)
# best: cr21
# best most simple: cr21

##### Cbest2 + re(ID, by Sex) + re(jDate) #####
cr24 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re")+
              s(jDate, bs="re"),
            # data
            data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cr24)


MuMIn::model.sel(cbest2, cr21, cr22, cr23, cr24)
# best: cr24
# best most simple: cr21

# Using REML
MuMIn::model.sel(cr21, cr22, cr23, cr24)
# same as above

### BEST forward selected model without Lat Long ######
cbest1 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, by=Sex, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest1)
appraise(cbest1, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(cbest1)
gam.check(cbest1)

### BEST forward selected model with Lat Long ######
cbest2 <- gam(formula = d13C ~ 
               # Fixed effects
               s(Agemean)+
               s(Agemean, by=Sex, k=8)+
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
draw(cbest2)


## Backward model selection without Lat Long#####
### C ~ Agemean + s(ID, Sex, bs="re") + s(LatLon, bs="re") #####
cb00 <- gam(formula = d13C ~ 
                # Fixed effects
                s(Agemean)+
                # Random effects 
                s(ID, by=Sex, bs="re")+
                s(Long, Lat, bs="re"),
              # data
              data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb00)

### remove 1 ####
### C ~ s(ID, Sex, bs="re")  + s(LatLon, bs="re") #####
cb01 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(VR, k=5)+
              # Random effects
              s(ID, by=Sex, bs="re")+
              s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb01)


### C ~ s(Agemean) + s(LatLon, bs="re") #####
cb02 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              #s(ID, Sex, bs="re")+
              s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb02)

### C ~ s(Agemean) + s(Id, by Sex, bs="re") #####
cb03 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, by=Sex, bs="re"),
              #s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb03)

### remove 2 ####
### C ~ SexF + SexM+ s(LatLon, bs="re") #####
cb04 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(VR, k=5)+
              # Random effects
              #s(ID, Sex, bs="re")+
              s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb04)

### C ~ s(ID, by Sex, bs="re") #####
cb05 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(VR, k=5)+
              # Random effects
              s(ID, by=Sex, bs="re"),
              #s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb05)

### C ~ s(Agemean) + s(LatLon, bs="re")  #####
cb06 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean, k=5),
              #SexM+
              #SexF+
              # Random effects
              #s(ID, Sex, bs="re")+
              #s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb06)

### C ~ s(Agemean) + s(ID, Sex, bs="re") #####
c00 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID,by= Sex, bs="re"),
              #s(Lat, Long, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(c00)

MuMIn::model.sel(cb00, cb01, cb02, cb03, cb04, cb05, cb06, c00)
# best: cb00
# best most simple: cb02

## Backward model selection with Lat Long#####
##### C ~ s(Agemean) + s(Agemean, by=Sex) + SexMSexF + s(Long, Lat) + re(ID) #####
cb20 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb20)

### remove 1 ####
##### C ~ s(Agemean, by=Sex) + SexMSexF + s(Long, Lat) + re(ID) #####
cb21 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb21)

##### ** C ~ s(Agemean) + SexMSexF + s(Long, Lat) + re(ID) #####
cb22 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb22)

##### C ~ s(Agemean) + s(Agemean, by=Sex) s(Long, Lat) + re(ID) #####
cb23 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb23)

##### C ~ s(Agemean) + s(Agemean, by=Sex) + SexMSexF + re(ID) #####
cb24 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb24)

##### C ~ s(Agemean) + s(Agemean, by=Sex) + SexMSexF + s(Long, Lat) #####
cb25 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb25)

MuMIn::model.sel(cb20, cb21, cb22, cb23, cb24, cb25)
#best: cb22

### remove 2 ####

##### C ~ SexMSexF + s(Long, Lat) + re(ID) #####
cb26 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb26)

##### C ~ s(Agemean, by=Sex) + s(Long, Lat) + re(ID) #####
cb27 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb27)

##### C ~ s(Agemean, by=Sex) + SexMSexF + re(ID) #####
cb28 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb28)

##### C ~ s(Agemean, by=Sex) + SexMSexF + s(Long, Lat) #####
cb29 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb29)

##### C ~ s(Agemean) + s(Long, Lat) + re(ID) #####
cb30 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb30)

##### C ~ s(Agemean) + SexMSexF + re(ID) #####
cb31 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb31)

##### C ~ s(Agemean) + SexMSexF + s(Long, Lat) #####
cb32 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb32)

##### C ~ s(Agemean) + s(Agemean, by=Sex) + re(ID) #####
cb33 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb33)

##### C ~ s(Agemean) + s(Agemean, by=Sex) + s(Long, Lat) #####
cb34 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb34)

##### C ~ s(Agemean) + s(Agemean, by=Sex) + SexMSexF  #####
cb35 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF,
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb35)

MuMIn::model.sel(cb22, cb26, cb27, cb28, cb29, cb30, cb31, cb32, cb33, cb34, cb35)
#best: cb22

### remove 3 ####
##### C ~ s(Long, Lat) + re(ID) #####
cb36 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb36)

##### C ~ SexMSexF + re(ID) #####
cb37 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb37)

##### C ~ SexMSexF + s(Long, Lat) #####
cb38 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb38)

##### C ~ s(Agemean, by=Sex) + re(ID) #####
cb39 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb39)

##### C ~ s(Agemean, by=Sex) + s(Long, Lat) #####
cb40 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb40)

##### C ~ s(Agemean, by=Sex) + SexMSexF + #####
cb41 <- gam(formula = d13C ~ 
              # Fixed effects
              #s(Agemean)+
              s(Agemean, by=Sex, k=8)+
              SexM + SexF,
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb41)

##### C ~ s(Agemean) + re(ID) #####
cb42 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb42)

##### C ~ s(Agemean) + s(Long, Lat) #####
cb43 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              #SexM + SexF+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb43)

##### C ~ s(Agemean) + SexMSexF + #####
cb44 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              #s(Agemean, by=Sex, k=8)+
              SexM + SexF,
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb44)

##### C ~ s(Agemean) + s(Agemean, by=Sex) #####
cb45 <- gam(formula = d13C ~ 
              # Fixed effects
              s(Agemean)+
              s(Agemean, by=Sex, k=8),
              #SexM + SexF+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=C_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cb45)

MuMIn::model.sel(cb22, cb36, cb37, cb38, cb39, cb40, cb41, cb42, cb43, cb44, cb45)
#best: cb22


## BEST overall selected model without Lat Long ######
cbest1 <- gam(formula = d13C ~ 
                # Fixed effects
                s(Agemean)+
                # Random effects 
                s(ID, by=Sex, bs="re")+
                s(Long, Lat, bs="re"),
              # data
              data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
appraise(cbest1, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(cbest1)
summary(cbest1)


## BEST overall selected model with Lat Long ######
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
draw(cbest2, select = "s(Long,Lat)")

cbest2_split <- gam(formula = d13C ~ 
                # Fixed effects
                s(Agemean)+
                SexM + SexF+
                s(Long)+
                s(Lat)+
                # Random effects 
                s(ID, bs="re"),
              # data
              data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest2_split)



## 3.3 Transormations #####

### BEST forward selected model ######
cbest_norm <- gam(formula = d13C ~ 
               # Fixed effects
               s(VR, k=5)+
               SexM+
               SexF+
               # Random effects
               s(ID, Sex, bs="re")+
               s(Lat, Long, bs="re"),
             # data
             data=C_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest_norm)
draw(cbest_norm)
appraise(cbest_norm,
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.99,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")

cbest_scale <- gam(formula = d13C ~ 
                     # Fixed effects
                     s(VR, k=5)+
                     SexM+
                     SexF+
                     # Random effects
                     s(ID, Sex, bs="re")+
                     s(Lat, Long, bs="re"),
                   # data
                   data=C_SIAM_scale, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest_scale)
appraise(cbest_scale,
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")

cbest_mm <- gam(formula = d13C ~ 
                  # Fixed effects
                  s(VR, k=5)+
                  SexM+
                  SexF+
                  # Random effects
                  s(ID, by=Sex, bs="re")+
                  s(Lat, Long, bs="re"),
                # data
                data=C_SIAM_mm, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(cbest_mm)
draw(cbest_mm)
appraise(cbest_mm,
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")

# --> Norm best QQ plot

# Visualisation ######

## VR #####
# basic gratia plot
draw(cbest, select="s(VR)")

# customizable ggplot
sm <- smooth_estimates(cbest, smooth = "s(VR)")

sm %>%
  add_confint() %>%
  ggplot(aes(y = est, x = VR)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(colour = "forestgreen", size = 1.5) +
  labs(y = "Partial effect",
       title = expression("Partial effect of" ~ f(VR)),
       x = expression("Vertebrae radius"))+
  pub_theme

## Sex #####
# Create a new data frame with a sequence of VR values
newdata <- data.frame(VR = seq(min(C_SIAM$VR), max(C_SIAM$VR), length.out = 100))

# Predict the fitted values for SexM=1 and SexF=0
newdata$SexM <- 1
newdata$SexF <- 0
newdata$ID <- (C_SIAM$ID)[1]
newdata$Sex <- levels(C_SIAM$Sex)[1]
newdata$Lat <- mean(C_SIAM$Lat)
newdata$Long <- mean(C_SIAM$Long)
newdata$predM <- predict(cbest, newdata, type="response")

# Predict the fitted values for SexM=0 and SexF=1
newdata$SexM <- 0
newdata$SexF <- 1
newdata$predF <- predict(cbest, newdata, type="response")

# Combine the predictions into a single data frame for plotting
plotdata <- rbind(
  data.frame(VR = newdata$VR, Sex = "M", Pred = newdata$predM),
  data.frame(VR = newdata$VR, Sex = "F", Pred = newdata$predF)
)
plotdata$Sex = ordered(plotdata$Sex, levels=c("M", "F"))

# Create the plot
ggplot(plotdata, aes(x=Sex, y=Pred)) +
  geom_boxplot(aes(fill=Sex), alpha=0.7) +
  labs(title="Predicted values for SexM and SexF",
       x="Sex",
       y="Predicted values")+
  scale_fill_viridis_d()+
  pub_theme

