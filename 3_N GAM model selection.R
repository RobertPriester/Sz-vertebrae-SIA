########## d15N GAMM code for Smooth hammerhead vertebrae SIA manuscript

# Input: 
# - either "VertebraeSIA_RAW" and "Vertebrae_age_reading" data for reproducing all calculations
# - or "VertebraeSIA_Aged" data and already compiled "SEVb" and "Overlap_summary" data to only reproduce manuscript plots with already processed data

# Output:
# - Processed and filtered SIA data with estimated corresponding ages
# - d15N model subset data
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

#over_list <- read.csv("Overlaps.csv", header=F)
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

SIA$Age <- SIA$growth+1.5 # centering time on first point

# Calculating tissue date estimates
#SIA$tissue_date <- SIA$Date - (months(SIA$FinalAge*12)-months(SIA$growth*12)) # calculating age of deposited tissue and subtracting it from capture date
#SIA$tissue_j <- julian(SIA$tissue_date)

SIA$Lat_pos <- SIA$Lat
SIA$Lat_pos[(SIA$Lat_pos<=0)] <- SIA$Lat[(SIA$Lat<=0)]*-1

rel_sum$Lifestage <- ordered(rel_sum$Lifestage, levels=c("Pre-birth", "Coastal Juv", "Pelagic Juv", "Adult"))



mycolors <- c("royalblue4", "red2")
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
  dplyr::select(ID, d13C, d15N, d34S, VR, Agemean, growth, CalcTL, growth, Sex, SexM, SexF, jDate, Lat, Long)

## 2.2 Outlier removal ####

# remove one outlier

N_Sex <- ggplot(SIA, aes(Sex, d15N))+
  geom_boxplot()+
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  pub_theme
N_Sex
# remove one outlier

### remove outliers
### N
# Calculate z-scores for specified columns
z_scores <- scale(SIAM[, "d15N"])
# Define a threshold for outliers (e.g., 3 standard deviations)
threshold <- 3
# Create a boolean mask for outliers
outliers_mask <- apply(z_scores, 2, function(x) abs(x) > threshold)
# Remove rows with outliers
N_SIAM <- SIAM[!rowSums(outliers_mask), ]

## 2.3 Scaling ####
# creating scaled subset
N_SIAM_scale <- N_SIAM %>%
  mutate(across(where(~!is.integer(.) & is.numeric(.)), scale))

N_SIAM_mm <- N_SIAM %>%
  mutate(across(where(~!is.integer(.) & is.numeric(.)), scale_min_max))



# 3 Model construction ####
## 3.1 Assessing distribution of response variable ####

shapiro.test(N_SIAM$d15N)
# p < 0.05 - not normally distributed

hist(N_SIAM$d15N)
lines(density(N_SIAM$d15N))

descdist(N_SIAM$d15N, discrete = FALSE)

N_dist <- fitdist(N_SIAM$d15N, "weibull")
plot(N_dist)


## 3.2 Basis Models #####

###  Zero model ####
n00 <- gam(formula = d15N ~ 
             # Fixed effects
             1,
           # data
           data=N_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Null space
summary(n00)

### Complete model####
nxx <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexF+
             SexM+
             # Fixed effects
             s(ID,k=8, bs="re")+
             s(Long, Lat, k=8, bs="re")+
             s(jDate, bs="re"),
           #s(Long, Lat, bs="re"),
           # data
           data=N_SIAM, method="REML", family=scat, select=F)
summary(nxx)
gam.check(nxx)
appraise(nxx, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.99,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(nxx)


## 3.3 Forward model selection #####
#### Fixed Effects #####
##### ** N ~ Age #####
n01 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n01)

##### *** N ~ VR #####
n02 <- gam(formula = d15N ~ 
             # Fixed effects
             s(VR, k=5),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n02)
gam.check(n02)

##### ** N ~ growth #####
n03 <- gam(formula = d15N ~ 
             # Fixed effects
             s(growth),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n03)

##### N ~ TL #####
n04 <- gam(formula = d15N ~ 
             # Fixed effects
             s(CalcTL),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n04)

##### N ~ Sex #####
n05 <- gam(formula = d15N ~ 
             # Fixed effects
             SexM + SexF,
           # data
           data=N_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Null space
summary(n05)

##### N ~ Date #####
n06 <- gam(formula = d15N ~ 
             # Fixed effects
             s(jDate),
           # data
           data=N_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n06)

##### N ~ Lat,Lon #####
n07 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n07)

##### N ~ ID #####
n08 <- gam(formula = d15N ~ 
             # Fixed effects
             s(ID, bs="fs"),
           # data
           data=N_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n08)

MuMIn::model.sel(n00, n01, n02, n03, n04, n05, n06, n07, n08)

##### <<<<<< #####

##### * N ~ Agemean + Sex #####
n09 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF,
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n09)
draw(n09)

##### * N ~ Agemean + s(Agemean by Sex) #####
n10 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF,
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n10)

##### N ~ te(Agemean by Sex) #####
n11 <- gam(formula = d15N ~ 
             # Fixed effects
             te(Agemean, by=Sex),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n11)

##### N ~ ti(Agemean by Sex) #####
n12 <- gam(formula = d15N ~ 
             # Fixed effects
             ti(Agemean, by = Sex),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n12)

#####  N ~ ti(Agemean, ID) #####
n13 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF,
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n13)

MuMIn::model.sel(n01, n09, n10, n11, n12, n13)


##### <<<<<< #####
##### Potentially remove - Capture location should be random factor only #####

##### *** C ~ Agemean + LatLon #####
n23 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n23)

##### C ~ Agemean + Sex + LatLon #####
n24 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=6)+
             SexM+
             SexF+
             s(Long, Lat, k=6, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n24)

##### C ~ Agemean + s(Agemean by Sex) + LatLon #####
n25 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n25)


##### C ~ s(Agemean) + ti(Agemean, ID) + LatLong #####
n26 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n26)


MuMIn::model.sel(n01, n23, n24, n25, n26)

##### <<<<<< ######
##### C ~ Agemean + LatLon + ti(Agemean, Long, Lat) #####
n27 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n27)

##### C ~ Agemean + Sex + LatLon + ti(Agemean, Long, Lat) #####
n28 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=6)+
             SexM+
             SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n24)

##### C ~ Agemean + s(Agemean by Sex) + LatLon + ti(Agemean, Long, Lat) #####
n29 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n29)


##### C ~ s(Agemean) + ti(Agemean, ID) + LatLong + ti(Agemean, Long, Lat) #####
n30 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n30)

MuMIn::model.sel(n23, n27, n28, n29, n30)

##### <<<<<< ######
##### C ~ Agemean + s(Agemean by Sex) + LatLon  + condition#####
n31 <- gam(formula = d15N ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
             s(condition),
           # data
           data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(n31)

MuMIn::model.sel(c01, n25, n27)

#### BEST fixed effect models ####
# without lat: C ~ s(Agemean)
nbest1 <- n01

# with lat long: C ~ s(Agemean) + s(Agemean by Sex) + SexM + SexF + s(LatLon)
nbest2 <- n23


### Random effects ######
## Testing best previous models (cbest1, cbest2) with different remaining random effects
#####  N00 #####
nr00 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean),
              # Random effects 
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr00)

#####  N01 + re(ID) #####
nr11 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr11)

##### N01 + re(ID, by Sex) #####
nr12 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, by=Sex, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr12)

##### N01 + re(jDate) #####
nr13 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(jDate, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr13)

##### N01 + re(Long Lat) #####
nr131 <- gam(formula = d15N ~ 
               # Fixed effects
               s(Agemean)+
               # Random effects 
               s(Long, Lat, bs="re"),
             # data
             data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr131)

MuMIn::model.sel(nr00, nr11, nr12, nr13, nr131)
# best: nr11
# best most simple: nr11

##### ** N01 + re(ID, by Sex) + re(jDate) #####
nr14 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, bs="re")+
              s(jDate, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr14)

##### N01 + re(ID, by Sex) + re(Lat, Lon) #####
nr15 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr15)

#####  N01 + re(ID, by Sex) re(jDate) + re(Lat, Lon) #####
nr16 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, bs="re")+
              s(jDate, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr16)

MuMIn::model.sel(nr00, nr11, nr12, nr13, nr131, nr14, nr15, nr16)
# best: nr14
# best most simple: nr14


##### <<<<<<<< #######
nr20 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean, k=5)+
              s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr20)

##### ** Nbest2 + re(ID) #####
nr21 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr21)

##### Nbest2 + re(ID, by Sex) #####
nr22 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, by=Sex, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr22)

##### Nbest2 + re(jDate) #####
nr23 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(jDate, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr23)

##### Nbest2 + re(jDate) #####
nr24 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re")+
              s(jDate, bs="re")
              ,
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nr24)


MuMIn::model.sel(nr20, nr21, nr22, nr23)
# best: nr21
# best most simple: nr21


### BEST forward selected model without Long, Lat ######
nbest1 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              # Random effects 
              s(ID, bs="re")+
              s(jDate, bs="re"),
            # data
            data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nbest1)
appraise(nbest1, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(nbest1)
gam.check(nbest1)


### BEST forward selected model with Long, Lat ######
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
draw(nbest2)
gam.check(nbest2)

## Backward model selection without Long, Lat #####

###  N ~ s(VR, k=5) + s(ID, Sex, bs="re") + s(jDate, bs="re") #####
nb00 <- nr11
summary(nb00)

###  N ~ s(ID, bs="re") #####
nb01 <- gam(formula = d15N ~ 
              # Fixed effects
              #s(Agemean)+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb01)

###  N ~ s(Agemean) #####
nb02 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb02)

MuMIn::model.sel(nb00, nb01, nb02)

# best: nb00

## Backward model selection with Long, Lat #####
nb20 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space

##### N ~ s(Long,Lat) + re(ID) #####
nb21 <- gam(formula = d15N ~ 
              # Fixed effects
              #s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb21)


##### N ~ s(Agemean) + re(ID) #####
nb22 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb22)

##### N ~ s(Agemean) + s(Long,Lat) #####
nb23 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb23)

##### N ~ s(Long,Lat) #####
nb24 <- gam(formula = d15N ~ 
              # Fixed effects
              #s(Agemean)+
              s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb24)

##### N ~ s(Agemean) #####
nb25 <- gam(formula = d15N ~ 
              # Fixed effects
              s(Agemean),
              #s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5))+
              # Random effects 
              #s(ID, bs="re"),
            # data
            data=N_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nb25)

MuMIn::model.sel(nb20, nb21, nb22, nb23, nb24, nb25)

# best: nb20
# best most simple: nb20


## BEST overall selected model without Long,Lat#####
nbest1 <- gam(formula = d15N ~ 
                # Fixed effects
                s(Agemean)+
                # Random effects 
                s(ID, bs="re"),
              # data
              data=N_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(nbest1)
appraise(nbest1, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(nbest1)
gam.check(nbest1)

## BEST overall selected model with Long, Lat ######
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
draw(nbest2)
gam.check(nbest2)


# Visualisation ######

## VR #####
# basic gratia plot
draw(nbest, select="s(VR)")

# customizable ggplot
sm <- smooth_estimates(nbest, smooth = "s(VR)")

sm %>%
  add_confint() %>%
  ggplot(aes(y = est, x = VR)) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci),
              alpha = 0.2, fill = "forestgreen") +
  geom_line(colour = "forestgreen", size = 1.5) +
  labs(y = "Partial effect",
       title = expression("Partial effect of" ~ f(VR)),
       x = expression("VR"))+
  pub_theme
