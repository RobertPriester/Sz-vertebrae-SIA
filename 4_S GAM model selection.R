########## d34S GAMM code for Smooth hammerhead vertebrae SIA manuscript

# Input: 
# - either "VertebraeSIA_RAW" and "Vertebrae_age_reading" data for reproducing all calculations
# - or "VertebraeSIA_Aged" data and already compiled "SEVb" and "Overlap_summary" data to only reproduce manuscript plots with already processed data

# Output:
# - Processed and filtered SIA data with estimated corresponding ages
# - d34S model subset data
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

SIA$Age <- SIA$Agemean+1.5 # centering time on first point

# Calculating tissue date estimates
SIA$tissue_date <- SIA$Date - (months(SIA$FinalAge*12)-months(SIA$Agemean*12)) # calculating age of deposited tissue and subtracting it from capture date
SIA$tissue_j <- julian(SIA$tissue_date)

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
  dplyr::select(ID, d13C, d15N, d34S, VR, Agemean, CalcTL, growth, Sex, SexM, SexF, jDate, Lat, Long)
SIAM$Cmm <- scale_min_max(SIAM$d13C)
SIAM$Nmm <- scale_min_max(SIAM$d15N)
SIAM$Smm <- scale_min_max(SIAM$d34S)

## 2.2 Outlier removal ####

S_Sex <- ggplot(SIA, aes(Sex, d34S))+
  geom_boxplot()+
  ylab(expression(paste(delta^{34}, "S (\u2030)")))+
  pub_theme
S_Sex
# no outliers

### S
S_SIAM <- SIAM


# 3 Model construction ####
## 3.1 Assessing distribution of response variable ####

shapiro.test(S_SIAM$d34S)
# p < 0.05 - not normally distributed

hist(S_SIAM$d34S)

descdist(S_SIAM$d34S, discrete = FALSE)

## 3.2 Basis Models #####

###  Zero model ####
s00 <- gam(formula = d34S ~ 
             # Fixed effects
             1,
           # data
           data=S_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Null space
summary(s00)

### Complete model####
sxx <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             te(Agemean, by=Sex)+
             # Fixed effects
             s(ID, bs="re")+
             s(jDate, bs="re")+
             s(Lat, Long, bs="re"),
           # data
           data=S_SIAM, method="REML", family=scat, select=F)
summary(sxx)
gam.check(sxx)
appraise(sxx, 
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.99,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")
draw(sxx)

## 3.3 Forward model selection #####
#### Fixed Effects #####
##### ** S ~ Age #####
s01 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s01)

##### *** S ~ VR #####
s02 <- gam(formula = d34S ~ 
             # Fixed effects
             s(VR, k=5),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s02)

##### * S ~ growth #####
s03 <- gam(formula = d34S ~ 
             # Fixed effects
             s(growth),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s03)

##### S ~ TL #####
s04 <- gam(formula = d34S ~ 
             # Fixed effects
             s(CalcTL),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s04)

##### S ~ Sex #####
s05 <- gam(formula = d34S ~ 
             # Fixed effects
             SexM + SexF,
           # data
           data=S_SIAM, method="ML", select=F) # select includes double penalty approach - penalizes Sull space
summary(s05)

##### S ~ Date #####
s06 <- gam(formula = d34S ~ 
             # Fixed effects
             s(jDate),
           # data
           data=S_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s06)

##### S ~ Lat,Lon #####
s07 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s07)

##### S ~ ID #####
s08 <- gam(formula = d34S ~ 
             # Fixed effects
             s(ID, bs="fs"),
           # data
           data=S_SIAM, method="ML", scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s08)

MuMIn::model.sel(s00, s01, s02, s03, s04, s05, s06, s07, s08)

##### <<<<< #####

##### * S ~ Agemean + Sex #####
s14 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF,
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s14)

##### S ~ Agemean + s(Agemean by Sex) #####
s15 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF,
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s15)

##### S ~ te(Agemean by Sex) #####
s16 <- gam(formula = d34S ~ 
             # Fixed effects
             te(Agemean, by=Sex),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s16)

##### S ~ ti(Agemean by Sex) #####
s17 <- gam(formula = d34S ~ 
             # Fixed effects
             ti(Agemean, by = Sex),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s17)

##### S ~ ti(Agemean, ID) #####
s18 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF,
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(s18)

MuMIn::model.sel(s01, s14, s15, s16, s17, s18)

##### <<<<<< #######
##### Potentially remove - Capture location should be random factor only #####

##### S ~ Agemean + LatLon #####
s23 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s23)

##### *** S ~ Agemean + Sex + LatLon #####
s24 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s24)

##### S ~ Agemean + s(Agemean by Sex) + LatLon #####
s25 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             s(Long,Lat,  k=10, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s25)

##### S ~ ti(Agemean, ID) #####
s26 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF+
             s(Lat, Long, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s26)


MuMIn::model.sel(s01, s23, s24, s25, s26)

##### <<<<<< #######
##### S ~ Agemean + LatLon #####
s27 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s27)

##### S ~ Agemean + Sex + LatLon #####
s28 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             SexM+
             SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s28)

##### *** S ~ Agemean + s(Agemean by Sex) + LatLon #####
s29 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             s(Agemean, by=Sex, k=8)+
             SexM + SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s29)

##### S ~ ti(Agemean, ID) #####
s30 <- gam(formula = d34S ~ 
             # Fixed effects
             s(Agemean, k=5)+
             t2(Agemean, by=Sex)+
             SexM+SexF+
             ti(Agemean, Long, Lat)+
             s(Long, Lat, k=5, bs="ds", m=c(1, 0.5)),
           # data
           data=S_SIAM, method="ML", family = scat, select=F) # select includes double penalty approach - penalizes Null space
summary(s30)


MuMIn::model.sel(s01, s24, s27, s28, s29, s30)

# best fixed effect models
# with and without lat: S ~ s(Agemean)
sbest <- s01

#### Random effects ######

## Testing best previous models (S01) with different remaining random effects
sr00 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space

##### S01 + re(ID) #####
sr11 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr11)

##### S01 + re(ID, Sex) #####
sr111 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, Sex, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr111)

##### S02 + re(jDate) #####
sr12 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(jDate, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr12)

##### S02 + re(ID, Sex) + re(jDate) #####
sr13 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, bs="re")+
              s(jDate, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr13)

##### S02 + re(Lat, Long) #####
sr14 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(Long, Lat, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr14)

##### *** S02 + re(ID) + re(Lat, Long) #####
sr15 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr15)

##### S02 + re(Sex) + re(Lat, Long) #####
sr16 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(jDate, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr16)

##### S02 + re(ID) + re(Lat, Long) #####
sr17 <- gam(formula = d34S ~ 
              # Fixed effects
              s(Agemean, k=5)+
              # Random effects
              s(ID, bs="re")+
              s(jDate, bs="re")+
              s(Long, Lat, bs="re"),
            # data
            data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sr17)

MuMIn::model.sel(sr00, sr11, sr111, sr12, sr13, sr14, sr15, sr16, sr17)


## BEST overall selected model #####
sbest <- gam(formula = d34S ~ 
               # Fixed effects
               s(Agemean),
             # data
             data=S_SIAM, method="REML", family = scat, select=F) # select includes double penalty approach - penalizes Sull space
summary(sbest)
appraise(sbest,
         method = "simulate",
         n_simulate = 1000,
         type = "deviance",
         level = 0.95,
         point_col = "#003366",
         point_alpha = 0.7,
         ci_alpha = 0.1,
         line_col = "darkred")  # --> no transformation necessary, qq plot good and resid. distributon reasonable
gam.check(sbest)
draw(sbest)

# Visualisation ######

## VR #####
# basic gratia plot
draw(sbest, select="s(VR)")

# customizable ggplot
sm <- smooth_estimates(sbest, smooth = "s(VR)")

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
newdata <- data.frame(VR = seq(min(S_SIAM$VR), max(S_SIAM$VR), length.out = 100))

# Predict the fitted values for SexM=1 and SexF=0
newdata$SexM <- 1
newdata$SexF <- 0
newdata$ID <- (S_SIAM$ID)[1]
newdata$Sex <- levels(S_SIAM$Sex)[1]
newdata$Lat <- mean(S_SIAM$Lat)
newdata$Long <- mean(S_SIAM$Long)
newdata$predM <- predict(sbest, newdata, type="response")

# Predict the fitted values for SexM=0 and SexF=1
newdata$SexM <- 0
newdata$SexF <- 1
newdata$predF <- predict(sbest, newdata, type="response")

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
