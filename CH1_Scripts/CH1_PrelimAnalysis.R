#### TL vs GV ####

library(stats)
library(MuMIn)
library(AICcmodavg)
library(interactions)
library(sandwich)
library(broom)
library(broom.mixed)
library(ggstance)
library(jtools)
library(lme4)
library(dplyr)
library(tidyverse)
library(lubridate)
library(sp)
library(geoR)
library(nlme)
library(gstat)
library(data.table)
library(ggplot2)
library(modelr)
library(rcompanion)
library(kimisc)
library(elevatr)


setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Research/PHD_Bio_CSU/CH1_YEWA_telomeres/Analysis/Data")

YEWA=read.csv("YEWA_TL_data/YEWA_TL_GV_NoHY_elev_2021.csv")
YEWA.coord =read.csv("YEWA_TL_data/YEWA_TL_GV_NoHY_2021_coord.csv")

# Make separate "Year," "Month," and "Day" into one "Date" variable
YEWA$Date <- paste(YEWA$Year, YEWA$Month, YEWA$Day, sep=",") %>% ymd() %>% as.Date()

# Make "Date" variable into Julian Date ("Jdate") variable
YEWA$Jday <- yday(YEWA$Date)


# Set specified variable as factors/numeric
as.factor(YEWA$Age)
as.factor(YEWA$Sex)
as.factor(YEWA$Loc_Num)
as.factor(YEWA$Plate)
as.numeric(YEWA$TS_Mean)
as.numeric(YEWA$Elevation_1)
as.numeric(YEWA$Tarsus)
as.numeric(YEWA$Jday)
as.numeric(YEWA$Lat)


### Variable distributions and transformation ####

hist(YEWA$TS_Mean) # right-skewed

plotNormalHistogram(YEWA$TS_Mean) #Right skewed

qqnorm(YEWA$TS_Mean)
qqline(YEWA$TS_Mean, 
       col="red")

T_sqrt = sqrt(YEWA$TS_Mean) # right-skewed
plotNormalHistogram(T_sqrt) # right-skewed
qqnorm(T_sqrt)
qqline(T_sqrt, 
       col="red")

T_log = log(YEWA$TS_Mean) 
plotNormalHistogram(T_log) # normal distribution
qqnorm(T_log)
qqline(T_log, 
       col="red")

YEWA$LogTS <- log(YEWA$TS_Mean) # will use log(TL) as response variable

#### Data visualization for TF85 ####

# LogTS vs GV (TF85)
ggplot(YEWA, aes(x=TF85_1, y=LogTS)) +
  geom_point(mapping = aes(x = TF85_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)

# LogTS vs GV (TF45)
ggplot(data=YEWA, aes(x=TF45_1, y=LogTS)) +
  geom_point(mapping = aes(x = TF45_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)

# LogTS vs GV (TF25)
ggplot(data=YEWA, aes(x=TF25_1, y=LogTS)) +
  geom_point(mapping = aes(x = TF25_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)

# LogTS vs GV (TF85) by Sex
ggplot(data = YEWA, aes(x=TF85_1, y=LogTS, color = Sex)) +
  geom_point(mapping = aes(x = TF85_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)

# LogTS vs GV (TF45) by Sex
ggplot(data = YEWA, aes(x=TF45_1, y=LogTS, color = Sex)) +
  geom_point(mapping = aes(x = TF45_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)

# LogTS vs GV (TF25) by Sex
ggplot(data = YEWA, aes(x=TF25_1, y=LogTS, color = Sex)) +
  geom_point(mapping = aes(x = TF25_1, y = LogTS)) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2)


### TF85 
#### Candidate model set using all combinations of variables up to 4 variables ####
# explanatory variables: TF85, Age, Sex, Lat, Loc_Num, Elevatrion, Tarsus, and Jday
# Response variable is LogTS

TLTF85Mod1 <- lmer(LogTS ~ TF85_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod2 <- lmer(LogTS ~ Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod3 <- lmer(LogTS ~ Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod4 <- lmer(LogTS ~ Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod5 <- lmer(LogTS ~ Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod6 <- lmer(LogTS ~ Elevation_1 + (1|Loc_Num), data=YEWA, REML=T)
TLTF85Mod7 <- lmer(LogTS ~ Tarsus + (1|Loc_Num), data=YEWA, REML=T)
TLTF85Mod8 <- lmer(LogTS ~ TF85_1 * Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod9 <- lmer(LogTS ~ TF85_1 * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod10 <- lmer(LogTS ~ TF85_1 * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod11 <- lmer(LogTS ~ TF85_1 * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod12 <- lmer(LogTS ~ TF85_1 * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod13 <- lmer(LogTS ~ TF85_1 * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod14 <- lmer(LogTS ~ Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod15 <- lmer(LogTS ~ Sex * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod16 <- lmer(LogTS ~ Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod17 <- lmer(LogTS ~ Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod18 <- lmer(LogTS ~ Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod19 <- lmer(LogTS ~ Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod20 <- lmer(LogTS ~ Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod21 <- lmer(LogTS ~ Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod22 <- lmer(LogTS ~ Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod23 <- lmer(LogTS ~ Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod24 <- lmer(LogTS ~ TF85_1 + Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod25 <- lmer(LogTS ~ TF85_1 + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod26 <- lmer(LogTS ~ TF85_1 + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod27 <- lmer(LogTS ~ TF85_1 + Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod28 <- lmer(LogTS ~ TF85_1 + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod29 <- lmer(LogTS ~ TF85_1 + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod30 <- lmer(LogTS ~ Sex + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod31 <- lmer(LogTS ~ Sex + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod32 <- lmer(LogTS ~ Sex + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod33 <- lmer(LogTS ~ Age + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod34 <- lmer(LogTS ~ Age + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod35 <- lmer(LogTS ~ Age + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod36 <- lmer(LogTS ~ Jday + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod37 <- lmer(LogTS ~ Jday + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod38 <- lmer(LogTS ~ Lat + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod39 <- lmer(LogTS ~ Lat + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod40 <- lmer(LogTS ~ TF85_1 * Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod41 <- lmer(LogTS ~ TF85_1 * Sex * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod42 <- lmer(LogTS ~ TF85_1 * Sex * Jday +  (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod43 <- lmer(LogTS ~ TF85_1 * Sex * Elevation_1 +(1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod44 <- lmer(LogTS ~ TF85_1 * Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod45 <- lmer(LogTS ~ TF85_1 * Age * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod46 <- lmer(LogTS ~ TF85_1 * Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod47 <- lmer(LogTS ~ TF85_1 * Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod48 <- lmer(LogTS ~ TF85_1 * Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod49 <- lmer(LogTS ~ TF85_1 * Lat * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod50 <- lmer(LogTS ~ TF85_1 * Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod51 <- lmer(LogTS ~ TF85_1 * Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod52 <- lmer(LogTS ~ TF85_1 * Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod53 <- lmer(LogTS ~ TF85_1 * Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod54 <- lmer(LogTS ~ TF85_1 * Elevation_1 * Tarsus+ (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod50 <- lmer(LogTS ~ TF85_1 + Sex + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod51 <- lmer(LogTS ~ TF85_1 + Sex + Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod52 <- lmer(LogTS ~ TF85_1 + Sex + Jday +  (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod53 <- lmer(LogTS ~ TF85_1 + Sex + Elevation_1 +(1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod54 <- lmer(LogTS ~ TF85_1 + Sex + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod55 <- lmer(LogTS ~ TF85_1 + Age + Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod56 <- lmer(LogTS ~ TF85_1 + Age + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod57 <- lmer(LogTS ~ TF85_1 + Age + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod58 <- lmer(LogTS ~ TF85_1 + Age + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod59 <- lmer(LogTS ~ TF85_1 + Lat + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod60 <- lmer(LogTS ~ TF85_1 + Lat + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod61 <- lmer(LogTS ~ TF85_1 + Lat + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod62 <- lmer(LogTS ~ TF85_1 + Jday + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod63 <- lmer(LogTS ~ TF85_1 + Jday + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod64 <- lmer(LogTS ~ TF85_1 + Elevation_1 + Tarsus+ (1|Loc_Num), data = YEWA, REML=T)
TLTF85Mod65 <- lmer(LogTS ~ (1|Loc_Num), data = YEWA, REML=T)



TLTF85model_list <- list(TLTF85Mod1, TLTF85Mod2, TLTF85Mod3, TLTF85Mod4, TLTF85Mod5, TLTF85Mod6, TLTF85Mod7, TLTF85Mod8, TLTF85Mod9, TLTF85Mod10,
                   TLTF85Mod11, TLTF85Mod12, TLTF85Mod13, TLTF85Mod14, TLTF85Mod15, TLTF85Mod16, TLTF85Mod17, TLTF85Mod18, TLTF85Mod19, TLTF85Mod20,
                   TLTF85Mod21, TLTF85Mod22, TLTF85Mod23, TLTF85Mod24, TLTF85Mod25, TLTF85Mod26, TLTF85Mod27, TLTF85Mod28, TLTF85Mod29,TLTF85Mod30, 
                   TLTF85Mod31, TLTF85Mod32, TLTF85Mod33, TLTF85Mod34, TLTF85Mod35, TLTF85Mod36, TLTF85Mod37, TLTF85Mod38, TLTF85Mod39,
                   TLTF85Mod40, TLTF85Mod41, TLTF85Mod42, TLTF85Mod43, TLTF85Mod44, TLTF85Mod45, TLTF85Mod46, TLTF85Mod47, TLTF85Mod48, TLTF85Mod49,
                   TLTF85Mod50, TLTF85Mod51, TLTF85Mod52, TLTF85Mod53, TLTF85Mod54, TLTF85Mod55, TLTF85Mod56, TLTF85Mod57, TLTF85Mod58,
                   TLTF85Mod59, TLTF85Mod60, TLTF85Mod61, TLTF85Mod62, TLTF85Mod63, TLTF85Mod64, TLTF85Mod65)

## Create AIC model selection table
TL_85_AIC <- aictab(TLTF85model_list, sort = TRUE, second.ord = TRUE)
TL_85_AIC # Model 13 is AICw = 0.81

summary(TLTF85Mod44) 

# LogTL vs GV (TF85) 
ggplot(YEWA, aes(x=TF85_1, y=LogTS)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF85)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# LogTL vs GV (TF85) by Sex
ggplot(YEWA, aes(x=TF45_1, y=LogTS, color = Sex)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", aes(fill=Sex), fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF25)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Log TL vs Tarsus
ggplot(YEWA, aes(x=Tarsus, y=LogTS)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", fullrange=F, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Tarsus Length (mm)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

# Log TL vs Tarsus by Sex
ggplot(YEWA, aes(x=Tarsus, y=LogTS, color = Sex)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", aes(fill = Sex), fullrange=F, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Tarsus Length (mm)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Tarsus by sex
ggplot(YEWA, aes(x=Sex, y=Tarsus)) + 
  geom_boxplot() +
  ylim(18, 23)


#### Results from top model: TLTF85Mod44 ####

sjPlot::tab_model(TLTF85Mod44,
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)", "TF85", "Tarsus", "TF85*Tarsus"),
                  dv.labels= "Genomic Vulnerability TF85 vs Telomere Length") # Conditional R2 = 0.386

# Plot of model TLTF85Mod44
effects_TF85 <- effects::effect(term= "TF85_1", mod= TLTF85Mod44)
summary(effects_TF85)
x_TF85 <- as.data.frame(effects_TF85)
GVTF85_plot <- ggplot() + 
  geom_point(data=YEWA, aes(TF85_1, LogTS)) + 
  geom_point(data=x_TF85, aes(x=TF85_1, y=fit), color="blue") +
  geom_line(data=x_TF85, aes(x=TF85_1, y=fit), color="blue") +
  geom_ribbon(data= x_TF85, aes(x=TF85_1, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="GV (TF85)", y="LogTL")
GVTF85_plot





### TF45
#### Candidate model set using all combinations of variables up to 4 variables ####
# explanatory variables: TF45, Age, Sex, Lat, Loc_Num, and Jday
# Response variable is LogTS

TLTF45Mod1 <- lmer(LogTS ~ TF45_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod2 <- lmer(LogTS ~ Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod3 <- lmer(LogTS ~ Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod4 <- lmer(LogTS ~ Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod5 <- lmer(LogTS ~ Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod6 <- lmer(LogTS ~ Elevation_1 + (1|Loc_Num), data=YEWA, REML=T)
TLTF45Mod7 <- lmer(LogTS ~ Tarsus + (1|Loc_Num), data=YEWA, REML=T)
TLTF45Mod8 <- lmer(LogTS ~ TF45_1 * Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod9 <- lmer(LogTS ~ TF45_1 * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod10 <- lmer(LogTS ~ TF45_1 * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod11 <- lmer(LogTS ~ TF45_1 * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod12 <- lmer(LogTS ~ TF45_1 * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod13 <- lmer(LogTS ~ TF45_1 * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod14 <- lmer(LogTS ~ Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod15 <- lmer(LogTS ~ Sex * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod16 <- lmer(LogTS ~ Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod17 <- lmer(LogTS ~ Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod18 <- lmer(LogTS ~ Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod19 <- lmer(LogTS ~ Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod20 <- lmer(LogTS ~ Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod21 <- lmer(LogTS ~ Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod22 <- lmer(LogTS ~ Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod23 <- lmer(LogTS ~ Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod24 <- lmer(LogTS ~ TF45_1 + Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod25 <- lmer(LogTS ~ TF45_1 + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod26 <- lmer(LogTS ~ TF45_1 + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod27 <- lmer(LogTS ~ TF45_1 + Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod28 <- lmer(LogTS ~ TF45_1 + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod29 <- lmer(LogTS ~ TF45_1 + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod30 <- lmer(LogTS ~ Sex + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod31 <- lmer(LogTS ~ Sex + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod32 <- lmer(LogTS ~ Sex + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod33 <- lmer(LogTS ~ Age + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod34 <- lmer(LogTS ~ Age + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod35 <- lmer(LogTS ~ Age + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod36 <- lmer(LogTS ~ Jday + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod37 <- lmer(LogTS ~ Jday + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod38 <- lmer(LogTS ~ Lat + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod39 <- lmer(LogTS ~ Lat + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod40 <- lmer(LogTS ~ TF45_1 * Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod41 <- lmer(LogTS ~ TF45_1 * Sex * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod42 <- lmer(LogTS ~ TF45_1 * Sex * Jday +  (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod43 <- lmer(LogTS ~ TF45_1 * Sex * Elevation_1 +(1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod44 <- lmer(LogTS ~ TF45_1 * Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod45 <- lmer(LogTS ~ TF45_1 * Age * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod46 <- lmer(LogTS ~ TF45_1 * Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod47 <- lmer(LogTS ~ TF45_1 * Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod48 <- lmer(LogTS ~ TF45_1 * Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod49 <- lmer(LogTS ~ TF45_1 * Lat * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod50 <- lmer(LogTS ~ TF45_1 * Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod51 <- lmer(LogTS ~ TF45_1 * Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod52 <- lmer(LogTS ~ TF45_1 * Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod53 <- lmer(LogTS ~ TF45_1 * Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF45Mod54 <- lmer(LogTS ~ TF45_1 * Elevation_1 * Tarsus+ (1|Loc_Num), data = YEWA, REML=T)

TLTF45model_list <- list(TLTF45Mod1, TLTF45Mod2, TLTF45Mod3, TLTF45Mod4, TLTF45Mod5, TLTF45Mod6, TLTF45Mod7, TLTF45Mod8, TLTF45Mod9, TLTF45Mod10,
                         TLTF45Mod11, TLTF45Mod12, TLTF45Mod13, TLTF45Mod14, TLTF45Mod15, TLTF45Mod16, TLTF45Mod17, TLTF45Mod18, TLTF45Mod19, TLTF45Mod20,
                         TLTF45Mod21, TLTF45Mod22, TLTF45Mod23, TLTF45Mod24, TLTF45Mod25, TLTF45Mod26, TLTF45Mod27, TLTF45Mod28, TLTF45Mod29, TLTF45Mod30, 
                         TLTF45Mod31, TLTF45Mod32, TLTF45Mod33, TLTF45Mod34, TLTF45Mod35, TLTF45Mod36, TLTF45Mod37, TLTF45Mod38, TLTF45Mod39, TLTF45Mod40, 
                         TLTF45Mod41, TLTF45Mod42, TLTF45Mod43, TLTF45Mod44, TLTF45Mod45, TLTF45Mod46, TLTF45Mod47, TLTF45Mod48,  TLTF45Mod49, TLTF45Mod50, 
                         TLTF45Mod51, TLTF45Mod52, TLTF45Mod53, TLTF45Mod54)

## Create AIC model selection table
TL_45_AIC <- aictab(TLTF45model_list, sort = TRUE, second.ord = TRUE)
TL_45_AIC # Model 12 is AICw = 0.79

summary(TLTF45Mod44) 

# LogTL vs GV (TF85) 
ggplot(YEWA, aes(x=TF45_1, y=LogTS)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF45)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.004, 0.010)

# LogTL vs GV (TF85)  by Sex
ggplot(YEWA, aes(x=TF45_1, y=LogTS, color = Sex)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", aes(fill=Sex), fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF45)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.004, 0.010)



# Plot of model 13: GV (TF45) vs LogTL by Sex
effects_TF45 <- effects::effect(term= "TF45_1", mod= TLTF45Mod13)
summary(effects_TF45)
x_TF45 <- as.data.frame(effects_TF45)
GVTF45_plot <- ggplot() + 
  geom_point(data=YEWA, aes(TF45_1, LogTS)) + 
  geom_point(data=x_TF45, aes(x=TF45_1, y=fit), color="blue") +
  geom_line(data=x_TF45, aes(x=TF45_1, y=fit), color="blue") +
  geom_ribbon(data= x_TF45, aes(x=TF45_1, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="GV (TF45)", y="LogTL")
GVTF45_plot



### TF25
#### Candidate model set using all combinations of variables up to 4 variables ####
# explanatory variables: TF25, Age, Sex, Lat, Loc_Num, and Jday
# Response variable is LogTS

TLTF25Mod1 <- lmer(LogTS ~ TF25_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod2 <- lmer(LogTS ~ Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod3 <- lmer(LogTS ~ Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod4 <- lmer(LogTS ~ Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod5 <- lmer(LogTS ~ Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod6 <- lmer(LogTS ~ Elevation_1 + (1|Loc_Num), data=YEWA, REML=T)
TLTF25Mod7 <- lmer(LogTS ~ Tarsus + (1|Loc_Num), data=YEWA, REML=T)
TLTF25Mod8 <- lmer(LogTS ~ TF25_1 * Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod9 <- lmer(LogTS ~ TF25_1 * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod10 <- lmer(LogTS ~ TF25_1 * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod11 <- lmer(LogTS ~ TF25_1 * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod12 <- lmer(LogTS ~ TF25_1 * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod13 <- lmer(LogTS ~ TF25_1 * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod14 <- lmer(LogTS ~ Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod15 <- lmer(LogTS ~ Sex * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod16 <- lmer(LogTS ~ Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod17 <- lmer(LogTS ~ Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod18 <- lmer(LogTS ~ Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod19 <- lmer(LogTS ~ Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod20 <- lmer(LogTS ~ Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod21 <- lmer(LogTS ~ Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod22 <- lmer(LogTS ~ Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod23 <- lmer(LogTS ~ Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod24 <- lmer(LogTS ~ TF25_1 + Sex + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod25 <- lmer(LogTS ~ TF25_1 + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod26 <- lmer(LogTS ~ TF25_1 + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod27 <- lmer(LogTS ~ TF25_1 + Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod28 <- lmer(LogTS ~ TF25_1 + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod29 <- lmer(LogTS ~ TF25_1 + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod30 <- lmer(LogTS ~ Sex + Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod31 <- lmer(LogTS ~ Sex + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod32 <- lmer(LogTS ~ Sex + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod33 <- lmer(LogTS ~ Age + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod34 <- lmer(LogTS ~ Age + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod35 <- lmer(LogTS ~ Age + Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod36 <- lmer(LogTS ~ Jday + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod37 <- lmer(LogTS ~ Jday + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod38 <- lmer(LogTS ~ Lat + Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod39 <- lmer(LogTS ~ Lat + Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod40 <- lmer(LogTS ~ TF25_1 * Sex * Age + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod41 <- lmer(LogTS ~ TF25_1 * Sex * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod42 <- lmer(LogTS ~ TF25_1 * Sex * Jday +  (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod43 <- lmer(LogTS ~ TF25_1 * Sex * Elevation_1 +(1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod44 <- lmer(LogTS ~ TF25_1 * Sex * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod45 <- lmer(LogTS ~ TF25_1 * Age * Lat + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod46 <- lmer(LogTS ~ TF25_1 * Age * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod47 <- lmer(LogTS ~ TF25_1 * Age * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod48 <- lmer(LogTS ~ TF25_1 * Age * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod49 <- lmer(LogTS ~ TF25_1 * Lat * Jday + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod50 <- lmer(LogTS ~ TF25_1 * Lat * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod51 <- lmer(LogTS ~ TF25_1 * Lat * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod52 <- lmer(LogTS ~ TF25_1 * Jday * Elevation_1 + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod53 <- lmer(LogTS ~ TF25_1 * Jday * Tarsus + (1|Loc_Num), data = YEWA, REML=T)
TLTF25Mod54 <- lmer(LogTS ~ TF25_1 * Elevation_1 * Tarsus+ (1|Loc_Num), data = YEWA, REML=T)

TLTF25model_list <- list(TLTF25Mod1, TLTF25Mod2, TLTF25Mod3, TLTF25Mod4, TLTF25Mod5, TLTF25Mod6, TLTF25Mod7, TLTF25Mod8, TLTF25Mod9, TLTF25Mod10,
                         TLTF25Mod11, TLTF25Mod12, TLTF25Mod13, TLTF25Mod14, TLTF25Mod15, TLTF25Mod16, TLTF25Mod17, TLTF25Mod18, TLTF25Mod19, TLTF25Mod20,
                         TLTF25Mod21, TLTF25Mod22, TLTF25Mod23, TLTF25Mod24, TLTF25Mod25, TLTF25Mod26, TLTF25Mod27, TLTF25Mod28, TLTF25Mod29, TLTF25Mod30,
                         TLTF25Mod31, TLTF25Mod32, TLTF25Mod33, TLTF25Mod34, TLTF25Mod35, TLTF25Mod36, TLTF25Mod37, TLTF25Mod38, TLTF25Mod39, TLTF25Mod40, 
                         TLTF25Mod41, TLTF25Mod42, TLTF25Mod43, TLTF25Mod44, TLTF25Mod45, TLTF25Mod46, TLTF25Mod47, TLTF25Mod48, TLTF25Mod49, TLTF25Mod50, 
                         TLTF25Mod51, TLTF25Mod52, TLTF25Mod53, TLTF25Mod54)

## Create AIC model selection table
TL_25_AIC <- aictab(TLTF25model_list, sort = TRUE, second.ord = TRUE)
TL_25_AIC # Model 12 is AICw = 0.22

summary(TLTF25Mod44) 

# GV (TF25) vs LogTL
ggplot(YEWA, aes(x=TF25_1, y=LogTS)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF25)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# GV (TF25) vs LogTL by Sex
ggplot(YEWA, aes(x=TF25_1, y=LogTS, color = Sex)) +
  geom_point(size = 0.9,
             alpha = 0.3) +
  geom_smooth(method = "lm", aes(fill=Sex), fullrange=T, alpha = 0.2) + 
  theme_bw(base_size = 25) +
  scale_color_brewer(type="qual", palette = 6) + labs(x="Genomic Vulnerability (TF25)", y="Log Telomere Length")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


#### Results from top model: TLTF25Mod12 ####

sjPlot::tab_model(TLTF25Mod44,
                  show.re.var= TRUE, 
                  pred.labels =c("(Intercept)", "TF25", "Sex", "TF25*Sex"),
                  dv.labels= "Genomic Vulnerability TF25 vs Telomere Length") # Conditional R2 = 0.319

#Plot of model 44: GV(TF25) vs LogTL by Sex
effects_TF25 <- effects::effect(term= "TF25_1", mod= TLTF25Mod12)
summary(effects_TF25)
x_TF25 <- as.data.frame(effects_TF25)
GVTF25_plot <- ggplot() + 
  geom_point(data=YEWA, aes(TF25_1, LogTS)) + 
  geom_point(data=x_TF25, aes(x=TF25_1, y=fit), color="blue") +
  geom_line(data=x_TF25, aes(x=TF25_1, y=fit), color="blue") +
  geom_ribbon(data= x_TF25, aes(x=TF25_1, ymin=lower, ymax=upper), alpha= 0.3, fill="blue") +
  labs(x="GV (TF25)", y="LogTL")
GVTF25_plot



#### Check for spatial autocorrelation ####

# Correlograms

library(geoR)
library(spdep)
library(gstat)
library(pgirmess)
library(ncf)

#Plot of LogTL across space
TL.matrix <- cbind(YEWA$Long, YEWA$Lat, YEWA$LogTS)
colnames(TL.matrix) <- c("Long", "Lat", "LogTS")
head(TL.matrix, 3)
plot(YEWA$Lat ~ YEWA$Long, pch = 21, bg = gray.colors(12))

#Plot of residuals across space
dat<-data.frame(YEWA$Long,YEWA$Lat,resids=resid(TLTF85Mod12))
coordinates(dat)<-c('YEWA.Long','YEWA.Lat')
bubble(dat,zcol='resids')

# distance matrix
coords <- cbind(YEWA$Long, YEWA$Lat)
colnames(coords) <- c("x", "y")
distmat <- as.matrix(dist(coords))
#maximum distance to consider in correlogram/variogram 
maxdist <- 2/3 * max(distmat) # maxdist = 32.765 decimal degrees


#correlograms with Monte Carlo test
LogTS <- YEWA$LogTS

#Using ncf
correlog.ncf <- ncf::correlog(x = YEWA$Long, y = YEWA$Lat, z = YEWA$LogTS, increment = 5, resamp = 99) 
plot(correlog.ncf) #Variable, but low correlation across distances, with no clear pattern
abline(h = 0)

#spline correlagram 
spline.corr <- spline.correlog(x = YEWA$Long, y = YEWA$Lat, z =
                                 YEWA$LogTS, xmax = maxdist, resamp = 100, type = "boot")
plot (spline.corr) # Very low correlation hovering around 0 across distances

# Monte-Carlo simulation of Moran's I
neigh <- dnearneigh(x = coords, d1 = 0, d2 = 50, longlat = F)
plot(neigh, coordinates(coords))
wts <- nb2listw(neighbours = neigh, style = "W", zero.policy = T)
mor.mc <- moran.mc(x = YEWA$LogTS, listw = wts, nsim = 999, zero.policy = T)
mor.norm <- moran.test(x = YEWA$LogTS, listw = wts, randomisation = F, zero.policy = T)
mor.mc # no correlation
mor.norm # no correlation

# Moran's I
correlog.sp <- data.frame(dist = seq(5, 0.5*max(distmat), by = 5), MoransI = NA, Null.UCL = NA, Pvalue = NA)
for(i in 1: nrow(correlog.sp)) {
  d.start <- correlog.sp[i, "dist"] - 5
  d.end <- correlog.sp[i, "dist"]
  neigh <- dnearneigh(x = coords, d1 = d.start, d2 = d.end, longlat = F)
  wts <- nb2listw(neighbours = neigh, style = 'W', zero.policy = T)
  mor.i <- moran.mc(x = YEWA$LogTS, listw = wts, nsim = 99, zero.policy = T)
  correlog.sp[i, "dist"] <- (d.end + d.start)/2
  correlog.sp[i, "MoransI"] <- mor.i$statistic
  correlog.sp[i, "Null.LCL"] <- quantile(mor.i$res, p=0.025)
  correlog.sp[i, "Null.UCL"] <- quantile(mor.i$res, p=0.975)
  correlog.sp[i, "Pvalue"] <- mor.i$p.value
}
plot(y=correlog.sp$MoransI, x = correlog.sp$dist)
abline(h=0)
lines(correlog.sp$dist, correlog.sp$Null.LCL, col = "red")
lines(correlog.sp$dist, correlog.sp$Null.UCL, col = "red") # no strong correlation



# Variograms

geo.TL <- as.geodata(TL.matrix)
dup.coords(TL.matrix)
geo.TL <- jitterDupCoords(geo.TL, 81)


# empirical variogram using geoR
TL.emp.geoR <- variog(geo.TL, max.dist = maxdist)
plot(TL.emp.geoR) # no change in semivariance with distance

TL.emp.geoR <- variog(geo.TL, max.dist = maxdist,
                   breaks = c(seq(0, maxdist, by = 3))) 
plot(TL.emp.geoR) # no change in semivariance with distance


# variogram using gstat
gstat.TL <- TL.matrix
gstat.TL <- as.data.frame(TL.matrix)
coordinates(gstat.TL) <- c("Long", "Lat")
emp.gstat <- variogram(LogTS ~ 1, cutoff=maxdist, width = 2, data = gstat.TL)
plot(emp.gstat)

# directional variogram using geoR
emp4.geoR <- variog4(geo.TL, maxdist= maxdist)
plot(emp4.geoR)

#directional variogram using gstat
emp4.gstat <- variogram(LogTS ~ 1, cutoff = maxdist, alpha = c(0, 45, 90, 135), gstat.TL)
plot(emp4.gstat)


