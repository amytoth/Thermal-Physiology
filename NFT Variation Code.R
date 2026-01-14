#BUMBLE BEE NFT VARIATION CODE
#Kurtt et al. 2026
#Created on R Studio Version 2023.03.0+386

##Load Packages
  library(ggplot2)
  library(ggforce)
  library(nlraa)
  library(mgcv)
  library(emmeans)
  library(nlme)
  library(readxl)
  
  setwd("/Users/tiger/OneDrive/Documents/Research/Bee Physiology")


#-------------------------------------------------------
# 1.Camera Validation 

#-- Load Data --
  dat<-read.csv("Camera Validation.csv")

#-- Model correlation using a Linear Model --
  ggplot(dat, aes(x= ThermoTemp, y= CameraTemp))+
    geom_point() + 
    geom_smooth(method = "lm")+ 
    labs(title= "Thermal Camera Validation", y= "Camera Temperature (°C)", 
         x= "Internal Temperature (°C)") + 
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  
  model <- lm(CameraTemp ~ ThermoTemp, data = dat)
  summary (model)

#-- Pearson correlation test --
  cor.test(dat$ThermoTemp, dat$CameraTemp, method = "pearson", 
           alternative = "two.sided", conf.level = 0.95)


#------------------------------------------------------------------------------
#2.Caste Comparison Model

#-- Load Data --
  ### CSV input
  qData <- read.csv("queen_nft_data_tidy.csv", stringsAsFactors = TRUE)
  dData <- read.csv("drone_nft_data_tidy.csv", stringsAsFactors = TRUE)
  wData <- read.csv("worker_nft_data_tidy.csv", stringsAsFactors = TRUE)

  ### Reformat data so each csv file corresponds to a given caste--
  dat <- rbind(data.frame(caste = "queen", qData),
               data.frame(caste = "drone", dData),
               data.frame(caste = "worker", wData))
  dat$caste <- as.factor(dat$caste)
  dat <- subset(dat, elapse_min < 10)

  ### Plot to ensure data input--
  ggplot(data = dat, aes(x = elapse_min, y = thoraxtemp, color = caste)) + 
    geom_point() + 
    geom_smooth(se = FALSE)

#-- Using Generalized Additive Models --
  ### Create GAM model
  fmm2 <- gam(thoraxtemp ~ caste + s(elapse_min, by = caste, bs = "cr") + 
                s(beeid, by = caste, bs = "re"),
              data = dat)
  
  #-- 
  prds <- predict_gam(fmm2, newdata = dat,interval = "conf", level = 0.95,
                      exclude = c("s(beeid):castedrone", "s(beeid):castequeen", 
                                  "s(beeid):casteworker"))
  datA <- cbind(dat, prds)
  
  anova(fmm2)

  ### Visualize individual curves 
  ggplot(data = datA, aes(x = elapse_min, y = thoraxtemp, color = caste)) + 
    facet_wrap(~caste) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = caste, color = NULL), 
                alpha = 1/4)

  ### Plot caste curves together --
  ggplot(data = datA, aes(x = elapse_min, y = thoraxtemp, color = caste)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = caste, color = NULL), 
                alpha = 1/4)

#-- Manuscript Figure 2.B -- 
  ggplot(data = datA, aes(x = elapse_min, y = thoraxtemp, color = caste)) + 
    geom_point(size=1.25, alpha=0.5) + 
    geom_line(linewidth= 1,aes(y = Estimate))+
    geom_hline(yintercept = 30, color = "grey30", linetype = "dashed", 
               linewidth=1)+ 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = caste, color = NULL), 
                alpha = 1/4) + 
    scale_color_manual(values = c("drone" = "dodgerblue3", "queen" = "yellow3", 
                                  "worker" = "limegreen"),
                       labels = c("drone" = "Drones", "queen" = "Queens", 
                                  "worker" = "Workers")) +
    scale_fill_manual(values = c("drone" = "dodgerblue3", "queen" = "yellow3", 
                                 "worker" = "limegreen"),
                      labels = c("drone" = "Drones", "queen" = "Queens", 
                                 "worker" = "Workers")) +
    theme_bw() + 
    theme(axis.title.y = element_text(size = 13,margin = margin(r=10)),
          axis.title.x = element_text(size = 13,margin = margin(r=10)))+
    xlab("Time (min)") + ylab(expression(T["est"]~(degree*C))) + 
    labs(color = "") + labs(fill = "")


#-- Run emmeans at different time points for direct comparisons -- 

  ### Time points are intervals of 0.5 minutes from 0 to 10
  ### Following example is for time point 2.5 (2 and a half minutes)
  emmeans(fmm2, ~caste, at = list(elapse_min = 2.5))
  pairs(emmeans(fmm2, ~caste, at = list(elapse_min = 2.5)))
  plot(pairs(emmeans(fmm2, ~caste, at = list(elapse_min = 2.5)))) + 
    geom_vline(xintercept = 2.5)

#-- When do castes reach a given temperature -- 
  timeAt30d <- subset(datA, Estimate > 28 & Estimate < 30.02)
  nrow(timeAt30d)
  aggregate(elapse_min ~ caste, FUN = mean, data = timeAt30d)

  
#------------------------------------------------------------------------------
#3.Incorperating Mass Effects

#-- Load data and format --
  ### Mass
  mass <- read_excel("1. Meta Data.xlsx")
  names(mass) <- c("beeid", "species", "caste", "NFT_date",
                   "collection_date", "site", "wing_score",
                   "ITD", "mass", "NFT_images", "processed",
                   "column1")
  mass$beeid <- as.factor(mass$beeid)
  mass2 <- na.omit(mass[, c(1:6, 8:11)])
  
  ### ITD
  ITD <- read_excel("1. Meta Data.xlsx")
  names(ITD) <- c("beeid", "species", "caste", "NFT_date",
                  "collection_date", "site", "wing_score",
                  "ITD", "mass", "NFT_images", "processed",
                  "column1")
  ITD$beeid <- as.factor(ITD$beeid)
  ITD2 <- na.omit(ITD[, c(1:6, 8:11)])
  
  ### Plot to check data
  ggplot(data = mass2, aes(x = caste, y = mass)) + 
    geom_point()
  ggplot(data = ITD2, aes(x = caste, y = ITD)) + 
    geom_point()
  
  ## Data cleaning
  mass2$caste <- tolower(mass2$caste)
  mass2$beeid <- tolower(gsub("Bimp", "", mass2$beeid))
  
  datAM <- merge(datA, mass2, by = "beeid")
  
  head(datAM)
  summary(mass2$beeid)
  summary(datA$beeid)  

#-- Anova to test differences in mass between castes --  
  maov<- aov(mass2$mass ~ factor(mass2$caste))
  summary(maov)
  TukeyHSD(maov)

#-- Manuscript Figure 2.A --
  ggplot(mass2, aes(x= caste, y=mass, fill= caste)) +
    geom_violin(alpha = 0.5)+
    geom_sina(alpha=0.5)+
    scale_fill_manual(values=c("drone"="dodgerblue3", "queen"="yellow3",
                               "worker"="limegreen"))+
    guides(fill=guide_legend(title="Legend"))+
    labs( x= "", y= "Mass (g)")+
    scale_x_discrete(labels= c("drone" = "Drones", "queen" = "Queens", 
                               "worker" = "Workers"))+
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "pointrange", color = "black")+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 13,margin = margin(r=10)))
  
#-- Comparing mass to intertegular distance --
  ### Load data
  ds<- read.csv("Drone Size.csv")
  ws<- read.csv("Worker Size.csv")
  qs<- read.csv("Queen Size.csv")
  
  ### Pearson tests
  cor.test(ds$ITD, ds$Mass, method = "pearson", alternative = "two.sided", 
           conf.level = 0.95)
  cor.test(ws$ITD, ws$Mass, method = "pearson", alternative = "two.sided", 
           conf.level = 0.95)
  cor.test(qs$ITD, qs$Mass, method = "pearson", alternative = "two.sided", 
           conf.level = 0.95)
  
  ### Plots for workers, drones, and queens respectively
  ggplot(ws, aes(x=ITD, y=Mass))+
    geom_point()+
    geom_smooth(method = "lm")+
    labs(title= "Worker Size Correlation", x= "Intertegular Distance (mm)", 
         y= "Mass (g)")+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  wmod <- lm(Mass~ITD, data = ws)
  summary (wmod)
  
  ggplot(ds, aes(x=ITD, y=Mass))+
    geom_point()+
    geom_smooth(method = "lm")+
    labs(title= "Drone Size Correlation", x= "Intertegular Distance (mm)", 
         y= "Mass (g)")+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  dmod <- lm(Mass~ITD, data = ds)
  summary (dmod)
  
  ggplot(qs, aes(x=ITD, y=Mass))+
    geom_point()+
    geom_smooth(method = "lm")+
    labs(title= "Queen Size Correlation", x= "Intertegular Distance (mm)", 
         y= "Mass (g)")+
    theme_minimal()+
    theme(plot.title = element_text(hjust=0.5))
  qmod <- lm(Mass~ITD, data = qs)
  summary (qmod)
  
#-- Fitting GAM model with mass --
  fmm3 <- gam(thoraxtemp ~ s(elapse_min, by = caste.x) + 
                s(beeid, bs = "re") + 
                s(mass, by = caste.x),
              data = datAM)
  
  anova(fmm3)
  plot(fmm3)
  
#-- Model with interaction terms --
  fmm4 <- gam(thoraxtemp ~ te(elapse_min, mass, by = caste.x) + 
                s(beeid, bs = "re"),
              data = datAM)
  
  anova(fmm4)
  plot(fmm4)
  
#-- Comparing mass size classes of workers and drones --
  ### We used 0.15 grams as it was approximately the lower quartile
  ### for both castes. 
  
  ggplot(data = subset(datAM, caste.x == "drone"),
         aes(x = mass)) + 
    geom_histogram()
  
  ## Create size classes
  datAM$size <- ifelse(datAM$mass < 0.15, "Small", "Large")
  
  ggplot(data = datAM, aes(x = elapse_min, y = thoraxtemp, color = size)) + 
    facet_wrap(~ caste.x) + 
    geom_point()
  
  ## Queens can be omitted as they all weigh more than 0.15 g
  datAM.nq <- subset(datAM[, c(1:11, 16:25)], caste.x != "queen")
  datAM.nq$size <- as.factor(datAM.nq$size)
  datAM.nq$caste_size <- as.factor(with(datAM.nq, paste(caste.x, size, 
                                                        sep = "_")))
  
  ## Create GAM model for drones and workers to compare size class differences
  fmm6 <- gam(thoraxtemp ~ size*caste.x + 
                s(elapse_min, by = caste_size) + 
                s(beeid, bs = "re"), 
              data = datAM.nq)
  anova(fmm6)
  
  prds <- predict_gam(fmm6, interval = "conf", exclude = "s(beeid)", 
                      level = 0.95)
  datAMA.nq <- cbind(datAM.nq, prds)

#-- Manuscript Figure 3 --
  ggplot(data = datAMA.nq, aes(x = elapse_min, y = thoraxtemp, color= size)) + 
    facet_wrap(~ caste.x) + 
    geom_point(size=1, alpha=0.8) + 
    geom_line(linewidth= 0.8,aes(y = Estimate))+
    facet_wrap(~ caste.x, labeller = as_labeller(c("worker"="Workers", 
                                                   "drone"="Drones"))) +
    geom_hline(yintercept = 30, color = "grey30", linetype = "dashed",
               linewidth=1)+ 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = size, color = NULL),
                alpha = 0.25) +
    scale_color_manual(values = c("Small" = "grey60", "Large" = "black")) + 
    scale_fill_manual(values = c("Small" = "grey50", "Large"= "black")) + 
    scale_shape_manual(values = c("Large" = 1, "Small" = 16))+
    theme_bw() + 
    xlab("Time (min)") + ylab(expression(T["est"]~(degree*C))) + 
    labs(color = "") + labs(fill = "")

  
#-- Comparing intertegular distance size classes of workers and drones --
  ### We used 3.16 grams as it was approximately the lower quartile
  ### for both castes.  
  ggplot(data = subset(datAI, caste.x == "drone"),
         aes(x = ITD)) + 
    geom_histogram()
  
  ## Data cleaning
  ITD2$caste <- tolower(ITD2$caste)
  ITD2$beeid <- tolower(gsub("Bimp", "", ITD2$beeid))
  
  datAI <- merge(datA, ITD2, by = "beeid")
  
  head(datAI)
  summary(ITD2$beeid)
  summary(datA$beeid)
  
#-- Fitting GAM model with ITD --  
  fmm7 <- gam(thoraxtemp ~ s(elapse_min, by = caste.x) + 
                s(beeid, bs = "re") + 
                s(ITD, by = caste.x),
              data = datAI)
  
  anova(fmm7)
  
  plot(fmm3)
  
  ## Create size classes
  datAI$size <- ifelse(datAI$ITD < 3.16, "Small", "Large")
  
  ggplot(data = datAI, aes(x = elapse_min, y = thoraxtemp, color = size)) + 
    facet_wrap(~ caste.x) + 
    geom_point()
  
  ## Queens can be omitted as they all larger than 3.16 mm
  datAI.nq <- subset(datAI[, c(1:11, 16:25)], caste.x != "queen")
  datAI.nq$size <- as.factor(datAI.nq$size)
  datAI.nq$caste_size <- as.factor(with(datAI.nq, paste(caste.x, size, sep = "_")))
  
  ## Create GAM model for drones and workers to compare size class differences
  fmm8 <- gam(thoraxtemp ~ size*caste.x + 
                s(elapse_min, by = caste_size) + 
                s(beeid, bs = "re"), 
              data = datAI.nq)
  anova(fmm8)
  
  prds <- predict_gam(fmm7, interval = "conf", exclude = "s(beeid)", level = 0.95)
  datAIA.nq <- cbind(datAI.nq, prds)

#-- Visualization of GAM model for ITD -- 
  ggplot(data = datAIA.nq, aes(x = elapse_min, y = thoraxtemp, color= size)) + 
    facet_wrap(~ caste.x) + 
    geom_point(size=1.25, alpha=0.5) + 
    geom_line(linewidth= 1,aes(y = Estimate))+
    facet_wrap(~ caste.x, labeller = as_labeller(c("worker"="Workers", "drone"="Drones"))) +
    geom_hline(yintercept = 30, color = "grey30", linetype = "dashed", linewidth=1)+ 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = size, color = NULL),alpha = 0.25) +
    scale_color_manual(values = c("Small" = "grey60", "Large" = "black")) + 
    scale_fill_manual(values = c("Small" = "grey50", "Large"= "black")) + 
    scale_shape_manual(values = c("Large" = 1, "Small" = 16))+
    theme_bw() + 
    xlab("Time (min)") + ylab(expression(T["est"]~(degree*C))) + 
    labs(color = "") + labs(fill = "")  
  
  
#------------------------------------------------------------------------------
#4.Seasonal Comparison Model  

#-- Load Data --  
  Season_Data <- read.csv("spring__fall_queen_data.csv", stringsAsFactors = TRUE)
  
  ### Format elapsed minutes to decimal values
  head(Season_Data)
  class(Season_Data$elapse_min)
  Season_Data$elapse_min <- as.character(Season_Data$elapse_min)
  Season_Data$elapse_seconds <- as.numeric(as.difftime(Season_Data$elapse_min, 
                                                       format = "%H:%M:%S", 
                                                       units = "secs"))
  Season_Data$elapse_min_dec <- Season_Data$elapse_seconds / 60
  
#-- Create GAM Model  
  seasongam <- gam(thoraxtemp ~ season + s(elapse_min_dec, by = season) 
                   + s(beeid, bs = "re"), data = Season_Data)
  
  anova(seasongam)
  
  prdz <- predict_gam(seasongam, interval = "conf", level = 0.95, 
                      exclude = "s(beeid)")
  names(prdz)
  datB <- cbind(Season_Data, prdz)
  
  datB <- datB[order(datB$season, datB$elapse_min_dec), ]
  
  head(datB)
  
  "yellow3" %in% colors()
  
  "orchid2" %in% colors()
  
  datB$season <- factor(datB$season, levels = c("fall","spring"), 
                        labels = c("Fall", "Spring"))
  
#-- Manuscript Figure 4B --
  ggplot(data = datB, aes(x = elapse_min_dec, y = thoraxtemp, color = season)) +  
    geom_point(size=1.25, alpha = 0.5) +
    geom_line(linewidth= 1, aes(group = season, y = Estimate)) +
    geom_hline(yintercept = 30, color = "grey30", linetype = "dashed", 
               linewidth=1)+
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = season, color = NULL),
                alpha = 0.25) + 
    scale_color_manual(values = c("Fall" = "yellow3", "Spring" = "orchid2")) +
    scale_fill_manual(values = c("Fall" = "yellow3", "Spring"= "orchid2")) +
    theme_bw() +
    theme(axis.title.y = element_text(size = 13,margin = margin(r=10)),
          axis.title.x = element_text(size = 13,margin = margin(r=10)))+
    xlab("Time (min)") + 
    ylab(expression(T["est"]~(degree*C))) +
    labs(color = "") + labs(fill = "")

#-- Run emmeans at different time points for direct comparisons -- 
  
  ### Time points are intervals of 1 minute from 0 to 10
  pairs(emmeans(seasongam, ~ season, at = list(elapse_min_dec = 0))) 
  
#-- Effects of Mass --
  ### Load data file
  SMass <- read.csv("Season mass.csv")
  
  ### Visualize to ensure input
  ggplot(data = SMass, aes(x = season, y = weight_g)) + 
    geom_point()
  
  ### Test for difference between groups
  t.test(weight_g~season, data=SMass, var.equal = TRUE)
  
  ### Create GAM model
  seasongam2 <- gam(thoraxtemp ~ s(timestamp, by = season) + s(weight_g) + s(beeid, bs = "re"),
                    data = Season_Data, method = "REML")
  anova(seasongam2)

#-- Manuscript Figure 4A --
  ggplot(data=SMass, aes(x= season, y=weight_g, fill= season)) +
    geom_violin(alpha = 0.5)+
    geom_sina(alpha=0.5)+
    scale_fill_manual(values=c("Fall"="yellow3", "Spring"="orchid2"))+
    labs( x= "", y= "Mass (g)")+
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), 
                 geom = "pointrange", color = "black")+
    theme_bw()+
    theme(legend.position = "none",
          axis.title.y = element_text(size = 13,margin = margin(r=10)))

  
  #### -- END -- ####