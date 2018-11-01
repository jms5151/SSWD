# analyze sea star data ---------------------------------------------------------
rm(list = ls()) # remove everything stored in environment

# load library
library(ordinal)

# load data 
sswd.surveys <- read.csv("Sunflower_data_complete_abundance_score_lower_limit_with_temp.csv", head=T, stringsAsFactors = F)

# format date, month, and year of surveys
sswd.surveys$Date <- as.Date(sswd.surveys$date, "%m/%d/%y")
sswd.surveys$Month <- format(sswd.surveys$Date, "%b")
sswd.surveys$Year <- format(sswd.surveys$Date, "%Y")

# change abundance bins and latitude to factors and scale temperature values (already normally distributed)
sswd.surveys$Abundance <- as.factor(sswd.surveys$abundance)
sswd.surveys$f_lat2 <- as.factor(sswd.surveys$lat2)

# run models with different SST and anomaly metrics
tempMetric <- c("maxsst30", "maxsst60", "maxsst90", "maxsst180", "maxsst360", "maxanom30", "maxanom60", "maxanom90", "maxanom180", "maxanom360")
daysToTemp <- c("maxsst30_days", "maxsst60_days", "maxsst90_days", "maxsst180_days", "maxsst360_days", "maxanom30_days", "maxanom60_days", "maxanom90_days", "maxanom180_days", "maxanom360_days")

for (i in 2:length(tempMetric)){
  sswd.surveys2 <- sswd.surveys
  sswd.surveys2$scaled_stemp2 <- scale(sswd.surveys2[,tempMetric[i]])
  sswd.surveys2$days <- scale(sswd.surveys2[,daysToTemp[i]])
  sswd_clmm_model <- clmm(Abundance ~ scaled_stemp2 + days + Year + (1|f_lat2) + (1|Month), data=sswd.surveys2, link = "probit")
  modName <- paste0("sswd_clmm_", tempMetric[i], ".rds")
  saveRDS(sswd_clmm_model, file = modName) # save model
  cat("finished running", tempMetric[i], "/")
}

# create table of model results for all models
SSWDaicDF <- c()

for (j in 1:length(tempMetric)){
  modName <- paste0("sswd_clmm_", tempMetric[j], ".rds")
  tempMod <- readRDS(modName)
  stat <- as.data.frame(summary(tempMod)$info)
  call <- as.data.frame(as.character(tempMod$call)[2])
  newRow <- cbind(stat, call, tempMetric[j])
  SSWDaicDF <- rbind(SSWDaicDF, newRow)
  assign(tempMetric[j], modName)
}

# Remove days metric for model with lowest AIC (maxanom60) and compare AIC values with and without days 
sswd.surveys2 <- sswd.surveys
sswd.surveys2$scaled_stemp2 <- scale(sswd.surveys2$maxanom60)
sswd_clmm_maxanom60_nodays <- clmm(Abundance ~ scaled_stemp2 + Year + (1|f_lat2) + (1|Month), data=sswd.surveys2, link = "probit")
saveRDS(sswd_clmm_maxanom60_nodays, "sswd_clmm_maxanom60_nodays.rds")
summary(sswd_clmm_maxanom60_nodays)

# Calculate Pseudo-R2
library(rcompanion)
null_model_no_temp <- clmm(Abundance ~ Year + (1|f_lat2) + (1|Month), data = sswd.surveys2)
nagelkerke(fit  = sswd_clmm_maxanom60_nodays, null = null_model_no_temp)