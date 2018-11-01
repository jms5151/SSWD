# Test proportional odds assumption for sea star data -----------------------------------------

# load libraries
library(Hmisc)
library(lme4)

# load data
sswd.surveys <- read.csv("Sunflower_data_complete_abundance_score_lower_limit_with_temp.csv", head=T, stringsAsFactors = F)

# format date, month, and year of surveys
sswd.surveys$Date <- as.Date(sswd.surveys$date, "%m/%d/%y")
sswd.surveys$Month <- format(sswd.surveys$Date, "%b")
sswd.surveys$Year <- format(sswd.surveys$Date, "%Y")

# change abundance bins and latitude to factors and scale temperature values (already normally distributed)
sswd.surveys$Abundance <- as.factor(sswd.surveys$abundance)
sswd.surveys$f_lat2 <- as.factor(sswd.surveys$lat2)

# compare logistic regression across categorical variables
sf <- function(y) {
  c('Y>=0' = qlogis(mean(y >= 0)),
    'Y>=1' = qlogis(mean(y >= 1)),
    'Y>=2' = qlogis(mean(y >= 2)),
    'Y>=3' = qlogis(mean(y >= 3)))
}

s <- with(sswd.surveys, summary(as.numeric(apply) ~ stemp2 + Year + (1|f_lat2) + (1|Month), fun=sf))

ab01 <- subset(sswd.surveys, Abundance == "0"| Abundance == "1")
abund01.logit <- glmer(Abundance ~ scaled_stemp2 + Year + (1|f_lat2) + (1|Month), family = binomial, data = ab01, control = glmerControl(optimizer = "bobyqa"))
summary(abund01.logit)

ab12 <- subset(sswd.surveys, Abundance == "1"| Abundance == "2")
abund12.logit <- glmer(Abundance ~ scaled_stemp2 + Year + (1|f_lat2) + (1|Month), family = binomial, data = ab12, control = glmerControl(optimizer = "bobyqa"))
summary(abund12.logit)

ab23 <- subset(sswd.surveys, Abundance == "2"| Abundance == "3")
abund23.logit <- glmer(Abundance ~ scaled_stemp2 + Year + (1|f_lat2) + (1|Month), family = binomial, data = ab23, control = glmerControl(optimizer = "bobyqa"))
summary(abund23.logit)

ab34 <- subset(sswd.surveys, Abundance == "3"| Abundance == "4")
abund34.logit <- glmer(Abundance ~ scaled_stemp2 + Year + (1|f_lat2) + (1|Month), family = binomial, data = ab34, control = glmerControl(optimizer = "bobyqa"))
summary(abund34.logit)