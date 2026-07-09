#### Validating bee heat budgets with photogrammetry ####

library(ggplot2)
library(dplyr)
library(tidyr)
library(car)
library(readr)
library(tidyverse)
library(emmeans)

df <- read.csv("3D_heat_budget_measurement_data_wide_format.csv")
df <- subset(df, caste=="worker")

#### Do estimated body volumes and surface areas differ from empirical measurements? workers only ####

# Assess normality 
qqPlot(df$head_volume_diff) #normal
qqPlot(df$head_surface_area_diff) # normal
shapiro.test(df$head_surface_area_diff) # normal
qqPlot(df$thorax_volume_diff) # normal
qqPlot(df$thorax_surface_area_diff) # normal
qqPlot(df$abdomen_volume_diff) # normal
shapiro.test(df$abdomen_volume_diff) #normal
qqPlot(df$abdomen_surface_area_diff) # normal
shapiro.test(df$abdomen_surface_area_diff) # normal
shapiro.test(df$abdomen_surface_area_diff)
qqPlot(df$all_tagmata_volume_diff) # normal
shapiro.test(df$all_tagmata_volume_diff)
qqPlot(df$all_tagmata_surface_area_diff) 
shapiro.test(df$all_tagmata_surface_area_diff) # normal

# paired t-tests and wilcoxon tests
t.test(df$head_volume_geom_estimate, df$head_volume_model, paired=TRUE)
t.test(df$head_surface_area_geom_estimate, df$head_surface_area_model, paired=TRUE)
t.test(df$thorax_volume_geom_estimate, df$thorax_volume_model, paired=TRUE)
t.test(df$thorax_surface_area_geom_estimate, df$thorax_surface_area_model, paired=TRUE)
t.test(df$abdomen_volume_geom_estimate, df$abdomen_volume_model, paired=TRUE) #not significant!
t.test(df$abdomen_surface_area_geom_estimate, df$abdomen_surface_area_model, paired=TRUE) #not significant!
t.test(df$all_tagmata_volume_geom_estimate, df$all_tagmata_volume_model, paired=TRUE)
t.test(df$all_tagmata_surface_area_geom_estimate, df$all_tagmata_surface_area_model, paired=TRUE)


#### What is the percent error in the measurements? ####

# Calculate the percent error in means
df <- df %>%
  mutate(
    head_volume_percent_diff = (abs(head_volume_geom_estimate - head_volume_model) / head_volume_model) * 100
  )
df <- df %>%
  mutate(
    head_surface_area_percent_diff = (abs(head_surface_area_geom_estimate - head_surface_area_model) / head_surface_area_model) * 100
  )
df <- df %>%
  mutate(
    thorax_volume_percent_diff = (abs(thorax_volume_geom_estimate - thorax_volume_model) / thorax_volume_model) * 100
  )
df <- df %>%
  mutate(
    thorax_surface_area_percent_diff = (abs(thorax_surface_area_geom_estimate - thorax_surface_area_model) / thorax_surface_area_model) * 100
  )
df <- df %>%
  mutate(
    abdomen_volume_percent_diff = (abs(abdomen_volume_geom_estimate - abdomen_volume_model) / abdomen_volume_model) * 100
  )
df <- df %>%
  mutate(
    abdomen_surface_area_percent_diff = (abs(abdomen_surface_area_geom_estimate - abdomen_surface_area_model) / abdomen_surface_area_model) * 100
  )
df <- df %>%
  mutate(
    all_tagmata_volume_percent_diff = (abs(all_tagmata_volume_geom_estimate - all_tagmata_volume_model) / all_tagmata_volume_model) * 100
  )
df <- df %>%
  mutate(
    all_tagmata_surface_area_percent_diff = (abs(all_tagmata_surface_area_geom_estimate - all_tagmata_surface_area_model) / all_tagmata_surface_area_model) * 100
  )

mean(df$head_volume_percent_diff)
sd(df$head_volume_percent_diff/sqrt(11))

mean(df$thorax_volume_percent_diff)
sd(df$thorax_volume_percent_diff/sqrt(11))

mean(df$all_tagmata_volume_percent_diff)
sd(df$all_tagmata_volume_percent_diff/sqrt(11))

mean(df$head_surface_area_percent_diff)
sd(df$head_surface_area_percent_diff/sqrt(11))

mean(df$thorax_surface_area_percent_diff)
sd(df$thorax_surface_area_percent_diff/sqrt(11))

mean(df$all_tagmata_surface_area_percent_diff)
sd(df$all_tagmata_surface_area_percent_diff/sqrt(11))

#### Are geometric estimates over- or under-estimate? Calculate mean signed error ####

mean(df$head_volume_geom_estimate-df$head_volume_model)
mean(df$thorax_volume_geom_estimate-df$thorax_volume_model)
mean(df$all_tagmata_volume_geom_estimate-df$all_tagmata_volume_model)

mean(df$head_surface_area_geom_estimate-df$head_surface_area_model)
mean(df$thorax_surface_area_geom_estimate-df$thorax_surface_area_model)
mean(df$all_tagmata_surface_area_geom_estimate-df$all_tagmata_surface_area_model)


#### How much does adding legs and wings increase measurements? ####

df$whole_body_volume_diff <- df$whole_body_volume_model-df$all_tagmata_volume_model

df$whole_body_surface_area_diff <- df$whole_body_surface_area_model-df$all_tagmata_surface_area_model

df$whole_body_volume_percent_increase_with_legs_wings <- (df$whole_body_volume_diff/df$whole_body_volume_model)*100

df$whole_body_surface_area_percent_increase_with_legs_wings <- (df$whole_body_surface_area_diff/df$whole_body_surface_area_model)*100

mean(df$whole_body_volume_percent_increase_with_legs_wings)
sd(df$whole_body_volume_percent_increase_with_legs_wings/sqrt(11))

mean(df$whole_body_surface_area_percent_increase_with_legs_wings)
sd(df$whole_body_surface_area_percent_increase_with_legs_wings/sqrt(11))


#### Allometric scaling ####

data <- read.csv("3D_heat_budget_measurement_data_new.csv")

#scaling stats for geometric estimate
#linear model and confidence intervals 
#scaling exponent is 0.574; 95% CI of the slope is (0.40, 0.75) -- hypometric scaling 
#interesting because we "expect" an exponent of 0.66; i.e., for things to scale to the ~2/3 
data_workers <- data %>%
  filter(sex == "W") %>%
  filter(measurement_type == "geom_estimate")

lm_geom <- lm(log10(all_tagmata_surface_area) ~ log10(all_tagmata_volume), data = data_workers)

summary(lm_geom)
confint(lm_geom)

##next up is the same thing for the 3D measures:
#exponent is 0.519; 95% CI slope (0.14, 0.90), also hypometric
data_workers_model <- data %>%
  filter(sex == "W") %>%
  filter(measurement_type == "model")

lm_model <- lm(log10(all_tagmata_surface_area) ~ log10(all_tagmata_volume), data = data_workers_model)

summary(lm_model)
confint(lm_model)

#whole body including legs, wings
#slope is 0.72 (0.27, 1.18)
lm_whole <- lm(log10(whole_body_surface_area_model) ~ log10(whole_body_volume_model), data = data_workers_model)

summary(lm_whole)
confint(lm_whole)

#are the models different from each other? Geometric assumption vs 3D vs whole body
#bind geom and model and whole body model
data_2 <- read.csv("sa_v_clean.csv")

#not significantly different - The scaling exponent is approximately the same for all methods
interaction_model <- lm(log_surface_area ~ log_volume * method, data = data_2)
summary(interaction_model)

additive_model <- lm(log_surface_area ~ log_volume + method, data = data_2)
anova(additive_model, interaction_model)

#testing is slope =2/3 
#do not reject the null hypothesis - so our geometric data is no diff than expectation of 0.66
#same with the other two
linearHypothesis(lm_geom, "log10(all_tagmata_volume_geom_estimate) = 0.6667")
linearHypothesis(lm_model, "log10(all_tagmata_volume_model) = 0.6667")
linearHypothesis(lm_whole, "log10(whole_body_volume_model) = 0.6667")


#### Heat Budget Estimation ####

#This heat budget uses Roberts and Harrison 1999 honey bee heat budget data
#we took the linear models produced for thorax, head, abdomen temperatures
#on air temperature to approximate segments temperatures from 20 to 50 degrees C

#we similarly used the linear models for metabolic heat production (mhp) and evaporative heat loss (ehl)

#Surface area is used in the calculations of long wave radiation (rgain - rloss) (lwr)

#convective heat loss is calculated by summing mhp, ehl, and lwr

#we used the surface areas reported in the Roberts and Harrison 1999 manuscript
#and calculated the low and high values (of lwr, and as a result of chl) based on the percent error of each segment's surface area measure
#that our manuscript reports in the table

data <- read_csv("roberts_model_long.csv")

#separate heat budgets
hb_geom <- data %>%
  filter(route %in% c("mhp", "ehl", "lwr_geom", "chl_geom"))
hb_model <- data %>%
  filter(route %in% c("mhp", "ehl", "lwr_3d", "chl_3d"))

##################################################
# are lines different for LWR? slope and intercept are significantly differ

hb_lwr <- data %>%
  filter(route %in% c("lwr_geom", "lwr_3d"))

model_lwr <- lm(value ~ air_temp * route, data = hb_lwr)
anova(model_lwr)
summary(model_lwr)
#post hoc test - compare slope between routes 
emtrends(model_lwr, pairwise ~ route, var = "air_temp")

#########################################################
# are the lines different for chl?

hb_chl <- data %>%
  filter(route %in% c("chl_geom", "chl_3d"))

model_chl <- lm(value ~ air_temp * route, data = hb_chl)
anova(model_chl) # no significant interaction between air temp and route. Slopes don't differ.
summary(model_chl)
