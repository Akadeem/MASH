# script_path <- rstudioapi::getActiveDocumentContext()$path
# setwd(dirname(script_path)) # my current Working Directory
rm(list = ls())

re_calib <- F
# To control:
# whether to run the calibration function (T),
# or to use the previously calibrated value directly (F).

# Use a loop to check whether packages are installed in R before loading them ---- # nolint
packages <- c("readxl", "tidyverse")
for (package in packages) {
  if (!requireNamespace(package, quietly = T)) {
    install.packages(package)
  }
}
library(readxl)
library(tidyverse)

# Transition probabilities between states ----
df <- read.csv("data/TPs.csv")
list2env(df, .GlobalEnv)
rm(df)


# Lifetable 2019 ----
age_n <- 1:101

tpDn_male <- read_excel("data/Table02.xlsx", range = "B5:B104", col_names = "prob")
tpDn_male <- rbind(tpDn_male, 1)
tpDn_female <- read_excel("data/Table03.xlsx", range = "B5:B104", col_names = "prob")
tpDn_female <- rbind(tpDn_female, 1)

mat_tpDn <- data.frame("age" = age_n, "male" = tpDn_male$prob, "female" = tpDn_female$prob)

rm(tpDn_male, tpDn_female)

# Risk ratios by age ----
## Table S5 and Table S6
df <- read.csv("data/prevObesity.csv")
prev_obes_lookup <- setNames(df$prevalence, df$agegroup)
age_grp1 <- cut(age_n, breaks = c(0, 1, 5, 11, 19, 39, 59, 101))
prev_obes <- prev_obes_lookup[age_grp1]

df <- read.csv("data/prevDiabetes.csv")
prev_diab_lookup <- setNames(df$prevalence, df$agegroup)
age_grp1 <- cut(age_n, breaks = c(0, 4, 9, 14, 17, 44, 64, 101))
prev_diab <- prev_diab_lookup[age_grp1]



# A multivariable analysis found that obesity and diabetes have relative risk ratios of NASH of 10.07 and 1.78, respectively
# Eskridge W, Vierling JM, Gosbee W, Wan GA, Hyunh ML, Chang HE. Screening for undiagnosed non-alcoholic fatty liver disease (NAFLD) and non-alcoholic steatohepatitis (NASH): A population-based risk factor assessment using vibration controlled transient elastography (VCTE). PLoS One. 2021;16(11):e0260320.
RR_obes <- 10.07
RR_diab <- 1.78

incRR <- (prev_obes * RR_obes + (1 - prev_obes) * 1) *
  (prev_diab * RR_diab + (1 - prev_diab) * 1)


age_grp1 <- cut(age_n, breaks = c(0, 1, 4, 5, 9, 11, 14, 17, 19, 39, 44, 59, 64, 101))
incRR <- setNames(incRR, age_grp1)


# Utility decrements----
df <- read.csv("data/Utilities.csv")
list2env(df, .GlobalEnv)
rm(df)


# Utilities by age ----
mat_uNN_lookup <- read.csv("data/uNN_lookup.csv")

age_grp2 <- cut(age_n, breaks = c(0, 24, 34, 44, 54, 64, 74, 101))
mat_utility <- mat_uNN_lookup[age_grp2, ]


# US population ----
mat_uspop <- read.csv("data/USPopulation.csv")

mat_uspop_10 <- data.frame(
  agegroup = c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
  male = c(
    sum(mat_uspop$male[1:2]), # 0-9
    sum(mat_uspop$male[3:4]), # 10-19
    sum(mat_uspop$male[5:6]), # 20-29
    sum(mat_uspop$male[7:8]), # 30-39
    sum(mat_uspop$male[9:10]), # 40-49
    sum(mat_uspop$male[11:12]), # 50-59
    sum(mat_uspop$male[13:14]), # 60-69
    sum(mat_uspop$male[15:16]), # 70-79
    sum(mat_uspop$male[17:18]) # 80+
  ),
  female = c(
    sum(mat_uspop$female[1:2]),
    sum(mat_uspop$female[3:4]),
    sum(mat_uspop$female[5:6]),
    sum(mat_uspop$female[7:8]),
    sum(mat_uspop$female[9:10]),
    sum(mat_uspop$female[11:12]),
    sum(mat_uspop$female[13:14]),
    sum(mat_uspop$female[15:16]),
    sum(mat_uspop$female[17:18])
  )
)

# Calibration ----
## load the previously calibrated incidence
inc_inuse <- read.csv("data/calibrated_inc.csv")$incidence

## if a re-calibration is wanted, then load the function and run it through
## the re-run result will replace the previously calibrated one
## and be saved as a csv file ("data/calibrated_inc.csv")

if (re_calib) {
  source("function/Calibration2.R")

  calibrated_incs <- NA
  opt_obj <- list()

  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 40973831,
    fn = CalibrateInc2,
    F0 = F,
    F1 = F,
    F2 = F,
    F3 = F,
    F4 = T,
    HCC = F,
    DCC = T,
    LT = F,
    PLT = F,
    AdultOnly = F,
    ExactNum = T,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par
  unique_names <- names(incRR[!duplicated(names(incRR))])
  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[1] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 6447,
    fn = CalibrateInc2,
    F0 = F,
    F1 = F,
    F2 = F,
    F3 = F,
    F4 = F,
    HCC = T,
    DCC = F,
    LT = F,
    PLT = F,
    AdultOnly = F,
    ExactNum = T,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[2] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.015,
    fn = CalibrateInc2,
    F0 = T,
    F1 = T,
    F2 = T,
    F3 = T,
    F4 = T,
    HCC = T,
    DCC = T,
    LT = T,
    PLT = T,
    AdultOnly = F,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[3] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.0645,
    fn = CalibrateInc2,
    F0 = T,
    F1 = T,
    F2 = T,
    F3 = T,
    F4 = T,
    HCC = T,
    DCC = T,
    LT = T,
    PLT = T,
    AdultOnly = F,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[4] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.0267,
    fn = CalibrateInc2,
    F0 = T,
    F1 = T,
    F2 = T,
    F3 = T,
    F4 = T,
    HCC = T,
    DCC = T,
    LT = T,
    PLT = T,
    AdultOnly = F,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[5] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.00178,
    fn = CalibrateInc2,
    F0 = F,
    F1 = F,
    F2 = F,
    F3 = F,
    F4 = T,
    HCC = F,
    DCC = T,
    LT = F,
    PLT = F,
    AdultOnly = T,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[6] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.0012,
    fn = CalibrateInc2,
    F0 = F,
    F1 = F,
    F2 = F,
    F3 = F,
    F4 = T,
    HCC = F,
    DCC = T,
    LT = F,
    PLT = F,
    AdultOnly = T,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[7] <- opt_obj$par


  opt_obj <- try(optim(
    par = 0.0004, # dummy staring parameter
    targetpre = 0.0279,
    fn = CalibrateInc2,
    F0 = T,
    F1 = T,
    F2 = T,
    F3 = T,
    F4 = T,
    HCC = T,
    DCC = T,
    LT = T,
    PLT = T,
    AdultOnly = F,
    ExactNum = F,
    method = "Brent",
    hessian = F,
    lower = 0, upper = 1
  ))

  print(paste0("Calibrated Baseline Incidence = ", signif(opt_obj$par * 100, 3), "%"))
  inc_by_age <- unique(incRR) * opt_obj$par

  inc_by_age <- setNames(sprintf("%.4f%%", inc_by_age * 100), unique_names)
  print(inc_by_age)
  calibrated_incs[8] <- opt_obj$par

  mean(calibrated_incs)
  # plot(x = c(1:8), y = calibrated_incs)

  mean(calibrated_incs[-1])
  # plot(x = c(1:7), y = calibrated_incs[-1])

  mean(calibrated_incs[-(1:2)])
  # plot(x = c(1:7), y = calibrated_incs[-(1:2)])

  mean(calibrated_incs[3:7])
  # plot(x = c(1:7), y = calibrated_incs[3:7])

  inc_inuse <- mean(calibrated_incs[3:7])

  write.csv(data.frame("incidence" = inc_inuse),
    file = "data/calibrated_inc.csv", row.names = FALSE
  )
}

print(paste0("Calibrated Baseline Incidence = ", signif(inc_inuse * 100, 3), "%"))
inc_by_age <- unique(incRR) * inc_inuse
unique_names <- names(incRR[!duplicated(names(incRR))])
inc_by_age <- setNames(sprintf("%.3f%%", inc_by_age * 100), unique_names)
print(inc_by_age)

Total_Inc <- 0
age_grp1 <- cut(age_n, breaks = c(
  0, 4, 9, 14, 19, 24, 29, 34, 39, 44, 49,
  54, 59, 64, 69, 74, 79, 84, 101
))
uspop_1 <- mat_uspop
for (gender in c("male", "female")) {
  uspop_1[1, gender] <- mat_uspop[1, gender] / 4
  uspop_1[2:17, gender] <- mat_uspop[2:17, gender] / 5
  uspop_1[18, gender] <- mat_uspop[18, gender] / 17
  Total_Inc <- Total_Inc + sum(uspop_1[age_grp1, gender] * incRR * inc_inuse)
}
print(paste0("Population Incident Cases = ", round(Total_Inc)))
print(paste0("Population Overall Incidence = ", signif(Total_Inc / sum(mat_uspop[, 2:3]) * 100, 4), "%"))

# QALY Function body ----
source("function/QALY_function_master.R")


# Results of all combinations
scenarios <- c("F0", "No_incidence", "No_NASH")
genders <- c("male", "female")
for (scenario in scenarios) {
  for (gender in genders) {
    name <- paste("Results", scenario, gender, sep = ".")
    assign(name, QALY_function_master(calibrated_inc = inc_inuse, scenario = scenario, gender = gender))
  }
}

# For Fig4 ----

# Reshape the population data to long format
mat_uspop_10_long <- pivot_longer(mat_uspop_10, cols = c(male, female), names_to = "gender", values_to = "population")

# Convert population to millions
mat_uspop_10_long$population <- mat_uspop_10_long$population / 1000000

# Plot the data
Fig4 <-
  ggplot(mat_uspop_10_long, aes(x = agegroup, y = population, fill = gender)) +
  geom_col(data = subset(mat_uspop_10_long, gender == "male"), aes(y = -population), width = 0.62) +
  geom_col(data = subset(mat_uspop_10_long, gender == "female"), width = 0.62) +
  scale_y_continuous(
    labels = function(x) format(abs(x), big.mark = ",", scientific = FALSE),
    breaks = seq(-25, 25, by = 5), limits = c(-25, 25)
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("male" = "#000000", "female" = "#7F7F7F"),
    labels = c("male" = "Male", "female" = "Female")
  ) +
  labs(x = "Age Group", y = "US population (millions)", fill = "") +
  theme_minimal() +
  theme(
    legend.position = c(.98, .98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "black"),
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )


# General function to draw the pyramid plots ----

source("function/pyramid.R")


# For Fig2 ----
mat_F0_with_NASH_male <- Results.F0.male$return.mat
mat_F0_with_NASH_female <- Results.F0.female$return.mat

mat_F0_without_NASH_male <- Results.No_incidence.male$return.mat
mat_F0_without_NASH_female <- Results.No_incidence.female$return.mat

male.morbidity2 <- -(((mat_F0_without_NASH_male$QALYbyAge / mat_F0_without_NASH_male$LEbyAge) - (mat_F0_with_NASH_male$QALYbyAge / mat_F0_with_NASH_male$LEbyAge)) * mat_F0_with_NASH_male$LEbyAge) / 1000
female.morbidity2 <- (((mat_F0_without_NASH_female$QALYbyAge / mat_F0_without_NASH_female$LEbyAge) - (mat_F0_with_NASH_female$QALYbyAge / mat_F0_with_NASH_female$LEbyAge)) * mat_F0_with_NASH_female$LEbyAge) / 1000

male.mortality2 <- -(mat_F0_without_NASH_male$LEbyAge - mat_F0_with_NASH_male$LEbyAge) * (mat_F0_with_NASH_male$QALYbyAge / mat_F0_with_NASH_male$LEbyAge) / 1000
female.mortality2 <- (mat_F0_without_NASH_female$LEbyAge - mat_F0_with_NASH_female$LEbyAge) * (mat_F0_with_NASH_female$QALYbyAge / mat_F0_with_NASH_female$LEbyAge) / 1000

Fig2 <- pyramid(
  mat_F0_with_NASH_male$age_range,
  male.morbidity2,
  female.morbidity2,
  male.mortality2,
  female.mortality2,
  14,
  "Total QALE losses associated with per case NASH",
  2,
  c(.96, .96),
  c("right", "top")
)

# For Fig3 ----

## Fig3.1 Age 0-9 ----
with_incidence_male <- Results.F0.male$return.cycle5
with_incidence_female <- Results.F0.female$return.cycle5

without_incidence_male <- Results.No_incidence.male$return.cycle5
without_incidence_female <- Results.No_incidence.female$return.cycle5

male.morbidity3.1 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.1 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.1 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.1 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.1 <- pyramid(
  Results.F0.male$return.cycle5$rownames.get.cycle_QALY..,
  male.morbidity3.1,
  female.morbidity3.1,
  male.mortality3.1,
  female.mortality3.1,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 0-9",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.2 Age 10-19 ----
with_incidence_male <- Results.F0.male$return.cycle15
with_incidence_female <- Results.F0.female$return.cycle15

without_incidence_male <- Results.No_incidence.male$return.cycle15
without_incidence_female <- Results.No_incidence.female$return.cycle15

male.morbidity3.2 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.2 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.2 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.2 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.2 <- pyramid(
  Results.F0.male$return.cycle15$rownames.get.cycle_QALY..,
  male.morbidity3.2,
  female.morbidity3.2,
  male.mortality3.2,
  female.mortality3.2,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 10-19",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.3 Age 20-29 ----
with_incidence_male <- Results.F0.male$return.cycle25
with_incidence_female <- Results.F0.female$return.cycle25

without_incidence_male <- Results.No_incidence.male$return.cycle25
without_incidence_female <- Results.No_incidence.female$return.cycle25

male.morbidity3.3 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.3 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.3 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.3 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.3 <- pyramid(
  Results.F0.male$return.cycle25$rownames.get.cycle_QALY..,
  male.morbidity3.3,
  female.morbidity3.3,
  male.mortality3.3,
  female.mortality3.3,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 20-29",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.4 Age 30-39 ----
with_incidence_male <- Results.F0.male$return.cycle35
with_incidence_female <- Results.F0.female$return.cycle35

without_incidence_male <- Results.No_incidence.male$return.cycle35
without_incidence_female <- Results.No_incidence.female$return.cycle35

male.morbidity3.4 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.4 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.4 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.4 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.4 <- pyramid(
  Results.F0.male$return.cycle35$rownames.get.cycle_QALY..,
  male.morbidity3.4,
  female.morbidity3.4,
  male.mortality3.4,
  female.mortality3.4,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 30-39",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.5 Age 40-49 ----
with_incidence_male <- Results.F0.male$return.cycle45
with_incidence_female <- Results.F0.female$return.cycle45

without_incidence_male <- Results.No_incidence.male$return.cycle45
without_incidence_female <- Results.No_incidence.female$return.cycle45

male.morbidity3.5 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.5 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.5 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.5 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.5 <- pyramid(
  Results.F0.male$return.cycle45$rownames.get.cycle_QALY..,
  male.morbidity3.5,
  female.morbidity3.5,
  male.mortality3.5,
  female.mortality3.5,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 40-49",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.6 Age 50-59 ----
with_incidence_male <- Results.F0.male$return.cycle55
with_incidence_female <- Results.F0.female$return.cycle55

without_incidence_male <- Results.No_incidence.male$return.cycle55
without_incidence_female <- Results.No_incidence.female$return.cycle55

male.morbidity3.6 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.6 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.6 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.6 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.6 <- pyramid(
  Results.F0.male$return.cycle55$rownames.get.cycle_QALY..,
  male.morbidity3.6,
  female.morbidity3.6,
  male.mortality3.6,
  female.mortality3.6,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 50-59",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.7 Age 60-69 ----
with_incidence_male <- Results.F0.male$return.cycle65
with_incidence_female <- Results.F0.female$return.cycle65

without_incidence_male <- Results.No_incidence.male$return.cycle65
without_incidence_female <- Results.No_incidence.female$return.cycle65

male.morbidity3.7 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.7 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.7 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.7 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.7 <- pyramid(
  Results.F0.male$return.cycle65$rownames.get.cycle_QALY..,
  male.morbidity3.7,
  female.morbidity3.7,
  male.mortality3.7,
  female.mortality3.7,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 60-69",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.8 Age 70-79 ----
with_incidence_male <- Results.F0.male$return.cycle75
with_incidence_female <- Results.F0.female$return.cycle75

without_incidence_male <- Results.No_incidence.male$return.cycle75
without_incidence_female <- Results.No_incidence.female$return.cycle75

male.morbidity3.8 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.8 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.8 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.8 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.8 <- pyramid(
  Results.F0.male$return.cycle75$rownames.get.cycle_QALY..,
  male.morbidity3.8,
  female.morbidity3.8,
  male.mortality3.8,
  female.mortality3.8,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 70-79",
  2,
  c(.96, .04),
  c("right", "bottom")
)

## Fig3.9 Age 80+ ----
with_incidence_male <- Results.F0.male$return.cycle90
with_incidence_female <- Results.F0.female$return.cycle90

without_incidence_male <- Results.No_incidence.male$return.cycle90
without_incidence_female <- Results.No_incidence.female$return.cycle90

male.morbidity3.9 <- -(((without_incidence_male$get.cycle_QALY. / without_incidence_male$get.cycle_LE.) - (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.)) * with_incidence_male$get.cycle_LE.) / 1000
female.morbidity3.9 <- (((without_incidence_female$get.cycle_QALY. / without_incidence_female$get.cycle_LE.) - (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.)) * with_incidence_female$get.cycle_LE.) / 1000

male.mortality3.9 <- -(without_incidence_male$get.cycle_LE. - with_incidence_male$get.cycle_LE.) * (with_incidence_male$get.cycle_QALY. / with_incidence_male$get.cycle_LE.) / 1000
female.mortality3.9 <- (without_incidence_female$get.cycle_LE. - with_incidence_female$get.cycle_LE.) * (with_incidence_female$get.cycle_QALY. / with_incidence_female$get.cycle_LE.) / 1000

Fig3.9 <- pyramid(
  Results.F0.male$return.cycle90$rownames.get.cycle_QALY..,
  male.morbidity3.9,
  female.morbidity3.9,
  male.mortality3.9,
  female.mortality3.9,
  4,
  "Lifetime QALE losses of subjects associated with NASH | Age 80+",
  2,
  c(.96, .04),
  c("right", "bottom")
)


# For Fig5 ----
inc_pat_10 <- matrix(
  data = rep(0, 18),
  nrow = 9,
  dimnames = list(
    c("0-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"),
    c("male", "female")
  )
)

for (gender in c("male", "female")) {
  inc_pat <- uspop_1[age_grp1, gender] * incRR * inc_inuse
  inc_pat_10[1, gender] <- sum(inc_pat[1:9])
  inc_pat_10[2, gender] <- sum(inc_pat[10:19])
  inc_pat_10[3, gender] <- sum(inc_pat[20:29])
  inc_pat_10[4, gender] <- sum(inc_pat[30:39])
  inc_pat_10[5, gender] <- sum(inc_pat[40:49])
  inc_pat_10[6, gender] <- sum(inc_pat[50:59])
  inc_pat_10[7, gender] <- sum(inc_pat[60:69])
  inc_pat_10[8, gender] <- sum(inc_pat[70:79])
  inc_pat_10[9, gender] <- sum(inc_pat[80:101])
}

male.morbidity5 <- male.morbidity2 * inc_pat_10[, "male"]
female.morbidity5 <- female.morbidity2 * inc_pat_10[, "female"]

male.mortality5 <- male.mortality2 * inc_pat_10[, "female"]
female.mortality5 <- female.mortality2 * inc_pat_10[, "female"]

Fig5 <- pyramid(
  mat_F0_with_NASH_male$age_range,
  male.morbidity5,
  female.morbidity5,
  male.mortality5,
  female.mortality5,
  150000,
  "Total QALE losses associated with all incident NASH cases per year",
  50000,
  c(.96, .96),
  c("right", "top")
)


# For Fig6 ----
with_incidence_male <- Results.No_NASH.male$return.pop
with_incidence_female <- Results.No_NASH.female$return.pop

without_incidence_male <- Results.No_incidence.male$return.pop
without_incidence_female <- Results.No_incidence.female$return.pop

male.morbidity6 <- -(((without_incidence_male$QALY_pop / without_incidence_male$LE_pop) - (with_incidence_male$QALY_pop / with_incidence_male$LE_pop)) * with_incidence_male$LE_pop)
female.morbidity6 <- (((without_incidence_female$QALY_pop / without_incidence_female$LE_pop) - (with_incidence_female$QALY_pop / with_incidence_female$LE_pop)) * with_incidence_female$LE_pop)

male.mortality6 <- -(without_incidence_male$LE_pop - with_incidence_male$LE_pop) * (with_incidence_male$QALY_pop / with_incidence_male$LE_pop)
female.mortality6 <- (without_incidence_female$LE_pop - with_incidence_female$LE_pop) * (with_incidence_female$QALY_pop / with_incidence_female$LE_pop)

Fig6 <- pyramid(
  with_incidence_male$age_range,
  male.morbidity6,
  female.morbidity6,
  male.mortality6,
  female.mortality6,
  3000000,
  "Total population QALE losses associated with NASH",
  1000000,
  c(.96, .04),
  c("right", "bottom")
)

# print out the total QALE losses
male_pop_loss <- -sum(male.morbidity6, male.mortality6)
female_pop_loss <- sum(female.morbidity6, female.mortality6)
print(paste0(
  "The total QALE losses among males are ",
  format(male_pop_loss, big.mark = ",")
))
print(paste0(
  "The total QALE losses among females are ",
  format(female_pop_loss, big.mark = ",")
))
print(paste0(
  "The total QALE losses among all US population are ",
  format(male_pop_loss + female_pop_loss, big.mark = ",")
))

# Open a PDF file
pdf("www/my_plots.pdf", width = 6, height = 4)

# Draw the plots

suppressWarnings({
  plot(Fig2)

  plot(Fig3.1)
  plot(Fig3.2)
  plot(Fig3.3)
  plot(Fig3.4)
  plot(Fig3.5)
  plot(Fig3.6)
  plot(Fig3.7)
  plot(Fig3.8)

  plot(Fig3.9)

  plot(Fig4)

  plot(Fig5)

  plot(Fig6)
})

# Close the PDF file
dev.off()

# TIFF format

# Save Fig2, Fig4, Fig5, and Fig6 as separate files
tiff("www/Fig2.tiff", width = 6, height = 4, units = "in", res = 300)
plot(Fig2)
dev.off()

tiff("www/Fig4.tiff", width = 6, height = 4, units = "in", res = 300)
plot(Fig4)
dev.off()

tiff("www/Fig5.tiff", width = 6, height = 4, units = "in", res = 300)
plot(Fig5)
dev.off()

tiff("www/Fig6.tiff", width = 6, height = 4, units = "in", res = 300)
plot(Fig6)
dev.off()

# Load grid related libraries
library(grid)
library(gridExtra)

compositefig3 <- grid.arrange(
  Fig3.1,
  Fig3.2,
  Fig3.3,
  Fig3.4,
  Fig3.5,
  Fig3.6,
  Fig3.7,
  Fig3.8,
  ncol = 2, nrow = 4
)

# Save the output to a TIFF file
ggsave("www/Fig3.tiff", compositefig3, width = 12, height = 16, dpi = 300)
