rm(list = ls())
library(truncnorm)
library(readxl)
library(ggplot2)
library(tidyr)
library(MASS)
library(mvtnorm)
library(invgamma)
<<<<<<< HEAD

don <- read_xlsx("E:/Phd_Paris Saclay/10Sep2024_Cours Calibration/Ball_drops_data.xlsx", sheet = 2)
=======
library(here)

#don <- read_xlsx("/Users/negar/Documents/phd/estimate_mixture_models-main/data/Ball_drops_data.xlsx", sheet = 2)
don <- read_xlsx(here("data", "Ball_drops_data.xlsx"), sheet = 2)

>>>>>>> 8450abc4c47461f1601519d1d05674d6497dd7ed
names(don) <- c("drop", "time", "Height", "Velocity")
don$drop <- as.factor(don$drop)
don <- don[don$drop == 1, ]

t <- don$time
y <- don$Height
length(t)
t_min <- min(t)
t_range <- max(t) - min(t)
t <- (t - t_min) / t_range
<<<<<<< HEAD
n <- length(y)
=======
n <- length(y)
a <- t_range



>>>>>>> 8450abc4c47461f1601519d1d05674d6497dd7ed
