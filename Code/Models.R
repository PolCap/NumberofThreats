# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Models.R         
# - DATE:        20/07/2021
# - DESCRIPTION: Fit bayesian models predicting number of stressors. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(tidyverse)
library(brms)
library(data.table)
library(ape)
library(geiger)
library(zoo)
library(dplyr)
library(expss)

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(DataPath, "/LPIData.RData"))

# Create a new variable for species names and order the threats

lpi <- lpi %>% mutate(animal=SpeciesName)

lpi$n.threat <- ordered(lpi$n.threat)

# Separate the data for Class and system for the models 

mam_ter <- lpi %>% filter(Class=="Mammals"&System=="Terrestrial")
mam_mar <- lpi %>% filter(Class=="Mammals"&System=="Marine")
mam_fres <- lpi %>% filter(Class=="Mammals"&System=="Freshwater")
bird_ter <- lpi %>% filter(Class=="Birds"&System=="Terrestrial")
bird_mar <- lpi %>% filter(Class=="Birds"&System=="Marine")
bird_fres <- lpi %>% filter(Class=="Birds"&System=="Freshwater")
amp_ter <- lpi %>% filter(Class=="Amphibians"&System=="Terrestrial")
amp_fres <- lpi %>% filter(Class=="Amphibians"&System=="Freshwater")
rep_ter <- lpi %>% filter(Class=="Reptiles"&System=="Terrestrial")
rep_mar <- lpi %>% filter(Class=="Reptiles"&System=="Marine")
rep_fres <- lpi %>% filter(Class=="Reptiles"&System=="Freshwater")
fish_mar <- lpi %>% filter(Class=="Bony Fish"&System=="Marine")
fish_fres <- lpi %>% filter(Class=="Bony Fish"&System=="Freshwater")
car_mar <- lpi %>% filter(Class=="Cartilaginous Fish"&System=="Marine")

# Set modelling parameters #####################################################

iter <- 10000
thin <- 0.0005*iter
warmup <- 0.1*iter

# Set weakly informed prior

prior <- c(prior(normal(-0.674, 1), class = Intercept, coef = 1),
           prior(normal(0, 1), class = Intercept, coef = 2),
           prior(normal( 0.674, 1), class = Intercept, coef = 3),
           prior(normal(0, 1), class = "b"),
           prior(exponential(1), class = sd))

# Body mass ####################################################################

bm_1 <- brm(brms::bf(n.threat~scale(log10(bm_g+1))+(1|SpeciesName)), 
            iter = iter, thin = thin, warmup = warmup,           
            prior= prior, family= cumulative("probit"),
            control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
            data = lpi, cores=20)

# Latitude #####################################################################

lat_1 <- brm(n.threat~scale(abs(Latitude)) + (1|SpeciesName), 
             iter = iter, thin = thin, warmup = warmup,           
             prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
            data = lpi, cores=20)

# General model System #########################################################

sys_1 <- brm(n.threat~System + (1|SpeciesName), 
             iter = iter, thin = thin, warmup = warmup,           
             prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
             data = lpi,
             cores=20)

# General model Class #########################################################

clas_1 <- brm(n.threat~Class + (1|SpeciesName), 
             iter = iter, thin = thin, warmup = warmup,           
             prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
             data = lpi,
             cores=20)

# Habitat specificity ############################################

hab_1 <- brm(n.threat~ scale(log(Habitat_breadth_IUCN+1))+ (1|SpeciesName), 
             iter = iter, thin = thin, warmup = warmup,
             prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
             data = lpi,
             cores=20)

# Trophic level ################################################################

tl_1 <- brm(n.threat~ Trophic_level + (1|SpeciesName), 
             iter = iter, thin = thin, warmup = warmup,
             prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2, 
             data = lpi,
             cores=20)

# Population density ###########################################################

pop_dens_1 <- brm(n.threat~ scale(log(pop_dens+1)) + (1|SpeciesName), 
                  iter = iter, thin = thin, warmup = warmup,           
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = lpi,
                  cores=20)

# Body size Interactions, also with system and Class ###########################

bm_mam_ter <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_ter, cores=20)
bm_mam_mar <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_mar, cores=20)
bm_mam_fre <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_fres, cores=20)
bm_bird_ter <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = bird_ter, cores=20)
bm_bird_mar <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = bird_mar, cores=20)
bm_bird_fre <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = bird_fres, cores=20)
bm_rep_ter <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_ter, cores=20)
bm_rep_mar <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_mar, cores=20)
bm_rep_fre <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_fres, cores=20)
bm_amp_ter <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = amp_ter, cores=20)
bm_amp_fre <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = amp_fres, cores=20)
bm_fish_mar <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = fish_mar, cores=20)
bm_fish_fre <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = fish_fres, cores=20)
bm_car_mar <- brm(n.threat~scale(log10(bm_g+1))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = car_mar, cores=20)


# Latitude interactions, system and Class ######################################

lat_mam_ter <- brm(n.threat~scale(abs(Latitude)) +(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_ter, cores=20)
lat_mam_mar <- brm(n.threat~scale(abs(Latitude)) +(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_mar,  cores=20)
lat_mam_fre <- brm(n.threat~scale(abs(Latitude)) +(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = mam_fres, cores=20)
lat_bird_ter <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                   data = bird_ter,cores=20)
lat_bird_mar <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                   data = bird_mar, cores=20)
lat_bird_fre <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                   data = bird_fres, cores=20)
lat_rep_ter <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_ter, cores=20)
lat_rep_mar <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_mar, cores=20)
lat_rep_fre <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = rep_fres, cores=20)
lat_amp_ter <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = amp_ter, cores=20)
lat_amp_fre <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = amp_fres, cores=20)
lat_fish_mar <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                   data = fish_mar, cores=20)
lat_fish_fre <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                   data = fish_fres, cores=20)
lat_car_mar <- brm(n.threat~scale(abs(Latitude))+(1|SpeciesName),
                  iter = iter, thin = thin, warmup = warmup,
                  prior= prior, family= cumulative("probit"), control = list(adapt_delta = .975, max_treedepth = 20), init_r = 0.2,
                  data = car_mar, cores=20)

# Save everything ##############################################################

setwd(ResultPath)
save(list = ls(pattern="_"),file = "Models.RData") 


