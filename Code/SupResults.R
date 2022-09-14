# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupResults.R         
# - DATE:        20/07/2021
# - DESCRIPTION: Represent the supplementary results. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(dplyr)
library(tidybayes)
library(bayestestR)
library(ggthemes)
library(cowplot)
library(magrittr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(rstan)
library(ggsci)
library(wesanderson)
library(brms)
library(brmstools)
library(rphylopic)
library(RCurl)
library(png)
library(patchwork)
library(performance)
library(ggridges)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")))

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(ResultPath, "/Models.RData"))
load(paste0(DataPath, "/LPIData.RData"))

# Figure S4: Habitat breadth among groups ######################################

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  ungroup() %>% 
  drop_na(Habitat_breadth_IUCN) %>% 
  group_by(Class) %>% 
  summarise(n=n(),
            median=median(log(Habitat_breadth_IUCN+1)))

(figureS2 <- lpi %>% 
   ggplot(aes(log(Habitat_breadth_IUCN+1), 
              Class)) +
    stat_halfeye(aes(fill=Class, color=Class),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) + 
    geom_text(data = sample_size, 
              aes(y=Class,x=median, label=paste("n=",n)),
              nudge_y = 0.2, nudge_x = 1, colour="grey45")+
   scale_fill_jco(name="Class", alpha = 0.5)+
    scale_colour_jco(name="Class")+
   theme(legend.position = "none") +
   labs(x="log(Habitat breadth)", y=""))

ggsave(filename = "Fig. S2.pdf", figureS2,path = ResultPath,
       height = 6, width = 6)

# Figure S3: Size distribution among groups ####################################

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  ungroup() %>%  
  drop_na(bm_g) %>% 
  group_by(Class) %>% 
  summarise(n=n(),
            mean=median(log(bm_g+1)))

(figureS3 <- lpi %>% 
  ggplot(aes(log(bm_g+1), Class, group=Class)) +
    stat_halfeye(aes(fill=Class, color=Class),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) + 
   scale_fill_jco(name="Class",alpha = 0.5)+
   scale_color_jco(name="Class")+
    geom_text(data = sample_size, 
              aes(y=Class,x=mean, 
                  label=paste("n=",n)),
              nudge_y = 0.2, nudge_x = 5.2, colour="grey45")+
    theme(legend.position = "none") +
  labs(x="log(Body mass)", y=""))

ggsave(filename = "Fig. S3.pdf", figureS3,
       path = ResultPath,
       height = 6, width = 6)

# Table S2: General models #####################################################

tableS2 <- rbind(describe_posterior(sys_1),
                 describe_posterior(clas_1),
                 describe_posterior(tl_1),
                 describe_posterior(bm_1),
                 describe_posterior(lat_1),
                 describe_posterior(pop_dens_1),
                 describe_posterior(hab_1))

#Correct the table

tableS2 <- tableS2 %>% 
  mutate(Parameter=gsub("b_System","", Parameter),
         Parameter=gsub("b_Class","", Parameter),
         Parameter=gsub("b_Trophic_level","", Parameter),
         Parameter=gsub("b_scalelogbm_gP1","Body mass", Parameter),
         Parameter=gsub("b_scaleabs","", Parameter),
         Parameter=gsub("b_","", Parameter),
         Parameter=gsub("b_scalepop_dens","Population density", Parameter),
         Parameter=gsub("b_scaleHabitat_breadth_IUCN",
                        "Habitat breadth", Parameter))

# Save it

setwd(ResultPath)
write.csv2(tableS2, "TableS2.csv",sep = ",")

# Table S3: Specific body mass models ##########################################

tableS3 <- rbind(describe_posterior(bm_amp_fre), 
                 describe_posterior(bm_amp_ter),
                 describe_posterior(bm_bird_fre),
                 describe_posterior(bm_bird_ter),
                 describe_posterior(bm_bird_mar),
                 describe_posterior(bm_fish_fre),
                 describe_posterior(bm_fish_mar),
                 describe_posterior(bm_car_mar),
                 describe_posterior(bm_mam_fre),
                 describe_posterior(bm_mam_ter),
                 describe_posterior(bm_mam_mar),
                 describe_posterior(bm_rep_ter),
                 describe_posterior(bm_rep_fre),
                 describe_posterior(bm_rep_mar))

#Correct the table

tableS3 <- tableS3 %>% 
  mutate(Parameter=gsub("b_scalelogbm_gP1","Body mass", Parameter),
         Parameter=gsub("b_","", Parameter))

# Save it

setwd(ResultPath)
write.csv2(tableS3, "TableS3.csv")

# Table S4: Specific latitude models ##########################################

tableS4 <- rbind(describe_posterior(lat_amp_fre), 
                 describe_posterior(lat_amp_ter),
                 describe_posterior(lat_bird_fre),
                 describe_posterior(lat_bird_ter),
                 describe_posterior(lat_bird_mar),
                 describe_posterior(lat_fish_fre),
                 describe_posterior(lat_fish_mar),
                 describe_posterior(lat_car_mar),
                 describe_posterior(lat_mam_fre),
                 describe_posterior(lat_mam_ter),
                 describe_posterior(lat_mam_mar),
                 describe_posterior(lat_rep_fre),
                 describe_posterior(lat_rep_ter),
                 describe_posterior(lat_rep_mar))

#Correct the table

tableS4 <- tableS4 %>% 
  mutate(Parameter=gsub("scaleabs","", Parameter),
         Parameter=gsub("b_","", Parameter))

# Save it

setwd(ResultPath)
write.csv2(tableS4, "TableS4.csv")

# Figure S1: Phylogenetic signal ###############################################

load(paste0(ResultPath, "/PSModels.RData"))

# Create a common data frame

psign_res <- rbind(data.frame(taxon= "Amhibians", system="Terrestrial", 
                              dist=psign_amp_ter$samples$H1), 
                   # data.frame(taxon= "Amhibians", system="Freshwater", 
                   #            dist=psign_amp_fre$samples$H1), 
                   data.frame(taxon= "Birds", system="Terrestrial", 
                              dist=psign_bird_ter$samples$H1),
                   data.frame(taxon= "Birds", system="Marine", 
                              dist=psign_bird_mar$samples$H1),
                   data.frame(taxon= "Birds", system="Freshwater", 
                              dist=psign_bird_fre$samples$H1),
                   data.frame(taxon= "Mammals", system="Terrestrial", 
                              dist=psign_mam_ter$samples$H1),
                   data.frame(taxon= "Mammals", system="Marine", 
                              dist=psign_mam_mar$samples$H1),
                   data.frame(taxon= "Mammals", system="Freshwater", 
                              dist=psign_mam_fre$samples$H1),
                   data.frame(taxon= "Reptiles", system="Terrestrial", 
                              dist=psign_rep_ter$samples$H1),
                   data.frame(taxon= "Reptiles", system="Freshwater", 
                              dist=psign_rep_fre$samples$H1),
                   data.frame(taxon= "Bony Fishes", system="Freshwater", 
                              dist=psign_fish_fre$samples$H1),
                   data.frame(taxon= "Bony Fishes", system="Marine", 
                              dist=psign_fish_mar$samples$H1),
                   data.frame(taxon= "Cartilaginous Fishes", system="Marine", 
                              dist=psign_car_mar$samples$H1))

# Calculate the median 

res_median <- psign_res %>%
  group_by(taxon, system) %>% 
  summarise(median=median(dist))

# Plot 

(ga <- psign_res %>% 
    ggplot(aes(y=taxon, x=dist, color=taxon)) +
    stat_halfeye(aes(color = taxon,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) +
    geom_text(data=res_median,
              aes(x = median, label = format(round(median, 2), nsmall = 2)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    scale_color_manual(values = wesanderson::wes_palette("Cavalcanti1", 
                                                         type= "continuous",
                                                         n = 6))+
    facet_wrap(.~system)+
    labs(x="Phylogenetic signal", y="")+
    theme(legend.position = "none")+
    xlim(-0.05,0.3))

# Save it

ggsave("Fig. S5.pdf", ga, 
       width = 8, height = 6,
       path = ResultPath)
