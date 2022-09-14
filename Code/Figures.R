# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        20/07/2021
# - DESCRIPTION: Represent the results. 
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
library(ggrepel)
library(ggsci)
library(wesanderson)
library(RColorBrewer)
library(posterior)
library(brms)
library(rphylopic)
library(RCurl)
library(png)
library(patchwork)

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

## Add an alpha value to a colour

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(ResultPath, "/Models.RData"))
load(paste0(DataPath, "/LPIData.RData"))

# Figure 1: Global map #########################################################

# First we summarise the Number of threats by country #

lpi_map <- lpi %>% 
  group_by(Country) %>% 
  summarise(mean_threats=mean(n.threat, na.rm=T)) %>% 
  rename(region = Country)  # Rename columns

# Now by each 5 of latitude

lpi_lat <- lpi %>% 
  mutate(Lat=cut(Latitude, 
                 breaks = seq(-80,85, by=5), include.lowest = T, 
                 labels = as.character(seq(-80,80, by=5)))) %>% 
  group_by(Lat) %>% 
  summarise(n=n(), 
            mean_threats=mean(n.threat, na.rm=T),
            se_threats=plotrix::std.error(n.threat, na.rm = T))

# World map data

world_map <- map_data("world")

# Data on countries 

lpi_map <- lpi_map %>%
  mutate(region = ifelse(region=="United Kingdom", "UK",region),
         region = ifelse(region=="French Southern Territories", 
                         "French Southern and Antarctic Lands", region),
         region = ifelse(region=="Antigua and Barbuda", "Antigua", region),
         region = ifelse(region=="Bahamas, The", "Bahamas", region),
         region = ifelse(region=="Cote d'Ivoire", "Ivory Coast", region),
         region = ifelse(region=="Congo, Dem. Rep.", 
                         "Democratic Republic of the Congo", region),
         region = ifelse(region=="Congo, Rep.", "Republic of Congo", region),
         region = ifelse(region=="Cabo Verde", "Cape Verde", region),
         region = ifelse(region=="Egypt, Arab Rep.", "Egypt", region),
         region = ifelse(region=="Falkland Islands (Malvinas)", 
                         "Falkland Islands", region),
         region = ifelse(region=="Heard Island And McDonald Islands", 
                         "Heard Island", region),
         region = ifelse(region=="Iran, Islamic Rep.", "Iran", region),
         region = ifelse(region=="St. Kitts and Nevis", "Nevis", region),
         region = ifelse(region=="Korea, Rep.", "South Korea", region),
         region = ifelse(region=="Lao PDR", "Laos", region),
         region = ifelse(region=="St. Lucia", "Saint Lucia", region),
         region = ifelse(region=="North Macedoni", "Macedonia", region),
         region = ifelse(region=="Russian Federation", "Russia", region),
         region = ifelse(region=="South Georgia And The South Sandwich Islands", 
                         "South Georgia", region),
         region = ifelse(region=="Saint Helena, Ascension And Tristan Da Cunha", 
                         "Saint Helena", region),
         region = ifelse(region=="Slovak Republic", "Slovakia", region),
         region = ifelse(region=="Trinidad and Tobago", "Trinidad", region),
         region = ifelse(region=="Taiwan, Province Of China", "Taiwan", region),
         region = ifelse(region=="United States", "USA", region),
         region = ifelse(region=="Venezuela, RB", "Venezuela", region),
         region = ifelse(region=="Virgin Islands (U.S.)", "Virgin Islands", region)) 

# Change name of some countries to match the ones in world map
# Join with global data 

lpi_map <- world_map %>% left_join(lpi_map, by = "region")

# Plot it

(g1 <- ggplot(lpi_map, aes(long, lat, group = group))+
  geom_polygon(aes(fill = mean_threats), 
               color = "black")+
  scale_fill_viridis_c("Number of threats", 
                       option = "A",direction = -1) +
  labs(y="Latitude", x="Longitude") +
    scale_y_continuous(expand = expansion(mult = c(0, 0)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.position = "none"))

(g2 <- lpi_lat %>% ggplot(aes(y=mean_threats, x=Lat, 
                              fill=mean_threats))+
    geom_bar(stat = "Identity") + 
    geom_errorbar(aes(ymin= mean_threats-se_threats, 
                      ymax = mean_threats+se_threats), width = 0.5)+
    geom_text(aes(y = 2.7, 
                  label=paste("n=", n)),
              hjust=0.8)+
    scale_fill_viridis_c(option = "A",
                         direction = -1, breaks=c(0,1,2,3),
                         limits=c(0,3)) +
    scale_y_continuous(limits=c(0,3), 
                       expand = expansion(mult = c(0, 0)))+
    scale_x_discrete(breaks=c(80,50,0,-50,-80))+
    labs(x="", y="Number of threats")+
    coord_flip()+
    theme(legend.position = "none"))

legend <- get_legend(g1+theme(legend.position = "bottom",
                              legend.direction = "horizontal"))

row1 <- plot_grid(g1, g2,
          align = "hv",
          labels = c('a', 'b'),
          label_size = 12, rel_widths = c(1,0.5))

(figure1 <- plot_grid(row1, legend, 
                     ncol = 1, 
                     rel_heights = c(1,0.1)))

ggsave("Figure1.pdf", figure1,
       width = 16, height = 8,
       path = ResultPath)


# Figure 2: Missing values ##################################################### 

## Panel a: Proportion of missing values ---------------------------------------- 

(g2 <- lpi %>% ungroup() %>%
    distinct(ID, .keep_all=T) %>% 
    select(bm_g,pop_dens,
           Habitat_breadth_IUCN,
           Trophic_level, 
           Latitude) %>%  # replace to your needs
    rename(BodyMass=bm_g,
           HabitatBreadth=Habitat_breadth_IUCN,
           TrophicLevel=Trophic_level,
           PopulationDensity=pop_dens) %>% 
    gather(key = "key", value = "val") %>%
    mutate(val = na_if(val, NaN),
           isna = is.na(val),
           key=factor(key,levels=c("Latitude", 
                                   "BodyMass","PopulationDensity", 
                                   "HabitatBreadth", "TrophicLevel"))) %>%
    group_by(key) %>%
    mutate(total = n()) %>%
    group_by(key, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100) %>% 
    ggplot() +
    geom_bar(aes(x = key, 
                 y = pct, fill=isna), 
             stat = 'identity', alpha=0.8) +
    scale_fill_manual(name = "", 
                      values = c('steelblue', 'tomato3'), 
                      labels = c("Present", "Missing")) +
    coord_flip() +
    theme(legend.position = "none")+
    scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    labs(x ='', y = "% of missing values"))

## Panel b: Raster with combinations of missigness ------------------------------ 

(g3 <-  lpi %>%ungroup() %>% 
    distinct(ID, .keep_all=T) %>% 
    mutate(id=row_number()) %>% 
    select(id, bm_g,pop_dens,
           Habitat_breadth_IUCN,
           Trophic_level, Latitude) %>% 
    rename(BodyMass=bm_g,
           HabitatBreadth=Habitat_breadth_IUCN,
           TrophicLevel=Trophic_level,
           PopulationDensity=pop_dens) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(val = na_if(val, NaN),
           isna = is.na(val),
           key=factor(key,levels=c("Latitude",
                                   "BodyMass", "PopulationDensity", 
                                   "HabitatBreadth", "TrophicLevel"))) %>% 
    ggplot(aes(key, reorder(id, desc(isna)), fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
                      values = c('steelblue', 'tomato3'),
                      labels = c("Present", "Missing")) +
    labs(x = "",
         y = "Population") +
    coord_flip()+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))

# Combine figures 

legend <- get_legend(g3+theme(legend.text = element_text(size = 14)))

(figure2 <- plot_grid(g2+theme(plot.margin = margin(5,0,0,0)), 
                      g3+theme(legend.position = "none", 
                               plot.margin = margin(5,0,0,0)), 
                      nrow = 1, align = "h", axis = "b", 
                      labels = "auto"))

(figure2 <- plot_grid(figure2, legend, rel_widths = c(1, 0.1)))

ggsave(figure2,filename = "Figure 2.pdf", 
       width = 12, height = 4,
       path = ResultPath)

# Figure 3: Factors ##################################################
## Panel a: System -----------------------------------------------------

# Determine sample size

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  filter(!is.na(System)) %>% 
  group_by(System) %>% 
  summarise(n=n()) %>% 
  left_join(bind_rows(as_tibble(conditional_effects(sys_1)$System)
                      %>% mutate(.width=0.95)))

# Plot

(g3a <- 
    bind_rows(as_tibble(conditional_effects(sys_1)$System)
              %>% mutate(.width=0.95),
              as_tibble(conditional_effects(sys_1, prob = 0.8)$System)
              %>% mutate(.width=0.8),
              as_tibble(conditional_effects(sys_1,prob = 0.5)$System)
              %>% mutate(.width=0.5)) %>% 
  rename(estimate=estimate__,
         .lower=lower__,
         .upper=upper__) %>% 
  mutate(.point="median",.interval="qi") %>% 
  select(System,estimate, .lower,.upper, .width, .point, .interval) %>% 
  ggplot(aes(x=estimate, 
             y=reorder(System, -estimate), 
             colour=System))+
  geom_pointinterval(aes(xmin=.lower, xmax=.upper)) +
  scale_x_continuous(breaks = c(1,2),                        
                     limits = c(1,2.1))+
  scale_colour_manual(values=c("#A1D6E2","#336B87", "#CC954E"))+
  labs(x="Number of threats", y = "System") +
  geom_text(data = sample_size, 
            aes(x=1, y=System, label = paste0("n=", n)),
            hjust=-.1)+
  theme(legend.position = "none"))


## Panel b: Class -----------------------------------------------------

# Determine sample size

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  filter(!is.na(Class)) %>% 
  group_by(Class) %>% 
  summarise(n=n()) %>% 
  left_join(bind_rows(as_tibble(conditional_effects(clas_1)$Class)
                      %>% mutate(.width=0.95)))

# Plot

(g3b <- 
    bind_rows(as_tibble(conditional_effects(clas_1)$Class)
              %>% mutate(.width=0.95),
              as_tibble(conditional_effects(clas_1, prob = 0.8)$Class)
              %>% mutate(.width=0.8),
              as_tibble(conditional_effects(clas_1,prob = 0.5)$Class)
              %>% mutate(.width=0.5)) %>% 
    rename(estimate=estimate__,
           .lower=lower__,
           .upper=upper__) %>% 
    mutate(.point="median",.interval="qi") %>% 
    select(Class,estimate, .lower,.upper, .width, .point, .interval) %>% 
    ggplot(aes(x=estimate, 
               y=reorder(Class, -estimate), 
               colour=Class))+
    geom_pointinterval(aes(xmin=.lower, xmax=.upper)) +
    scale_x_continuous(breaks = c(1,2,3),                        
                       limits = c(0.5,3))+
    scale_colour_jco(name="Class")+
    labs(x="Number of threats", y = "Class") +
    geom_text(data = sample_size, 
              aes(x=0.5, y=Class, label = paste0("n=", n)),
              hjust=-.1)+
    theme(legend.position = "none"))

## Panel c: Trophic level -------------------------------------------------------

# Determine sample size

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  filter(!is.na(Trophic_level)) %>% 
  group_by(Trophic_level) %>% 
  summarise(n=n()) %>% 
  left_join(bind_rows(as_tibble(conditional_effects(tl_1)$Trophic_level)
                      %>% mutate(.width=0.95)))

# Plot

(g3c <- 
    bind_rows(as_tibble(conditional_effects(tl_1)$Trophic_level)
              %>% mutate(.width=0.95),
              as_tibble(conditional_effects(tl_1, prob = 0.8)$Trophic_level)
              %>% mutate(.width=0.8),
              as_tibble(conditional_effects(tl_1,prob = 0.5)$Trophic_level)
              %>% mutate(.width=0.5)) %>% 
    rename(estimate=estimate__,
           .lower=lower__,
           .upper=upper__) %>% 
    mutate(.point="median",.interval="qi") %>% 
    select(Trophic_level,
           estimate, .lower,.upper, .width, .point, .interval) %>% 
    ggplot(aes(x=estimate, 
               y=reorder(Trophic_level, -estimate), 
               colour=Trophic_level))+
    geom_pointinterval(aes(xmin=.lower, xmax=.upper)) +
    scale_x_continuous(breaks = c(1,2),                        
                       limits = c(0.5,2.1))+
    scale_colour_manual(values = wes_palette("Cavalcanti1"))+
    labs(x="Number of threats", y = "Trophic level") +
    geom_text(data = sample_size, 
              aes(x=0.5, y=Trophic_level, label = paste0("n=", n)),
              hjust=-.1)+
    theme(legend.position = "none"))


## Panel d: Body mass ------------------------------------------------------
# Use the conditional effects function to predict the 

p <- conditional_effects(bm_1, spaghetti = T,
                         ndraws=500)

# Separate into a different data set

data <- p$bm_g

# Generate the individual draws 

draws <- attributes(p$bm_g)$spaghetti

# Generate the plot 

(g3d <- draws %>% 
    ggplot(aes(x = log10(bm_g+1), 
               y = estimate__)) +
    geom_line(color="#216894", alpha=0.02, aes(group=sample__))+
     geom_line(data = data, aes(y=estimate__, 
                                x=log10(bm_g+1)),
               colour="#216894", size=1)+
    scale_y_continuous(breaks = c(0,1,2,3),                        
                       limits = c(0,3))+
    labs(x = "Log10(Body mass+1)",
         y = "Number of threats")) 

## Panel e: Latitude ------------------------------------------------------

# Use the conditional effects function to predict the 

p <- conditional_effects(lat_1, spaghetti = T,
                         ndraws=500)

# Separate into a different data set

data <- p$Latitude

# Generate the individual draws 

draws <- attributes(p$Latitude)$spaghetti

# Generate the plot 

(g3e <- draws %>% 
    ggplot(aes(x = abs(Latitude), 
               y = estimate__)) +
    geom_line(color="#949228", alpha=0.02, aes(group=sample__))+
    geom_line(data = data, aes(y=estimate__, 
                               x=abs(Latitude)),
              colour="#949228", size=1)+
    scale_y_continuous(breaks = c(0,1,2,3),                        
                       limits = c(0,3))+
    labs(x = "Absolute latitude",
         y = "Number of threats")) 

## Panel f: Population density --------------------------------------------------

# Use the conditional effects function to predict the 

p <- conditional_effects(pop_dens_1, spaghetti = T,
                         ndraws=500)

# Separate into a different data set

data <- p$pop_dens

# Generate the individual draws 

draws <- attributes(p$pop_dens)$spaghetti

# Generate the plot 

(g3f <- draws %>% 
    ggplot(aes(x = pop_dens, 
               y = estimate__)) +
    geom_line(color="#0396A6", alpha=0.02, aes(group=sample__))+
    geom_line(data = data, aes(y=estimate__, 
                               x=pop_dens),
              colour="#0396A6", size=1)+
    scale_y_continuous(breaks = c(0,1,2,3),                        
                       limits = c(0,3))+
    labs(x = "Population density",
         y = "Number of threats")) 

## Panel g: Habitat breath ------------------------------------------------------

# Use the conditional effects function to predict the 

p <- conditional_effects(hab_1, spaghetti = T,
                         ndraws=500)

# Separate into a different data set

data <- p$Habitat_breadth_IUCN

# Generate the individual draws 

draws <- attributes(p$Habitat_breadth_IUCN)$spaghetti

# Generate the plot 

(g3g <- draws %>% 
    ggplot(aes(x = Habitat_breadth_IUCN, 
               y = estimate__)) +
    geom_line(color="#A6323B", alpha=0.02, aes(group=sample__))+
    geom_line(data = data, aes(y=estimate__, 
                               x=Habitat_breadth_IUCN),
              colour="#A6323B", size=1)+
    scale_y_continuous(breaks = c(0,1,2,3),                        
                       limits = c(0,3))+
    labs(x = "Habitat breadth",
         y = "Number of threats")) 

## Combine figure 3 -------------------------------------------------------------

(row1 <- g3a+g3b+g3c+plot_layout(nrow = 1)+plot_annotation(tag_levels = "a")& 
   theme(plot.tag = element_text(face = "bold"),
         plot.margin = margin(0,0,0,0)))
(row2 <- g3d+g3e+ylab("")+g3f+ylab("")+g3g+ylab("")&
    plot_annotation(tag_levels = list(c("d", "e", "f", "g")))&
    plot_layout(nrow=1)& 
    theme(plot.tag = element_text(face = "bold"),
          plot.margin = margin(0,0,0,0)))

 (figure3 <- plot_grid(row1,
                        row2, 
                        nrow = 2))

ggsave(figure3, filename = "Figure3.pdf",
       height = 8, width = 12,path = ResultPath)

# Figure 4: Body mass ##########################################################

# Sample size

sample_size <- lpi %>% 
  filter(!is.na(bm_g)) %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(System, Class) %>% 
  summarise(n=n()) %>% 
  filter(n>1) %>% 
  rename(taxon=Class,
         system=System)

# Do summaries for each model

mu <- as_draws_df(bm_amp_fre) %>% 
  pivot_longer(ends_with("_gP1")) %>% 
  mutate(name=factor(name),
         system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(bm_amp_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(bm_bird_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(bm_bird_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Birds"))%>% 
  bind_rows(bm_bird_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Birds"))%>% 
  bind_rows(bm_car_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Cartilaginous Fish"))%>% 
  bind_rows(bm_fish_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Bony Fish"))%>% 
  bind_rows(bm_fish_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Bony Fish"))%>% 
  bind_rows(bm_mam_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Mammals")) %>% 
  bind_rows(bm_mam_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Mammals")) %>% 
  bind_rows(bm_mam_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(bm_rep_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_rep_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_rep_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("_gP1")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Reptiles"))  

(figure4 <- mu %>% 
  ggplot(aes(x = value, y = taxon, colour=taxon)) +
  stat_halfeye(point_interval = mode_hdi, 
               .width = c(0.5, .8, 0.95),
               normalize = "panels") +
  facet_wrap(system~.)+
  geom_vline(xintercept = 0, 
             linetype = "dashed", colour="grey50") +
  labs(x="Size effects (body mass)", y = "Taxon") +
  scale_colour_jco(name="Class")+
  geom_text(data=sample_size,
            aes(y=taxon, x=3, label=paste("n =", n)),
            vjust=-.1)+
  coord_cartesian(clip = "off")+
    xlim(-3, 3)+
  theme(legend.position = "none"))

ggsave(figure4, file="Figure4.pdf",
       width = 10, height = 6,
       path = ResultPath)


# Figure 5: Latitude ###########################################################

# Sample size

sample_size <- lpi %>% 
  filter(!is.na(Latitude)) %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(System, Class) %>% 
  summarise(n=n()) %>% 
  filter(n>1) %>% 
  rename(taxon=Class,
         system=System)

# Do summaries for each model

mu <- as_draws_df(lat_amp_fre) %>% 
  pivot_longer(ends_with("Latitude")) %>% 
  mutate(name=factor(name),
         system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(lat_amp_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(lat_bird_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(lat_bird_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Birds"))%>% 
  bind_rows(lat_bird_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Birds"))%>% 
  bind_rows(lat_car_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Cartilaginous Fish"))%>% 
  bind_rows(lat_fish_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Bony Fish"))%>% 
  bind_rows(lat_fish_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Bony Fish"))%>% 
  bind_rows(lat_mam_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Mammals")) %>% 
  bind_rows(lat_mam_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Mammals")) %>% 
  bind_rows(lat_mam_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(lat_rep_fre %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Freshwater",
                     taxon="Reptiles")) %>% 
  bind_rows(lat_rep_mar %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Marine",
                     taxon="Reptiles")) %>% 
  bind_rows(lat_rep_ter %>% 
              as_draws_df() %>%
              pivot_longer(ends_with("Latitude")) %>% 
              mutate(name=factor(name),
                     system="Terrestrial",
                     taxon="Reptiles"))  

(figure5 <- mu %>% 
    ggplot(aes(x = value, y = taxon, colour=taxon)) +
    stat_halfeye(point_interval = mode_hdi, 
                 .width = c(0.5, .8, 0.95),
                 normalize = "panels") +
    facet_wrap(system~.)+
    geom_vline(xintercept = 0, 
               linetype = "dashed", colour="grey50") +
    labs(x="Size effects (absolute latitude)", y = "Taxon") +
    scale_colour_jco(name="Class")+
    geom_text(data=sample_size,
              aes(y=taxon, x=3, label=paste("n =", n)),
              vjust=-.1)+
    coord_cartesian(clip = "off")+
    xlim(-3, 3)+
    theme(legend.position = "none"))

ggsave(figure5, file="Figure5.pdf",
       width = 10, height = 6,
       path = ResultPath)
