# Body mass and latitude as global predictors of vertebrate populations exposure to multiple threats

Pol Capdevila<sup>1</sup>*, Nicola Noviello<sup>1</sup>, Louise McRae<sup>2</sup>, Robin Freeman<sup>2</sup>, Christopher F Clements<sup>1</sup>

<sup>1</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK. 

<sup>2</sup>Institute of Zoology, Zoological Society of London, Regent’s Park, London NW1 4RY, UK.

#### Contact: pcapdevila.pc[at]gmail.com

---

## Abstract

_The interactive effects of multiple threats are one of the main causes of biodiversity loss, yet our understanding of what predisposes species to be impacted by multiple threats remains limited. Here we analyse a global dataset of over 7000 marine, freshwater, and terrestrial vertebrate populations, alongside trait, threat, and geographical data, to identify the factors influencing the number of threats a species is subjected to at the population level. Out of a suite of predictors tested, we find that body mass and latitude both are broadly available for vertebrate species and influence the number of threats a population is subjected to. Larger bodied species and those nearer the equator are typically affected by a higher number of threats. However, whilst this pattern broadly holds across ecosystems for most taxa, amphibians and reptiles show opposing trends. We suggest that latitude and body mass should be considered as key predictors to identify which vertebrate populations are likely to be impacted by multiple threats. These general predictors can help to better understand the impacts of the Anthropocene on global vertebrate biodiversity and design effective conservation policies._

---

## Data

- __`LPIData.RData`__: data set used for the analyses. Contains information about the ID of the population time series in the living planet database. The species scientific name. The number of threats at which the population is exposed. The body mass of the species. The body size of the species. The habitat breadth. The trophic level. The taxonomic realm, class, order, family, genus, species and subspecies. IUCN red list category. Latitude and longitude of the population. System.   
- __`Body_Mass_means.csv`__: contains the mean values of the adult body weight (g) and size (cm) of the studied species.
- __`Amphibians.csv`__: contains trait information for amphibians from [Etard et al. 2020](https://doi.org/10.1111/geb.13184).
- __`Birds.csv`__: contains trait information for birds from [Etard et al. 2020](https://doi.org/10.1111/geb.13184).
- __`Mammals.csv`__: contains trait information for mammals from [Etard et al. 2020](https://doi.org/10.1111/geb.13184).
- __`Reptiles.csv`__: contains trait information for reptiles from [Etard et al. 2020](https://doi.org/10.1111/geb.13184).
- __`amphibians.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).  
- __`birds.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).  
- __`mammals.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).
- __`reptiles.nex`__: phylogeny of amphibians obtained from [vertlife.org](https://vertlife.org/).
- __`fishtree.tre`__: phylogeny of fishes obtained from [Rabosky et al. 2014](https://www.nature.com/articles/s41586-018-0273-1).
- __`PopDensity.csv`__: human density data obtained from [Hurtt et al. 2011](https://link.springer.com/article/10.1007/s10584-011-0153-2).

---

# Code

To run the statistical analyses we used different R scripts: 

- __`Models.R`__: code to analyse the factors influencing the number of threats at which a population is exposed to.
- __`Figures.R`__: code to create the figures and tables of the study. 
- __`SupFigs.R`__: code to represent the supplementary results and figures. 
- __`PhyloSignal.R`__: code to analyse the phylogenetic signal of the number of threats. 

---

# Software

_R version 4.0.2 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .
