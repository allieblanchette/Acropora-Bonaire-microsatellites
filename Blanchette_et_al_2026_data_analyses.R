########################################################################################################################################################################################################################################################################################################################################################################
# MICROSATELLITE ALLELE DATA ANALYSES
# Purpose: analyzing (pre-cleaned/pre-wrangled) microsatellite data of Acropora cervicornis and Acropora palmata from Reef Renewal Bonaire's gene bank
# Data: Cleaned Acropora palmata and Acropora cervicornis microsatellite allele tables from Reef Renewal Bonaire's gene bank
# Last updated on: 2026 February 27
# Last updated by: Allie Blanchette
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 
########################################################################################################################################################################################################################################################################################################################################################################

##### Publication citation:
# Title: Genetic diversity and connectivity of endangered Acropora spp. corals in Bonaire, Caribbean Netherlands: implications for restoration and gene banking
# Authors: Allison Blanchette, Francesca Virdis, Sanne Tuijten, Pearl Rivers Key, Don Levitan, Sarah E. Lester, Andrew Rassweiler





##### Clear working environment
rm(list=ls())

##### Load in relevant packages 

# For data wrangling
library(dplyr)
library(tidyverse)

# For analyses
library(adegenet)     # find.clusters()
library(poppr)        # for bulk of population genetics analyses
# library(pegas)      # for determining if the populations are in Hardy-Weinberg equilibrium using the function hw.test(). Masks mst and amova from poppr so only use as needed
library(genepop)      # for test_LD() and nulls()
library(vegan)        # for Mantel test 
library(related)      # for coancestry()
library(iNEXT)        # for rarefaction at species level
#library(pegas)       # masks multiple other functions from other relevant packages. Only use as needed. Used for rarefactions per locus
library(PopGenReport) # for null.all 
library(graph4lg)     # for genind_to_genepop
# library(gdistance)  # for calculating distances between genotype harvesting locations. Masks many functions, only use as needed
library(hierfstat)    # for locus descriptive summary stats (Na, He, Ho)
library(sf)           # for geographic distance calcs in Mantel test
library(reshape2)     # for melt()

# For spatial anlayses and mapping
library(sf)
library(terra)
library(tmap)

# For plotting 
library(scales)     # for number_format()
library(cowplot)    # for ggdraw()
library(patchwork)  # for plot_layout(), plot_annotation()
library(magick)     # for draw_image()
library(ggfortify)  # for plotting PCA
library(ggrepel)    # for plotting PCA
#library(ggspatial) # for turning spatial object into ggplot2 layer
library(viridis)    # for colorblind-friendly palettes




####################################################################################################################################################################################
## Loading in data - df2genind 
####################################################################################################################################################################################

##### set working directory
#setwd("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis")
setwd("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Final code and data\\Github\\Data analyses")
#setwd("D:\\Science\\Allie")


##### Load in data
# make row names the lineage ID (first column in csv)
Apal_wG32_df <- read.csv("Apal32_251209.csv", row.names = 1)         
Acer_triploid_df <- read.csv("Acer_Triploid_formatted_251209.csv", row.names = 1)
Acer_df <- read.csv("Acer_251209.csv", row.names = 1)


##### Create genind objects! 
Apal_g0_f0 <- 
df2genind(
  X = Apal_wG32_df[, -which(names(Apal_wG32_df) == "Pop")],    # locus allele data is everything in allele_table except the "Pop" column
  sep = "/",                                                   # allele values per locus separated by "/" symbol
  ncode = 3,                                                   # number of digits per allele value
  pop = Apal_wG32_df$Pop,                                      # Population data is in Pop column  
  NA.char = "NA",                                              # missing data is represented by NA
  ploidy = 2,                                                  # For diploid-formatted data 
  type = "codom",                                              # Microsatellites are codominant (as opposed to 'PA' presence/absence markers like AFLP)
)

Acer_T_g0_f0 <- 
  df2genind(
    X = Acer_triploid_df[, -which(names(Acer_triploid_df ) == "Pop")],   # locus allele data is everything in allele_table except the "Pop" column
    sep = "/",                                                           # allele values per locus separated by "/" symbol
    ncode = 3,                                                           # number of digits per allele value
    pop = Acer_triploid_df$Pop,                                          # Population data is in Pop column  
    NA.char = "NA",                                                      # missing data is represented by NA
    ploidy = 3,                                                          # For diploid-formatted data 
    type = "codom",                                                      # Microsatellites are codominant (as opposed to 'PA' presence/absence markers like AFLP)
  )

Acer_g0_f0 <- 
  df2genind(
    X = Acer_df[, -which(names(Acer_df ) == "Pop")],        # locus allele data is everything in allele_table except the "Pop" column
    sep = "/",                                              # allele values per locus separated by "/" symbol
    ncode = 3,                                              # number of digits per allele value
    pop = Acer_df$Pop,                                      # Population data is in Pop column  
    NA.char = "NA",                                         # missing data is represented by NA
    ploidy = 2,                                             # For diploid-formatted data 
    type = "codom",                                         # Microsatellites are codominant (as opposed to 'PA' presence/absence markers like AFLP)
  )


##### Review
(geneious_genind <- Apal_g0_f0)
(geneious_genind <- Acer_T_g0_f0)
(geneious_genind <- Acer_g0_f0)

summary(geneious_genind)
pop(geneious_genind)
locNames(geneious_genind)



####################################################################################################################################################################################
## Further genind data set-up 
####################################################################################################################################################################################


##### load in subpopulation geographic data and assign to genind

# Acropora palmata
Subpop_centroids_AP <- read.csv("Subpops_centroids.csv") %>% 
    filter(Species == "AP") %>% 
    dplyr::select(-Species) 

rownames(Subpop_centroids_AP) <- Subpop_centroids_AP$Pop # important to move the Pop ID over to row names so that these character values don't mess with dist() (it warns about NAs otherwise)

Subpop_centroids_AP <- Subpop_centroids_AP[2:3]

Apal_g0_f0@other$xy <- Subpop_centroids_AP


# Acropora cervicornis
Subpop_centroids_AC <- read.csv("Subpops_centroids.csv") %>% 
    filter(Species == "AC") %>% 
    dplyr::select(-Species) 

rownames(Subpop_centroids_AC) <- Subpop_centroids_AC$Pop # important to move the Pop ID over to row names so that these character values don't mess with dist() (it warns about NAs otherwise)

Subpop_centroids_AC <- Subpop_centroids_AC[2:3]

Acer_g0_f0@other$xy <- Subpop_centroids_AC




##### lineage filtering
# Apal_g1 = AP_G32 excluded                  -- G32 = identical MLG identity with G02
# Acer_g1 = AC_G01 excluded, AC_G21 retained -- G01 = triploid, 
# Acer_g2 = both AC_G01 and AC_G21 excluded  -- G21 = collection location unknown, inconclusive w/G25

# Acropora palmata
Apal_g1_f0 <- Apal_g0_f0[!(indNames(Apal_g0_f0) == "AP_G32")]                  
  
# Acropora cervicornis
Acer_g1_f0 <- Acer_g0_f0[!(indNames(Acer_g0_f0) == "AC_G01")]                  
Acer_g2_f0 <- Acer_g0_f0[!(indNames(Acer_g0_f0) %in% c("AC_G01", "AC_G21"))]   



##### f1 locus filtering
# f1 = loci with >15% missing data removed. This 15% value is important in order to remove worst loci, while still retaining enough to distinguish MLGs

# Acropora palmata
Apal_g0_f1 <- missingno(Apal_g0_f0, type = "loci", cutoff = 0.15)     # No loci with missing values above 15% found
Apal_g1_f1 <- missingno(Apal_g1_f0, type = "loci", cutoff = 0.15)     # No loci with missing values above 15% found    

# Acropora cervicornis
Acer_g0_f1 <- missingno(Acer_g0_f0, type = "loci", cutoff = 0.15)     # 5 loci removed: 1195, 207, 6212, 0513, 2637 
Acer_g1_f1 <- missingno(Acer_g1_f0, type = "loci", cutoff = 0.15)     # 5 loci removed: 1195, 207, 6212, 0513, 2637 
Acer_g2_f1 <- missingno(Acer_g2_f0, type = "loci", cutoff = 0.15)     # 5 loci removed: 1195, 207, 6212, 0513, 2637 

# Check if any loci are uninformative after first pass of filtering
# informloci() function removes uninformative loci. The default is 2/N (N = individuals)
informloci(Apal_g1_f1)     
informloci(Acer_g1_f1)
informloci(Acer_g2_f1)

nLoc(Apal_g1_f1)           # 11 = no loci were removed 
nLoc(Acer_g1_f1)           # 6 = no loci were removed
nLoc(Acer_g2_f1)           # 6 = no loci were removed



##### Create version of f1 data in which all lineages are in 1 Bonaire population (ie no sub-pops) - used for null allele analyses and calculating summary stats 
# Using Apal_g1 because don't want duplicate MLG G32 influencing allele frequency calculations
# Using Acer_g1 because don't want triploid MLG G01 influencing allele frequency calculations
Apal_g1_f1_singlePop <- Apal_g1_f1
pop(Apal_g1_f1_singlePop) <- factor(rep("Bonaire", nInd(Apal_g1_f1_singlePop)))

Acer_g1_f1_singlePop <- Acer_g1_f1
pop(Acer_g1_f1_singlePop) <- factor(rep("Bonaire", nInd(Acer_g1_f1_singlePop)))



##### f2 locus filter
# f2 = loci significantly deviating from Hardy-Weinberg Equilibrium removed   

# Acropora palmata
remove_loci_AP <- c("X166_PET", "X207_VIC") 
keep_loci_AP <- setdiff(locNames(Apal_g0_f1), remove_loci_AP)
Apal_g0_f2 <- Apal_g0_f1[loc = keep_loci_AP]
Apal_g1_f2 <- Apal_g1_f1[loc = keep_loci_AP]

locNames(Apal_g0_f2)   # should be 9 remaining loci
locNames(Apal_g1_f2)   # should be 9 remaining loci



# Acropora cervicornis -- all f1 loci for A. cervicornis are in HWE, so none actually need to be removed here. I'm just creating the f2 version for consistency in naming
Acer_g0_f2 <- Acer_g0_f1
Acer_g1_f2 <- Acer_g1_f1
Acer_g2_f2 <- Acer_g2_f1

length(indNames(Acer_g0_f2))    # should be 25 lineages (1-27, -G02 out of commission, -G25 inconclusive)
length(indNames(Acer_g1_f2))    # should be 24 lineages (1-27, -G01 triploid, -G02 out of commission, -G25 inconclusive)
length(indNames(Acer_g2_f2))    # should be 23 lineages (1-27, -G01 triploid, -G02 out of commission, -G21 inconclusive, -G25 inconclusive)



##### f3 locus filter
# f3 = narrowing down Apal loci to the same 6 in Acer
# Note that while 207 (not in HWE) gets removed, 181 and 166 remain
AC_6loci <- locNames(Acer_g1_f1)
Apal_g1_f3 <- Apal_g1_f0[loc = AC_6loci]
locNames(Apal_g1_f3)   # should be 6

Apal_g1_f3_singlePop <- Apal_g1_f3
pop(Apal_g1_f3_singlePop) <- factor(rep("Bonaire", nInd(Apal_g1_f3_singlePop)))




##### Testing if removing certain loci makes a difference (i.e., sensitivity tests)

# Set genind object
set_coral_sens <- Apal_g1_f1

# Set locus/loci to remove
#remove_locus_test <- "X207_VIC" 
#remove_locus_test <- "X181_NED" 
remove_locus_test <- "X166_PET" 
#remove_locus_test <- c("X166_PET", "X207_VIC") 

# Set loci to keep
keep_loci_sens <- setdiff(locNames(set_coral_sens), remove_locus_test)

# Create genind with certain locus/loci removed
sens_test <- set_coral_sens[loc = keep_loci_sens]




####################################################################################################################################################################################
## Convert to genepop file 
####################################################################################################################################################################################

# genepop objects used in test_LD() and nulls() functions

##### General genepop
# genind_to_genepop converts an object of class genind into a text file compatible with GENEPOP software, as opposed to adegenet which converts an object of class genind into an object of class genpop
# Intended to be used with diploid msat data with alleles coded with 2 or 3 digits
# "There is only one population in your dataset" warning is fine (yes, it's all just pop = "Bonaire")
genind_to_genepop(Apal_g1_f1_singlePop, output = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\Other programs\\genepop\\Apal_g1_f1_singlePop_genepop.txt")
genind_to_genepop(Acer_g1_f1_singlePop, output = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\Other programs\\genepop\\Acer_g1_f1_singlePop_genepop.txt")


##### Targeted for use in micro-checker

# Create separate genind object
Apal_g1_f1_singlePop_microchecker <- Apal_g1_f1_singlePop
Acer_g1_f1_singlePop_microchecker <- Acer_g1_f1_singlePop

# Replace individual names with a constant
indNames(Apal_g1_f1_singlePop_microchecker) <- rep("Ind", nInd(Apal_g1_f1_singlePop_microchecker))
indNames(Acer_g1_f1_singlePop_microchecker) <- rep("Ind", nInd(Acer_g1_f1_singlePop_microchecker))

# Write Genepop file
# "There is only one population in your dataset" warning is fine (yes, it's all just pop = "Bonaire")
genind_to_genepop(Apal_g1_f1_singlePop_microchecker, output = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\null alleles\\micro-checker\\Apal_g1_f1_singlePop_microchecker.txt")
genind_to_genepop(Acer_g1_f1_singlePop_microchecker, output = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\null alleles\\micro-checker\\Acer_g1_f1_singlePop_microchecker.txt")




####################################################################################################################################################################################
## Set up plotting aesthetics
####################################################################################################################################################################################

# Sub-pop colors: 1. Central, 2. Klein, 3. Lac Bay, 4. North, 5. South
scales::show_col(viridis::viridis_pal(option = "D")(5))
pop.colors <- c("#440154FF", "#3B528BFF", "#21908CFF", "#5DC863FF", "#FDE725FF")

# Sub-pop order
pop.order <- c("North", "Central", "South", "Klein", "Lac Bay")

# Plot text sizes
axis.title.sz = 12
axis.txt.sz = 10
leg.txt.sz = 8
annotate.sz = 5
title.sz = 12
panel.txt.sz = 13

# Plot feature sizes
point.sz = 3
regln.sz = 0.75




####################################################################################################################################################################################
## Linkage Disequilibrium
####################################################################################################################################################################################


##### Linkage Disequilibrium:
# LD Definition: certain alleles are statistically linked with other alleles on other loci (ie loci are located near each other on the same chromosome, negating the assumption of independent assortment)
# Processes that produce LD: asexual reproduction, non-random mating, selection, population differentiation
# Ia = Index of Association. It ests to what extent individuals are the same or different at one locus are also the some or different at other loci
# rBarD = Standardized Index of Association (Agapapow & Burt 2001) to account for number of loci


########## poppr::ia ##########

##### Set coral dataset
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 
#set_coral <- Apal_g1_f1
set_coral <- Acer_g1_f1


##### Set myTheme overide in pair.ia() to be able to plot heatmap (myTheme no longer compatible with ggplot >=3.4)

# Get the actual function
f <- getFromNamespace("pair.ia", "poppr")

# Get its internal environment
fun_env <- environment(f)

# Create a valid theme and bind it inside that environment
if (!exists("myTheme", envir = fun_env, inherits = FALSE)) {
  assign("myTheme", ggplot2::theme_minimal(base_size = 12), envir = fun_env)
} else {
  # Unlock and replace if already there
  unlockBinding("myTheme", fun_env)
  assign("myTheme", ggplot2::theme_minimal(base_size = 12), envir = fun_env)
  lockBinding("myTheme", fun_env)
}


##### Across all loci. If p > 0.05 then accept null of Linkage Equilibrium
# p > 0.05 --> LE 
# p < 0.05 --> LD  (reject null, supports that some loci are correlated)  
set.seed(123)
ia(set_coral, sample = 999, quiet=TRUE)


##### Pairwise comparisons among all loci
# values in heatmap = p-value
# colors in heatmap = r-bar
# the redder the color and lower the p-value, the more tightly correlated the 2 loci are
set.seed(123)
pair.ia(set_coral, sample = 999)



########## genepop::test_LD ##########

##### set working directory (necessary for nulls to pull the txt file properly)
# genepop files created in section 'Convert to genepop file' above
setwd("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\Other programs\\genepop")


##### Run
test_LD(
  "genepop_Apal.txt",
  outputFile = "testLD_Apal.csv",
  settingsFile = "",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)

test_LD(
  "genepop_Acer.txt",
  outputFile = "testLD_Acer.csv",
  settingsFile = "",
  dememorization = 10000,
  batches = 100,
  iterations = 5000,
  verbose = interactive()
)







####################################################################################################################################################################################
## Check for null alleles
####################################################################################################################################################################################

# Null alleles: alleles that do not appear in a PCR product, for reasons such as mutations in flanking regions of the msat that prevent primer annealing, preferential amplification of short alleles, or slippage during aplification (Chapuis and Estoup 2006, Dabrowski et al 2015)

########## PopGenReport::null.all ##########

# The function null.all determines the frequency of null alleles at each locus of a genind object
# Step 1. makes a bootstrap estimate (based on the observed allele frequencies) of the probability of seeing the number of homozygotes observed for each allele. If there are a large number of null alleles present at a locus, it would result in multiple alleles at a locus having an excess of homozygotes
# homozygotes$observed = observed number of homozygotes for allele at each locus 
# homozygotes$bootstrap = distribution of the expected number of homozygotes for each allele at each locus based on observed allele frequencies 
# homozygotes$probability.obs = summary table given the probability of observing the number of homozygotes (homozygotes$probability.obs)
# Step 2. estimates the frequency of null alleles and a bootstrap confidence interval for each locus using the methods of Chakraborty et al. (1994) and Brookfield (1996). 
# null.allele.freq = null allele frequency estimates determined from the observed heterozygosity and homozygosity at a locus. 
# summary1 (Chakraborty et al. 1994) - use if there is some missing allele data
# summary2 (Brookfield (1996)        - use if all individuals have data
# median, 2.5th, and 97.5th percentiles are from bootstrap estimates of null allele frequencies obtained by resampling the individual genotypes from the original genind object
# ***If the 95% confidence interval includes zero(eg lower CI <0), then frequency of null alleles at a locus does not significantly differ from zero***

# Run
null_Apal <- null.all(Apal_g1_f1)
# null_AcerT <- null.all(AcerT_f1) # Error in null.all(AcerT_f1) : One or more populations has a ploidy other than 2! Script stopped!
null_Acer <- null.all(Acer_g1_f1)

# Call summary result info (this can take a minute or so)
# Check 2.5th percentile. If all values <0, then null allele frequency is not significantly different from 0
# Results stored in QAQC analyses and data, tab Final locus screening
null_Apal$null.allele.freq$summary1 
null_Acer$null.allele.freq$summary1





########## genepop::nulls ##########

##### set working directory (necessary for nulls to pull the txt file properly)
# genepop files created in section 'Convert to genepop file' above
setwd("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\Other programs\\genepop")


# Convert to genepop file
# genind_to_genepop converts an object of class genind into a text file compatible with GENEPOP software, as opposed to adegenet which converts an object of class genind into an object of class genpop
# Intended to be used with diploid msat data with alleles coded with 2 or 3 digits
# "There is only one population in your dataset" warning is fine (yes, it's all just pop = "Bonaire")
genind_to_genepop(Apal_g1_f1_singlePop, output = "genepop_Apal.txt")
genind_to_genepop(Acer_g1_f1_singlePop, output = "genepop_Acer.txt")


# Run
# NOTE - nulls cannot overwrite previous files, so would need to delete or rename previous version when re-running this function
nulls("genepop_Apal.txt",
      
      # Default method - Maximum likelihood using EM algorithm of Dempster et al 1977 which assumes that there IS a null allele.
      # output = "nulls_Apal_default.csv",  
      
      # # Apparent nulls method - maximum likelihood assuming that apparent nulls are technical failures independent of genotype 
      # output = "nulls_Apal_ApparentNulls.csv",
      # nullAlleleMethod = "ApparentNulls",
      
      # Brookfield (1996) estimator - recommended by PopGenReport to use when there is NOT missing data
      output = "nulls_Apal_B96.csv",
      nullAlleleMethod = "B96",
      
      CIcoverage = 0.95)

nulls("genepop_Acer.txt",
      
      # Default method - Maximum likelihood using EM algorithm of Dempster et al 1977 which assumes that there IS a null allele.
      # output = "nulls_Acer_default.csv",  
      
      # # Apparent nulls method - maximum likelihood assuming that apparent nulls are technical failures independent of genotype
      # output = "nulls_Acer_ApparentNulls.csv",
      # nullAlleleMethod = "ApparentNulls",
      
      # Brookfield (1996) estimator - recommended by PopGenReport to use when there is NOT missing data
      output = "nulls_Acer_B96.csv",
      nullAlleleMethod = "B96",
      
      CIcoverage = 0.95)

# Output interpretation
# https://genepop.curtin.edu.au/Option8.html
# Homoz. = number of observed homzoygotes that are actually homozygote (true homozygote)
# Null Heter. = number observed homozygotes that are predicted to actually be heterozygotes with 1 null allele (false homozygote)
# 'No info for CI' in confidence interval summary table means sample size was likely too small to generate CIs


####################################################################################################################################################################################
## Genotype confidence & MLLs
####################################################################################################################################################################################
########## Genotype accumulation curves ##########

##### Set coral settings
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 

# Apal
set_coral <- Apal_g0_f0
title <- expression(italic("Acropora palmata"))
set_ylims <- c(0, 34)
set_ybreaks <- c(0, 5, 10, 15, 20, 25, 30, 34)

# Acer
set_coral <- Acer_T_g0_f0
title <- expression(italic("Acropora cervicornis"))
set_ylims <- c(0, nInd(set_coral))
set_ybreaks <- seq(0, nInd(set_coral), b = 5)


##### Run genotype_curve()
# columns = number of loci 1-10
# row = resampling iteration (so for each number of loci, you have 1000 replicate curves)
# values = number of unique MLGs detected
# this is 'out' from within the function (see gentAnywhere(genotype_curve))
# the labs() function within genotype_curve() is no longer compatible with more recent ggplot2 version. So set plot = FALSE and re-create the plot myself
gc <- poppr::genotype_curve(set_coral,
               sample = 1000,   
               quiet = FALSE,
               thresh = 1,
               plot = FALSE,
               drop = TRUE,
               dropna = TRUE
)


##### Wrangle gc to plot
# Code pulled from within genotype_curve() function, found with getAnywhere()
outmelt <- as.data.frame.table(gc, responseName = "MLG", stringsAsFactors = FALSE) # melts gc into long-format. Note that sample and NumLoci are the dimension names of gc - see dimnames(gc)
outmelt$sample <- as.integer(outmelt$sample)    # make sample and NumLoci integers
outmelt$NumLoci <- as.integer(outmelt$NumLoci)  # make sample and NumLoci integers


##### Create individual plots per species

#genocurve_apal <-
genocurve_acer <-
ggplot(outmelt, aes(x = NumLoci, y = MLG)) +
  geom_boxplot(aes(group = factor(NumLoci))) + 
  ggtitle(title)+
  scale_x_continuous(breaks = seq(min(outmelt$NumLoci), max(outmelt$NumLoci)),
                     expand = c(0, 0.125))+
  scale_y_continuous(limits = set_ylims, breaks = set_ybreaks)+
  xlab("")+
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, size = title.sz),
    axis.title.x = element_text(size = axis.title.sz),
    axis.text.x = element_text(size = axis.txt.sz),
    axis.title.y = element_blank(), 
    axis.text.y = element_text(size = axis.txt.sz),
  )


##### Final plotting and save

# Combine and label plots as panels
genocurve_both <- 
  genocurve_apal + genocurve_acer +
  plot_layout(ncol = 2,      # side-by-side
              widths = c(1, 1)) +  # equal widths
  plot_annotation(
    tag_levels = "A",         # automatically label plots A and B
    tag_prefix = "",          # optional: removes "Figure" prefix
    tag_sep = ""              # just "A" and "B" as labels
  )

  
genocurve_final <-
  ggdraw()+
  draw_plot(genocurve_both, x = 0.02, y = 0.02, width = 0.96, height = 0.96)+                        # Add the plots, shrink them slightly in y/x
  draw_label("Number of MLGs", x = 0.02, y = 0.5, angle = 90, size = axis.title.sz, vjust = 0.5)+  # Y-axis label
  draw_label("Number of loci", x = 0.5, y = 0.02, vjust = 0, size = axis.title.sz)                 # X-axis label
  
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\geno_accumulation_curves.tiff", genocurve_final, dpi=600, units="cm", width = 18, height = 8)




########## Psex: probability of encountering a genotype more than once by chance ##########


##### Some notes:
# Psex is the probability of encountering a given genotype more than once by chance.
# The basic equation from Parks and Werth (1993) is
# $$p_{sex} = 1 - (1 - p_{gen})^{G})$$
# ### G is the number of multilocus genotypes
# ### pgen is the probability of a given genotype (see pgen for its calculation)
# For a given value of alpha (e.g. alpha = 0.05), genotypes with psex < alpha can be thought of as a single genet whereas genotypes with psex > alpha do not have strong evidence that members belong to the same genet (Parks and Werth, 1993).



##### Set up genclone objects fresh (to ensure not the collapsed MLL version from above)
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 

# Apal full dataset (including G32 clone, n = 9 loci in HWE)
(set_gc_obj <- as.genclone(Apal_g0_f2))
set_title <- "Apal"

# Acer full dataset (no filtering of loci, n = 6 loci in HWE, diploid AC_G01 retained)
(set_gc_obj <- as.genclone(Acer_g0_f2))
set_title <- "Acer"


##### Run psex

# run psex
Acrop_psex <- psex(set_gc_obj,
                 pop = NULL,         # Ignore sub-population strata because Bonaire in panmictic
                 by_pop = FALSE,     # Ignore sub-population strata because Bonaire in panmictic
                 freq = NULL,        # NULL indicates that the allele frequencies will be determined via round-robin approach in rraf. If this matrix or vector is not provided, zero-value allele frequencies will automatically be corrected
                 G = NULL,           # an integer vector specifying the number of observed genets. If NULL, this will be the number of original multilocus genotypes for method = "single"
                 method = "single"   # Using method = "single" (default) indicates that the calculation for psex should reflect the probability of encountering a second genotype. Using method = "multiple" gives the probability of encountering multiple samples of the same genotype 
)

# make psex output dataframe for inspection
Acrop_psex_df <- as.data.frame(Acrop_psex)

# Visual aid for psex: Index = lineages, psex = probability of obtaining that MLG by chance
plot(Acrop_psex,
     col = ifelse(Acrop_psex > 0.05, "red", "blue"),   # plot any probabilities greater than 0.05 (5%) in red
     main = title)






########## Multi-Locus Lineages (MLLs) ##########
# MLG = Multi-Locus Genotype
# MLL = Multi-Locus Lineage (collapsed MLGs that are very similar into same MLL)
# Conducted using tutorial: https://grunwaldlab.github.io/poppr/articles/mlg.html#filtered-contracted 
# Apal note - both Apal datasets stay at 34 MLGs and 34 MLLs
# Acer note - both Acer datasets collapse 26 MLGs into 24 MLLs



##### Acropora palmata

# Set up original data
(Apal.gc <- as.genclone(Apal_g0_f0))                # Convert from genind object to genclone and check # MLGs, # individuals, # loci
mlg.table(Apal.gc)                                  # Plot abundances of individuals per geno per population (make sure to manually check if MLGs are repeated BETWEEN pops): MLG.26 (NOT AC_G26) duplicated in Klein and South
xt <- apply(tab(Apal.gc), 1, paste, collapse = "")  # Determine MLG identity of each individual
rank(xt, ties.method = "first")                     # Identify ranked MLG assignment from poppr for RRB's individuals

# Choose threshold and visualize performance of algorithms
# threshold = minimum genetic distance at which two individuals would be considered from different clonal lineages
# ideally choose a threshold at which all algorithms perform the same (ie same #MLLs output)
Apal_alg <- filter_stats(Apal.gc, distance = bruvo.dist, replen = c(rep(3, 11)), plot = TRUE)    # Threshold should be value inbetween first small peak and taller main 2nd peak (eg ~0.1)
print(farthest_thresh <- cutoff_predictor(Apal_alg$farthest$THRESHOLDS))                         # Threshold for "farthest" algorithm = 0.20803
print(average_thresh <- cutoff_predictor(Apal_alg$average$THRESHOLDS))                           # Threshold for "average" algorithm =  0.20803
print(nearest_thresh <- cutoff_predictor(Apal_alg$nearest$THRESHOLDS))                           # Threshold for "nearest" algorithm =  0.201638

# Collapse MLGs into MLLs 
# use t = 0.1 because that is where the dip in the filter_stats() plot is at and while the cutoff_predictor() recommends 0.2, 0.2 inaccurately collapses AP_G06 and AP_G29, but t = 0.1 does not
mlg.filter(Apal.gc,                                      # Use genclone object
           missing = "asis",                             # "mean" = impute allele frequencies on missing data, "asis" = ignore missing data 
           algorithm = "farthest_neighbor",              # set algorithm to determine type of clustering (method of collapsing MLGs); "farthest neighbor" merges clusters based on the maximum distance between points in either cluster, it is the most conservative and the default for mlg.filter              
           distance = bruvo.dist,                        # set distance-type for comparing MLGs
           replen = c(rep(3, 11))                        # set repeat length of each locus
) <- 0.1                                                 # set MLGs at threshold (IMPORTANT- this must be done here, the argument "threshold = " does something different, and it still elludes me a bit)
Apal.gc                                                  # View updated collapsed MLLs. Note that under "contracted MLGs" the metadata of threshold, distance, and algorithm should show up

# View which MLGs (according to original RRB naming) got collapsed into MLLs
mlg.id(Apal.gc)                                          # Identify which genos got collapsed into the same MLL                                            
mlg.table(Apal.gc)                                       # Plot abundances of individuals per geno (see where the MLLs came from)








##### Acropora cervicornis

# Set up original data
(Acer.gc <- as.genclone(Acer_T_g0_f0))              # Convert from genind object to genclone and check # MLGs, # individuals, # loci
mlg.table(Acer.gc)                                # Plot abundances of individuals per geno per population (make sure to manually check if MLGs are repeated BETWEEN pops): MLG.26 (NOT AC_G26) duplicated in Klein and South
xt <- apply(tab(Acer.gc), 1, paste, collapse = "")  # Determine MLG identity of each individual
rank(xt, ties.method = "first")                     # Identify ranked MLG assignment from poppr for RRB's individuals

# Choose threshold and visualize performance of algorithms
# threshold = minimum genetic distance at which two individuals would be considered from different clonal lineages
# ideally choose a threshold at which all algorithms perform the same (ie same #MLLs output)
# warnings are all just about NAs induced by coercion, probably having to do with all the 0s in triploid format, or high amounts of missing data from 5 loci
Acer_alg <- filter_stats(Acer.gc, distance = bruvo.dist, replen = c(rep(3, 11)), plot = TRUE)    # Threshold should be value inbetween first small peak and taller main 2nd peak (eg ~0.1)
print(farthest_thresh <- cutoff_predictor(Acer_alg$farthest$THRESHOLDS))                         # Threshold for "farthest" algorithm = 0.2211062
print(average_thresh <- cutoff_predictor(Acer_alg$average$THRESHOLDS))                           # Threshold for "average" algorithm =  0.2211062
print(nearest_thresh <- cutoff_predictor(Acer_alg$nearest$THRESHOLDS))                           # Threshold for "nearest" algorithm =  0.2247721

# Collapse MLGs into MLLs 
# use t = 0.1 because while the cutoff_predictor() and plot indicates 0.22, 0.22 inaccurately collapses AC_G08 and AC_G23, but t = 0.1 does not
mlg.filter(Acer.gc,                                       # Use genclone object
           missing = "asis",                               # "mean" = impute allele frequencies on missing data, "asis" = ignore missing data 
           algorithm = "farthest_neighbor",                # set algorithm to determine type of clustering (method of collapsing MLGs); "farthest neighbor" merges clusters based on the maximum distance between points in either cluster, it is the most conservative and the default for mlg.filter              
           distance = bruvo.dist,                          # set distance-type for comparing MLGs
           replen = c(rep(3, 11))                          # set repeat length of each locus
           ) <- 0.1                                       # set MLGs at threshold (IMPORTANT- this must be done here, the argument "threshold = " does something different, and it still elludes me a bit)
Acer.gc                                                   # View updated collapsed MLLs. Note that under "contracted MLGs" the metadata of threshold, distance, and algorithm should show up

# View which MLGs (according to original RRB naming) got collapsed into MLLs
mlg.id(Acer.gc)                                          # Identify which genos got collapsed into the same MLL                                            
mlg.table(Acer.gc)                                       # Plot abundances of individuals per geno (see where the MLLs came from)






####################################################################################################################################################################################
## Descriptive statistics and literature comparisons
####################################################################################################################################################################################

##### Read in locus name key
locus_name_key <- read.csv("locus_name_key.csv") %>% 
  mutate(Locus_simple = Locus) %>% 
  mutate(Locus = as.character(Locus),
         Locus = recode(Locus, "585" = "0585"),
         Locus = recode(Locus, "513" = "0513"))


########## Descriptive stats from our study ##########


##### Set coral species
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 
set_coral <- Apal_g1_f1_singlePop
#set_coral <- Acer_g1_f1_singlePop


##### Number of alleles 
loctab <- locus_table(set_coral, lev = "allele")
Num_alleles <- data.frame(Locus_tech = row.names(loctab), loctab) %>% 
  select(Locus_tech, allele) %>% 
  rename(Na = allele)


##### Heterozygosities
basic <- basic.stats(set_coral, diploid = TRUE) # pull basic stats (including heterozygosities). Notice that Hs = Ht (because Hs = He averaged across sub-pops, Ht = total He across all individuals, and we only using singlepop data here)

# pull Ho and He per locus
Hobsv_loci <- basic$perloc[1] 
Hexp_loci <- basic$perloc[2] 

# pull overall Ho and He
Hobsv_overall <- basic$overall[1, drop = FALSE]  
Hexp_overall <- basic$overall[2, drop = FALSE]

# combine per locus with Overall 
Hobsv <- rbind(Hobsv_loci, Hobsv_overall) 
Hobsv <- as.data.frame(rownames_to_column(Hobsv, var = "Locus_tech"))  %>%  # make row names (loci IDs) into a Locus column
  mutate(Locus_tech = recode(Locus_tech, "12" = "mean")) %>% # row 12 in Apal becomes mean
  mutate(Locus_tech = recode(Locus_tech, "7" = "mean")) # row 6 in Acer becomes mean

Hexp <- rbind(Hexp_loci, Hexp_overall)
Hexp <- as.data.frame(rownames_to_column(Hexp, var = "Locus_tech"))  %>%  # make row names (loci IDs) into a Locus column
  rename(He = Hs) %>% 
  mutate(Locus_tech = recode(Locus_tech, "12" = "mean")) %>% # row 12 in Apal becomes mean
  mutate(Locus_tech = recode(Locus_tech, "7" = "mean")) # row 6 in Acer becomes mean


##### Wrangle together table of observed and expected heterozygosities, and Hardy Weinberg Equilibrium deviations
# hw.test outputs:
# chi2 deviance value
# Pr(chi^2 >) = analytical p-value
# Pr. exact = p-value based on a permutation procedure (Pr. exact)
# *** Both sets of p values indicate if loci deviate significantly from the null expectation of HWE
hwe_matrix <- as.data.frame(pegas::hw.test(set_coral, B = 1000))  # Performs 1000 permutations. Note that pegas masks some functions in poppr so only calling it here as needed
hwe <- data.frame(Locus_tech = row.names(hwe_matrix), hwe_matrix) %>% 
  dplyr::rename(p = Pr.exact) %>% 
  dplyr::select(Locus_tech, p)


##### Combine heterozygosity and HWE calculations
summ_stats <- Num_alleles %>% 
  left_join(Hobsv) %>% 
  left_join(Hexp) %>% 
  left_join(hwe) %>% 
  left_join(locus_name_key) %>% 
  relocate(Locus, .before = Na)

###### Save csv!

# Apal
#write.csv(summ_stats, "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\summ_stats_AP.csv", row.names = F)

# Acer
#write.csv(summ_stats, "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\summ_stats_AC.csv", row.names = F)

# Apal f3
#write.csv(summ_stats, "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\summ_stats_AP_f3.csv", row.names = F)




##### Wrangle our descriptive statistics and prepare for joining with other studies

# Double check number of genotypes in dataset (to be incorporated into descriptives_RRB_AX below)
mlg(Apal_g1_f1_singlePop) # 34
mlg(Acer_g1_f1_singlePop) # 24

# Double check number of loci
nLoc(Apal_g1_f1_singlePop) # 11
nLoc(Acer_g1_f1_singlePop) # 6

# Acropora palmata
# NOTE - make sure this summ_stats df was created using A. palmata data above
descriptives_RRB_AP <- summ_stats %>%
    select(Locus, Na, He) %>%
    mutate(Author = "Present study",
           Year = "Present study",
           Title = "Present study",
           Species = "Apal",
           Population = "Bonaire",
           Method = "microsats",
           Genos = as.numeric(34),
           Num_loci = 11,
           Alleles_per_geno = Na/Genos) %>% # standardize number of alleles by number of unique multi-locus genotypes assessed
    left_join(locus_name_key)

# Acropora cervicornis
# NOTE - make sure this summ_stats df was created using A. cervicronis data above
descriptives_RRB_AC <- summ_stats %>%
  select(Locus, Na, He) %>% 
  mutate(Author = "Present study",
         Year = "Present study",
         Title = "Present study",
         Species = "Acer",
         Population = "Bonaire",
         Method = "microsats",
         Genos = as.numeric(24),
         Num_loci = 6,
         Alleles_per_geno = Na/Genos) %>% # standardize number of alleles by number of unique multi-locus genotypes assessed
  left_join(locus_name_key)

# Combine
descriptives_RRB <- descriptives_RRB_AP %>% 
  full_join(descriptives_RRB_AC)




########## Literature comparisons of descriptive stats ##########

##### Load in and wrangle microsatellite descriptive statistics from other studies
descriptives_others_raw <- read.csv("litcomp_others_perlocus_260203.csv")
  
descriptives_others <- descriptives_others_raw %>% 
  mutate(Year = as.character(Year),
         Genos = as.numeric(Genos),
         He = ifelse(!is.na(Hexp_report), Hexp_report, Hexp_calc),
         Locus_simple = as.character(Locus_simple),
         Alleles_per_geno = allele/Genos) %>%          # standardize number of alleles by number of unique multi-locus genotypes assessed)
  filter(Author != "Porto-Hannes",                     # removing because it's the only study missing locus 181 for Apal
         Author != "Dominguez-Maldonado et al.",       # removing because I don't 100% trust the accuracy of their Na methods/reuslts (used gels to estimate allele values?)
         Author != "Calle-Trivino et al.",             # removing because it's the only study missing locus 0585 for Acer
         Author != "Yetsko et al.") %>%                # removing because they used a subset of available genos from nursery and it's not clear how they chose genos (possibly biased towards variation?)   
  rename(Na = allele) %>% 
  left_join(locus_name_key)

loci_to_use <- descriptives_others %>% 
  group_by(Species) %>% 
  mutate(tot_studies = n_distinct(Title)) %>% 
  group_by(Species, Locus_simple, tot_studies) %>% 
  summarise(studies_used = n_distinct(Title)) %>% 
  mutate(prop_studies_used = studies_used/tot_studies) %>% 
  mutate(Use = ifelse(prop_studies_used == 1, "Yes", "No"),
         Use = ifelse(Locus_simple %in% c(192), "No", Use),  # I did not use locus 192 at all, but it was used in 100% of Apal studies. Set 192 to never use
         Use = ifelse((Species == "Acer" & Locus_simple %in% c(1195, 207, 6212, 0513, 2637)), "No", Use)  # These 5 loci were the ones that got filtered out in Acer f1. Set them to "no" use when in Acer
) %>%  
  left_join(locus_name_key)


  

##### General Hexp summary, using full locus datasets
Hexp_others_avgs <- descriptives_others %>% 
  group_by(Author, Year, Title, Species, Population) %>% 
  summarize(Hexp_mean = mean(He))


##### Join together descriptive stats and summarize

# Big join
litcomp_all <- descriptives_RRB %>%
    full_join(descriptives_others) %>%
    left_join(loci_to_use) %>%
    filter(Use == "Yes",
           !(Population %in% c("Big Cay, Honduras", "Cordelia Shoal, Honduras", "Smith Bank, Honduras")))  # just using "Honduras, all sites" rather than the 3 sub-sites (even though their reporting for number of MLGs total was a little unclear) to to keep the population scaling more consistent and make the table a little cleaner

# Load in Florida population key (necessary for next step)
FL_update_key <- read.csv("FL_update_key.csv") %>% 
  mutate(Year = as.character(Year))

# Summarize and clean up for publication-readiness
litcomp_summ_stats <- litcomp_all %>%
    group_by(Author, Year, Title, Species, Population) %>%
    summarise(mean_Na = mean(Na),
              mean_MLGs = round(mean(Genos), digits = 0), # Note: there were 2 cases in which 1 less MLG was analyzed for locus 0585: Honduras (52, instead of 53) and USVI (26, instead of 27) in Baums 2010. You can see this by adding Geno in the group_by line above and commenting out this line
              mean_He = mean(He, na.rm = T),
              loci_used = paste(unique(Locus), collapse = ", ")) %>% 
  left_join(FL_update_key) %>% 
  mutate(mean_Na_per_MLG = mean_Na/mean_MLGs,
         Reference = paste(Author, Year, sep = ", "),
         Reference = recode(Reference, "Present study, Present study" = "Present study"),
         Species = recode(Species, "Apal" = "A. palmata"),
         Species = recode(Species, "Acer" = "A. cervicornis"),
         Population = recode(Population, "Mona" = "Mona Island"),
         Population = recode(Population, "Honduras, all sites" = "Honduras"),
         Population = ifelse(!is.na(FL_pop_update), FL_pop_update, Population)) %>% 
  select(Species, Population, mean_Na_per_MLG, mean_Na, mean_MLGs, mean_He, Reference, loci_used, Title, Author, Year) %>%  # re-ordering in desired format for publication
  arrange(desc(Species), desc(mean_Na_per_MLG))

write.csv(litcomp_summ_stats, "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\litcomp_summ_stats.csv", row.names = F)


##### Calculate percent of populations that RRB has greater allelic richness than

# RRB allelic richnes: Acropora palmata
(ar_rrb <- litcomp_summ_stats$mean_Na_per_MLG[litcomp_summ_stats$Species == "A. palmata" & litcomp_summ_stats$Reference == "Present study"]) # pull allelic richness for present study      
coral_name <- "A. palmata"

# RRB allelic richnes: Acropora cervicornis
(ar_rrb <- litcomp_summ_stats$mean_Na_per_MLG[litcomp_summ_stats$Species == "A. cervicornis" & litcomp_summ_stats$Reference == "Present study"]) # pull allelic richness for present study 
coral_name <- "A. cervicornis"

# General calcs
(tot_pops <- length(litcomp_summ_stats$mean_Na_per_MLG[litcomp_summ_stats$Species == coral_name])-1)  # total number of populations being compared (subtracting out present study): 14 Apal, 7 Acer
(pops_less <- length(litcomp_summ_stats$mean_Na_per_MLG[litcomp_summ_stats$Species == coral_name & litcomp_summ_stats$mean_Na_per_MLG < ar_rrb])) # number of populations with allelic richness less than Present study: 9 apal, 3 acer
(pops_more <- length(litcomp_summ_stats$mean_Na_per_MLG[litcomp_summ_stats$Species == coral_name & litcomp_summ_stats$mean_Na_per_MLG > ar_rrb])) # number of populations with allelic richness greater than Present study: 5 apal, 4 acer
pops_less/tot_pops # proportion of populations with allelich richness lower than RRB's in present study



##### Plot relationship between sample size (Ng) and genetic diversity (Na, allelic richness, He)
ggplot(data = litcomp_summ_stats, aes(x = mean_MLGs, y = mean_Na))+
  geom_point(size = 2, aes(color = Population))+    
  geom_smooth(method = "loess")+
  theme_classic() 

ggplot(data = litcomp_summ_stats, aes(x = mean_MLGs, y = mean_Na_per_MLG))+
  geom_point(size = 2, aes(color = Population))+    
  geom_smooth(method = "loess")+
  theme_classic() 

ggplot(data = litcomp_summ_stats, aes(x = mean_MLGs, y = mean_He))+
  geom_point(size = 2, aes(color = Population))+    
  geom_smooth(method = "loess")+
  theme_classic() 



####################################################################################################################################################################################
## Rarefactions of allelic richness
####################################################################################################################################################################################
########## iNEXT ##########

##### For further info on iNEXT
# further detailed information: http://chao.stat.nthu.edu.tw/wordpress/software_download/.
#'A Quick Introductin to iNEXT via Examples' : Appendix S1 in Hsieh et al. 2016, and now an R vignette
vignette("Introduction", package="iNEXT")

##### Other iNEXT extensions/adaptations
# iNEXT.4steps is Chao's proposed protocol for comparing diversity across multiple assemblages (Chao et al. 2020)
# iNEXT.3D is Chao's proposed method for comparing 3 dimensions of diversity across/among assemblages: taxonomic, phylogenetic, functional (Chao et al. 2021)
# iNEXT.beta3D is Chao's proposed method for standardizing beta diversity for comparisons across datasets  (Chao et al. 2023)
# iNEXT.link extends 3D, beta.#D, and 4steps to ecological networks

##### Citation:
# "If you publish your work based on the results from the iNEXT package, you should make references to the following methodology paper (Chao et al. 2014) and the application paper (Hsieh, Ma & Chao, 2016)"

# Data format:
# matrix, data.frame, or list
# in data frame: rows = species (ie alleles), column = abundances




##### Load in allele abundances formatted for iNEXT
# Data format:
# matrix, data.frame, or list
# in data frame: rows = species (ie alleles), column = abundances

Apal_g1_f1_inext <- read.csv("Apal_g1_f1_inext.csv") %>% 
  column_to_rownames(var = "unique_allele_ID") %>% 
  select(allele_abundance)

Acer_g1_f1_inext <- read.csv("Acer_g1_f1_inext.csv") %>%
  column_to_rownames(var = "unique_allele_ID") %>%
  select(allele_abundance)

# sens_test_inext <- read.csv("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\sens_test_inext.csv") %>% 
#   column_to_rownames(var = "unique_allele_ID") %>% 
#   select(allele_abundance)




##### Data check
# Apal_g1_f1:  742 total alleles documented in Apal  = 34 genos * 22 allele locations (11 markers * diploid) - 6 NAs
# Acer_g1_f1:  278 total alleles documented in Acer  = 24 genos * 12 allele location (6 markers * diploid) - 10 NAs
sum(Apal_g1_f1_inext$allele_abundance)
sum(Acer_g1_f1_inext$allele_abundance)
# sum(sens_test_inext$allele_abundance)



##### Set coral and number of loci used for coral
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 

# Acropora palmata
# set_coral <- Apal_g1_f1_inext
# set_alleles <- 22   # number of alleles (diploid)
# set_loci <- 11      # number of loci
# set_coral_full <- 'Acropora palmata'

# Acropora cervicornis
set_coral <- Acer_g1_f1_inext
set_alleles <- 12  # number of alleles (diploid)
set_loci <- 6      # number of loci
set_coral_full <- 'Acropora cervicornis'


# sensitivity test
# set_coral <- sens_test_inext
# set_alleles <- 20   # number of alleles (diploid)
# set_loci <- 10      # number of loci
# set_coral_full <- 'sensitivity test Apal 166'



##### Run iNEXT
# endpoint = an integer specifying the sample size that is the endpoint for R/E calculation; If NULL, then endpoint=double the reference sample size
# knots = an integer specifying the number of equally-spaced knots (40, by default) between size 1 and the endpoint
# nboot = an integer specifying the number of bootstrap replications; default is 50
output <- iNEXT(set_coral, q = 0, datatype = "abundance", nboot = 1000) # note this takes several minutes to run with nboot = 1000
output
ChaoRichness(set_coral)




##### Calculate what % of estimated asymptote has already been sampled: 88.91% for Apal | 77.84% for Acer!

# Observed allelic diversity
ChaoRichness(set_coral)[1,1] # total alleles
ChaoRichness(set_coral)[1,1]/set_loci # average alleles (divided by number of marker loci)

# Estimated asymptotic alleleic diversity
ChaoRichness(set_coral)[1,2] # total alleles
ChaoRichness(set_coral)[1,2]/set_loci # average alleles (divided by number of marker loci)

# Percent of asympotic already observed
ChaoRichness(set_coral)[1,1]/ChaoRichness(set_coral)[1,2]





##### Set up for ggiNEXT
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 

# set asymptote
asymptote <- ChaoRichness(set_coral)[1,2]
asymptote

# set relevant breaks for x-axis (must multiply by number of loci to convert from number MLGs sampled to number alleles documented, which is what iNEXT wants to plot)
break_x_max <- (output$iNextEst$size_based[nrow(output$iNextEst$size_based),2]) # default endpoint is 2x sample size
break_x_max
set_breaks_x <- (seq(0, break_x_max*1.1, by = 10*set_alleles))
set_breaks_x
set_breaks_y <- seq(0, asymptote*1.2, by = 2*set_loci)
set_breaks_y

# set axis location of species name annotation
annot_x <- 0.7*(break_x_max)
annot_x
annot_y <- 0.1*(asymptote)  
annot_y





##### plot with ggiNEXT + modifications

#apal_inext_plot <-
acer_inext_plot <-
ggiNEXT(output, type = 1, color.var = "Order.q", grey = T)+
  
  # Re-scale y-axis to be: Average unique alleles per locus (divide by # msats loci used)
  scale_y_continuous(labels = function(y) number_format(accuracy = 1)(y/set_loci),
                     name = "",
                     breaks = set_breaks_y)+         
  
  # RE-scale x-axis: Number MLGs (divide by # number alleles documented per mlg = loci*2)
  scale_x_continuous(labels = function(x) number_format(accuracy = 1)(x/set_alleles),
                     name = "",
                     breaks = set_breaks_x)+    
  
  # Set up remaining aesthetics
  geom_hline(yintercept = asymptote, linetype="solid", color = "black", linewidth = 1)+ 
  theme_classic()+
  theme(legend.position = "none")+
  annotate("text", x= annot_x, y=annot_y, 
           label = paste0("italic('", set_coral_full, "')"),
           parse =TRUE)


# IMPORTANT: save each panel separately. This is necessary because iNEXT was re-applying one species' re-scaling (eg set_breaks) to both plots, such that one panel always had incorrect axes
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig2\\msats_fig2a.tiff", apal_inext_plot, dpi=600, units="cm", width = 8.5, height = 8.5)
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig2\\msats_fig2b.tiff", acer_inext_plot, dpi=600, units="cm", width = 8.5, height = 8.5)



fig2<-
  ggdraw()+
  draw_image("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig2\\msats_fig2a.tiff",
             x=0.025, y=.5, width=.45, height=.45)+
  draw_image("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig2\\msats_fig2b.tiff",
             x=0.025, y=0.025, width=.45, height=.45)+
  draw_plot_label(label= c("A","B"), 
                  size=panel.txt.sz,
                  x= c(0.075, 0.075), 
                  y = c(1, 0.525))+
  draw_label("Number of MLGs", x = 0.25, y = 0.02, vjust = 0, size = axis.title.sz) +
  draw_label("Average number of alleles", x = 0.02, y = 0.55, vjust = 0.5, angle = 90, size = axis.title.sz)


ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig2\\msats_fig2.tiff", fig2, dpi=600, units="cm", width = 17, height = 17)






########## per locus: pegas ##########


##### Open tiff device and set coral settings
# NOTE - in this script, each analysis is written out once, and you change whether you are analyzing the A. palmata or A. cervicornis data using "set_coral <- " 

# Close any open devices
dev.off()
dev.off()

# Apal
# tiff("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\figSX_Apal_rares.tiff", width = 18, height = 15, units = "cm", res = 600) # important for plot to get exported correctly
# set_coral <- Apal_g1_f1
# set_par <- par(mfrow = c(3, 4), mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))
# set_title <- "Acropora palmata"

# Acer
tiff("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\figSX_Acer_rares.tiff", width = 18, height = 11, units = "cm", res = 600) # important for plot to get exported correctly
set_coral <- Acer_g1_f1
set_par <- par(mfrow = c(2, 4), mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))
set_title <- "Acropora cervicornis"


##### run pegas:: functions

# Convert to loci format
locidata <- pegas::as.loci(set_coral)

# Loop over loci and plot rarefactions (suppress individual labels/axes)
for (locus in attr(locidata, "locicol")) {
  pegas::rarefactionplot(locidata[, locus, drop = FALSE],
                         xlab = "", ylab = "")
}


# Add shared axes
mtext("Sample size (MLGs)", side = 1, line = 2.5, outer = TRUE, cex = 1.2)
mtext("Allelic richness", side = 2, line = 2.5, outer = TRUE, cex = 1.2)


# Add overall italic title
mtext(bquote(italic(.(set_title))),
      side = 3, line = 1.5, outer = TRUE, cex = 1.5)


##### Close device
# important for plot to get exported correctly
dev.off()



####################################################################################################################################################################################
## Mapping & spatial analyses
####################################################################################################################################################################################
########## Load in spatial data ##########

##### Load in RRFB donor-colony coordinate data
Sourcepatches_Acp_gx <- read.csv("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genotype source patches\\Sourcepatches_Acp_coords.csv") %>% 
  mutate_at(c("Long", "Lat"), as.numeric)

remove_list_AC <- c("G02", "G21", "G25")

Source_AC <- Sourcepatches_Acp_gx %>% 
  filter(Species == "AC",
         !Genotype %in% remove_list_AC)

remove_list_AP <- c("G03", "G37", "G38", "G39", "G40", "G41", "G42", "G43", "G44")

Source_AP <- Sourcepatches_Acp_gx %>% 
  filter(Species == "AP",
         !Genotype %in% remove_list_AP)


##### Load in sub-population/genotype designation key and extra meta-data from RRB

# Subpop key
Subpops_key <- read.csv("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genotype source patches\\Subpops_key.csv") %>% 
  separate(ID, into = c("Species", "Genotype"), sep = "_")

# RRB metadata
Sourcepatches_Acropora <- read.csv("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genotype source patches\\Sourcepatches_Acropora.csv") 


##### Add additional metadata and subpop key to rrb harvesting coordinates
Source_AP_meta <- Source_AP %>% 
  mutate(Species = "AP") %>% 
  left_join(Sourcepatches_Acropora) %>% 
  left_join(Subpops_key) %>% 
  mutate(ID = paste(Species, Genotype, sep = "_"))

write.csv(Source_AP_meta, "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\Source_AP_meta.csv", row.names = F)

Source_AC_meta <- Source_AC %>% 
  mutate(Species = "AC") %>% 
  left_join(Sourcepatches_Acropora) %>% 
  left_join(Subpops_key) %>% 
  mutate(ID = paste(Species, Genotype, sep = "_"))

write.csv(Source_AC_meta, "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\Source_AC_meta.csv", row.names = F)


##### Set crs for donor-colony coordinates
harvest_points_AP <- st_as_sf(Source_AP_meta, coords= c("Long","Lat"), crs= 4326) %>%
  st_transform(crs = 32619)

harvest_points_AP_g1 <- st_as_sf(Source_AP_meta, coords= c("Long","Lat"), crs= 4326) %>%
  filter(Genotype != "G32") %>% 
  st_transform(crs = 32619)

harvest_points_AC <- st_as_sf(Source_AC_meta, coords= c("Long","Lat"), crs= 4326) %>%
  st_transform(crs = 32619)

harvest_points_AC_g2 <- st_as_sf(Source_AC_meta, coords= c("Long","Lat"), crs= 4326) %>%
  filter(Genotype != "G01") %>%  # note, G21 is already excluded because we set its harvest coordinates as unknown
  st_transform(crs = 32619)


##### Load in Bonaire shapefile 

# Load in Bonaire shapefile
bonaire.map <- sf::st_read("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Bonaire shapefiles\\BES_admin\\BES_admin\\BES_adm1.shp") %>%
  dplyr::filter(NAME_1 == "Bonaire") %>% 
  st_transform(crs= 4326) %>%                                 
  st_transform(crs = 32619)



########## Rasterize and transition matrix ##########
# Necessary for calculating shortest *marine* distance (avoiding land)

##### Create buffer around island (current shapefile stops right at island edges)
bonaire_buffered <- st_buffer(bonaire.map, dist = 2000) # 2km buffer around the island


##### Define raster extent and resolution
r <- raster::raster(raster::extent(bonaire_buffered), res = 100, crs = sp::CRS(st_crs(bonaire_buffered)$proj4string))  # resolution = 100 meters, Inf values produced with lower resolutions (eg 25, 50)


##### Assign 0 to land, NA to background (sea)
land_raster <- raster::rasterize(bonaire.map, r, field=1)
land_raster[!is.na(land_raster)] <- 0   # land = 0
land_raster[is.na(land_raster)] <- 1    # water = 1
raster::crs(land_raster) # double-check units are in meters:  +proj=utm +zone=19 +datum=WGS84 +units=m +no_defs 


##### Create transition layer; function = mean or min of adjacent cells
tr <- gdistance::transition(land_raster, transitionFunction = mean, directions = 4)


##### Correct for geographic distances
# transition layer's distances are in steps, need to be transformed back into meters. Resolution = 100 (set above in raster()), therefore...
# up/down and left/right steps = 1 step = 100m
# diagonal steps = sqrt(2) steps = 100*sqrt(2) (ie pyhtagorean theorem)
tr_corrected <- gdistance::geoCorrection(tr, type="c")



########## Calculate least-cost marine distances between donor-colonies ##########

##### Convert to SpatialPoints
# Make sure to use g1 and g2 versions here (reduced to just to MLGs being analyzed in Mantel test)
pts_AP <- as(harvest_points_AP_g1, "Spatial")
pts_AC <- as(harvest_points_AC_g2, "Spatial")
pts_AC_wG01 <- as(harvest_points_AC, "Spatial")


##### Calculate least-cost distance matrices

# Calculate
leastcost_AP <- gdistance::costDistance(tr_corrected, pts_AP)
leastcost_AC <- gdistance::costDistance(tr_corrected, pts_AC)
leastcost_AC_wG01 <- gdistance::costDistance(tr_corrected, pts_AC_wG01)

# Plot for verification it worked
#pts <- pts_AP
#pts <- pts_AC
plot(land_raster, col=c("white","blue"))
points(pts, pch=16, col="black")
lc_paths <- gdistance::shortestPath(tr, pts[1,], pts[18,], output="SpatialLines") # To determine which 2 points you're measuring between, change the pts[X,] values
plot(lc_paths, add=TRUE, col="orange", lwd=2)
lc_paths <- gdistance::shortestPath(tr, pts[18,], pts[22,], output="SpatialLines") # To determine which 2 points you're measuring between, change the pts[X,] values
plot(lc_paths, add=TRUE, col="green", lwd=2)
lc_paths <- gdistance::shortestPath(tr, pts[8,], pts[20,], output="SpatialLines") # To determine which 2 points you're measuring between, change the pts[X,] values
plot(lc_paths, add=TRUE, col="red", lwd=2)


##### Re-assign geno IDs to dist_matrix

# Apal
ids_AP <- pts_AP$ID  # or any vector of individual IDs, length = 34
rownames(leastcost_AP) <- ids_AP
colnames(leastcost_AP) <- ids_AP
leastcost_AP

# Acer
ids_AC <- pts_AC$ID  # or any vector of individual IDs, length = 23
rownames(leastcost_AC) <- ids_AC
colnames(leastcost_AC) <- ids_AC
leastcost_AC




########## Figure 1 maps ##########

##### Set map feature aesthetics
bbox <- st_bbox(bonaire.map) # set bounding box
border.lwd = 0.5
scalebar.lwd = 0.1
scale.txt = 2.25
sites.sz <- 2
coord.txt.sz <- 8


##### Acropora palmata
apal_harvest <-
  ggplot() +
  geom_sf(data = bonaire.map, fill = "gray", color = "black", lwd = border.lwd) +
  geom_sf(data = harvest_points_AP, color = "black", aes(fill = Pop), size = sites.sz, shape = 21) +
  coord_sf(
    xlim = c(bbox$xmin - 12000, bbox$xmax + 3000),  # <-- extra 1000 meters on the left
    ylim = c(bbox$ymin - 3000, bbox$ymax),
    expand = FALSE)+
  scale_fill_manual(values = pop.colors, breaks = pop.order) +
  annotation_scale(style = "bar", location = "tr", bar_cols = c("black", "black"), text_cex = 0.7, width_hint = 0.35, height = unit(0.1, "cm")) +
  theme_classic()+                
  theme(axis.text.x = element_text(color="black",            
                                   size = coord.txt.sz,
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1), 
        axis.text.y = element_text(color = "black",
                                   size = coord.txt.sz),
        legend.title = element_blank(),      
        legend.text = element_blank(),
        legend.position = "none")

ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig1\\fig1a_apal.tiff", apal_harvest, dpi=600, units="cm", width = 8.4, height = 7)




##### Acer 
acer_harvest <- 
ggplot() +
  geom_sf(data = bonaire.map, fill = "gray", color = "black", lwd = border.lwd) +
  geom_sf(data = harvest_points_AC, color = "black", aes(fill = Pop), size = sites.sz, shape = 21) +
  coord_sf(
    xlim = c(bbox$xmin - 12000, bbox$xmax + 3000),  # <-- extra 1000 meters on the left
    ylim = c(bbox$ymin - 3000, bbox$ymax),
    expand = FALSE)+
  scale_fill_manual(values = pop.colors, breaks = pop.order) +
  annotation_scale(style = "bar", location = "tr", bar_cols = c("black", "black"), text_cex = 0.7, width_hint = 0.35, height = unit(0.1, "cm")) +
  theme_classic()+                
  theme(axis.text.x = element_text(color="black",            
                                   size = coord.txt.sz,
                                   angle = 45,
                                   vjust = 1,
                                   hjust = 1), 
        axis.text.y = element_text(color = "black",
                                   size = coord.txt.sz),
        legend.title = element_blank(),      
        legend.text = element_text(color="black", size = leg.txt.sz),
        legend.position = "none")
  
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig1\\fig1a_acer.tiff", acer_harvest, dpi=600, units="cm", width = 8.4, height = 7)
  


##### Extract the legend 
fig1_legend <- get_legend(
  acer_harvest + 
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))+
    theme(legend.position = "right"))

ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig1\\fig1_legend.tiff", fig1_legend, dpi=600, units="cm", width = 2, height = 3)





########## Figure 1 histograms of nearest neighbor distances ##########

########## Convert to spatial format and calculate distances nearest neighboring donor-colony

# NOTE:
# Using simpler st_distance() method rather than least-cost distance calculations because...
# least-cost only really works here with res = 100m, which would not work for Apal which has 13 donor colonies with neast neighbor <100m away
# st_distance() is WAY less labor-intensive than a high-resolution raster
# nearest neighbors are all close enough that moving around land is not a major obstacle

##### Acropora palmata

# Convert to spatial format and calculate distances among all donor colonies
Source_AP_sf <- Source_AP %>%
  st_as_sf(coords = c("Long", "Lat"), crs = "EPSG: 4326") %>%   # Convert to spatial format
  st_distance()                                                 # Calculate all distnaces to one another
diag(Source_AP_sf) <- NA                                        # Set diagonal to NA (want to remove distance to self)

# Reduce dataset to nearest neighbors
ap_nn <- apply(Source_AP_sf, 1, FUN = min, na.rm = T)           # Find the minimum distance to each site (so, neaerest neighbor)
ap_nn_df <- as.data.frame(ap_nn)              # Convert list to data frame
ap_nn_df <- data.frame(                                     # Add back in meta data of each genotype
  Genotype = Source_AP$Genotype,
  NearestNeighborDistance = ap_nn
) %>%
  left_join(Source_AP_meta)

# Some basic calcs
sum(ap_nn_df$NearestNeighborDistance>300)                       # Number of lineages with nearest neighbor >300m away: 22                 
sum(ap_nn_df$NearestNeighborDistance>300)/35                    # Percent of lineages with nearest neighbor >300m away: 62.9%
sum(ap_nn_df$NearestNeighborDistance<300)                       # Number of lineages with nearest neighbor <300m away: 13
sum(ap_nn_df$NearestNeighborDistance<300)/35                    # Percent of lineages with nearest neighbor <300m away: 37.1%




##### Acropora cervicornis

# Pull least-cost distance between AC_G01 to AC_G03 (since the euclidean distnace calculated below is an 'as the crow flies' under-estimate)
# land_raster created in "Rasterize and transition matrix"
# tr_corrected created in "Rasterize and transition matrix"
# pts_AC_wG01 created in "Calculate least-cost..."
pts <- pts_AC_wG01
plot(land_raster, col=c("white","blue"))
points(pts, pch=16, col="black")
lc_paths <- gdistance::shortestPath(tr_corrected, pts[1,], pts[4,], output="SpatialLines") # To determine which 2 points you're measuring between, change the pts[X,] values
plot(lc_paths, add=TRUE, col="green", lwd=3)
(AC_G01_lcd_m <- gdistance::costDistance(tr_corrected, pts[1,], pts[4,]))
(AC_G01_lcd_km <- AC_G01_lcd_m/1000)  # manually checked in Google Earth - 11.8km is accurate


# Convert to spatial format and calculate distances among all donor colonies
Source_AC_sf <- Source_AC %>%
  st_as_sf(coords = c("Long", "Lat"), crs = "EPSG: 4326") %>%   # Convert to spatial format
  st_distance()                                                 # Calculate all distnaces to one another
diag(Source_AC_sf) <- NA                                        # Set diagonal to NA (want to remove distance to self)

# Reduce dataset to nearest neighbors
ac_nn <- apply(Source_AC_sf, 1, FUN = min, na.rm = T)           # Find the minimum distance to each site (so, nearest neighbor)
ac_nn_df <- as.data.frame(ac_nn)                                # Convert list to data frame
ac_nn_df <- data.frame(                                         # Add back in meta data of each genotype
  Genotype = Source_AC$Genotype,
  NearestNeighborDistance = ac_nn
) %>%
  left_join(Source_AC_meta) %>% 
  mutate(NearestNeighborDistance = ifelse(Genotype == "G01", AC_G01_lcd_m, NearestNeighborDistance)) # replacing G01's distance with least-cost marine distance around the island (rather than as the crow flies *over* the island). This is the only nn case that requires least-cost calculation (because AC_G01 is in Lac Bay by itself)

# Some basic calcs
(sum(ac_nn_df$NearestNeighborDistance>300)/24)*100              # Percent of lineages with nearest neighbor >300m away: 100%



##### Histograms!

# close any tiff devices that may have been previously open
dev.off()
dev.off()

# open tiff device
tiff('C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\fig1\\fig1b_matchX.tiff', units="cm", width = 16, height = 5, res=600)

# set panels and margins and tiff environment
par(mfrow = c(1,2), mar = c(1.5, 1.5, 0.5, 0.5), oma = c(1.5, 1.5, 0.5, 0.5)) # c(bottom, left, top, right) # oma = outer margins

# Apal
#set_title <- expression(italic("Acropora palmata"))
hist(ap_nn,
     xlab = "",
     ylab = "",
     main = "",
     xlim = c(0,12000),
     #cex.lab = 1,
     cex.axis = 0.75,
     #breaks = 20
     )

# Acer
#set_title <- expression(italic("Acropora cervicornis"))
hist(ac_nn_df$NearestNeighborDistance,
     xlab = "",
     ylab = "",
     main = "",
     #cex.lab = 1,
     cex.axis = 0.75,
     breaks = 20) 

# Add shared axes
# side: 1=bottom, 2=left, 3=top, 4=right
mtext("Distance to nearest neighbor (m)", side = 1, line = 0.5, outer = TRUE, cex = 0.75)  
mtext("Number of donor colonies", side = 2, line = 0.5, outer = TRUE, cex = 0.75)


# dev.off
dev.off()

#



########## Percent patches harvested from ##########

##### Load in stinapa survey areas
survey.areas <- sf::st_read("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Acropora patch mapping\\STINAPA_Mapped_Areas.shp") %>%
  st_transform(crs= 3857) %>%                                 # QGIS layer was EPSG: 3857
  st_transform(crs = 32619)                                   # transform to 32619 (UTM 19N) because it enables st_intersects(): 4326 is coordinates, so unites are in degrees, but 32619 is a projected crs (flattened to 2D) so its units are in meters



##### Calculate total shelf area covered in 19 Mapped Areas
survey.areas$area_calculated <- st_area(survey.areas)
shelf_area <- sum(survey.areas$area_calculated)
shelf_area
(shelf_area_km2 <- shelf_area/(10^6))



##### Calculate percent of coastline captured by survey area polygons

# Extract coastline as linestring
coastline <- st_boundary(bonaire.map)

# Add buffer to survey area polygons because survey areas dont fully overlap with coastline, given that they're offshore
survey.areas.buf <- st_buffer(survey.areas, dist = 20) # 20 meter buffer

# Double-check by visualising survey areas along Bonaire
tmap_mode("view")
tm_shape(bonaire.map)+
  tm_polygons(fill = "tan")+
  tm_shape(survey.areas.buf)+
  tm_polygons(fill = "blue")

# Get intersections of survey polygons and coastline: rows are segments of coastline that are overlapped by stinapa survey polygons
coast_captured <- st_intersection(coastline, survey.areas.buf)

# Coastline length
total_coast_length <- sum(st_length(coastline))     
total_coast_length                                  # 146198.4 meters
total_coast_length/1000                             # 146.1984 kilometers (seems right, cursory perimeter of Bonaire + Klein = 110km in Google Earth)

# Overlap length
captured_length <- sum(st_length(coast_captured))
captured_length                                     # 14151.63 meters  
captured_length/1000                                # 14.15163 kilometers

# Percent overlap
percent_captured <- as.numeric(captured_length / total_coast_length * 100)
percent_captured  




##### Set coral species and convert RRB harvesting source patch coordinates to sf object using Lat/Long columns
# Load in as original crs (4326, units = degrees) then convert to projected crs in UTM zone 19n (EPSG: 32619, units = meters)

# Acropora palmata
# species <- "APAL"
# rrb.points <- st_as_sf(Source_AP_meta, coords= c("Long","Lat"), crs= 4326) %>%
#   st_transform(crs = 32619)
# set_title <- "Acropora palmata"

# Acropora cervicornis
species <- "ACER"
rrb.points <- st_as_sf(Source_AC_meta, coords= c("Long","Lat"), crs= 4326) %>%
  st_transform(crs = 32619)
set_title <- "Acropora cervicornis"




##### Load in and tidy up Acropora patch polygons shapefile 
stinapa.map <- sf::st_read("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Acropora patch mapping\\STINAPA_mapping_patches.shp") %>% 
  st_transform(crs= 4326) %>%                                 # QGIS layer was EPSG: 4326
  st_transform(crs = 32619) %>%                               # transform to 32619 (UTM 19N) because it enables st_intersects(): 4326 is coordinates, so unites are in degrees, but 32619 is a projected crs (flattened to 2D) so its units are in meters
  mutate(Species = recode(Species, "ACEr" = "ACER")) %>%      # clean up typo
  rename("Long" = "Lon") %>%                                  # clean up typo (necessary?)
  filter(Species != "NA",                                     # remove entries in which Species is unknown
         #Area_m2 > 5,                                       # remove very small polygons, which are more likely to be singular colonies than actual patches, and less therefore likely to have been fully established                               
         LiveTissue > 1,                                      # remove polygons that were mostly dead 
         Density > 1                                         # remove polygons that were very sparse                                                             
  ) %>% 
  filter(Species == species)                                  # filter for only Acer or Apal, depending on "species <- "_" " designation above




##### Edit polygons from stinapa map

# Step 1. Create a buffer around polygons (dist is in meters)
buffered <- st_buffer(stinapa.map, dist = 40) 

# Step 2. Combine overlapping polygons to get a single polygon anywhere individual polygons overlapped. You can plot to check what distance is working best for your data and then refine
combined <- st_union(buffered) %>%
  st_cast("POLYGON") 

# Step 3. Convert the combined polygons from sfc format back to sf
combined_sf <- st_sf(polygon_id = 1:length(combined),
                     geometry = combined) 

# Visualize 3 steps to double check they worked properly
# tmap_mode("view")
# tm_shape(stinapa.map)+tm_polygons()   # original polygons
# tm_shape(buffered)+tm_polygons()      # polygons with 20m buffer (overlapping)
# tm_shape(combined_sf)+tm_polygons()   # overlapped polygons merged
# tm_shape(combined_sf)+tm_polygons()+  # merged polygons with harvesting points overlaid
#   tm_shape(rrb.points)+
#   tm_dots(fill = "Active", size = 0.5)+
#   tm_labels("Active", options = opt_tm_text(point.label = T, point.label.gap = 2))




##### Create a TRUE/FALSE matrix for all intersections: rows are mapped polygons and columns are RRB points (35 Apal lineages and 24 Acer lineages with origin location data)
# TRUE = intersection occurs = 1
# FALSE = no intersection = 0
poly_with_points <- st_intersects(combined_sf, rrb.points, sparse = F) 

# Reduce across polygons (columns) to single TRUE/FALSE per point
point_overlaps <- colSums(poly_with_points) > 0  

# Add as new column to rrb.points
rrb.points.overlap <- rrb.points
rrb.points.overlap$overlap <- point_overlaps




##### Calculate percent mapped stands that RRB sampled from
# For each polygon (row) add up the total number of 'TRUE's (intersections) using rowSums(). Then, add up the number of polygons that had at least 1 point intersection (rowSums > 0) and divide by the total number of polygons (nrow()) --> giving you percent of polygons that have been harvest from at least once

# raw number of mapped stands (ie merged polygons) that RRB harvested from at least once
rrb.stands.intersect <- sum(rowSums(poly_with_points) > 0)       
rrb.stands.intersect

# raw number of mapped stands (ie merged polygons) that RRB never harvested from 
rrb.stands.never.intersect <- sum(rowSums(poly_with_points) == 0)       
rrb.stands.never.intersect

# total number of mapped stands (ie merged polygons)
mapped.stands <- nrow(poly_with_points)                     
mapped.stands 
rrb.stands.never.intersect + rrb.stands.intersect # double-check: should = mapped.stands

# percent of mapped stands harvested from
rrb.stands.intersect/mapped.stands * 100   

# percent of mapped stands never harvested from
rrb.stands.never.intersect/mapped.stands * 100 




##### Create a TRUE/FALSE matrix for all intersections: rows are RRB points (35 Apal lineages, 26 Acer lineages) and columns are survey area boxes (19)
# TRUE = intersection occurs = 1
# FALSE = no intersection = 0
points_with_surveyareas <- st_intersects(rrb.points, survey.areas, sparse = F) 




##### Calculate raw number of RRB harvesting locations inside and outside of survey area boxes

# Number of harvesting locations OUTSIDE of survey areas
rrb.stands.outside <- sum(rowSums(points_with_surveyareas) == 0)                     
rrb.stands.outside

# Total number of stands harvested from INSIDE survey areas (whether intersected mapped polygon or not)
rrb.stands.inside <- sum(rowSums(points_with_surveyareas) > 0)                     
rrb.stands.inside

# Double-check:  outside +inside = should equal total # lineages for that species
rrb.stands.outside+rrb.stands.inside

# Number of stands harvested from INSIDE survey areas BUT do NOT INTERSECT mapped polygons
inside.non.intersect <- rrb.stands.inside - rrb.stands.intersect
inside.non.intersect

# Number of stands harvested from that were not in STINAPA's maps (both inside and outside survey areas)
rrb.stands.outside + inside.non.intersect

# NOTE - the Acer calcs are off by 1 because G24 was just outside of a survey area, but got included in a polygon with added buffer
# So it got included in the poly_with_points as an intersection but not in the surveyareas_with_points
# For now, I just adjusted the test be 14 outside and 5 inside (instead of 15 out, 4 in), but a long term solution is to update the st_intersects code for inside/outside (or fudge it by buffer survey area polygon as well)




##### Visualizing maps

# set view (play) vs plot (real) modes
tmap_mode("view")
#tmap_mode("plot")

# base play map
tm_shape(survey.areas)+tm_polygons(fill = "lightblue")+                      # 19 surveyed areas
  tm_shape(combined_sf)+tm_polygons(fill = "orange")+                                # merged polygons
  tm_shape(rrb.points)+                                                              # rrb harvesting locations
  tm_dots(size = 0.75,
          fill = "black"
          #fill = "Active",
          #fill.scale = tm_scale(value = c(Active = "black", Dead = "red"))
  )#+
#tm_labels("Genotype", options = opt_tm_text(point.label = T, point.label.gap = 2)) 


##### overlap map - Figure S1
coral_title <- "Acropora palmata"
#coral_title <- "Acropora cervicornis"

tmap_mode("plot")


##### Overlap map (Figure S3)

# Set map feature aesthetics
map.title.sz = 2
border.lwd = 3
scalebar.lwd = 10
comp.sz = 1.5
comp.txt = 1.5
sites.sz <- 0.75
harvest.cols <- c("navy", "orange2")

# Map 
#overlapmap_AP <-
  #overlapmap_AC <-
  tm_shape(bonaire.map)+
  tm_polygons(lwd = border.lwd)+
  tm_shape(rrb.points.overlap)+
  tm_symbols(size = sites.sz,
             fill = "overlap",
             fill.scale = tm_scale(values = harvest.cols))+
  #tm_scalebar(position = c("left", "bottom"), breaks = c(0,10), text.color = "white", text.size = 0, color.light = "black", lwd = scalebar.lwd)+
  #tm_compass(position = c("right", "top"), size = comp.sz, text.size = comp.txt)+
  tm_layout(frame = FALSE, 
            legend.show = F)


tmap_save(filename = "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\overlapmap_AP.tiff", overlapmap_AP, dpi=600, units="cm", width = 8.5)
#tmap_save(filename = "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\overlapmap_AC.tiff", overlapmap_AC, dpi=600, units="cm", width = 8.5)





####################################################################################################################################################################################
## PCA
####################################################################################################################################################################################
########## Impute missing data ##########


##### Set coral
#set_coral <- Apal_g1_f1
set_coral <- Acer_g2_f1
#set_coral <- sens_test

##### Check missing data first

# Use info_table to summarize and visualize missing data
info_table(
  set_coral,
  type = "missing",
  percent = FALSE,
  plot = TRUE,
  df = FALSE,
  returnplot = FALSE,
  low = "blue",
  high = "red",
  plotlab = FALSE,
  scaled = TRUE
)

# Count NAs in genind object manually
# NOTE - there will be more missing values here than actual missing alleles because of how the data is organized in genind
sum(is.na(set_coral@tab))

# Look at actual data to verify  sum of NAs accurately represents missing alleles
options(max.print = 1e6)  # must set max print larger in order to see full allele frequency table
set_coral@tab
# 18 NAs in Apal because 1195 (n = 6 alleles) is missing data from AP_G02, AP_G30, AP_G19, so 6 possible alleles * 3 individuals missing any allele data = 18
# 31 NAs in Acer because 182 (n = 8 alleles) in missing data from 2 individuals, 1490 (n = 5 alleles) is missing data from 2 individuals, and 166 (n = 5 alleles) is missing data from 1 individual. Altogether, (2*8)+(2*5)+(1*5)=31 

##### Impute missing data with means (because PCA cannot handle missing data) and convert to dataframe
# Number of missing values replaced should = sum(is.na(set_coral@tab))
set_coral@tab                                             # Check original allele frequency data. Note NAs and that data are all integers ranging 0-2 (ie rates of occurrence of each allele per individual)
impute_coral <- missingno(set_coral, type = "mean")       # Impute missing allele values (NAs in @tab). You'll get a warning message about integers  
impute_coral@tab                                          # Check imputed allele frequencies. Note that the data is all numeric decimals, no longer integers                   

# Reset max.print to a reasonable number of elements
options(max.print = 1000)



########## Run PCA ##########


##### Run PCA via ade4

# Exploratory run to see screeplot
# pca.exp <- dudi.pca(impute_coral, center=TRUE, scale=FALSE) 

# Choose number of axes to retain
# dudi.pca() above returned a barplot of eigenvalues (screeplot) and asks for number of retained principal components
# eigenvalues represent the amount of genetic diversity represented by each PC
# A sharp decrease in the eigenvalues is usually indicative of the boundaries between relevant structures and random noise." - Jombart 2016 Introduction...

# Final runs
#pca.apal <- dudi.pca(df = impute_coral, center = TRUE, scale = FALSE, scannf = FALSE, nf = 4) # 4 PCs retained for Apal
pca.acer <- dudi.pca(df = impute_coral, center = TRUE, scale = FALSE, scannf = FALSE, nf = 4)  # 4 PCs retained for Acer

# Exploratory plotting 

# Apal
#s.label(pca.apal$li) # with geno labels
#s.class(pca.apal$li, fac=pop(impute_coral), col=funky(15)) # with population labels

# Acer
# s.label(pca.acer$li) # with geno labels
# s.class(pca.acer$li, fac=pop(impute_coral), col=funky(15)) # with population labels




##### Set up for plotting PCA from ade4 via ggplot

# Turn scores from PCA into dataframe
#scores_df <- as.data.frame(pca.apal$li)
scores_df <- as.data.frame(pca.acer$li)

# Add metadata to scores
scores_df$ind <- rownames(scores_df)
scores_df$pop <- pop(set_coral) 
scores_df$pop <- recode(
  pop(set_coral),
  "LacBay" = "Lac Bay"
)

# Pull polygon values
hull_df <- scores_df %>%
  group_by(pop) %>%
  slice(chull(Axis1, Axis2))

# Pull PC1 and PC2 values - Apal
# eig <- pca.apal$eig                       # eigenvalues
# var_exp <- round(100 * eig / sum(eig), 1) # % variance explained
# pc1_pct <- var_exp[1]                     # PC1 percentage
# pc2_pct <- var_exp[2]                     # PC2 percentage

# Pull PC1 and PC2 values - Acer
eig <- pca.acer$eig                       # eigenvalues
var_exp <- round(100 * eig / sum(eig), 1) # % variance explained
pc1_pct <- var_exp[1]                     # PC1 percentage
pc2_pct <- var_exp[2]                     # PC2 percentage

# Apal settings
set_title <- expression(italic("Acropora palmata"))
set.legpos <- "none"

# Acer settings
set_title <- expression(italic("Acropora cervicornis"))
set.legpos <- "none"




##### Plot ade4 PCA with ggplot

# Plot using ggplot
#apal_pca <-
acer_pca <-
ggplot(scores_df, aes(x = Axis1, y = Axis2, color = pop)) +
  geom_hline(yintercept=0, linetype="dashed", color = "gray")+
  geom_vline(xintercept=0, linetype="dashed", color = "gray")+
  geom_point(size = point.sz) +
  scale_colour_manual(values = pop.colors,
                      breaks = pop.order)+
  scale_fill_manual(values = pop.colors,
                    breaks = pop.order)+
  geom_polygon(data = hull_df,
               aes(x = Axis1, y = Axis2, group = pop, fill = pop, color = pop),
               alpha = 0.25,
               inherit.aes = FALSE) +
  theme_classic()+
  ggtitle(set_title)+
  labs(
    x = paste0("PC1 (", pc1_pct, "%)"),
    y = paste0("PC2 (", pc2_pct, "%)")) +
  theme(axis.title.x = element_text(color="black", size = axis.title.sz), 
        axis.title.y=element_text(color="black", size = axis.title.sz),
        axis.text.x = element_text(color="black", size = axis.txt.sz), 
        axis.text.y = element_text(color="black", size = axis.txt.sz),
        legend.title = element_blank(),
        legend.text = element_text(size = leg.txt.sz),
        plot.title = element_text(size = title.sz, hjust = 0.5),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"))



# Extract the legend fro apal_pca
legend <- get_legend(
  apal_pca + 
    guides(color = guide_legend(override.aes = list(shape = 15)),
           fill = guide_legend(override.aes = list(shape = 15))) +
    theme(legend.position = "bottom")         # position at bottom
)


# Combine top and bottom plots
combined_plots <- plot_grid(
  apal_pca,
  acer_pca,
  ncol = 1,
  align = "v",
  labels = c("A", "B"),  # optional labels
  label_size = 13,
  rel_heights = c(1, 1)
)


# Add legend at the very bottom
fig3 <- plot_grid(
  combined_plots,
  legend,
  ncol = 1,
  rel_heights = c(1, 0.1)  # legend takes 10% of height
)


# Save!
ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\msats_fig3.tiff", fig3, dpi=600, units="cm", width = 8.5, height = 17)





####################################################################################################################################################################################
## find.clusters() via adegenet
####################################################################################################################################################################################

##### find.clusters() was intended to be used to inform K in DAPC (adegenet, Jombart)
# For more detail on methods, see tutorials_popgen.R script


##### Set up data

# Set coral
set_coral <- Apal_g1_f1   
#set_coral <- Acer_g2_f1  
#set_coral <- sens_test


# Set max number clusters depending on data set (must be 1 less than the totla number of individuals)
max_clusters <- length(indNames(set_coral)) - 1    # 33 for Apal, 22 for Acer


##### Use find.clusters() to identify clusters

# Run find.clusters()
# Note - do NOT use the n.clust argument here. If you use n.clust, then only that single k value is assessed. Instead, just enter the number interactively after BIC plot created
grp <- find.clusters(set_coral, max.n.clust= max_clusters)
grp <- find.clusters(set_coral)

# Choose the number PCs to retain (>= 1): From Jombart & Collins 2022 dapc tutorial: "Apart from computational time, there is no reason for keeping a small number of components" (they used ~double the #PCs on returned plot, in order to ensure all actual PCs were retained)
# Choose the number PCs to retain (>= 1): 50 Apal
# Choose the number PCs to retain (>= 1): 50 Acer

# END: we infer admixture here because there should be an elbow of minimum BIC value for a biologically realistic k. Note, k cannot = 1, so given how wonky these BIC plots look, I conclude admixture. Given the admixture, DAPC provides no benefit over PCA, so I will pursue the more commonly-understood PCA for multi-variate visualization of allele frequencies






####################################################################################################################################################################################
## F-statistics & analogues (G''st, Jost's D, Gst, theta, Fis) pairwise comparisons among sub-populations and sig. diff. from 0
####################################################################################################################################################################################


########## Fst - Comparing different methods ##########

##### set coral species settings

# Acropora palmata - 11 loci
# lab_order <- c("North", "Central", "South", "Klein", "LacBay")  
# set_coral_fst <- Apal_g1_f1
# set_name_sp <- "AP.csv"

# Acropora palmata - 6 loci (reduced to same 6 loci in Acer for comparability)
# lab_order <- c("North", "Central", "South", "Klein", "LacBay")  
# set_coral_fst <- Apal_g1_f3
# set_name_sp <- "AP6.csv"

# Acropora cervicornis - excluding Lac Bay from Acer sub-pop analyses
lab_order <- c("North", "Central", "South", "Klein")             
set_coral_fst <- Acer_g2_f1
set_name_sp <- "AC.csv"

# Sensitivity test for removing certain loci
#set_coral_fst <- sens_test



##### Run Gst Nei - mmod
# This function calculates Nei's Gst between all combinations of populaitons in a genind objec (following Nei's method and using Nei and Chesser's estimators for Hs and Ht).
# Nei 1973, Nei and Chesser 1983
(gst_Nei <- pairwise_Gst_Nei(set_coral_fst, linearized = FALSE) %>%  
    round(digits = 3))

##### Run Nei 1987 - hierftstat
# note - This is misleadingly NOT the same thing as Nei's Gst (just running for exploratory purposes)
# (genet.dist(set_coral_fst, method = "Nei87") %>%  
#     round(digits = 3))

##### Run theta Weir & Cockerham 1984 - hierfstat
# note - only works on diploid or haploid data
(theta_WC84 <- genet.dist(set_coral_fst, method = "WC84") %>%  
    round(digits = 3))

##### Run G''st Hedrick - mmod
# This function calculates G''st, a measure of genetic differentiation, between all combinations of populaitons in a genind object
# G''st is Nei's Gst adjusted for # alleles in subpop (because Nei's Gst cannot reach 1 when mutation rate is high) and sample size of subpopulations
# Hedrick 2005, Meirmans and Hedrick 2011
(gpst_Hedrick <- pairwise_Gst_Hedrick(set_coral_fst, linearized = FALSE) %>%  
    round(digits = 3))

##### Run Jost's D - mmod
(JostD <- pairwise_D(set_coral_fst, linearized = FALSE, hsht_mean = "arithmetic") %>%  
    round(digits = 3))





########## Fst - Visualise pairwise  ##########

##### For each coral data type (AP 11 loci, AP 6 loci, AC), set input matrix 

# set_matrix <- gst_Nei
# set_name_mat <- "gst_"

# set_matrix <- theta_WC84
# set_name_mat <- "theta_"

# set_matrix <- gpst_Hedrick
# set_name_mat <- "gpst_"

set_matrix <- JostD
set_name_mat <- "Jostd_"


##### Create table of Fst and other methods
# Note - just using the generic term fst below, but the code is actually for theta, G'st, and Jost's D

# Finalize formatting
fstat.mat <- as.matrix(set_matrix)     # Convert to matrix
fstat.mat1 <- fstat.mat[lab_order,]      # Change order of rows and cols
fstat.mat2 <- fstat.mat1[, lab_order]    # Change order of rows and cols
fstat.mat2[upper.tri(fstat.mat2)] <- NA  # Convert upper triangle values (above diagonal) to NAs
fstat_table <- as.data.frame(fstat.mat2) # Convert to data frame 
fstat_table

# Save csv
set_filepath <- "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\"
write.csv(fstat_table, paste(set_filepath, set_name_mat, set_name_sp, sep = ""), row.names = F)



########## Fst - Create exploratory heat map as visual aid (Jenkins tutorial) ##########

# Create a data.frame
ind <- which(upper.tri(fst.mat2), arr.ind = TRUE)
fst.df <- data.frame(Site1 = dimnames(fst.mat2)[[2]][ind[,2]],
                    Site2 = dimnames(fst.mat2)[[1]][ind[,1]],
                    Fst = fst.mat2[ ind ])

# Keep the order of the levels in the data.frame for plotting 
fst.df$Site1 = factor(fst.df$Site1, levels = unique(fst.df$Site1))
fst.df$Site2 = factor(fst.df$Site2, levels = unique(fst.df$Site2))

# Convert minus values to zero
fst.df$Fst[fst.df$Fst < 0] = 0

# Fst italic label
fst.label = expression(italic("F")[ST])

# Extract middle Fst value for gradient argument
mid = max(fst.df$Fst)/2

# Plot heatmap
ggplot(data = fst.df, aes(x = Site1, y = Site2, fill = Fst))+
  geom_tile(colour = "black")+
  geom_text(aes(label = Fst), color="black", size = 3)+
  scale_fill_gradient2(low = "blue", mid = "pink", high = "red", midpoint = mid, name = fst.label, limits = c(0, max(fst.df$Fst)), breaks = c(0, 0.05, 0.10))+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0), position = "right")+
  theme(axis.text = element_text(colour = "black", size = 10, face = "bold"),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 10)
  )




########## Fst - Is it significantly different from 0? ##########

##### Perform bootstrapping across pairwise F-statistic values to see if Global Estimate is significantly different form zero

##### Set coral and label order
# set coral
#set_coral_boot <- Apal_g1_f1
#set_coral_boot <- Apal_g1_f3
#set_coral_boot <- Acer_g2_f1
#set_coral_boot <- sens_test


##### mmod - G''st Hedrick
# Chao bootstrap
set.seed(123) # to make the results reproducible, otherwise the values will be slightly different each time you run the bootstrapping.
bs_G <- chao_bootstrap(set_coral_boot)
summarise_bootstrap(bs_G, Gst_Hedrick)
# Global Estimates based on average heterozygosity:
# Apal (11 loci):      -0.0099	(-0.083 - 0.063, 95% CI)
# Apal (6 Acer loci):  -0.0294	(-0.122 - 0.064, 95% CI)
# Acer:                 0.0198	(-0.137 - 0.177, 95% CI)


##### mmod - Jost's D
# Chao bootstrap
# Warning message: harmonic mean is undefined - this is fine, it's just because some Jost's D estimates were negatives (ie 0). Use Global estimate based on average heterozygosity, not harmonic mean
set.seed(123) # to make the results reproducible, otherwise the values will be slightly different each time you run the bootstrapping.
bs_D <- chao_bootstrap(set_coral_boot)
summarise_bootstrap(bs_D, D_Jost)
# Global Estimates based on average heterozygosity:
# Apal (11 loci):     -0.0063	(-0.055 - 0.043, 95% CI)    
# Apal (6 Acer loci): -0.0183	(-0.080 - 0.043, 95% CI)
# Acer:               0.0142	(-0.105-0.134, 95% CI)        



##### hierfstat - Weir and Cockerham 1984
# hierfstat::pp.fst() does pairwise fst comparisons following Weir and Cockerham 1984 n(https://rdrr.io/cran/hierfstat/man/ppfst.html), so I would assume boot.ppst does as well. Plus WC84 seems to be the primary method of estimating Fst for hierfstat
set.seed(123) # to make the results reproducible, otherwise the values will be slightly different each time you run the bootstrapping. 
boot.ppfst(set_coral_boot, nboot=1000, quant=c(0.025,0.975), diploid = T)



##### Other potential methods that I was getting error messages for. Unresolved.

# Jacknife populations
# NOTE - this is not working, returns 'Error: unable to find an inherited method for function 'pop' for signature 'x = "NULL"' 
# In https://github.com/dwinter/mmod/issues Winter hasn't responded since 2018 and multiple people have had issues relating to bootstrap
# ChatGPT insists the error is in summarise_bootstrap. My guess is the error is in jacknife_populations
# jacknife doesn't even work on nancycats, so the problem is beyond my dataset
# data(nancycats)
# bs_jk <- mmod::jacknife_populations(nancycats, sample_frac = 0.8, nreps = 1000)
# mmod::summarise_bootstrap(bs_jk, Gst_Hedrick)


# StAMPP - G'st Hedrick
# NOTE - Was getting error messages, not resolved
# set.seed(123) # to make the results reproducible, otherwise the values will be slightly different each time you run the bootstrapping.
# stamppGst(set_coral_boot,
#   nboots = 1000,
#   percent = 95
# )







########## Fis - Comparing different methods ##########


##### hierftstat - Nei 1987

# Use in final paper because it subsets by subpopulation

# set coral species settings
#set_coral_fis <- Apal_g1_f1
#set_coral_fis <- Apal_g1_f3
set_coral_fis <- Acer_g2_f1
#set_coral_fis <- sens_test
  
fis_Nei87 <- basic.stats(set_coral_fis, diploid = TRUE)
fis_Nei87$Fis
apply(fis_Nei87$Fis, MARGIN = 2, FUN = mean, na.rm = TRUE) %>%
  round(digits = 3)
fis_Nei87$overall[9]






##### hierfstat - Weir & Cockerham 1984

# For exploratory purposes only, not using in final paper. Can get Fis per locus and overall, but not divided by subpops. Results similar enough to Nei87

# # set coral species settings
# # NOTE - must be in data frame form for wc()
# #set_coral_WCfis <- Apal_wG32_df
# #set_coral_WCfis <- Apal_wG32_df %>% select(X182_FAM, X0585_NED, X5047_FAM, X181_NED, X1490_NED, X166_PET)
# set_coral_WCfis <- Acer_df %>%  rownames_to_column() %>%   filter(rowname != "AC_G01", rowname != "AC_G21") %>% select(X182_FAM, X0585_NED, X5047_FAM, X181_NED, X1490_NED, X166_PET, -rowname)
#   
# 
# # set data in correct format for wc(): must be data frame, first column is Pop, rest of columns are loci with 1 column per locus, no separator punctuation between alleles, as.numeric()
# # NAs induced by coercion is okay (turning the NANA character into NA numeric)
# set_coral_WCfis[, 2:ncol(set_coral_WCfis)] <- lapply(set_coral_WCfis[, 2:ncol(set_coral_WCfis)], function(x) {as.numeric(gsub("/", "", x))})
# 
# # run wc()
# fis_WC84 <- wc(set_coral_WCfis, diploid = T, pol = 0.0)
# fis_WC84$FIS






####################################################################################################################################################################################
## Mantel tests
####################################################################################################################################################################################

##### Make sure geographic data is loaded in
# If not, run  Mapping & spatial analyses > sections 1 [Load...]  --> 3 [Calculate...]
leastcost_AP
leastcost_AC

##### Convert to dist object for Mantel test

# Apal
dist_obj <- as.dist(leastcost_AP)
dist_obj

# Acer
dist_obj <- as.dist(leastcost_AC)
dist_obj

# Convert from meters to kilometers
dist_obj_km <- dist_obj/1000
dist_obj_km # make sure to double check in Google Earth that these are correct


##### Calculate genetic distances among individuals

# set coral
#set_coral <- Apal_g1_f1
set_coral <- Acer_g2_f1
#set_coral <- sens_test
# dist_matrix <- dist_matrix_AP
# dist_matrix <- dist_matrix_AC

# Bruvo's distance among individuals
# mutstep = 1 for stepwise mutations (default)
(dist_bruvo <- bruvo.dist(set_coral, replen = c(rep(3, 11))))



##### Reorder Bruvo distnace matrix to match geographic distance  
# !! a crucial step to make sure Mantel test is run properly !!

# check genetic distance matrix
# rownames(dist_bruvo)
# colnames(dist_bruvo)

# check geo distance matrix
# rownames(dist_matrix)
# colnames(dist_matrix)

# Match
dist_bruvo2 <- as.dist(as.matrix(dist_bruvo)[labels(dist_obj), labels(dist_obj)])
rownames(dist_bruvo2)
colnames(dist_bruvo2)




##### Mantel test to compare genetic and geographic distances using both ade4 and vegan 

# Choose your [genetic distance] weapon
genedist <- dist_bruvo2

# Choose your [physical distance] weapon
physdist <- dist_obj_km


# ade4 Mantel
# 'complete enumeration' warning means there are <999 possible unique permutations (because small sample size), so instead ade4 just works with all possible permutations
# ade4::mantel.randtest(genedist, physdist, 999) # For APAL: correlation = 0.01776595, p-value = 0.346 | For ACER: correlation =  0.03660642, p = 0.299 | Overall: no, genetic and geo distances are not more correlated than random expectation

# vegan Mantel (vegan supposedly more flexible because it does not mind the negative genetic distance values)
m <- vegan::mantel(genedist, physdist)                 
m # For APAL: Mantel statistic r = 0.01629, p = 0.368 | For ACER: correlation =  0.04469,   p = 0.245 | Overall: same pattern as above  







##### Plot

# Set title
set_title <- expression(italic("Acropora palmata"))
#set_title <- expression(italic("Acropora cervicornis"))

# Combine genetic and physical distance data
genedist_vec <- as.vector(genedist)
physdist_vec <- as.vector(physdist)
mantel_vecs <- cbind(genedist_vec, physdist_vec)


# Plot genetic distance as a function of physical distance

apal_mantel <-
#acer_mantel <- 
  ggplot(data = mantel_vecs, aes(x = physdist_vec, y = genedist_vec))+
  geom_point(size = point.sz-1, color = "black")+
  scale_y_continuous(limits = c(0.15, 0.75), breaks = seq(0.2, 0.8, by = 0.1))+   # Apal range 0.19-0.65; Acer range: 0.17-0.73
  #scale_x_continuous(limits = c(0, 45), breaks = seq(0, 45, by = 10))+           # Apal range 0-44;      Acer range: 0.3-30m
  theme_classic()+
  ggtitle(set_title)+
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        axis.text.x = element_text(color="black",size = axis.txt.sz), 
        axis.text.y = element_text(color="black",size = axis.txt.sz),
        plot.title = element_text(size = title.sz, hjust = 0.5))


fig4 <- plot_grid(
  apal_mantel,
  acer_mantel,
  ncol = 1,
  align = "v",
  labels = c("A", "B"),  # optional labels
  label_size = 13,
  rel_heights = c(1, 1)
)

fig4_labeled <- ggdraw() +
  # Add the plots, shrink them slightly in y/x
  draw_plot(fig4, x = 0.1, y = 0.1, width = 0.85, height = 0.85) +
  # X-axis label
  draw_label("Physical distance (km)", x = 0.55, y = 0.02, vjust = 0, size = 14) +
  # Y-axis label
  draw_label("Bruvo genetic distance", x = 0.02, y = 0.55, vjust = 0.5, angle = 90, size = 14)


ggsave("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Figures\\msats_fig4.2.tiff", fig4_labeled, dpi=600, units="cm", width = 8.5, height = 17)














####################################################################################################################################################################################
## Relatedness
####################################################################################################################################################################################

########## Data formatting guidelines ##########

# 1. It should be a text file (not an Excel file);
# 2. It should be space- or tab-delimited;
# 3. Missing data must be represented by zeros (0) or NA; and
# 4. There should not be a row of column names in the genotype file.
# 5. If you want to analyze relatedness values based on pre-defined groups (e.g. compare relatedness within versus among groups), then the first two (2) characters of each individual ID should represent the labels for each group

########## Summary output ##########


##### Load in data

# Allie's computer - final
# Note - for Acer, w read in Acer_g2_f1 to retain G21 in individual pairwise comparisons, but 'Un' pop (just G21) gets dropped from analyses at sub-population level
#input <- readgenotypedata("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\Apal_g1_f1_related.txt")
input <- readgenotypedata("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\Acer_g2_f1_related.txt")

# Allie's computer - other 
#input <- readgenotypedata ("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\Acer_g1_f1_related.txt")
#input <-  readgenotypedata("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\sens_test_related.txt")
#input <- readgenotypedata("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Fragment Analysis Data!\\Allele Tables\\For analysis\\Apal_g1_f3_related.txt")


# Scott's computer
# input <- readgenotypedata("D:\\Science\\Allie\\Apal_g1_f1_related.txt")
#input <- readgenotypedata("D:\\Science\\Allie\\Acer_g2_f1_related.txt")
#input <- readgenotypedata("D:\\Science\\Allie\\sens_test_related.txt")
#input <- readgenotypedata("D:\\Science\\Allie\\Apal_g1_f3_related.txt")

##### Calculate all relatedness estimators

# Calculate all estimators at once (if want confidence intervals, should run dyad and trioml separately)
related_est <- coancestry(input$gdata, lynchli =1, lynchrd =1, quellergt =1, ritland =1, wang =1, dyadml = 1, trioml = 1)

# Create data frame out of relatedness output
related_est_df <- as.data.frame(related_est$relatedness)

# Create cleaned up data frame of relatedness outputs to save
related_est_df_cleaned <- related_est_df %>% 
    separate(ind1.id, into = c("Population_1", "Species_1", "Genotype_1"), sep = "_") %>% 
    separate(ind2.id, into = c("Population_2", "Species_2", "Genotype_2"), sep = "_")

#write.csv(related_est_df_cleaned, file = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\related_est_df_Apal_g1_f1.csv", row.names = F)
#write.csv(related_est_df_cleaned, file = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\related_est_df_Acer_g2_f1.csv", row.names = F)

#write.csv(related_est_df_cleaned, file = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\related_est_df_Acer_g1_f1.csv", row.names = F)
#write.csv(related_est_df_cleaned, file = "C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data files\\Genetics\\Downstream_created\\related_est_df_Apal_f3.csv", row.names = F)


# Calculate overall means of each estimator
related_est_overall <- related_est_df %>% 
    summarise(trioml_avg = mean(trioml),
              wang_avg = mean(wang),
              lynchli_avg = mean(lynchli),
              lynchrd_avg = mean(lynchrd),
              ritland_avg = mean(ritland),
              quellergt_avg = mean(quellergt),
              dyadml_avg = mean(dyadml)) %>% 
    mutate(group = "overall")

# Calculate means per subpopulation of each estimator and join in the overall averages calculated above
related_est_summary <- related_est_df %>% 
    filter(group %in% c("CeCe", "KlKl", "LaLa", "NoNo", "SoSo")) %>%   # Only want the within-group relatedness values to average. This also removes G21 from 'Un' (unknown) location
    group_by(group) %>% 
    summarise(trioml_avg = mean(trioml),
              wang_avg = mean(wang),
              lynchli_avg = mean(lynchli),
              lynchrd_avg = mean(lynchrd),
              ritland_avg = mean(ritland),
              quellergt_avg = mean(quellergt),
              dyadml_avg = mean(dyadml)) %>% 
    full_join(related_est_overall)



# Save for focal coral species

# Apal
#write.csv(related_est_summary, "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Apal_g1_f1.csv", row.names = F)

# Acer
#write.csv(related_est_summary, "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Acer_g1_f1.csv", row.names = F)
#write.csv(related_est_summary, "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Acer_g2_f1.csv", row.names = F)

# Apal f3
#write.csv(related_est_summary, "C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Apal_g1_f3.csv", row.names = F)


########## Relatedness heatmap ##########

##### Prepare data

# Pull all individual genotypes per species
all_inds <- mixedsort(unique(c(related_est_df_cleaned$Genotype_1, 
                               related_est_df_cleaned$Genotype_2)))

# wrangle data
relatedness_to_plot <- related_est_df_cleaned %>%
  select(Genotype_1, Genotype_2, trioml) %>%
  rbind(setNames(.[, c("Genotype_2", "Genotype_1", "trioml")],                     # make sure every pair has a value in both direction
                 c("Genotype_1", "Genotype_2", "trioml"))) %>%
  rbind(data.frame(Genotype_1 = all_inds, Genotype_2 = all_inds, trioml = 1)) %>%  # add diagonal
  filter(match(Genotype_2, all_inds) >= match(Genotype_1, all_inds)) %>%           # keep only lower triangle
  mutate(Genotype_1 = factor(Genotype_1, levels = all_inds),                       # mutate as factor to make sure ggplot orders them sequentially
         Genotype_2 = factor(Genotype_2, levels = all_inds))


##### Plot

# set title
#set_title <- "Acropora palmata"
set_title <- "Acropora cervicornis"

# Red-blue heatmap
ggplot(relatedness_to_plot , aes(x = Genotype_2, y = Genotype_1, fill = trioml)) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = round(trioml, 2)), size = 2.5) +
  geom_hline(yintercept = seq(0.5, length(all_inds) + 0.5, 1), color = "gray", linewidth = 0.5) +
  geom_vline(xintercept = seq(0.5, length(all_inds) + 0.5, 1), color = "gray", linewidth = 0.5) +
  scale_fill_gradient(low = "blue", high = "red", na.value = "transparent") +
  scale_x_discrete(limits = mixedsort(unique(relatedness_to_plot $Genotype_1))) +
  scale_y_discrete(limits = mixedsort(unique(relatedness_to_plot $Genotype_2))) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(fill = "Relatedness\n(trioml)")+
  ggtitle(set_title)


# Stoplight heatmap
ggplot(relatedness_to_plot , aes(x = Genotype_2, y = Genotype_1,
                                 fill = cut(trioml, breaks = c(-Inf, 0.125, 0.5, Inf), 
                                            labels = c("< 0.125", "0.125 - 0.5", "> 0.5")))) +
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = round(trioml, 2)), size = 2.5) +
  geom_hline(yintercept = seq(0.5, length(all_inds) + 0.5, 1), color = "gray", linewidth = 0.5) +
  geom_vline(xintercept = seq(0.5, length(all_inds) + 0.5, 1), color = "gray", linewidth = 0.5) +
  scale_fill_manual(values = c("< 0.125" = "darkgreen", 
                               "0.125 - 0.5" = "orange", 
                               "> 0.5" = "red"),
                    na.value = "transparent") +
  scale_x_discrete(limits = mixedsort(unique(relatedness_to_plot $Genotype_1))) +
  scale_y_discrete(limits = mixedsort(unique(relatedness_to_plot $Genotype_2))) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(fill = "Relatedness\n(trioml)")+
  ggtitle(set_title)





########## Compare expected vs observed relatedness values across subgroups ##########
# NOTE: this can take a long time. For trioml with 1000 iterations, it took about 1-1.5 hours on Scott's computer

##### Load in data
# input <- load-ins above

##### set working directory for R to save expected and observed reltedness values as csvs. Make sure to move and rename after so that they don't get re-written!
setwd("C:\\Users\\rassweiler\\Documents\\PhD Research!\\Bonaire\\Data analyses\\Genetics\\relatedness\\related and coancestry\\grouprel outputs")
#setwd("D:\\Science\\Allie\\grouprel outputs 260125")

##### Run grouprel()
# only 1 estimator can be used at a time
# this can take a long time because relatedness gets recalculated for each iteration
# Apal_g1_f1 duration: 5.5 seconds per iteration. Run 9:20am - 10:50
# Acer g2_f1 duration: 1.5 seconds per iteration. Run 10:50am - 11:20am
# Apal g1_f3 duration: 4 seconds per iteration. Run 11:20am - 12:30 (estimated)
par(mfrow = c(2, 3))
grouprel(genotypes = input$gdata, estimatorname="trioml", usedgroups = "all", iterations = 1000)


# CAUTION: grouprel() function is only appropriate if the possible combinations of data greatly outnumbers the number of iterations performed
# Use a website like: 
# https://www.statskingdom.com/combinations-calculator.html#google_vignette
# https://stattrek.com/online-calculator/combinations-permutations
# https://www.calculatorsoup.com/calculators/discretemathematics/combinations.php
# r = subpop sample size (smallest is conservative)
# n = total number genotypes
# output combinations/num groups = total number of non-matching combinations per group possible. 
# APAL r = 4, n = 34, groups = 5. Possible combos = 46376/5 = 9275.2
# ACER r = 2 (minimum), n = 23, groups = 4. Possible combos = 253/4 = 63.25
# ACER r = 5 (median), n = 23, groups = 4. Possible combos = 33649/4 = 8412.25



##### Set up grouprel() output data to make my own version of historgrams produced by grouprel()

# IMPORTANT DATA FORMATTING
# rows = iterations, columns = sub-pops + average across subpops (last column)
# in 'expected', the columns are unlabeled
# The order of groups (i.e., the order of columns) is the same order as they are encountered in the genotype file (and as printed in the \observed-r.csv" file). The last column is the average relatedness value across all groups for each iteration.

# read in Apal_g1_f1 data
# obsv <- read.csv("260125_1000iter_trioml\\observed-r_Apal_g1_f1.csv")
# expec <- read.csv("260125_1000iter_trioml\\expectedrel_Apal_g1_f1.csv")
# obsv_overall <- read.csv("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Apal.csv")[6,2]
# col_renames <- c(obsv$within, "Overall") # pull population names in the order they're listed in 'observed', plus add in "Overall" at end
# colnames(expec)[2:7] <- col_renames  # rename columns in expected based on subpop ordering in observed

# read in Acer_g2_f1 data
# obsv <- read.csv("260125_1000iter_trioml\\observed-r_Acer_g2_f1.csv")
# expec <- read.csv("260125_1000iter_trioml\\expectedrel_Acer_g2_f1.csv")
# obsv_overall <- read.csv("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Acer_g2_f1.csv")[5,2]
# col_renames <- c(obsv$within, "Overall") # pull population names in the order they're listed in 'observed', plus add in "Overall" at end
# colnames(expec)[2:6] <- col_renames  # rename columns in expected based on subpop ordering in observed

# read in Apal_g1_f3 data
obsv <- read.csv("260125_1000iter_trioml\\observed-r_Apal_g1_f3.csv")
expec <- read.csv("260125_1000iter_trioml\\expectedrel_Apal_g1_f3.csv")
obsv_overall <- read.csv("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Tables\\related_est_summary_Apal_g1_f3.csv")[6,2]
col_renames <- c(obsv$within, "Overall") # pull population names in the order they're listed in 'observed', plus add in "Overall" at end
colnames(expec)[2:7] <- col_renames  # rename columns in expected based on subpop ordering in observed



###### Histograms!


# FIRST: must open tiff before creating plot
# Also, set for-loop settings

# Apal initial set up 
# tiff("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\figSX_Apal_rels.tiff", width = 18, height = 11, units = "cm", res = 600) # important for plot to get exported correctly
# set_title <- "Acropora palmata"
# panel_pops <- c("NoNo", "CeCe", "SoSo", "KlKl", "LaLa", "Overall")
# panel_titles <- c("North", "Central", "South", "Klein", "Lac Bay", "Overall")

# Acer initial set up 
# tiff("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\figSX_Acer_rels.tiff", width = 18, height = 11, units = "cm", res = 600) # important for plot to get exported correctly
# set_title <- "Acropora cervicornis"
# panel_pops <- c("NoNo", "CeCe", "SoSo", "KlKl", "Overall")
# panel_titles <- c("North", "Central", "South", "Klein", "Overall")

# Apal f3 initial set up 
tiff("C:\\Users\\rassweiler\\Documents\\Publishing\\Ch1_Microsatellites RRB Acropora\\Supplement\\figSX_Apal_f3_rels.tiff", width = 18, height = 11, units = "cm", res = 600) # important for plot to get exported correctly
set_title <- "Acropora palmata (with 6 Acer loci)"
panel_pops <- c("NoNo", "CeCe", "SoSo", "KlKl", "LaLa", "Overall")
panel_titles <- c("North", "Central", "South", "Klein", "Lac Bay", "Overall")


# Set plot panel framing
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))



# For-loop for Apal
for (i in 1:6){ 
  
  exp_vals <- expec[,panel_pops[i]]
  
  if(i == 6){
    
    obsv_val <- obsv_overall
    
  } else{
    
    obsv_val <- obsv$relvalues[obsv$within == panel_pops[i]]
    
  }
  
  title <- panel_titles[i]
  
  hist(exp_vals,
       xlab = "",
       ylab = "",
       main = title)
  abline(v = obsv_val, col = "red", lwd = 7)
  
}




# For-loop for Acer
for (i in 1:5){ 
  
  exp_vals <- expec[,panel_pops[i]]
  
  if(i == 5){
    
    obsv_val <- obsv_overall
    
  } else{
    
    obsv_val <- obsv$relvalues[obsv$within == panel_pops[i]]
    
  }
  
  title <- panel_titles[i]
  
  hist(exp_vals,
       xlab = "",
       ylab = "",
       main = title)
  abline(v = obsv_val, col = "red", lwd = 7)
  
}





# Add shared axes
mtext("Relatedness", side = 1, line = 2.5, outer = TRUE, cex = 1.2)
mtext("Frequency", side = 2, line = 2.5, outer = TRUE, cex = 1.2)


# Add overall italic title
mtext(bquote(italic(.(set_title))),
      side = 3, line = 1.5, outer = TRUE, cex = 1.5)




# Close device: important for plot to get exported correctly
dev.off()
#












