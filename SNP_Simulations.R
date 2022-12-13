# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GENERATE SNP PARAMETERS AND RUN SIMULATIONS IN STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files for use in the 
# Hendrikse Hybridization project, which analyzes how well the population clustering software STRUCTURE
# does at detecting species level differences in complex hybridization scenarios.

# It generates 2 different scenarios, both using DNA (SNP) markers: one with 4 species, and one with 8 species.
# The Arlequin outputs from fastSimcoal2 are converted to genind files, and then files are cleaned up in order
# to match the required file/folder structure for the project.

# THIS SCRIPT SHOULD ONLY NEED TO BE RUN ONCE! After running, genind objects will be created, which can be 
# used for all downstream steps.

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# Set working directory to the RAID1 partition on the server, to allow for enough space for simulation outputs
# This folder is linked to the home folder
# sim.wd <- "~/Documents/hendrickseHybridization/SNP_Demo/"
# sim.wd <- "/home/chendrikse/Shared/HybridSimulation/"
gitHub.wd <- "/Users/clhen/Documents/HybridSimulation/"
setwd(gitHub.wd)

# ---- FUNCTIONS ----
# Function converting Arlequin output to a single genind object (through gtypes format)
strataG_arp2gen <- function(params, repNumber){
  # Extract marker type from params argument
  marker <- params$settings$genetics$fsc.type
  # Read in the Arlequin file, convert it to a gtype object, then to a genind object
  arp <- fscReadArp(params, sim=c(1,repNumber), marker = marker)
  gtype <- df2gtypes(arp, ploidy = 2)
  genind <- gtypes2genind(gtype)
  # In the "other" slot of the genind object, pass the name of the simulation scenario, and return
  genind@other <- list(params$label)
  return(genind)
}

# Function for converting all of the Arlequin files in a directory to genind, generating a list of genind objects
convertAllArp <- function(arp.path, params){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing simulation outputs
  setwd(arp.path)
  # Create an empty list object to receive list of genind.
  # The length of this list is the number of replicates, which is specified as a numeric vector
  genind.list <- vector("list",length=length(dir()[str_detect(dir(), pattern = ".arp")]))
  fscReps <- seq(1, length(genind.list))
  # Move up one directory, in order for the fscReadArp command (within strataG_arp2gen) to work
  setwd("..")
  # Convert all Arlequin files to a list of genind objects
  for(i in 1:length(genind.list)){
    genind.obj <- strataG_arp2gen(params, rep=i)
    genind.list[[i]] <- genind.obj
  }
  # Reset to original working directory, and return a list of genind objects
  setwd(original.wd)
  return(genind.list)
}

# ---- PARAMETER VALUES ----
# Specify number of simulation replicates
num_reps <- 5
fscVersion <- "fsc2709"

# DEMES
# Specify number of total individuals, for all simulations. 
# This value is 10,000 individuals, leading to 20,000 haplotypes
nInd <- 10000

# Specify number of demes (either 4 or 8)
deme <- fscDeme(deme.size = nInd, sample.size = 10)
demes4 <- fscSettingsDemes(deme, deme, deme, deme)
demes8 <- fscSettingsDemes(deme, deme, deme, deme, deme, deme, deme, deme)

# MIGRATION: No migration between demes, so migration parameters can just be excluded

# HISTORICAL EVENTS
# Number of historical events = number of demes - 1 (all other demes merging into origin deme)
hist.event1 <- fscEvent(event.time = 100000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 0)
hist.event2 <- fscEvent(event.time = 100000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 0)
hist.event3 <- fscEvent(event.time = 100000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 0)
histEvent_4sp <- fscSettingsEvents(hist.event1,hist.event2,hist.event3)
hist.event4 <- fscEvent(event.time = 100000, source = 4, sink = 0, prop.migrants = 1, migr.mat = 0)
hist.event5 <- fscEvent(event.time = 100000, source = 5, sink = 0, prop.migrants = 1, migr.mat = 0)
hist.event6 <- fscEvent(event.time = 100000, source = 6, sink = 0, prop.migrants = 1, migr.mat = 0)
hist.event7 <- fscEvent(event.time = 100000, source = 7, sink = 0, prop.migrants = 1, migr.mat = 0)
histEvent_8sp <- fscSettingsEvents(hist.event1,hist.event2,hist.event3,
                                   hist.event4,hist.event5,hist.event6,hist.event7)

# GENETIC PARAMETERS
# Medium low mutation rate, increased sequence length/blocks/chromosomes
dna_mutRate <- 1e-7
dna <- fscBlock_dna(sequence.length = 500, mut.rate = dna_mutRate)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, dna, dna, dna, dna, dna, dna, 
                                   num.chrom = 100)

# ---- WRITE PARAMETER FILES AND RUN SIMULATIONS ----
# Navigate to the folder containing simulation outputs
sim.wd <- paste0(gitHub.wd,"Scenarios/")
setwd(sim.wd)
# Declare scenario name variables, which will be used for downstream file processing
scenario_DNA_4Deme <- "DNA_4sp_10ns"
scenario_DNA_8Deme <- "DNA_8sp_10ns"
# Write parameter files
# 4 Demes
DNA_4sp_10ns.params <- fscWrite(demes = demes4, 
                                events = histEvent_4sp, genetics = DNAgenetics, 
                                label = scenario_DNA_4Deme, use.wd=TRUE)
# 8 Demes
DNA_8sp_10ns.params <- fscWrite(demes = demes8, 
                                events = histEvent_8sp, genetics = DNAgenetics, 
                                label = scenario_DNA_8Deme, use.wd=TRUE)
# Run simulations
# 4 Demes
DNA_4sp_10ns.params <- fscRun(DNA_4sp_10ns.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# 8 Demes
DNA_8sp_10ns.params <- fscRun(DNA_8sp_10ns.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# ---- CONVERT ARLEQUIN OUTPUTS TO STRATAG ----
# 4 Demes
# Make a list of genind objects
setwd(paste0(sim.wd, scenario_DNA_4Deme))
DNA_4sp_10ns.genList <- convertAllArp(arp.path = paste0(sim.wd, scenario_DNA_4Deme),
                                      params = DNA_4sp_10ns.params)
# Save genind objects to .Rdata objects, to be read in later
for(i in 1:length(DNA_4sp_10ns.genList)){
  saveRDS(object = DNA_4sp_10ns.genList[i], file = paste0(scenario_DNA_4Deme,"_1_",i,"_genind.Rdata"))
}

# 8 Demes
# Make a list of genind objects
setwd(paste0(sim.wd, scenario_DNA_8Deme))
DNA_8sp_10ns.genList <- convertAllArp(arp.path = paste0(sim.wd, scenario_DNA_8Deme),
                                      params = DNA_8sp_10ns.params)
# Save genind objects to .Rdata objects, to be read in later
for(i in 1:length(DNA_8sp_10ns.genList)){
  saveRDS(object = DNA_8sp_10ns.genList[i], file = paste0(scenario_DNA_4Deme,"_1_",i,"_genind.Rdata"))
}

# ---- MOVE FILES AND CLEANUP SIMULATION OUTPUTS ----
# The commands in this section are meant to remove unwanted fastSimcoal2 output files, 
# and move the required output files into folders that will allow them to be processed easily.
setwd(sim.wd)

# Remove .log files generated from fastSimcoal2 runs
file.remove("DNA_4sp_10ns.log", "DNA_8sp_10ns.log")

# Move .par files to the parameter files directory
file.rename(from="DNA_4sp_10ns.par", to=paste0(gitHub.wd,"SimParFiles/DNA/DNA_4sp_10ns.par"))
file.rename(from="DNA_8sp_10ns.par", to=paste0(gitHub.wd,"SimParFiles/DNA/DNA_8sp_10ns.par"))

# Cleanup files/folders
# 4 Demes
setwd(paste0(sim.wd, scenario_DNA_4Deme))
# Create subfolder for Arlequin files
dir.create("ArpFiles")
# Move Arlequin output files to newly created folder
file.copy(from=list.files(pattern = ".arp"), to="ArpFiles")
file.remove(list.files(pattern = ".arp"))
file.copy(from=list.files(pattern = ".Rdata"), to="ArpFiles")
file.remove(list.files(pattern = ".Rdata"))
# Remove .arb and .simparam files
file.remove(list.files(pattern = ".arb")) ; file.remove(list.files(pattern = ".simparam")) 
# Create (empty) subfolder for STRUCTURE input files
dir.create("StructIn")

# 8 Demes
setwd(paste0(sim.wd, scenario_DNA_8Deme))
# Create subfolder for Arlequin files
dir.create("ArpFiles")
# Move Arlequin output files to newly created folder
file.copy(from=list.files(pattern = ".arp"), to="ArpFiles")
file.remove(list.files(pattern = ".arp"))
file.copy(from=list.files(pattern = ".Rdata"), to="ArpFiles")
file.remove(list.files(pattern = ".Rdata"))
# Remove .arb and .simparam files
file.remove(list.files(pattern = ".arb")) ; file.remove(list.files(pattern = ".simparam")) 
# Create (empty) subfolder for STRUCTURE input files
dir.create("StructIn")
