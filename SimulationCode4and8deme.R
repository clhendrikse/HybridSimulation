

library(diveRsity)
library(adegenet)

#_____________________________4 Demes____________________________________________________________________________________
setwd("/Users/CHendrikse/Documents/fsc26_win64/4Demes/")
#4demes file n=1 --------------------------------------------------------------------------------------------------------

#converting .arp file to a .gen file
arp2gen("4Demes_1_1.arp")

#turning file into a genind object
deme4ind1 <- read.genepop("4Demes_1_1.gen", ncode = 3)
deme4ind1@tab

#checking values
nLoc(deme4ind1)
# ncol(deme4ind1@tab)

individuals4demes1 <- seppop(deme4ind1,res.type="genind")
names(individuals4demes1)

species1_1 <- individuals4demes1[["pop1"]]
species2_1 <- individuals4demes1[["pop2"]]
species3_1 <- individuals4demes1[["pop3"]]
species4_1 <- individuals4demes1[["pop4"]]


#hybridizing pop1 and pop2
hybridsplusparents4demes1 <- hybridize(species1_1, species2_1, pop = "hybrid", n=10)
hybridsplusparents4demes1@pop
#storing all species info into one genind object 4_1 referring to 4 demes n=1
  #note: hybridize does not appear to store parents
final_genind4_1 <- repool(species1_1, species2_1, hybridsplusparents4demes1, species3_1, species4_1)
final_genind4_1@pop
#making function to store to vector
final_genind_df4_1 <- genind2df(final_genind)
actual_values4_1 <- as.vector(final_genind_df$pop)
actual_values4_1

#writing values into text file

is.vector(actual_values4_1)

genind2structure(final_genind4_1, file="parentandhybrid", pops=FALSE)

#4demes n=2------------------------------------------------------------------------------------------------------------------

#converting .arp file to a .gen file
arp2gen("4Demes_1_2.arp")

#turning file into a genind object
deme4ind2 <- read.genepop("4Demes_1_2.gen", ncode = 3)
deme4ind2@tab

#checking values
nLoc(deme4ind2)
# ncol(deme4ind1@tab)

individuals4demes2 <- seppop(deme4ind2,res.type="genind")
names(individuals4demes2)

species1_2 <- individuals4demes2[["pop1"]]
species2_2 <- individuals4demes2[["pop2"]]
species3_2 <- individuals4demes2[["pop3"]]
species4_2 <- individuals4demes2[["pop4"]]


#hybridizing pop1 and pop2
hybridsplusparents4demes2 <- hybridize(species1_2, species2_2, pop = "hybrid", n=10)
hybridsplusparents4demes2@pop
#storing all species info into one genind object
#note: hybridize does not appear to store parents
final_genind4_2 <- repool(species1_2, species2_2, hybridsplusparents4demes2, species3_2, species4_2)
final_genind4_2@pop
#making function to store to vector
final_genind_df4_2 <- genind2df(final_genind4_2)
actual_values4_2 <- as.vector(final_genind_df4_2$pop)
actual_values4_2

#writing values into text file

is.vector(actual_values4_2)

genind2structure(final_genind4_2, file="parentandhybrid4_2", pops=FALSE)



#____________________________________________8 Demes___________________________________________________________________
setwd("/Users/CHendrikse/Documents/fsc26_win64/8Demes/")
#8Demes n=1------------------------------------------------------------------------------------------------------------

#converting .arp file to a .gen file
arp2gen("8Demes_1_1.arp")

#turning file into a genind object
deme8ind1 <- read.genepop("8Demes_1_1.gen", ncode = 3)
deme8ind1@tab

#checking values
nLoc(deme8ind1)
# ncol(deme4ind1@tab)

individuals8demes1 <- seppop(deme8ind1,res.type="genind")
names(individuals8demes1)

species8deme1_1 <- individuals8demes1[["pop1"]]
species8deme2_1 <- individuals8demes1[["pop2"]]
species8deme3_1 <- individuals8demes1[["pop3"]]
species8deme4_1 <- individuals8demes1[["pop4"]]
species8deme5_1 <- individuals8demes1[["pop5"]]
species8deme6_1 <- individuals8demes1[["pop6"]]
species8deme7_1 <- individuals8demes1[["pop7"]]
species8deme8_1 <- individuals8demes1[["pop8"]]

#hybridizing pop1 and pop2
hybridsplusparents8demes1 <- hybridize(species8deme1_1, species8deme2_1, pop = "hybrid", n=10)
hybridsplusparents8demes1@pop
#storing all species info into one genind object
#note: hybridize does not appear to store parents
final_genind8_1 <- repool(species8deme1_1, species8deme2_1, hybridsplusparents8demes1, 
                          species8deme3_1, species8deme4_1, species8deme5_1, species8deme6_1,
                          species8deme7_1, species8deme8_1)
final_genind8_1@pop
#making function to store to vector
final_genind_df8_1 <- genind2df(final_genind8_1)
actual_values8_1 <- as.vector(final_genind_df8_1$pop)
actual_values8_1
is.vector(actual_values8_1)

#convert to structure file
genind2structure(final_genind8_1, file="parentandhybrid8_1", pops=FALSE)



#8Demes n=2------------------------------------------------------------------------------------------------------------

#converting .arp file to a .gen file
arp2gen("8Demes_1_2.arp")

#turning file into a genind object
deme8ind2 <- read.genepop("8Demes_1_2.gen", ncode = 3)
deme8ind2@tab

#checking values
nLoc(deme8ind2)
# ncol(deme4ind1@tab)

individuals8demes2 <- seppop(deme8ind2,res.type="genind")
names(individuals8demes2)

species8deme1_2 <- individuals8demes2[["pop1"]]
species8deme2_2 <- individuals8demes2[["pop2"]]
species8deme3_2 <- individuals8demes2[["pop3"]]
species8deme4_2 <- individuals8demes2[["pop4"]]
species8deme5_2 <- individuals8demes2[["pop5"]]
species8deme6_2 <- individuals8demes2[["pop6"]]
species8deme7_2 <- individuals8demes2[["pop7"]]
species8deme8_2 <- individuals8demes2[["pop8"]]

#hybridizing pop1 and pop2
hybridsplusparents8demes2 <- hybridize(species8deme1_2, species8deme2_2, pop = "hybrid", n=10)
hybridsplusparents8demes2@pop
#storing all species info into one genind object
#note: hybridize does not appear to store parents
final_genind8_2 <- repool(species8deme1_2, species8deme2_2, hybridsplusparents8demes2, 
                          species8deme3_2, species8deme4_2, species8deme5_2, species8deme6_2,
                          species8deme7_2, species8deme8_2)
final_genind8_2@pop
#making function to store to vector
final_genind_df8_2 <- genind2df(final_genind8_2)
actual_values8_2 <- as.vector(final_genind_df8_2$pop)
actual_values8_2
is.vector(actual_values8_2)

#convert to structure file
genind2structure(final_genind8_2, file="parentandhybrid8_2", pops=FALSE)





#genind to structure file code (citation in my drive)------------------------------------
genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- adegenet::nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:adegenet::nPop(obj)
    names(popnums) <- as.character(unique(adegenet::pop(obj)))
    popcol <- rep(popnums[as.character(adegenet::pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- adegenet::locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=adegenet::nLoc(obj),
                           dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste("^", L, "\\.", sep=""), 
                              dimnames(obj@tab)[[2]]), 
                        drop = FALSE] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}

read.csv2("4deme/parentandhybridset/Results/parentandhybridset_run_1_f", skip = 55)
