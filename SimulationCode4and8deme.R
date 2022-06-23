

library(diveRsity)
library(adegenet)
library(pophelper)
library(xlsx)
library(sjmisc)

#_____________________________4 Demes____________________________________________________________________________________
setwd("/Users/CHendrikse/Documents/fsc26_win64/4Demes/")
#4demes file n=1 --------------------------------------------------------------------------------------------------------

#converting .arp file to a .gen file
arp2gen("4Demes_1_1.arp")

#turning file into a genind object
deme4ind1 <- read.genepop("4Demes_1_1.gen", ncode = 3)
deme4ind1@tab

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
final_genind_df4_1 <- genind2df(final_genind4_1)
actual_values4_1 <- as.vector(final_genind_df4_1$pop)
actual_values4_1

#writing values into text file

is.vector(actual_values4_1)

genind2structure(final_genind4_1, file="parentandhybrid", pops=FALSE)

#4demes n=2------------------------------------------------------------------------------------------------------------------

setwd("/Users/CHendrikse/Documents/fsc26_win64/4Demes/")
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
#_____________________________8 Demes_____________________________________________________________________________________
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
setwd("/Users/CHendrikse/Documents/fsc26_win64/8Demes/")
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
#----


#___________genind to structure file code (citation in my drive)__________________________________________________
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

  
demeQmat <- readQ("4deme/parentandhybridset/Results4deme1/results4deme1")  
str(demeQmat)  


#function that takes the Q values and creates a new column that labels STRUCTURE's decision changing the value in this function should allow the other functions to still run
structuredecision <- function(obj){
  
  #putting individuals in groups depending on what STRUCTURE decides (base deme=4)
      #pure
      pure1 <- rownames(obj[which(obj$Cluster1 > 0.9),])  
      pure2 <- rownames(obj[which(obj$Cluster2 > 0.9),])
      pure3 <- rownames(obj[which(obj$Cluster3 > 0.9),])
      pure4 <- rownames(obj[which(obj$Cluster4 > 0.9),])
      #hybrids
      hybrids1 <- rownames(obj[which(obj$Cluster1 > 0.4 & obj$Cluster1 < .6),])
      hybrids2 <- rownames(obj[which(obj$Cluster2 > 0.4 & obj$Cluster2 < .6),])
      hybrids3 <- rownames(obj[which(obj$Cluster3 > 0.4 & obj$Cluster3 < .6),])
      hybrids4 <- rownames(obj[which(obj$Cluster4 > 0.4 & obj$Cluster4 < .6),])
      hybridinds <- c(hybrids1, hybrids2, hybrids3, hybrids4)
      if(ncol(obj) == 8){
        hybridinds <- c(hybrids5, hybrids6, hybrids7, hybrids8)  
      }
      hybridinds <- unique(hybridinds)
      guesses <- c(pure1,pure2,pure3,pure4,hybridinds)
#for 8 demes
      
      #creating base value
      pure5<- 0
      pure6<- 0
      pure7<- 0
      pure8<- 0
    if(ncol(obj) == 8){
      pure5 <- rownames(obj[which(obj$Cluster5 > 0.9),])  
      pure6 <- rownames(obj[which(obj$Cluster6 > 0.9),])
      pure7 <- rownames(obj[which(obj$Cluster7 > 0.9),])
      pure8 <- rownames(obj[which(obj$Cluster8 > 0.9),])
      
      hybrids5 <- rownames(obj[which(obj$Cluster5 > 0.4 & obj$Cluster5 < .6),])
      hybrids6 <- rownames(obj[which(obj$Cluster6 > 0.4 & obj$Cluster6 < .6),])
      hybrids7 <- rownames(obj[which(obj$Cluster7 > 0.4 & obj$Cluster7 < .6),])
      hybrids8 <- rownames(obj[which(obj$Cluster8 > 0.4 & obj$Cluster8 < .6),])
      
      
      hybridinds <- c(hybridinds, hybrids5, hybrids6, hybrids7, hybrids8)
      guesses <- c(pure1,pure2,pure3,pure4,hybridinds, pure5, pure6, pure7, pure8)
    }
      v <- 1
      
      #creating initial decisions object
      if(row.names(obj)[1] %in% pure1) {
        decision <- "pop1"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure2) {
        decision <- "pop2"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure3) {
        decision <- "pop3"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure4) {
        decision <- "pop4"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure5) {
        decision <- "pop5"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure6) {
        decision <- "pop6"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure7) {
        decision <- "pop7"
        v <- v+1
      }
      if(row.names(obj)[1] %in% pure8) {
        decision <- "pop8"
        v <- v+1
      }
      if(row.names(obj)[1] %in% hybridinds) {
        decision <- "hybrid"
        v <- v+1
      }
      if(row.names(obj)[1] %nin% guesses) {
        decision <- "unknown"
        v <- v+1
      }
      
      
      
      while(v <= length(obj$Cluster4)){
        if(row.names(obj)[v] %in% pure1) {
          decision <- c(decision, "pop1")
          
        }
        if(row.names(obj)[v] %in% pure2) {
          decision <- c(decision,"pop2")
         
        }
        if(row.names(obj)[v] %in% pure3) {
          decision <- c(decision,"pop3")
         
        }
        if(row.names(obj)[v] %in% pure4) {
          decision <- c(decision,"pop4")
         
        }
        if(row.names(obj)[v] %in% pure5) {
          decision <- c(decision,"pop5")
          
        }
        if(row.names(obj)[v] %in% pure6) {
          decision <- c(decision,"pop6")
          
        }
        if(row.names(obj)[v] %in% pure7) {
          decision <- c(decision,"pop7")
          
        }
        if(row.names(obj)[v] %in% pure8) {
          decision <- c(decision,"pop8")
          
        }
        if(row.names(obj)[v] %in% hybridinds) {
          decision <- c(decision,"hybrid")
          
        }
        if(row.names(obj)[v] %nin% guesses) {
          decision <- c(decision, "unknown")
        
        }
         v <- v+1
      }
        
      length(decision)
        
        obj$Decision <- decision
      
      tablewithdecision <- obj
      return(tablewithdecision)

}
tabledecisions <- structuredecision(demeQmat$results4deme1)
tabledecisions

#function adding all true values
AddingTrueGroup <- function(tablewithdecisions, lablevector){
 tablewithdecisions$ActualVals <- lablevector
  return(tablewithdecisions)
}
tablewithboth <- AddingTrueGroup(tabledecisions,actual_values4_1)

#comparing groups and giving either TRUE or FALSE values depending on if STRUCTURE was correct in the individuals placement
CompareQ <- function(tabletocompare){
     bookmark <- 1
     #creating initial comparison object
     if(tabletocompare$Decision[1] == tabletocompare$ActualVals[1]){
       bookmark <- bookmark +1
       comparison <- c(TRUE)
     }else{
       comparison <- c(FALSE)
       bookmark <- bookmark + 1
     }
     
      #looping the table to make each comparison value
     while(bookmark <= length(tabletocompare$Cluster1)){
           if(tabletocompare$Decision[1] == tabletocompare$ActualVals[1]){
             bookmark <- bookmark +1
             comparison <- c(comparison, TRUE)
           }else{
             comparison <- c(comparison, FALSE)
             bookmark <- bookmark + 1
           }
     }
  
     tabletocompare$Comparison <- comparison
     newtable <- tabletocompare
     return(newtable)
}
completetable <- CompareQ(tablewithboth)

#Function calculating proportion of successes_________________________________________________________
ProportionOfSuccess <- function(finaltable){
  
  numsuccess <- sum(finaltable$Comparison)
  percentsuccess <- numsuccess/length(finaltable$Comparison)
  
  return(percentsuccess)
  
}
SuccessProp <- ProportionOfSuccess(completetable)
SuccessProp
#outputs proportion of successes