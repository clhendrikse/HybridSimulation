remotes::install_github('royfrancis/pophelper')
library(remotes)
library(diveRsity)
library(adegenet)
library(pophelper)
library(sjmisc)
library(dplyr)
library(xlsx)



#functions ---------------------------------------------------------------------
# Function converting .arp file to a .gen file 
makeStructure <- function(filename){
  
  addarp <- paste(arpfilename, ".arp", sep="")
  arp2gen(addarp)
  
  #turning file into a genind object
  addgen<- paste(arpfilename, ".gen", sep="")
  genindobj <- read.genepop(addgen, ncode = 3)
  
  individuals <- seppop(genindobj,res.type="genind")
  
  
  #hybridizing pop1 and pop2
  hybrids <- hybridize(individuals[["pop1"]], individuals[["pop2"]], pop = "hybrid", n=10)
  
  #storing all species info into one genind object 4_1 referring to 4 demes n=1
  #note: hybridize does not appear to store parents
  final_genind <- repool(individuals[["pop1"]], individuals[["pop2"]], hybrids)
  a <- 3
  while(a <= length(names(individuals))){
    popname <- paste("pop", a, sep = "")
    final_genind <- repool(final_genind, individuals[[popname]])
    
    a <- a+1
  }
  
  final_genind@pop
  #making function to store to vector
  final_genind_df <- genind2df(final_genind)
  actual_values <- as.vector(final_genind_df$pop)
  #writing values into structure file
  genind2structure(final_genind, file="parentandhybrid", pops=FALSE)
  
  return(actual_values)
}

#genind to structure file code (citation in my drive)
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
} #put pops = true need to call it differently


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

#function adding all true values
AddingTrueGroup <- function(tabledecisions){
  a <- 1
  b<-a*10
  c <- b-9
  #setting base value for the vector
  vect <- "unknown"
 
  while(a <= length(tabledecisions$Cluster1)){
     decisionvect <- tabledecisions$Decision[c:b]
     

    if(length(str_find(decisionvect,tabledecisions$Decision[a])) >= 5){
     vect[c:b] <- tabledecisions$Decision[a] 
    }
     
     
    if(a%%10 == 0){
      b <- b+10
      c<- c+10
    }
    a <- a+1
  }
  tabledecisions$ActualVals <- vect
  return(tabledecisions)
}

#comparing groups and giving either TRUE or FALSE values depending on if STRUCTURE was correct in the individuals placement
CompareQ <- function(tabletocompare, originalLabels){
  bookmark <- 1
  #creating initial comparison object
  if(tabletocompare$Decision[1] %in% tabletocompare$ActualVals[1]){
    bookmark <- bookmark +1
    comparison <- c(TRUE)
  }else{
    comparison <- c(FALSE)
    bookmark <- bookmark + 1
  }
  
  #looping the table to make each comparison value
  while(bookmark <= length(tabletocompare$Cluster1)){
    if(tabletocompare$Decision[bookmark] %in% tabletocompare$ActualVals[bookmark]){
      bookmark <- bookmark +1
      comparison <- c(comparison, TRUE)
    }else{
      comparison <- c(comparison, FALSE)
      bookmark <- bookmark + 1
    }
  }
  
  tabletocompare$Comparison <- comparison
  newtable <- tabletocompare
  newtable$Orignial <- originalLabels
  return(newtable)
  print(newtable)
}

#Function calculating proportion of successes
ProportionOfSuccess <- function(finaltable){
  
  
  hybridinds <- str_find(finaltable$ActualVals, "hybrid")
  pureinds <- str_find(finaltable$ActualVals, "pop")
  
  numsuccess <- sum(finaltable$Comparison)
  numhybridsuccess <- sum(finaltable$Comparison[hybridinds])
  numpuresuccess <-sum(finaltable$Comparison[pureinds])
  
  percentsuccess <- numsuccess/length(finaltable$Comparison)
  percenthybridsuccess <- numhybridsuccess/length(hybridinds)
  percentPureSuccess <- numpuresuccess/length(pureinds)
  
  propsuccessvect <- paste("Total Proportion of Success: ", percentsuccess, sep = "")
  propsuccesshybrid <- paste("Hybrid Proportion of Success: ", percenthybridsuccess, sep = "")
  propsuccesspure <- paste("Pure Proportion of Success: ", percentPureSuccess, sep = "")
  
  print(propsuccessvect)
  print(propsuccesspure)
  print(propsuccesshybrid)
  
  #write.xlsx(finaltable,file = "FinalTable")

}

maketable <- function(demeQmat, results, originallabels){
  #after genind2Structure
  str(demeQmat) 
  clusterdf <- demeQmat
  tabledecisions <- structuredecision(results)
  tablewithboth <- AddingTrueGroup(tabledecisions)
  completetable <- CompareQ(tablewithboth, originallabels)
  print(completetable)
  SuccessProp <- ProportionOfSuccess(completetable)
  SuccessProp
}


#Pipeline-----------------------------------------------------------------------



setwd("/Users/CHendrikse/Documents/fsc26_win64")
setwd("/Users/clhen/Documents/HybridSimulation/SimParFiles/8DemeFixed/")
arpfilename<- "8DemeFixed_1_1" #without .arp


labels <-makeStructure(arpfilename)

demeQmat <- readQ("/Users/clhen/Documents/HybridSimulation/Structure/Outputs/parentandhybrid8Deme_8_1_f")  
demeQmat <- readQ("/Users/clhen/Documents/Internship/results4deme1.txt") 
results <- demeQmat$parentandhybrid8Deme_8_1_f
hold <- maketable(demeQmat, results, labels)
