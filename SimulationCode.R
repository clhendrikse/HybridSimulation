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
  
  #takes the file name and adds .arp so it can be found in the file
  addarp <- paste(arpfilename, ".arp", sep="")
  arp2gen(addarp)
  
  #turning file into a genind object
  addgen<- paste(arpfilename, ".gen", sep="")
  genindobj <- read.genepop(addgen, ncode = 3)
  
  #taking the genpop object and separating it into genind objects
  individuals <- seppop(genindobj,res.type="genind")
  
  
  #hybridizing pop1 and pop2
  hybrids <- hybridize(individuals[["pop1"]], individuals[["pop2"]], pop = "hybrid", n=10)
  
  #storing all species info into one genind object 4_1 referring to 4 demes n=1
  #note: hybridize does not store parents
  #repooling the first two species and the hybrids
  final_genind <- repool(individuals[["pop1"]], individuals[["pop2"]], hybrids)
  
  #adding in the other species individuals after the hybrids (sp1,sp2,hybrids,sp3,sp4...)
  a <- 3
  while(a <= length(names(individuals))){
    popname <- paste("pop", a, sep = "")
    final_genind <- repool(final_genind, individuals[[popname]])
    
    a <- a+1
  }
  
  #making function to store to vector
  final_genind_df <- genind2df(final_genind)
  actual_values <- as.vector(final_genind_df$pop)
  
  
  #writing values into structure file
  genind2structure(final_genind, file="parentandhybrid", pops=FALSE)
  
  #return the labels for individuals in the order they are in the file
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
 #write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
} #put pops = true need to call it differently


#function that takes the Q values and creates a new column that labels STRUCTURE's decision changing the value in this function should allow the other functions to still run
structuredecision <- function(obj){
  
  #putting individuals in groups depending on what STRUCTURE decides (max is 10 clusters)
  #pure
  pure1 <- rownames(obj[which(obj$Cluster1 > 0.9),])  
  pure2 <- rownames(obj[which(obj$Cluster2 > 0.9),])
  pure3 <- rownames(obj[which(obj$Cluster3 > 0.9),])
  pure4 <- rownames(obj[which(obj$Cluster4 > 0.9),])
  pure5 <- rownames(obj[which(obj$Cluster5 > 0.9),])  
  pure6 <- rownames(obj[which(obj$Cluster6 > 0.9),])
  pure7 <- rownames(obj[which(obj$Cluster7 > 0.9),])
  pure8 <- rownames(obj[which(obj$Cluster8 > 0.9),])
  pure9 <- rownames(obj[which(obj$Cluster9 > 0.9),])
  pure10 <- rownames(obj[which(obj$Cluster10 > 0.9),])
  #hybrids
  hybrids1 <- rownames(obj[which(obj$Cluster1 > 0.4 & obj$Cluster1 < .6),])
  hybrids2 <- rownames(obj[which(obj$Cluster2 > 0.4 & obj$Cluster2 < .6),])
  hybrids3 <- rownames(obj[which(obj$Cluster3 > 0.4 & obj$Cluster3 < .6),])
  hybrids4 <- rownames(obj[which(obj$Cluster4 > 0.4 & obj$Cluster4 < .6),])
  hybrids5 <- rownames(obj[which(obj$Cluster5 > 0.4 & obj$Cluster5 < .6),])
  hybrids6 <- rownames(obj[which(obj$Cluster6 > 0.4 & obj$Cluster6 < .6),])
  hybrids7 <- rownames(obj[which(obj$Cluster7 > 0.4 & obj$Cluster7 < .6),])
  hybrids8 <- rownames(obj[which(obj$Cluster8 > 0.4 & obj$Cluster8 < .6),])
  hybrids9 <- rownames(obj[which(obj$Cluster7 > 0.4 & obj$Cluster7 < .6),])
  hybrids10 <- rownames(obj[which(obj$Cluster8 > 0.4 & obj$Cluster8 < .6),])
  hybridinds <- c(hybrids1,hybrids2,hybrids3,hybrids4,hybrids5,hybrids6,hybrids7,hybrids8,hybrids9,hybrids10)
  hybridinds <- unique(hybridinds)
  guesses <- c(pure1, pure2, pure3, pure4, hybridinds, pure5, pure6, pure7, pure8, pure9, pure10)
  #for 8 demes

  #creating initial decisions object of what cluster structure decided the individual belongs to
  v <- 1
  decision <- as.numeric(0)
  if(row.names(obj)[1] %in% pure1) {
    decision <- "cluster1"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure2) {
    decision <- "cluster2"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure3) {
    decision <- "cluster3"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure4) {
    decision <- "cluster4"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure5) {
    decision <- "cluster5"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure6) {
    decision <- "cluster6"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure7) {
    decision <- "cluster7"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure8) {
    decision <- "cluster8"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure9) {
    decision <- "cluster9"
    v <- v+1
  }
  if(row.names(obj)[1] %in% pure10) {
    decision <- "cluster10"
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
  

  #while loop to add in the entire string of decisions of what cluster STRUCTURE grouped them in
  while(v <= length(obj$Cluster1)){
    if(row.names(obj)[v] %in% pure1) {
      decision <- c(decision, "cluster1")
      
    }
    if(row.names(obj)[v] %in% pure2) {
      decision <- c(decision,"cluster2")
      
    }
    if(row.names(obj)[v] %in% pure3) {
      decision <- c(decision,"cluster3")
 
    }
    if(row.names(obj)[v] %in% pure4) {
      decision <- c(decision,"cluster4")

    }
    if(row.names(obj)[v] %in% pure5) {
      decision <- c(decision,"cluster5")

    }
    if(row.names(obj)[v] %in% pure6) {
      decision <- c(decision,"cluster6")

    }
    if(row.names(obj)[v] %in% pure7) {
      decision <- c(decision,"cluster7")
      values <- c(values,v)
      
    }
    if(row.names(obj)[v] %in% pure8) {
      decision <- c(decision,"cluster8")

    }
    if(row.names(obj)[v] %in% pure9) {
      decision <- c(decision,"cluster9")

    }
    if(row.names(obj)[v] %in% pure10) {
      decision <- c(decision,"cluster10")

    }
    if(row.names(obj)[v] %in% hybridinds) {
      decision <- c(decision,"hybrid")
 
    }
    if(row.names(obj)[v] %nin% guesses) {
      decision <- c(decision, "unknown")

    }
    v <- v+1
  }
  
  #adding in the Decision column
  obj$Decision <- decision
  
  #saving entire table into new object and returning it
  tablewithdecision <- obj
  return(tablewithdecision)
  
}

#function adding all true values
#clusters based on which cluster is most commonly used for that group in STRUCTURE, may be inaccurate if STRUCTURE does not run properly
AddingTrueGroup <- function(tabledecisions){
  
  #adding initial values for while loop a is the individual number and b/c are the threshold for the number of individuals from each species
  a <- 1
  b<-a*10
  c <- b-9
  #setting base value for the vector
  vect <- "unknown"
 
  #while loop looking for the number of times the individual matches clusters with other individuals in the same species
  #if the cluster is used 5 or more times, it is assumed to be the correct cluster
  
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
  #putting the vector of actual clusters into a new column
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
  
  #adding comparison to new column in the table
  tablewithboth$Comparison <- comparison
  #putting the table into a new object
  newtable <- tablewithboth
  #adding original species into the table 
  newtable$Orignial <- labels
  
  #printing and returning the table
  print(newtable)
  return(newtable)
}

#Function calculating total proportion of successes, proportion of success for pure species, and proportion of success for hybrids
#will not run properly if structure does not use the correct number of clusters
ProportionOfSuccess <- function(finaltable,numberIndsPerSpecies){
  
  #finding the hybrid and pure individuals
  hybridinds <- str_find(finaltable$ActualVals, "hybrid")
  pureinds <- str_find(finaltable$ActualVals, "cluster")
  
  #finding the number of times the decision matches the actual cluster
  numsuccess <- sum(finaltable$Comparison)
  numhybridsuccess <- sum(finaltable$Comparison[hybridinds])
  numpuresuccess <-sum(finaltable$Comparison[pureinds])
  
  #finding the proportion of success for hybrids, pure, and all
  percentsuccess <- numsuccess/length(finaltable$Comparison)
  percenthybridsuccess <- numhybridsuccess/length(hybridinds)
  percentPureSuccess <- numpuresuccess/length(pureinds)
  
  #vectorizing the outputs
  propsuccessvect <- paste("Total Proportion of Success: ", percentsuccess, sep = "")
  propsuccesshybrid <- paste("Hybrid Proportion of Success: ", percenthybridsuccess, sep = "")
  propsuccesspure <- paste("Pure Proportion of Success: ", percentPureSuccess, sep = "")
  
  #printing the results with a warning that the way structure runs may cause an inaccurate result
  print(propsuccessvect)
  print(propsuccesspure)
  print(propsuccesshybrid)
  warning("Correct values based on Structure")
  
  #taking the number of individuals per species and finding the result for each group if two or more groups are put into the same cluster, a warning will print
  a<- numberIndsPerSpecies
  clusters <- as.numeric(0)
  while(a <= length(finaltable$Cluster1)){
    clusters <- c(clusters, finaltable$Decision[a])
    a <- a+numberIndsPerSpecies
    uniqueclusters <- unique(clusters)
  }
  if(length(uniqueclusters) != length(clusters))
  {
    warning("Number of clusters in decision does not match number of species")
  }
  

}

#master function to create the table and calculate proportion of success
maketable <- function(demeQmat, results, originallabels,numberIndsPerSpecies){
  #after genind2Structure
  str(demeQmat) 
  clusterdf <- demeQmat
  tabledecisions <- structuredecision(results)
  #outputs a warning if there are more than 10 species
  if(ncol(tabledecisions) >11){
    warning("Functions only work with 10 species or less")
    return(0)
  }
  #adding the true clusters
  tablewithboth <- AddingTrueGroup(tabledecisions)
  #adding comparison (TRUE/FALSE) values and the original labels of the species they belong to
  completetable <- CompareQ(tablewithboth, labels)
  #prints the table
  print(completetable)
  #prints the proportion of success for the pure species, the hybrids, and all species
  SuccessProp <- ProportionOfSuccess(completetable, numberIndsPerSpecies)
  SuccessProp
}


#Pipeline-----------------------------------------------------------------------



#set the working directory and list the file name
setwd("/Users/CHendrikse/Documents/fsc26_win64/8DemesLessMigration/")
arpfilename<- "8DemesLessMigration_1_1" #without .arp

#make the structure file and store the original species groups
labels <-makeStructure(arpfilename)

#reads structure's results and makes the table and outputs the proportion of success
demeQmat <- readQ("/Users/CHendrikse/Documents/HybridSimulation/Structure/Outputs/parentandhybrid8DemeFixed_8_1_f")  
demeQmat <- readQ("/Users/clhen/Documents/Internship/results4deme1.txt") 
results <- demeQmat$parentandhybrid8DemeFixed_8_1_f
  #input the number of individuals for the species so the CompareQ function can compare the number of clusters found to the number of clusters there should be
numberIndsPerSpecies <- 10
hold <- maketable(demeQmat, results, labels, numberIndsPerSpecies)

