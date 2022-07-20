remotes::install_github('royfrancis/pophelper')

library(remotes)
library(diveRsity)
library(adegenet)
library(pophelper)
library(sjmisc)
library(dplyr)

#functions ---------------------------------------------------------------------
# Function converting .arp file to a .gen file 
makeStructure <- function(filename,a, numinds, ScenarioFolder){

  ArpFileLocation <- paste0(getwd(),"/", ScenarioFolder,"ArpFiles/", filename, ".arp")
  #takes the file name and adds .arp so it can be found in the file\
  arp2gen(ArpFileLocation)
  filename <- paste0(getwd(),"/", ScenarioFolder,"ArpFiles/", filename)
  #turning file into a genind object
  addgen<- paste(filename, ".gen", sep="")
  genindobj <- read.genepop(addgen, ncode = 3)
  
  #taking the genpop object and separating it into genind objects
  individuals <- seppop(genindobj,res.type="genind")
  
  #hybridizing pop1 and pop2
  hybrids <- hybridize(individuals[["pop1"]], individuals[["pop2"]], pop = "hybrid", n=numinds)
  
  #storing all species info into one genind object 4_1 referring to 4 demes n=1
  #note: hybridize does not store parents
  #repooling the first two species and the hybrids
  final_genind <- repool(individuals[["pop1"]], individuals[["pop2"]], hybrids)
  final_genind@pop
  #adding in the other species individuals after the hybrids (sp1,sp2,hybrids,sp3,sp4...)
  b <- 3
  while(b <= length(names(individuals))){
    popname <- paste("pop", b, sep = "")
    final_genind <- repool(final_genind, individuals[[popname]])
    
    b <- b+1
  }
  
  #making function to store to vector
  final_genind_df <- genind2df(final_genind)
  actual_values <- as.vector(final_genind_df$pop)
  
  #writing values into structure file

  # filenameforstruct <- paste("parentandhybrid", a, sep="")
  filenameforstruct <- paste0("parentandhybrid", a, ".str")
  genind2structure(final_genind, file= filenameforstruct, pops=FALSE)
  
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
  names(tab)[1] <- ""
 write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
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
CompareQ <- function(tabletocompare, originalLabels, numberIndsPerSpecies){
  

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
  tabletocompare$Comparison <- comparison
  return(tabletocompare)
  }

#Finds the proportion of times STRUCTURE incorrectly labeled individuals from two different species as the same
#Note: if each species has 10 inds, but STRUCTURE labeled 20 inds in the same cluster, all 20 inds are labeled as incorrect
FindMergeError <- function(tabletocompare){
  #taking the number of individuals per species and finding the result for each group if two or more groups are put into the same cluster, a warning will print
  
  #starting point for finding the population clusters
  x<- numberIndsPerSpecies -1
  #Creating initial object for the clusters and number of times the error occurres
  clusters <- as.numeric(0)
  numMergeError <- 0
  
  #finds all the species' cluster decisions and look for repeats using unique clusters
  while(x <= length(tabletocompare$Cluster1)){
    clusters <- c(clusters, tabletocompare$Decision[x])
    x <- x+numberIndsPerSpecies
    uniqueclusters <- unique(clusters)
  }

  #Removes initial "0" value
  clusters <- clusters[-1]
  uniqueclusters <- uniqueclusters[-1]
  #If there are repeats of the clusters find the duplicates and place the error in the comparison column
  if(length(uniqueclusters) != length(clusters)){
    
    #looks for duplicates
    duplicates <- duplicated(clusters)
    #finds position of duplicates
    positionduplicates <- str_find(duplicates, TRUE)
    #finds the name of the clusters that are repeated
    repeatedclusters <- clusters[duplicates]
    
    #initial counter value for while function
    d <- 1
    while (d < length(tabletocompare$Decision)){
      #looks to see if STRUCTURE's decision matches the name of a repeated cluster name
      findcluster <- tabletocompare$Decision[d] == repeatedclusters
      
      #if one of the repeated cluster names are found, the comparison value changes to merge error
      if(str_contains(findcluster, TRUE) == TRUE){
        tabletocompare$Comparison[d] <- "Merge Error"
      }
      
      d <- d+1
    }
    #outputs warning (might remove due to the code now being able to catch the error)
    warning("Number of clusters in decision does not match number of species")
    #counts how many times the merge error occurs
    positionofError <- str_find(tabletocompare$Comparison,"Merge")
    numMergeError <- length(positionofError)
  
  }  
  #finds proportion of individuals with the error and returns the amount
  PropMergeError <- numMergeError/length(tabletocompare$Comparison)
  return(PropMergeError)
    

}

FindP2HError <- function(completetable){
ClusterBookmark <- 1
numP2HError <- 0
while(ClusterBookmark <= length(completetable$Cluster1)){
  Calledhybrid <- str_contains(completetable$Decision[ClusterBookmark], "hybrid" )
  IndisHybrid <- str_contains(completetable$ActualVals[ClusterBookmark], "hybrid")
  if(Calledhybrid == TRUE && IndisHybrid == FALSE){
    completetable$Comparison[ClusterBookmark] <- "P2H"
    numP2HError <- numP2HError + 1
  }
  
  ClusterBookmark <- ClusterBookmark +1
}
  
  propP2HError <- numP2HError/length(completetable$Comparison)
  return(propP2HError)
   
}

FindH2PError <- function(completetable){
  ClusterBookmark <- 1
  numH2PError <- 0
  while(ClusterBookmark <= length(completetable$Cluster1)){
    Calledhybrid <- str_contains(completetable$Decision[ClusterBookmark], "hybrid" )
    IndisHybrid <- str_contains(completetable$ActualVals[ClusterBookmark], "hybrid")
    if(Calledhybrid == FALSE && IndisHybrid == TRUE){
      completetable$Comparison[ClusterBookmark] <- "H2P"
      numH2PError <- numP2HError + 1
    }
    
    ClusterBookmark <- ClusterBookmark +1
  }
  
  propH2PError <- numH2PError/length(completetable$Comparison)
  return(propH2PError)
  
}

FindPUnknownError <- function(completetable){
    ClusterBookmark <- 1
    numPUError <- 0
    while(ClusterBookmark <= length(completetable$Cluster1)){
      CalledUnknown <- str_contains(completetable$Decision[ClusterBookmark], "unknown" )
      IndisPure <- str_contains(completetable$ActualVals[ClusterBookmark], "cluster")
      if(CalledUnknown== TRUE && IndisPure == TRUE){
        completetable$Comparison[ClusterBookmark] <- "P Unknown"
        numPUError <- numPUError + 1
      }
      
      ClusterBookmark <- ClusterBookmark +1
    }
    
    propPUError <- numPUError/length(completetable$Comparison)
    return(propPUError)
    

}

FindHUnknownError <- function(completetable){
  ClusterBookmark <- 1
  numHUError <- 0
  while(ClusterBookmark <= length(completetable$Cluster1)){
    CalledUnknown <- str_contains(completetable$Decision[ClusterBookmark], "unknown" )
    IndisHybrid <- str_contains(completetable$ActualVals[ClusterBookmark], "Hybrid")
    if(CalledUnknown== TRUE && IndisHybrid == TRUE){
      completetable$Comparison[ClusterBookmark] <- "H Unknown"
      numHUError <- numHUError + 1
    }
    
    ClusterBookmark <- ClusterBookmark +1
  }
  
  propHUError <- numHUError/length(completetable$Comparison)
  return(propHUError)
  
  
}

#master function to create the table and calculate proportion of success
maketable <- function(results, originallabels,numberIndsPerSpecies,percenterror){
  #after genind2Structure
  tabledecisions <- structuredecision(results)
  #outputs a warning if there are more than 10 species
  if(ncol(tabledecisions) >11){
    warning("Functions only work with 10 species or less")
  }
  #adding the true clusters
  tabledecisions <- AddingTrueGroup(tabledecisions)
  #adding comparison (TRUE/FALSE) values and the original labels of the species they belong to
  completetable <- CompareQ(tabledecisions, labels, numberIndsPerSpecies)
  #prints the table
  return(completetable)
}


#put this in a while function
#creates array to store proportion of individuals in each run that have the merge error
CreateMergeErrorArray <- function(ScenarioLoc){
  #while loop to make multiple matrices into an array
  ArrayBookmark <- 1
  #Finds number of scenarios saved in wd for z value in the array
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #number of x slots in the array
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector() #finds all files in the Arp folder
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  #list of all the arp files
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles) #counts how many arp files are in the folder
  
  
  
  #number of y slots in the array
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  #creates initial array with the numbers found above
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))
  
  
  #while loop for the array
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    MergeErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    MergeErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propMergeError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindMergeError(allTables[[makeTablesCounter]])
        propMergeError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      MergeErrorMatrix[Matrixbookmark,] <- propMergeError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- MergeErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}

CreateP2HErrorArray <- function(ScenarioLoc){
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector()
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles)
  
  
  
  #y
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))
  

  
  
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    P2HErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    P2HErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propP2HError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindP2HError(allTables[[makeTablesCounter]])
        propP2HError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      P2HErrorMatrix[Matrixbookmark,] <- propP2HError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- P2HErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}

CreateH2PErrorArray <- function(ScenarioLoc){
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector()
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles)
  
  
  
  #y
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))

  
  
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    H2PErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    H2PErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propH2PError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindH2PError(allTables[[makeTablesCounter]])
        propH2PError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      H2PErrorMatrix[Matrixbookmark,] <- propH2PError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- H2PErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}

CreatePUnknownErrorArray <- function(ScenarioLoc){
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector()
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles)
  
  
  
  #y
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))
  
  
  
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    PUnknownErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    PUnknownErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propPUnknownError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindPUnknownError(allTables[[makeTablesCounter]])
        propPUnknownError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      PUnknownErrorMatrix[Matrixbookmark,] <- propPUnknownError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- PUnknownErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}

CreateHUnknownErrorArray <- function(ScenarioLoc){
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector()
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles)
  
  
  
  #y
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))
  
  
  
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    HUnknownErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    HUnknownErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propHUnknownError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindHUnknownError(allTables[[makeTablesCounter]])
        propHUnknownError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      HUnknownErrorMatrix[Matrixbookmark,] <- propHUnknownError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- HUnknownErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}

CreateCorrectPropArray <- function(ScenarioLoc){
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  totalFilesInFSCOut <- list.files(path = arpfilesloc)
  listArpFiles <- vector()
  arpFiles <- list.files(path = arpfilesloc)
  arpFileBookmark <- 1
  NumtotalArpFiles <- length(arpFiles)
  
  while(arpFileBookmark <= NumtotalArpFiles ){
    if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
      listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
    }
    arpFileBookmark <- arpFileBookmark +1
  }
  NumArpFiles <- length(listArpFiles)
  
  
  
  #y
  structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
  structOutFolders <- structOutFolders[-1]
  structOutFolders <- structOutFolders[-1]
  StructFolderX <- structOutFolders[ArrayBookmark] #calls singular OutFolder
  StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
  NumStructOutFiles <- length(StructOutFiles)
  arpnums <- vector()
  arpnumsbookmark <- 1
  while(arpnumsbookmark <= length(listArpFiles)){
    arpnums <- c(arpnums, arpnumsbookmark)
    arpnumsbookmark <- arpnumsbookmark +1
  }
  
  dimensionsforArray <- c(NumArpFiles, NumStructOutFiles, length(ScenarioList))
  ErrorArray <- array(dim=dimensionsforArray, dimnames = list(arpnums, StructOutFiles, ScenarioList))
  
  
  
  while(ArrayBookmark <= length(ScenarioList)){
    basepath <- paste0(ScenarioLoc,"/", ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    
    #finds the arp files to create labels
    a <- 1
    totalFilesInFSCOut <- list.files(path = arpfilesloc)
    b <- 1
    d <- 1
    ArpFilesInFSCOut <- as.numeric(0)
    totalFilesInFSCOut
    
    #only search for .arp files
    while(b <= length(totalFilesInFSCOut)){
      searchforarp <- str_contains(totalFilesInFSCOut[b], ".arp")
      if(searchforarp == TRUE){
        ArpFilesInFSCOut[d] <- totalFilesInFSCOut[b]
        d <- d+1
      }
      
      b <- b+1
    }
    ScenarioFolder <- paste0(ScenarioList[ArrayBookmark],"/")
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-1]
    structOutFolders <- structOutFolders[-1]
    CorrectPropErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    CorrectPropErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propCorrectPropError <- vector()
      
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
    
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        
        propHunknown <- FindHUnknownError(allTables[[makeTablesCounter]])
        propPUnknown <- FindPUnknownError(allTables[[makeTablesCounter]])
        propMergeErrror <- FindMergeError(allTables[[makeTablesCounter]])
        propP2H <- FindP2HError(allTables[[makeTablesCounter]])
        propH2P <- FindH2PError(allTables[[makeTablesCounter]])
        tempvect <- 1 - propHunknown - propPUnknown - propMergeErrror - propP2H - propH2P
        
        propCorrectPropError[makeTablesCounter] <- tempvect
        makeTablesCounter <- makeTablesCounter+1
      }
      CorrectPropErrorMatrix[Matrixbookmark,] <- propCorrectPropError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- CorrectPropErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
}
#Pipeline-----------------------------------------------------------------------

#Setting values the functions need
numberIndsPerSpecies <- 10
setwd("/Users/CHendrikse/Documents/REUHybridSimulation/Scenarios/")
ScenarioLoc <- getwd()
StructOutLoc <- "/Users/CHendrikse/Documents/REUHybridSimulation/Scenarios/4Deme10ind/"
CreateP2HErrorArray(ScenarioLoc)
CreateH2PErrorArray(ScenarioLoc)
CreatePUnknownErrorArray(ScenarioLoc)
CreateHUnknownErrorArray(ScenarioLoc)
CreateCorrectPropArray(ScenarioLoc)

#Before Running Structure----

arpfilesloc <- getwd()
listArpFiles <- vector()
arpFiles <- list.files(path = arpfilesloc)
arpFileBookmark <- 1
NumtotalArpFiles <- length(arpFiles)

while(arpFileBookmark <= NumtotalArpFiles ){
  if(str_contains(arpFiles[arpFileBookmark], ".arp") == TRUE){
    listArpFiles <- c(listArpFiles, arpFiles[arpFileBookmark]) 
  }
  arpFileBookmark <- arpFileBookmark +1
}
a <- 1
arpfilename <- as.numeric(0)
while(a <= length(listArpFiles)){
  arpfilename<- paste("4DemeFixed_1_", a, sep = "") #without .arp
filename <- arpfilename
#make the structure file and store the original species groups
labels <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder)
a<- a+1
}



#After Structure: reads structure's results and makes the table and outputs the proportion of success----
MergeErrorNumVect <- CreateMergeErrorArray(ScenarioLoc)