remotes::install_github('royfrancis/pophelper')

#--------------------------------------------
#HybridSimulation

#Script includes three parts
#The first is used before running STRUCTURE after making .arp files in FastSimCoal 2.6. the function makeStructure makes hybrid individuals pools all individuals into a genind
#    object and uses the function genind2structure (Clark 2017) to create a STRUCTURE input file
#The second part includes functions that create an array with the percentage of errors and "success" with the results matching the expected outcomes which is done by creating a
#    table for each run. The individuals are listed along with the correct placement then a column with STRUCTURE's decision is added based on Q values. Q vals beteen .4 and .6 are
#    listed as hybrids Q values above .9 are labeled as pure. the code then compares the correct value with the decision and labels if it is correct or if an error has been made
#The final part is listed as Archived Functions which were originally used to make arrays finding proportions of different errors

#Future change plans are noted in the comments and can be found by finding the word change

library(remotes)
library(diveRsity)
library(adegenet)
library(pophelper)
library(sjmisc)
library(dplyr)

#functions ---------------------------------------------------------------------
# Function uses the location of the arp files to make hybrids, pools them all together and makes a STRUCTURE input file using genind2structure 
makeStructure <- function(filename,a, numinds, ScenarioFolder, InFunction){
  #Finds the arp files location so the files can be looped to create the STRUCTURE input files for each arp file
  ArpFileLocation <- paste0(getwd(),"/", ScenarioFolder,"ArpFiles/", filename, ".arp")
  
  #Uses the name of the scenario folder to find the number of individuals per species (will need to be changed later in event name of scenarios need to be changed)
  splitScenario <- strsplit(ScenarioFolder, split = "_")
  numindstotal <- splitScenario[[1]][3]
  numinds <- gsub("ns/", "", numindstotal)
  as.numeric(numinds)
  
  #takes the file name and adds .arp so it can be found in the file
  arp2gen(ArpFileLocation)
  filename <- paste0(getwd(),"/", ScenarioFolder,"ArpFiles/", filename)

  
  #turning file into a genind object
  addgen<- paste(filename, ".gen", sep="")
  genindobj <- read.genepop(addgen, ncode = 3)
  
  #taking the genpop object and separating it into genind objects
  individuals <- seppop(genindobj,res.type="genind")
  
  #making hybrids from pop1 and pop2
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
  #change so parentandhybrid01, parentandhybrid02... parentandhybrid10
  filenameforstruct <- paste("parentandhybrid", a, sep="")
  filenameforstruct <- paste0("parentandhybrid", a, ".str")
  if(InFunction == FALSE){
  genind2structure(final_genind, file= filenameforstruct, pops=FALSE)
  }
  
  #return the labels for individuals in the order they are in the file
  return(actual_values)
}

#genind to structure file code (Clark 2017)
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
} 


#function that takes the Q values and creates a new column that labels STRUCTURE's decision changing the value in this function should allow the other functions to still run
#values between .4 and .6 are labeled as hybrid, above .9 are labeled as pure, and any other values are labeled as unknown
structuredecision <- function(obj){
  
  #putting individuals in groups depending on what STRUCTURE decides 
  obj$Decision<-0 #declare a new column, to hold Decisions
  for (i in 1:nrow(obj)){ #for all rows e.g. all individuals
    #first check if the individual is a pure species (>0.9) and if yes, assign it to the cluster for the column number which is >0.9
    if (any(obj[i,]>0.9)) obj$Decision[i]<-paste0("Cluster",(which(obj[i,]>0.9)))
    #check if its called as a hybrid
    if(any(obj[i,]>0.4 & obj[i,]<0.6)) obj$Decision[i]<-paste0("hybrid",(which(obj[i,]>0.4 & obj[i,]<0.6)))
    #Otherwise, e.g. it didn't get called a parent cluster or hybrid, call it unknown
    if(obj$Decision[i]==0)  obj$Decision[i]<-"unknown"
    
  }
  return(obj)
}

#function adding all true values
#clusters based on which cluster is most commonly used for that group in STRUCTURE, may be inaccurate if STRUCTURE does not run properly
AddingTrueGroup <- function(tabledecisions){
  
  #adding initial values for while loop a is the individual number and b/c are the threshold for the number of individuals from each species
  a <- 1
  b<-a*10
  c <- b-9
  #creating vector to store true values changed from vect <- "unknown" on 9/11/22
  vect <- vector()
 
  #while loop looking for the number of times the individual matches clusters with other individuals in the same species
  #if the cluster is used 5 or more times, it is assumed to be the correct cluster
  while(a <= length(tabledecisions$Cluster1)){
     decisionvect <- tabledecisions$Decision[c:b]
     
  #looks to see if a decision is decided five or more times within a cluster, if it is, that is determined to be the "true" value
    if(length(str_find(decisionvect,tabledecisions$Decision[a])) >= 5){
     vect[c:b] <- tabledecisions$Decision[a] 
    }
     
     #looks one cluster at a time when the cluster has been read through, it finds the next cluster (Change later so numinds value is used here rather than increasing by 10)
    if(a%%10 == 0){
      b <- b+10
      c<- c+10
    }
    a <- a+1
  }
  #putting the vector of actual clusters into a new column
  tabledecisions$ActualVals <- vect
  
  #returns the table with decisions
  return(tabledecisions)
}

#comparing groups and giving either TRUE or FALSE values depending on if STRUCTURE was correct in the individuals placement
CompareQ <- function(tabletocompare, originalLabels, numberIndsPerSpecies){
  
  #creating placeholder for loop
  bookmark <- 1
  
  #creating initial comparison object if the true value and decision are the same, the comparison is labeled as TRUE, if not, it is FALSE
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
  tabletocompare$Original <- originalLabels
  return(tabletocompare)
  }

#Finds the proportion of times STRUCTURE incorrectly labeled individuals from two different species as the same
#Note: if each species has 10 inds, but STRUCTURE labeled 20 inds in the same cluster, all 20 inds are labeled as incorrect

FindMergeError <- function(tabletocompare, numberIndsPerSpecies){
  #taking the number of individuals per species and finding the result for each group if two or more groups are put into the same cluster, a warning will print
#browser()
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

#Finds the P2H error by looking at if the decison is Hybrid when the true value is Pure 
FindP2HError <- function(completetable){
#loop placemarker
ClusterBookmark <- 1

#starts with initial error value of 0
numP2HError <- 0
#while loop if decision is labeled as hybrid, and the true value is not labeled as hybrid, the comparison label changes to label the error and the count of the error increases by 1
while(ClusterBookmark <= length(completetable$Cluster1)){
  #looks for the true and STRUCTUREs decision to be hybrid
  Calledhybrid <- str_contains(completetable$Decision[ClusterBookmark], "hybrid" )
  IndisHybrid <- str_contains(completetable$Original[ClusterBookmark], "hybrid")
  
  #if the decision is labeled hybrid, but it shouldn't the P2H label is added
  if(Calledhybrid == TRUE && IndisHybrid == FALSE){
    completetable$Comparison[ClusterBookmark] <- "P2H"
    numP2HError <- numP2HError + 1
  }
  #placeholder increases by 1
  ClusterBookmark <- ClusterBookmark +1
}
 
#find the proportion of the error by dividing the number of times it occurs by the total number of individuals 
  propP2HError <- numP2HError/length(completetable$Comparison)
  return(propP2HError)
   
}

#Finds the H2P error by looking at if the decison is Pure when the true value is Hybrid
FindH2PError <- function(completetable){
  
  #placeholder for loop
  ClusterBookmark <- 1
  
  #initial value of 0 for H2P error
  numH2PError <- 0
  
  #looks to see if the decision and true value are hybrid, if the decision is not labeled as hybrid and the true value is hybrid, the H2P label is added to the comparisson column
  while(ClusterBookmark <= length(completetable$Cluster1)){
    Calledhybrid <- str_contains(completetable$Decision[ClusterBookmark], "hybrid" )
    IndisHybrid <- str_contains(completetable$Original[ClusterBookmark], "hybrid")
    if(Calledhybrid == FALSE && IndisHybrid == TRUE){
      completetable$Comparison[ClusterBookmark] <- "H2P"
      numH2PError <- numH2PError + 1
    }
    
    ClusterBookmark <- ClusterBookmark +1
  }
  
  #finds the proportion of the error by dividing number of inds with the error by the total number of inds
  propH2PError <- numH2PError/length(completetable$Comparison)
  return(propH2PError)
  
}

#Finds the P2H error by looking at if the decison is unknown when the true value is Pure 
FindPUnknownError <- function(completetable){

    ClusterBookmark <- 1
    numPUError <- 0
    
    #while loop looking for the decision to be labeled as unknown and the true value is one of the pure species, if so the P Unknown label is added
    while(ClusterBookmark <= length(completetable$Cluster1)){
      #lookds for decision and true group
      CalledUnknown <- str_contains(completetable$Decision[ClusterBookmark], "unknown" )
      IndisPure <- str_contains(completetable$Original[ClusterBookmark], "pop")
      
      #adds label for inds with p unknown error
      if(CalledUnknown== TRUE && IndisPure == TRUE){
        completetable$Comparison[ClusterBookmark] <- "P Unknown"
        numPUError <- numPUError + 1
      }
      
      ClusterBookmark <- ClusterBookmark +1
    }
    
    #finds prop of the error by dividing number of inds with error by total number of inds
    propPUError <- numPUError/length(completetable$Comparison)
  
    return(propPUError)
    

}

#Finds the P2H error by looking at if the decison is unknown when the true value is hybrid 
FindHUnknownError <- function(completetable){
  ClusterBookmark <- 1
  numHUError <- 0
  
  #looks for decision to be labeled as unknown and the true value to be hybrid, if so, H unknown label is added
  while(ClusterBookmark <= length(completetable$Cluster1)){
    CalledUnknown <- str_contains(completetable$Decision[ClusterBookmark], "unknown" )
    IndisHybrid <- str_contains(completetable$Original[ClusterBookmark], "Hybrid")
    if(CalledUnknown== TRUE && IndisHybrid == TRUE){
      completetable$Comparison[ClusterBookmark] <- "H Unknown"
      numHUError <- numHUError + 1
    }
    
    ClusterBookmark <- ClusterBookmark +1
  }
  
  #finds prop by dividing number of inds with error by total number of inds
  propHUError <- numHUError/length(completetable$Comparison)
  return(propHUError)
  
  
}

#Calls the FindError functions to find proportion of error depending on what error type is specified. Options are Merge, P2H, H2P, PUnknown/P Unknown, HUnknown/H Unknown
FindErrors <- function(table, errortype){
  if(errortype == "Merge"){
    propError <- FindMergeError(table)
  }
  if(errortype == "P2H"){
    propError <- FindP2HError(table)
  }
  if(errortype == "H2P"){
    propError <- FindH2PError(table)
  }
  if(errortype == "PUnknown" || errortype == "P Unknown"){
    propError <- FindPUnknownError(table)
  }
  if(errortype == "HUnknown" || errortype == "H Unknown"){
    propError <- FindHUnknownError(table)
  }
  return(propError)
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
  completetable <- CompareQ(tabledecisions, originallabels, numberIndsPerSpecies)
  #prints the table
  print(completetable)
  return(completetable)
}

#Uses FindErrors() to decide which FindXError() function to use and creates array of the error. First specifies axes initial vals then inputs the information using the findError functions.
CreateErrorArray <- function(ScenarioLoc, errortype){

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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    
    #Finding number of inds change later
    totalstr <- ScenarioList[ArrayBookmark]
    totalstr[3]
    data <- strsplit(totalstr, split = "_")
    numindstotal <- data[[1]][3]
    numinds <- gsub("ns", "", numindstotal)
    numberIndsPerSpecies <- as.numeric(numinds)
    if (str_contains(ScenarioFolder, "20") == TRUE){
      numberIndsPerSpecies <- 20
    }
    if(str_contains(ScenarioFolder, "10") == TRUE){
      numberIndsPerSpecies <- 10
    }
    NumArpFiles <- length(ArpFilesInFSCOut)
    a <- 1
    labelsList <- list()
    #while a is less than or equal to the number of arp files, create labels for each file
    while(a <= NumArpFiles){
      arpfilename <- ArpFilesInFSCOut[a]
      arpfilename = unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      FileWithoutArp <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
      #make the structure file and store the original species groups
      InFunction <- TRUE
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder, InFunction)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-c(1,2)]
    ErrorMatrix <- matrix()
    Matrixbookmark <- 1
    StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
    StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
    NumStructOutFiles <- length(StructOutFiles)
    ErrorMatrix <- matrix(nrow = length(structOutFolders), ncol = NumStructOutFiles)
    
    while(Matrixbookmark <= length(structOutFolders)){
      StructFolderX <- structOutFolders[Matrixbookmark] #calls singular OutFolder
      StructOutFiles <- list.files(path = paste0(basepath ,"/", StructFolderX)) #finds all Structure output files within folder
      NumStructOutFiles <- length(StructOutFiles)
      makeTablesCounter <- 1
      propError <- vector()
      demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
      
      allTables <- vector(mode = "list", length = NumStructOutFiles) #creates initial object for all the tables for each run
      
      
      #creates vector of proportion of error for each structure file in a run
      while(makeTablesCounter <= length(StructOutFiles)){
        demeQmat <- readQ(paste(basepath,"/", StructFolderX,"/",StructOutFiles[makeTablesCounter],sep = ""))  
        
        results <- demeQmat[[StructOutFiles[makeTablesCounter]]]
        
        percenterror <- 0
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[1]], numberIndsPerSpecies, percenterror)
        tempvect <- FindErrors(allTables[[makeTablesCounter]], errortype)
        propError[makeTablesCounter] <- tempvect #check for error here
        makeTablesCounter <- makeTablesCounter+1
      }
      ErrorMatrix[Matrixbookmark,] <- propError
      Matrixbookmark <- Matrixbookmark + 1
    }
    
    ErrorArray[,,ArrayBookmark] <- ErrorMatrix #through here
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
  return(ErrorArray)
}


#creates an array for proportion of successful labels by taking 1 and subtracting the proportion of each error first specifies axes initial vals then inputs the information using the findError functions.
CreateCorrectPropArray <- function(ScenarioLoc){
#browser()
  ArrayBookmark <- 1
  ScenarioList <- list.dirs(path = ,full.names = FALSE, 
                            recursive = FALSE)
  #x
  basepath <- paste0(ScenarioLoc, ScenarioList[ArrayBookmark])
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    basepath <- paste0(ScenarioLoc, ScenarioList[ArrayBookmark])
    arpfilesloc <- paste0(basepath, "/ArpFiles" )
    scenarioName <- strsplit(basepath, "/")
    scenarioVariables <- strsplit(scenarioName[[1]][7], "_")
    numberIndsPerSpecies <- scenarioVariables[[1]][3]
  
    numberIndsPerSpecies <- as.numeric(gsub("ns", "", numberIndsPerSpecies))
    
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
      InFunction <- TRUE
      labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder, InFunction)
      a<- a+1
    }
    
    #input filename in ""
    #lists directories is specified structure path
    #Change path to ../ when arp file separate
    structOutFolders <- list.dirs(path= basepath, full.names = FALSE)[-1]
    structOutFolders <- structOutFolders[-c(1,2)]
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
        allTables[[makeTablesCounter]] <- maketable(results, labelsList[[makeTablesCounter]], numberIndsPerSpecies, percenterror)
        propHunknown <- FindHUnknownError(allTables[[makeTablesCounter]])
        propPUnknown <- FindPUnknownError(allTables[[makeTablesCounter]])
        propMergeError <- FindMergeError(allTables[[makeTablesCounter]], numberIndsPerSpecies)
        propP2H <- FindP2HError(allTables[[makeTablesCounter]])
        propH2P <- FindH2PError(allTables[[makeTablesCounter]])
        tempvect <- 1 - propHunknown - propPUnknown - propMergeError - propP2H - propH2P
        
        propCorrectPropError[makeTablesCounter] <- tempvect #look for error
        makeTablesCounter <- makeTablesCounter+1
      }
      CorrectPropErrorMatrix[Matrixbookmark,] <- propCorrectPropError
      Matrixbookmark <- Matrixbookmark + 1
    }

    ErrorArray[,,ArrayBookmark] <- CorrectPropErrorMatrix
    ArrayBookmark <- ArrayBookmark +1
  }
  print(ErrorArray)
  return(ErrorArray)
}




#Archived functions----
#Each function creates an array for a different kind of error
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    structOutFolders <- structOutFolders[-c(1,2)]
    
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
  ScenarioList <- list.dirs(path = ScenarioLoc,full.names = FALSE, 
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    if (str_contains(ScenarioFolder, "20") == TRUE){
      numberIndsPerSpecies <- 20
    }
    if(str_contains(ScenarioFolder, "10") == TRUE){
      numberIndsPerSpecies <- 10
    }
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
    structOutFolders <- structOutFolders[-c(1,2)]
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
  return(ErrorArray)
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    structOutFolders <- structOutFolders[-c(1,2)]
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    structOutFolders <- structOutFolders[-c(1,2)]
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
  structOutFolders <- structOutFolders[-c(1,2)]
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
    structOutFolders <- structOutFolders[-c(1,2)]
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
#Pipeline-----------------------------------------------------------------------

#Setting values the functions need ----

#input number of inds per species to decide how many hybrids to make (Should be changed later due to makestructure() using the file names rather than this input)
numberIndsPerSpecies <- 20
setwd("/Users/CHendrikse/Documents/HybridSimulation/Scenarios/4Deme20inds/StructOut7/")
setwd("/Users/CHendrikse/Documents/H")
setwd("/Users/clhen/Documents/HybridSimulation/TempDirectories/20ms") #for my laptop
ScenarioLoc <- getwd()
setwd("/Users/clhen/Documents/HybridSimulation/Scenarios/")
ScenarioLoc <- "/Users/clhen/Documents/HybridSimulation/TempDir15ms/"




#Before running structure ----
  ScenarioList <- list.dirs(path = , full.names = FALSE, recursive = FALSE)
  ScenarioNum <- ScenarioList[10]
  basepath <- paste0(ScenarioLoc,"/", ScenarioNum)
  arpfilesloc <- paste0(basepath, "/ArpFiles" )
  InFunction <- FALSE
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
  ScenarioFolder <- paste0(ScenarioNum,"/")
  NumArpFiles <- 10
  a <- 1
  labelsList <- list()
  #while a is less than or equal to the number of arp files, create labels for each file
  while(a <= NumArpFiles){
    arpfilename <- ArpFilesInFSCOut[a]
    arpfilename <- unlist(strsplit(arpfilename, split='.', fixed=TRUE))[1]
    #make the structure file and store the original species groups
    labelsList[[a]] <-makeStructure(arpfilename, a,numberIndsPerSpecies, ScenarioFolder, InFunction)
    a<- a+1
  }
  




#After Structure: reads structure's results and makes the table and outputs the proportion of success----
  ScenarioLoc <- "/Users/clhen/Documents/HybridSimulation/Scenarios/"
  correctarray <- CreateCorrectPropArray(ScenarioLoc)
  ErrorArray <- CreateErrorArray(ScenarioLoc, "Merge")
  
  #Creating boxplot for data
  boxplot(list(ErrorArray[,,1],ErrorArray[,,2],ErrorArray[,,3],ErrorArray[,,4],ErrorArray[,,5],ErrorArray[,,6],ErrorArray[,,7],ErrorArray[,,8],ErrorArray[,,9],ErrorArray[,,10]), ylim = c(0,1), pch = 19)
  boxplot(list(correctarray[,,1],correctarray[,,2],correctarray[,,3],correctarray[,,4],correctarray[,,5],correctarray[,,6],correctarray[,,7],correctarray[,,8], correctarray[,,9],correctarray[,,10],
               correctarray[,,11],correctarray[,,12],correctarray[,,13], correctarray[,,14],correctarray[,,15],
               correctarray[,,16],correctarray[,,17],correctarray[,,18], correctarray[,,19]), ylim = c(0,1), pch = 19)
  
  