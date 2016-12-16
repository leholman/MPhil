###Simulation of Breeding Program. 
##Luke E Holman, 2016

#RequiredLibraries----
library(plyr)
library(pedigree)

#***********GLOBAL EDITABLE PARAMETERS***********----

##Burn in generation parameters
#Number of burn in generations
burningens <- 5
#Number of Individuals Per Generation
Burngensize <- 1000
#Sex ratio 
sexratio <- 0.5
#Mating proportion
matingprop <- 1
Bdatarows <- Burngensize*burningens

##Current breeding scheme
#Number of tanks
selondaBreedingtanksSimulated <- 5
#Number of generations
breedinggens <- 3
#Number of broodstock fish
broodstocksize <- 50
#tank sex ratio
tanksexratio <- 0.4
BreedingRows <- breedinggens*broodstocksize

##Trait and selective breeding parameters
#Trait hereditability
traithered <- 0.4
#Trait mean
sizemean <- 400
#Trait standard devaition
sizesd <- 100
#Number of individuals in the grow out
growoutN <- 10000
#Number of individuals pre-selected
preselection <- 750
#BV EBV corralation
BLUPcorr <- 0.4
#Number of tanks for selected individuals
selectedgenerationTanks <- 4
#Number of selected individuals
selectedgenerationNumber <- 200
#Total number of selected generations
SelectionGenerationsSimulated <- 10


##COLONY/SIMPED
#number of loci simulated
numberofloci <-100
#proportion of error simulated
errorprop <- 0.02
#type of error simulated
errortype <- "none"


#***********Economic Modelling Paramters***********----

##Overall
yearsSimulated <- 20

##Growout
#Length of time in months to slaughter
TimeToSlaughter <- 1.5
#Feed conversion ratio
FeedConversionRatio <- 2.4
#Cost in euros per kilo of fish produced
TotalPerKiloCosts <- 3.5 
#Number of fish produced each year
numberofselectedfishproduced <- 100000000
NumberOfProductionFishAtSlaughter <- numberofselectedfishproduced
UnselectedProductionFishWeights <- rnorm(NumberOfProductionFishAtSlaughter,mean=0.4,sd=0.1)

##BreedingProgramCosts
#Cost in euros of molecular pedigree reconstruction
PerFishGenotypeCost <- 7
#Number of fish used in molecular pedigree reconstruction
fishGenotyped <- 750
#Size of growout
GrowOutSize <- 10000
#Monthly cost in euros of selective breeding program management
SBManagementFeeMonthly <- 7000
#Additional labour cost for selective breeding program per grow out
LabourForSelectionPerFish <- 1200 
#Other costs per grow out
OtherSelectiveBreedingCosts <- 0

##Sales
#Profit in euros per kilo of fish sold
perkiloprofit <- 0.9

##All Replicates----
#Total simulation number of replicates
totalreps <- 100

#***********Functions for simulation***********----

##Adding missing or errronous alleles to SNP data
AddMissingDataSNP <- function(input,proportion,type){
  if(!is.data.frame(input)){stop('Input does not exist or is not a Dataframe')}
  if(type=="missing"){
    workingdata <- as.matrix(input[,-1])
    prop <- length(workingdata)/2*proportion
    places <- data.frame("columns"=sample(seq(1,length(workingdata[1,]),by=2),prop,replace=TRUE),"rows"=sample(length(workingdata[,1]),prop,replace=TRUE))
    places$columns2 <- places$columns + 1
    for (i in 1:length(places$columns)){
      workingdata[places$rows[i],places$columns[i]] <- 0
      workingdata[places$rows[i],places$columns2[i]] <- 0
    }
    return(as.data.frame(cbind(input[,1],workingdata)))
    rm(workingdata,prop,places)
  }
  if(type=="error"){
    workingdata <- as.matrix(input[,-1])
    prop <- length(workingdata)/2*proportion
    places <- data.frame("columns"=sample(seq(1,length(workingdata[1,]),by=2),prop,replace=TRUE),"rows"=sample(length(workingdata[,1]),prop,replace=TRUE))
    places$columns2 <- places$columns + 1
    for (i in 1:length(places$columns)){
      workingdata[places$rows[i],places$columns[i]] <- sample(c(1,2),1)
      workingdata[places$rows[i],places$columns2[i]] <- sample(c(1,2),1)
    }
    return(cbind(input[,1],as.data.frame(workingdata)))
    rm(workingdata,prop,places)
  }
  if(type=="none"){
    return(input)
    rm(workingdata,prop,places)
  }
}


##Function that adds a random number between 10 and 80 to each value from a vector
Error_n_check = function(vec){
  perm_vec = vec+rnorm(length(vec),0,runif(1,10,80))
  cor_value =cor(perm_vec,vec)
  return(list(vec=perm_vec,cor=cor_value))
}


###Functions that take a dataframe of 200 inds with 'ID', 'Parent1' and 'Parent2'columns and returns alist of IDs minimising sibs in 4 tanks

##This function shuffles and scores the number of sib groups
ShuffleNScore <- function(frame){
  newframe <- frame[sample(1:200),]
  newframe$sibgroup <- paste(newframe$Parent1,newframe$Parent2,sep="_")
  T1 <- table(newframe$sibgroup[1:50])
  T2 <- table(newframe$sibgroup[51:100])
  T3 <- table(newframe$sibgroup[101:150])
  T4 <- table(newframe$sibgroup[151:200])
  Tall <-c(T1,T2,T3,T4)
  finalscore <-sum(Tall[Tall>1])
  indOrder <- newframe$ID
  return(list(order=indOrder,score=finalscore))
  rm(newframe,T1,T2,T3,T4,Tall,finalscore,indOrder)
}

##This function returns a character list
TankPickSibSplit <- function(ped){
  require('pedigree')
  frame <- SelectedParents
  n_iterations = lapply(1:1000, function(x) ShuffleNScore(SelectedParents))
  score_values = sapply(n_iterations, '[[', 'score')
  permutations <-n_iterations[[which.min(abs(score_values))]]
  return(permutations$order)  
}




#***********Simulation Loop Start***********----
#@#Create matrices to store results from whole simulation reps
replicateresults <- matrix(ncol=SelectionGenerationsSimulated+2,nrow=totalreps)
replicateInbrresults <- matrix(ncol=SelectionGenerationsSimulated+6,nrow=totalreps)
BioEconomicResults <- matrix(ncol=(yearsSimulated*12)+2,nrow=totalreps)
BioEconomicIncomings <- matrix(ncol=(yearsSimulated*12)+2,nrow=totalreps)
BioEconomicOutgoings <- matrix(ncol=(yearsSimulated*12)+2,nrow=totalreps)
colnames(BioEconomicResults) <- c("Repeat",seq(1,yearsSimulated*12+1,1))

#@#Start loop for complete replicates
for (allrep in 1:totalreps){
#@#Restore h2 from earlier simulations
traithered <- 0.4

#BurnIn of Wild Population----

#@#Create dataframe for wild population
originPOP <- data.frame("ID"=rep(NA,Bdatarows),"Parent1"=rep(NA,Bdatarows),"Parent2"=rep(NA,Bdatarows),"Gen"=rep(NA,Bdatarows),"Sex"=rep(NA,Bdatarows))
count <- 1

#@#start loop over wild population burn ingenerations
for (gen in 1:burningens){
  #@#if it's the first generation then start by setting zeros as parents for founders 
  if (gen==1){
    for (i in 1:Burngensize){
      genID <- paste("B",gen-1,sep="")
      originPOP[i,] <- c(paste(genID,i,sep="_"),0,0,"B0",sample(c("1","2"),1,prob=c(sexratio,(1-sexratio))))
      count <- count + 1
    }
  }
  else{
    #@#otherwise identify the current and parent generations and then randomly subsample parental allocation
    for (i in 1:Burngensize){
      genID <- paste("B",gen-1,sep="")
      parentgen <- paste("B",gen-2,sep="")
      potentialparents <- originPOP[!is.na(originPOP$ID),]
      potentialparents <- potentialparents[potentialparents$Gen==parentgen,]
      potentialparents <- potentialparents[potentialparents$ID %in% sample(potentialparents$ID,(Burngensize*matingprop)),]
      originPOP[count,] <- c(paste(genID,i,sep="_"),sample(potentialparents$ID[potentialparents$Gen==parentgen & potentialparents$Sex==1],1),sample(potentialparents$ID[potentialparents$Gen==parentgen & potentialparents$Sex==2],1),genID,sample(c("1","2"),1,prob=c(sexratio,(1-sexratio))))
      count <- count + 1
    } 
  }
}
#@#loop finished here, get rid of potential parent file and write burn in generations to the fish bank
rm(potentialparents)

runningFishBank <- cbind(originPOP,matrix(ncol=3,nrow=length(originPOP$ID)))

#Selonda Breeding Program as is for x generations----

##Non unifrom breeding spline setup

  #@#these data are form the parental contribution work provided by Xelect
  inds <- c(65,65,59,59,42,24,22,16,16,15,14,13,12,11,11,10,9,8,7,7,7,6,6,6,5,5,5,5,5,4,4,4,4,4,3,3,2,1,0)
  #@#here we set splines for both male and female contribution, some fish may have different male and female contirbutions 
  MSpline<-smooth.spline(1:length(inds),inds)
  FSpline<-smooth.spline(1:length(inds),inds) 

#@#note the final burn in generation and produce a blank variable for the tank names
finalburningen <- paste("B",burningens-1,sep="")
tank_ids <-c()

#@#loop over the number of tanks
for (tank in 1:selondaBreedingtanksSimulated){
  breedingresults <- data.frame("ID"=rep(NA,BreedingRows),"Parent1"=rep(NA,BreedingRows),"Parent2"=rep(NA,BreedingRows),"Gen"=rep(NA,BreedingRows),"Sex"=rep(NA,BreedingRows),"GeneticComponent"=rep(NA,BreedingRows),"EnvironComponent"=rep(NA,BreedingRows),"Phenotype"=rep(NA,BreedingRows))
  BrCount <- 1
  
  #@#loop over the generations within the tanks
  for (BrGen in 1:breedinggens){
    if (BrGen==1){
      #@#if it's the first generation we need to set the population from the wild burnin fish...
      parents <-sample(originPOP$ID[originPOP$Gen==finalburningen],broodstocksize)
      #@#...so wildcaughtF0 becomes our first parents for the simulation
      wildcaughtF0 <- originPOP[originPOP$ID %in% parents,]  
      BreedingGen <- paste("F",BrGen,sep="")
      #@#We call different parameters here as px, this is because we're adding each generation as blocks to the final datasheet, we do this becuase its more convenient than looping over a list. 
      p1 <- c(paste(BreedingGen,tank,1:broodstocksize,sep="_"))
      p2 <- sample(wildcaughtF0$ID[wildcaughtF0$Sex==1],broodstocksize,prob=predict(MSpline,seq(1,length(inds),length.out=length(wildcaughtF0$ID[wildcaughtF0$Sex==1])))[[2]],replace=TRUE)
      p3 <- sample(wildcaughtF0$ID[wildcaughtF0$Sex==2],broodstocksize,prob=predict(FSpline,seq(1,length(inds),length.out=length(wildcaughtF0$ID[wildcaughtF0$Sex==2])))[[2]],replace=TRUE)
      p4 <- rep(BreedingGen,broodstocksize)
      p5 <- sample(c(1,2),broodstocksize,prob=c(1-tanksexratio,tanksexratio),replace=TRUE)
      p6 <- rep(NA,broodstocksize)
      p7 <- rep(NA,broodstocksize)
      p8 <- rep(NA,broodstocksize)
      breedingresults[1:broodstocksize,] <-as.matrix(cbind(p1,p2,p3,p4,p5,p6,p7,p8))
      rm(p1,p2,p3,p4,p5,p6,p7,p8)
    }
    else{
      #@#the same as above but adding the block to the final results is more complex as we are defining the location in the block
      parentGenbr <- paste("F",BrGen-1,sep="")
      parents <- breedingresults[breedingresults$Gen==parentGenbr & !is.na(breedingresults$Gen),]
      BreedingGen <- paste("F",BrGen,sep="")
      p1 <- c(paste(BreedingGen,tank,1:broodstocksize,sep="_"))
      p2 <- sample(parents$ID[parents$Sex==1],broodstocksize,prob=predict(MSpline,seq(1,length(inds),length.out=length(parents$ID[parents$Sex==1])))[[2]],replace=TRUE)
      p3 <- sample(parents$ID[parents$Sex==2],broodstocksize,prob=predict(FSpline,seq(1,length(inds),length.out=length(parents$ID[parents$Sex==2])))[[2]],replace=TRUE)
      p4 <- rep(BreedingGen,broodstocksize)
      p5 <- sample(c(1,2),broodstocksize,prob=c(1-tanksexratio,tanksexratio),replace=TRUE)
      p6 <- rep(NA,broodstocksize)
      p7 <- rep(NA,broodstocksize)
      p8 <- rep(NA,broodstocksize)
      breedingresults[((broodstocksize*(BrGen-1))+1):((broodstocksize*(BrGen-1))+broodstocksize),] <-as.matrix(cbind(p1,p2,p3,p4,p5,p6,p7,p8)) 
    }
    
  }
  #@#here we are saving the individuals in the breeding scheme and calculating inbreeding for each tank
  wilds <-cbind(originPOP,matrix(ncol=3,nrow=length(originPOP$ID)))
  colnames(wilds) <- colnames(breedingresults)
  bigped <- rbind(wilds,breedingresults)
  colnames(runningFishBank) <- colnames(breedingresults)
  runningFishBank <- rbind(runningFishBank,breedingresults)
  bigped$Inbr <- calcInbreeding(bigped)

  #@#here we save the tank as a dataframe for later use
  assign(paste("CurrentSchemeTank",tank,sep=""),bigped)
  tank_ids <-c(tank_ids,paste("CurrentSchemeTank",tank,sep=""))
}

#Establish Phenotypes For Parent Fish ----

lastbreedinggen <- paste("F",breedinggens,sep="")
#@#here we assemble a list of all parents in tanks
listOfParents <- c()
for (tank in tank_ids){
  workingtank <- get(tank)
  for (eachparent in workingtank$ID[workingtank$Gen==lastbreedinggen])
    listOfParents <- c(listOfParents,eachparent)
}
#@#now we generate parental phenotypes according to the herditability and trait mean
parentalPhenotypes <- data.frame("ID"=listOfParents,"Genetic"=rnorm(length(listOfParents), mean=(sizemean*traithered),sd=sqrt((sizesd^2)*traithered)),"Environmental"=rnorm(length(listOfParents), mean=(sizemean*(1-traithered)),sd=sqrt((sizesd^2)*(1-traithered))))
parentalPhenotypes$Phenotype <- parentalPhenotypes$Genetic+parentalPhenotypes$Environmental

#Breed Growout Fish From Parents----

#@#first we set some parameters and an empty dataframe for breeding growout fish
pertankfish <- growoutN/selondaBreedingtanksSimulated
GrowOutFish <-as.data.frame(matrix(ncol=8,nrow=growoutN))
colnames(GrowOutFish)<- colnames(breedingresults)
count <- 1
#@#now we loop over the tank dataframes we created
for (tank in tank_ids){
  #@#again we are allocating a matrix to the final growout with each of px being used to form a column.
  workingtank <- get(tank)
  start <- 1+((count-1)*2000)
  end <- count*2000
  p1 <- paste("S0",start:end,sep="_")
  p2 <- sample(workingtank$ID[workingtank$Sex=="1" & workingtank$Gen==lastbreedinggen],pertankfish,prob=predict(MSpline,seq(1,length(inds),length.out=length(workingtank$ID[workingtank$Sex=="1" & workingtank$Gen==lastbreedinggen])))[[2]],replace=TRUE)
  p3 <- sample(workingtank$ID[workingtank$Sex=="2" & workingtank$Gen==lastbreedinggen],pertankfish,prob=predict(MSpline,seq(1,length(inds),length.out=length(workingtank$ID[workingtank$Sex=="2" & workingtank$Gen==lastbreedinggen])))[[2]],replace=TRUE)
  p4 <- rep("S0",pertankfish)
  p5 <- sample(c(1,2),pertankfish,prob=c(1-tanksexratio,tanksexratio),replace=TRUE)
  p6 <- rep(NA,pertankfish)
  p7 <- rep(NA,pertankfish)
  p8 <- rep(NA,pertankfish)
  GrowOutFish[start:end,] <-as.matrix(cbind(p1,p2,p3,p4,p5,p6,p7,p8)) 
  
  count <- count+1
}
#@#here we call up the ID of both parents and then calculate the genetic value for the offspring as a mean of the two parents
for (i in 1:growoutN){
  wSire <- GrowOutFish$Parent1[i]
  wDam <- GrowOutFish$Parent2[i]
  GrowOutFish$GeneticComponent[i] <- mean(c(parentalPhenotypes$Genetic[parentalPhenotypes$ID==wSire],parentalPhenotypes$Genetic[parentalPhenotypes$ID==wDam]))  
}
#@#We generate a random environmental effect here and then calculate phneotype for each growout fish
GrowOutFish$EnvironComponent <- rnorm(growoutN, mean=(sizemean*(1-traithered)),sd=sqrt((sizesd^2)*(1-traithered)))
GrowOutFish$Phenotype <- as.numeric(GrowOutFish$GeneticComponent)+ GrowOutFish$EnvironComponent

#Pre-SelectFish On Phenotype and calculateEBV----
#@#order fish by phenotype
orderedGrowOut <- GrowOutFish[order(-GrowOutFish$Phenotype),]
#@#subset by preselection number
Preselected <- orderedGrowOut[1:preselection,]    
#@#set genetic value as breeding value 
BV <- as.numeric(Preselected$GeneticComponent)


##We error function the function 1000 times and select the iterate that is closest to the value we desire
n_iterations = lapply(1:1000, function(x) Error_n_check(BV))
cor_values = sapply(n_iterations, '[[', 'cor')
permutations <-n_iterations[[which.min(abs(cor_values - BLUPcorr))]]
#@#this function stops everything if the random function above cannot get close to the desired BLUP
if(sqrt((BLUPcorr-permutations$cor)^2) > 0.2) {stop("BLUP correlation cannot be performed - change BLUP accuracy")}

##EBV is expressed as additative genetic value above the mean of the growout
Preselected$EBV <- permutations$vec - mean(as.numeric(GrowOutFish$GeneticComponent))

#Segregate Selected Fish Calculate Gain and Randomise Next Generation's Additative Effect----
EBVordered <- Preselected[order(-Preselected$EBV),]
#@#here we check that fish and tanks can be split up as designed
if ((selectedgenerationNumber/selectedgenerationTanks)%%1==0){print("Splitting Selected Fish Into Tanks")}else{stop("Number of Tanks and Fish do Not Fit Nicely STOP and change")}
#@#we shuffle the fish that are going to be selected
EBVordered[1:selectedgenerationNumber,] <- EBVordered[sample(1:selectedgenerationNumber),]
#@#then we pull them out already shuffled
SelectedFish <- EBVordered[1:selectedgenerationNumber,]
#@#they are allocated to a new variable
BreedingProgramParentalFish <- SelectedFish
#@#and the trait increase is recorded
traitincrease <- mean(as.numeric(SelectedFish$GeneticComponent)) - mean(as.numeric(GrowOutFish$GeneticComponent))
#@#we increase the next generations mean value
NewMeanTraitValue <- mean(GrowOutFish$Phenotype)+ traitincrease
##SetVa
Va <- (sizesd^2)*traithered
#@#calculate reduction in varaince for this generation
reduction <- (1-(traithered*0.843))
#@#new variance
Va <- reduction*Va
#@#and then adjust the phenotypes for the folowing generation - since the fish are recorded with the old value this is fine, we do this because it makes it easier to initiate further generations  
SelectedFish$GeneticComponent <- rnorm(length(SelectedFish$GeneticComponent), mean=(NewMeanTraitValue*traithered),sd=sqrt(Va))
SelectedFish$EnvironComponent <- rnorm(length(SelectedFish$GeneticComponent), mean=(NewMeanTraitValue*(1-traithered)),sd=sqrt(6000))
SelectedFish$Phenotype <- SelectedFish$GeneticComponent + SelectedFish$EnvironComponent
SelectedParents <- SelectedFish
#@#new trait hered
traithered <- Va/(Va+6000)
#SetDataframeForResultsChecking----
#@#now we are setting up dataframes to store the fish from the new selection program
newrows <- as.data.frame(matrix(ncol=8,nrow=((SelectionGenerationsSimulated)*growoutN)))
colnames(newrows) <- colnames(GrowOutFish)
BreedingProgramResults <- rbind(GrowOutFish,newrows)
newrows <- as.data.frame(matrix(ncol=8,nrow=((SelectionGenerationsSimulated)*selectedgenerationNumber)))
colnames(newrows) <- colnames(GrowOutFish)
BreedingProgramParentalFish <- rbind(BreedingProgramParentalFish[,1:8],newrows)
PreParentageFish <- c()

#StartGenerational Loops Here----



#@#Start with a loop over all generations
for (SelectedGeneration in 1:SelectionGenerationsSimulated){
  #@#record parents and current generation labels
  CurrentSelectedGen <- paste("S",SelectedGeneration,sep="")
  lastbreedinggen <- paste("S",SelectedGeneration-1,sep="")
  ##Set tanks for new fish
  tank_ids <- c()
  #@#then loop over tanks within generations
  for (newtank in 1:selectedgenerationTanks){
    #@#Reorder the selected parents to minimise sib-sib groupings
    order <- TankPickSibSplit(SelectedParents)
    SelectedParents <- SelectedParents[match(order,SelectedParents$ID),]
    #@#the below expresson calls out the tank based on order - its a bit clumsy like this but should scale as values change
    workingtank <- SelectedParents[(1+((newtank-1)*(selectedgenerationNumber / selectedgenerationTanks))):(newtank*(selectedgenerationNumber / selectedgenerationTanks)),]
    assign(paste("CurrentSelectedTank",newtank,sep=""),workingtank)
    tank_ids <- c(tank_ids,paste("CurrentSelectedTank",newtank,sep=""))
  }
  
  #Breed Growout Fish From Parents
  #@#now we're doing a growout for the seelcted fish
  pertankfish <- growoutN/selectedgenerationTanks
  GrowOutFish <-as.data.frame(matrix(ncol=8,nrow=growoutN))
  colnames(GrowOutFish)<- colnames(breedingresults)
  count <- 1
  for (tank in tank_ids){
    #@#as before using px to generate a matrix
    workingtank <- get(tank)
    start <- 1+((count-1)*pertankfish)
    end <- count*pertankfish
    p1 <- paste(CurrentSelectedGen,start:end,sep="_")
    p2 <- sample(workingtank$ID[workingtank$Sex=="1" & workingtank$Gen==lastbreedinggen],pertankfish,prob=predict(MSpline,seq(1,length(inds),length.out=length(workingtank$ID[workingtank$Sex=="1" & workingtank$Gen==lastbreedinggen])))[[2]],replace=TRUE)
    p3 <- sample(workingtank$ID[workingtank$Sex=="2" & workingtank$Gen==lastbreedinggen],pertankfish,prob=predict(MSpline,seq(1,length(inds),length.out=length(workingtank$ID[workingtank$Sex=="2" & workingtank$Gen==lastbreedinggen])))[[2]],replace=TRUE)
    p4 <- rep(CurrentSelectedGen,pertankfish)
    p5 <- sample(c(1,2),pertankfish,prob=c(1-tanksexratio,tanksexratio),replace=TRUE)
    p6 <- rep(NA,pertankfish)
    p7 <- rep(NA,pertankfish)
    p8 <- rep(NA,pertankfish)
    GrowOutFish[start:end,] <-as.matrix(cbind(p1,p2,p3,p4,p5,p6,p7,p8)) 
    
    count <- count+1
  }
  
  #@#here we are allocating genetic, enviornmental and phenotypic values to the offspring
  for (i in 1:growoutN){
    wSire <- GrowOutFish$Parent1[i]
    wDam <- GrowOutFish$Parent2[i]
    GrowOutFish$GeneticComponent[i] <- mean(c(SelectedParents$GeneticComponent[SelectedParents$ID==wSire],SelectedParents$GeneticComponent[SelectedParents$ID==wDam]))  
  }
  GrowOutFish$EnvironComponent <- rnorm(growoutN, mean=(NewMeanTraitValue*(1-traithered)),sd=sqrt(6000))
  GrowOutFish$Phenotype <- as.numeric(GrowOutFish$GeneticComponent)+ GrowOutFish$EnvironComponent
  
  #Pre-SelectFish On Phenotype and calculateEBV----
  #@#order the growout by phenotype
  orderedGrowOut <- GrowOutFish[order(-GrowOutFish$Phenotype),]
  #@#create a dataframe for preselection
  Preselected <- orderedGrowOut[1:preselection,]    
  #@#subset another dataframe for COLONY to use later
  PreParentageFish <- rbind(PreParentageFish,Preselected)
  #@#set the breeding values to a variable 
  BV <- as.numeric(Preselected$GeneticComponent)
  
  
  ##We then iterate the function 1000 times and select the iterate that is closest to the value we desire
  #@#as before we add error to the breeding values and then select the value for comparison
  n_iterations = lapply(1:1000, function(x) Error_n_check(BV))
  cor_values = sapply(n_iterations, '[[', 'cor')
  permutations <-n_iterations[[which.min(abs(cor_values - (BLUPcorr)))]]
  
  
  ##EBV is expressed as additative genetic value above the mean of the growout
  Preselected$EBV <- permutations$vec - mean(as.numeric(GrowOutFish$GeneticComponent))
  
  #Segregate Selected Fish Calculate Gain and Randomise Next Generation's Additative Effect----
  #@#order the preselected fish by EBV
  EBVordered <- Preselected[order(-Preselected$EBV),]
  #@#Check to make sure fish can be split into tanks
  if ((selectedgenerationNumber/selectedgenerationTanks)%%1==0){print("Splitting Selected Fish Into Tanks")}else{stop("Number of Tanks and Fish do Not Fit Nicely STOP and change")}
  #@#shuffle selected fish
  EBVordered[1:selectedgenerationNumber,] <- EBVordered[sample(1:selectedgenerationNumber),]
  #@#subset shuffled fish
  SelectedFish <- EBVordered[1:selectedgenerationNumber,]
  #@#add fish to parent list
  BreedingProgramParentalFish[(1+(SelectedGeneration*selectedgenerationNumber)):((SelectedGeneration+1)*selectedgenerationNumber),] <- SelectedFish[,1:8]
  #@#calculate trait increase
  traitincrease <- mean(as.numeric(SelectedFish$GeneticComponent)) - mean(as.numeric(GrowOutFish$GeneticComponent))
  #@#calculate new mean trait value
  NewMeanTraitValue <- mean(GrowOutFish$Phenotype)+ traitincrease
  #@#calculate reduction in varaince for this generation
  reduction <- (1-(traithered*0.843))
  #@#new variance
  Va <- reduction*Va
  #@#change selected fish genetic and environmental components to simulate increase of mean for offspring
  SelectedFish$GeneticComponent <- rnorm(length(SelectedFish$GeneticComponent), mean=(NewMeanTraitValue*traithered),sd=sqrt(Va))
  SelectedFish$EnvironComponent <- rnorm(length(SelectedFish$GeneticComponent), mean=(NewMeanTraitValue*(1-traithered)),sd=sqrt(6000))
  SelectedFish$Phenotype <- SelectedFish$GeneticComponent + SelectedFish$EnvironComponent
  SelectedParents <- SelectedFish
  #@#new hered
  traithered <- Va/(Va+6000)
  #RecordGrowout Fish
  BreedingProgramResults[(((SelectedGeneration)*growoutN)+1):((SelectedGeneration+1)*growoutN),] <- GrowOutFish
  mean(GrowOutFish$Phenotype)
  
  #@#turn the below boxplot on for phenotype over years
  # boxplot(BreedingProgramResults$Phenotype~BreedingProgramResults$Gen,xlab="Generation",ylab="Weight At Slaughter (g)")
  
  
}
#@#Put all the data into a large dataframe and calculate inbreeding
all <- rbind(runningFishBank,BreedingProgramResults)
all$inbr <- calcInbreeding(ped = all)

#@#here we calculate per generation inbreeding including the selected (S) base generation (B) and current scheme (F)
InbrB4 <- mean(all$inbr[all$Gen=="B4"])
InbrF1 <- mean(all$inbr[all$Gen=="F1"])
InbrF2 <- mean(all$inbr[all$Gen=="F2"])
InbrF3 <- mean(all$inbr[all$Gen=="F3"])
InbrS0 <- mean(all$inbr[all$Gen=="S0"])
InbrS1 <- mean(all$inbr[all$Gen=="S1"])
InbrS2 <- mean(all$inbr[all$Gen=="S2"])
InbrS3 <- mean(all$inbr[all$Gen=="S3"])
InbrS4 <- mean(all$inbr[all$Gen=="S4"])
InbrS5 <- mean(all$inbr[all$Gen=="S5"])
InbrS6 <- mean(all$inbr[all$Gen=="S6"])
InbrS7 <- mean(all$inbr[all$Gen=="S7"])
InbrS8 <- mean(all$inbr[all$Gen=="S8"])
InbrS9 <- mean(all$inbr[all$Gen=="S9"])
InbrS10 <- mean(all$inbr[all$Gen=="S10"])

#@#we put the results into the results matrix
replicateInbrresults[allrep,] <- c(allrep,InbrB4,InbrF1,InbrF2,InbrF3,InbrS0,InbrS1,InbrS2,InbrS3,InbrS4,InbrS5,InbrS6,InbrS7,InbrS8,InbrS9,InbrS10)


#@#we calcuate the mean trait value for each generation
S0 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S0"],na.rm=TRUE)
S1 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S1"],na.rm=TRUE)
S2 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S2"],na.rm=TRUE)
S3 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S3"],na.rm=TRUE)
S4 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S4"],na.rm=TRUE)
S5 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S5"],na.rm=TRUE)
S6 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S6"],na.rm=TRUE)
S7 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S7"],na.rm=TRUE)
S8 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S8"],na.rm=TRUE)
S9 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S9"],na.rm=TRUE)
S10 <- mean(BreedingProgramResults$Phenotype[BreedingProgramResults$Gen=="S10"],na.rm=TRUE)

#@#we add them tio the result sheet
replicateresults[allrep,] <- c(allrep,S0,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10)

####BIOECONOMIC SECTION

##SimulationSetup
#@#set up a month balacen sheet
balancesheet <- matrix(ncol=(yearsSimulated*12),nrow=3)
#@#the below number is the number of months until S1 fish are availiable for production
monthsuntilnewfish <-36
#@#The below number is the number of months until S2 fish are produced in the hatchery  (- 1 to allow new fish to be in the correct month)
monthsuntilnewbroodstockgen <-71
ProductionBroodstockGen <- 1
#@#Fish is a variable that contains a summary of fish in production
Fish <- matrix(ncol=2,nrow=0)
#@#normal harvest is an estimate of the profit from a unselected harvest  
normalharvest <- sum(rnorm(numberofselectedfishproduced,0.4,0.1)*perkiloprofit)
#@#Cost of genotpying fish
SelectiveBreedingGrowOutCost <- (fishGenotyped*PerFishGenotypeCost)
#@#The below values indicate how frequently the genotyping costs are paid
growoutcostMonths <- seq(1,yearsSimulated*12, by=36)





####Actual Sim----
#@#start monthly sim
for (month in 1:(yearsSimulated*12)){
  #@#Convert/constrain Fish to matrix
  Fish <- as.matrix(Fish,ncol=2)
  #@#take one month from time until fish ready
  Fish[,2]<- as.numeric(Fish[,2]) -1
  #@#select the current production broodstock
  BreedingGen <- paste("S",ProductionBroodstockGen,sep="")
  #@#set costs and sales at nil for the month
  sales <- 0
  costs <- 0
  costs <- SBManagementFeeMonthly
  #@# if the current month is a growout cost month the below expression adds the cost of a full growout including food and genotyping to the monthly costs
  if(month %in% growoutcostMonths){GrowOutCost <- sum((as.numeric(all$Phenotype[all$Gen==BreedingGen]))/1000)* TotalPerKiloCosts
                                   costs <- (costs + GrowOutCost + SelectiveBreedingGrowOutCost + LabourForSelectionPerFish)  } 
  #@# if there is a new broodstock fish available then the current breeding gen is selected and new fish are also selected for trait values
  if(monthsuntilnewbroodstockgen<1){ProductionBroodstockGen <- ProductionBroodstockGen+1
                                    monthsuntilnewbroodstockgen <- 36}
  #@# if new production fish can be started they are added to the growout fish
  if(monthsuntilnewfish<1){
    Fish <-rbind(Fish,c(BreedingGen,18))
    monthsuntilnewfish <- monthsuntilnewfish+12 }
  #@# if production fish are availiable for harvest they are harvested and the profit against the unselected fish is recorded
  if(0 %in% Fish[,2]){thisharvestgen <- Fish[which((Fish[,2] == "0")==TRUE),1]
                      harvest <- sum(sample(as.numeric(all$Phenotype[all$Gen==thisharvestgen]),numberofselectedfishproduced,replace=TRUE)/1000)*perkiloprofit
                      sales <- harvest - normalharvest
                      Fish <- matrix(Fish[-which(Fish[,2] == "0"),],ncol = 2)}
  #@#the below line prints some details about each month
  print(paste("Costs - ",costs,", Sales - ",sales,". Current ProductionBroodtock -",BreedingGen,sep=""))
  #@#the next two lines count the months down until new fish and new broodstock
  monthsuntilnewbroodstockgen <- monthsuntilnewbroodstockgen - 1
  monthsuntilnewfish <- monthsuntilnewfish -1
  #@#record all the balances
  balancesheet[1,month] <- costs
  balancesheet[2,month] <- sales
  balancesheet[3,month] <- balancesheet[2,month] - balancesheet[1,month] 
  }

#@#the below loop calcuates profit per month and then adds the results to the results sheet
  for (month in 1:(yearsSimulated*12)){
   
    BioEconomicResults[allrep,month] <- balancesheet[3,month]
    BioEconomicIncomings[allrep,month] <- balancesheet[2,month]
    BioEconomicOutgoings[allrep,month] <- balancesheet[1,month] 
  }
BioEconomicResults[allrep,1] <- allrep


##SIMPED & COLONY FILES----


##GenerateSimPedPedigree
#@#calculate the size of the dataframe       
size <- length(PreParentageFish$ID)
#@#generate a dataframe in the SimPed format for the entire loop without S0 individuals
SimPed <- data.frame("Ped"=rep(1,size),"ID"=PreParentageFish$ID,"SireID"=PreParentageFish$Parent1,"DamID"=PreParentageFish$Parent2,"Sex"=PreParentageFish$Sex,"Generate"=rep(1,size))
#@#set items as character
SimPed$SireID <-as.character(SimPed$SireID)
SimPed$DamID <-as.character(SimPed$DamID)
#@#set number of S0 individuals
numadded <- length(BreedingProgramParentalFish$ID[BreedingProgramParentalFish$Gen=="S0"])
#@#set S0 parents individually as founder individuals to simulate parentage panel being developed
S0Parents <- data.frame("Ped"=rep(1,numadded),"ID"=BreedingProgramParentalFish$ID[BreedingProgramParentalFish$Gen=="S0"],"SireID"=rep(0,numadded),"DamID"=rep(0,numadded),"Sex"=BreedingProgramParentalFish$Sex[BreedingProgramParentalFish$Gen=="S0"],"Generate"=rep(1,numadded))
#@#bind the SimPed file and the base generation
SimPed <- rbind(S0Parents,SimPed)
#@#write the SimPed file out
write.table(SimPed,file="/Users/Luke/Desktop/Simpeds/PED.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#@#get SimPed to produce genotypes
system("(cd /Users/Luke/Desktop/Simpeds; ./sim infile.txt )")

##InputData
#@#pull the genotype data back
loopGenetic <- read.table(file="/Users/Luke/Desktop/Simpeds/genoout.out",sep=" ")
#@#pull out information items
loopGenetic1 <- data.frame("ID"=loopGenetic$V3,"Parent1"=loopGenetic$V5,"Parent2"=loopGenetic$V7,"Sex"=loopGenetic$V9) 
#@#pull out genetic information
loopGenetic2 <- loopGenetic[,12:((numberofloci*2)+11)]
#@#set up clean dataset for COLONY files
loopGenetic <- cbind(loopGenetic1,loopGenetic2)



##COLONY file per generation

#@#loop over the selected generations 
for (SelectedGeneration in 1:SelectionGenerationsSimulated){
    #@# note down offspring and parent generation in the loop
    offGen <- paste("S",SelectedGeneration,sep="")
    parentGen <- paste("S",SelectedGeneration-1,sep="")
    #@# set run name
    runname <- paste("RUN_",allrep,"_",parentGen,offGen,sep="")
    #@#subset loop offspring and parents by generation
    loopoffspring <- loopGenetic[loopGenetic$ID %in% all$ID[all$Gen==offGen],]
    loopSires <- loopGenetic[loopGenetic$ID %in% all$ID[all$Gen==parentGen & all$Sex==1],]
    loopDams <- loopGenetic[loopGenetic$ID %in% all$ID[all$Gen==parentGen & all$Sex==2],]
    #@#write out pedigree for error checking
    write.csv(loopoffspring,file=paste("/Users/Luke/Desktop/Simpeds/colonyfiles/pedigree/",runname,"PED.txt",sep=""))
    
    #@#get rid of unwanted variables eg. sex
    loopoffspring <-loopoffspring[c(-2,-3,-4)]
    loopSires <-loopSires[c(-2,-3,-4)]
    loopDams <-loopDams[c(-2,-3,-4)]
    
    
    #@#set location for colony file
    location <- paste("/Users/Luke/Desktop/Simpeds/colonyfiles/",runname,".dat",sep="")
    
    #COLONYFILE
    #@#below script writes colony file
    write(runname, file = location)
    write(runname, file = location,append=TRUE)
    write(paste(length(loopoffspring$ID),"    !Off Number",sep=""), file = location, append=TRUE )
    write(paste(numberofloci,"          !Loci NUmber",sep=""), file = location, append=TRUE )
    write("1234      ! Seed for random number generator", file =location, append=TRUE )
    write("1         ! 0/1=Not updating/updating allele frequency", file =location, append=TRUE )
    write("2         ! 2/1=Dioecious/Monoecious species", file =location, append=TRUE )
    write("1         ! 0/1=Inbreeding absent/present", file =location, append=TRUE )
    write("0         ! 0/1=Diploid species/HaploDiploid species", file =location, append=TRUE )
    write("0  0      ! 0/1=Polygamy/Monogamy for males & females", file =location, append=TRUE )
    write("0           ! 0/1=Clone inference =No/Yes", file =location, append=TRUE )
    write("1           ! 0/1=Full sibship size scaling =No/Yes", file =location, append=TRUE )
    write("0         ! 0/1/2=no prior/sibship size prior/sibship complexity prior", file =location, append=TRUE )
    write("0         ! 0/1=Unknown/Known population allele frequency", file =location, append=TRUE )
    write("1         ! Number of runs", file =location , append=TRUE )
    write("2         ! 1/2/3/4 = Short/Medium/Long/VeryLong run", file =location, append=TRUE )
    write("1         ! 0/1=Monitor method by Iterate#/Time in second", file =location, append=TRUE )
    write("1         ! Monitor interval in Iterate# / in seconds", file =location, append=TRUE )
    write("0         ! 0/1=DOS/Windows version", file =location, append=TRUE )
    write("1         ! 0/1/2=Pair-Likelihood-Score(PLS)/Full-Likelihood(FL)/FL-PLS-combined(FPLS) method", file =location, append=TRUE )
    write("1         ! 0/1/2/3=Low/Medium/High/VeryHigh precision", file =location , append=TRUE )
    write("", file =location, append=TRUE )
    write(paste("M",1:numberofloci,sep="_",collapse = " "), file =location, append=TRUE )
    write((paste((seq(0,0, length.out=numberofloci)), collapse =" ")), file =location, append=TRUE )
    write((paste((seq(0,0, length.out=numberofloci)), collapse =" ")), file =location, append=TRUE )
    write((paste((seq(0.0,0.0, length.out=numberofloci)), collapse =" ")), file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write.table(loopoffspring, file=location, quote=FALSE, row.names=FALSE, col.names=FALSE, append =TRUE)
    write("", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("1.0  1.0", file =location, append=TRUE )
    write(paste(length(loopSires$ID),"  ",length(loopDams$ID),sep=""), file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write.table(loopSires, file=location, quote=FALSE, row.names=FALSE, col.names=FALSE, append =TRUE)
    write("", file =location, append=TRUE )
    write.table(loopDams, file=location, quote=FALSE, row.names=FALSE, col.names=FALSE, append =TRUE)
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )
    write("", file =location, append=TRUE )
    write("0", file =location, append=TRUE )   
}

print(paste(allrep,"% Complete",sep=""))
}

##Write Out results
write.csv(BioEconomicResults,file="~/Desktop/simulationresults/BioEconomic.csv")
write.csv(replicateInbrresults,file="~/Desktop/simulationresults/Inbreeding.csv")
write.csv(replicateresults,file="~/Desktop/simulationresults/GeneticGain.csv")
write.csv(BreedingProgramResults,file="~/Desktop/simulationresults/lastgenFish.csv")
write.csv(BioEconomicIncomings,file="~/Desktop/simulationresults/BioEcoIN.csv")
write.csv(BioEconomicOutgoings,file="~/Desktop/simulationresults/BioEcoOUT.csv")
