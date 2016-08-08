## R version of Ecopath balance by Sarah Gaichas and Kerim Aydin
## Version 0.04 3/11/13
#Modified by Sean Lucey 6/14
## Function ecopathR takes as input 3 csv files and optional
## output file name.  

ecopathR <- function(modfile, dietfile, pedfile, outname = FALSE){

  require(MASS); require(data.table)

  model <- as.data.table(read.csv(modfile))  # Basic parameters, detritus fate, catch, discards in that order
  diet  <- as.data.table(read.csv(dietfile)) # diet matrix
  ped   <- as.data.table(read.csv(pedfile))  # pedigree file

  #Remove first column of diet matrix (assumed to be names)
  diet.group <- names(diet)[1]
  diet <- diet[, which(!names(diet) == diet.group), with = F]   
  diet[is.na(diet)] <- 0
  
  #Remove first column of pedigree (assumed to be names)
  ped.group <- names(ped)[1]
  ped <- ped[, which(!names(ped) == ped.group), with = F]     

  #Get number of groups, living, dead, and gear.  Must be input in that order
  ngroups <- length(model[, Type])
  nliving <- length(model[Type <  2, Type])
  ndead   <- length(model[Type == 2, Type])
  ngear   <- length(model[Type == 3, Type])

  nodetrdiet <- diet[1:nliving,]
  model[is.na(DetInput), DetInput := 0]

  # fill in GE and QB from inputs
  model[, GE := as.double(ProdCons)]
  model[is.na(ProdCons), GE := PB / QB]
  model[is.na(QB),       QB := PB / GE]

  # define catch, discards, necessary sums
    catchmat<-model[,(10+ndead+1):(10+ndead+ngear)]
    discardmat<-model[,(10+ndead+1+ngear):(10+ndead+(2*ngear))]
    totcatchmat<-catchmat+discardmat
    
    # KYA 1/16/14 Need if statement here because rowSums fail if only one 
    # fishery (catch is vector instead of matrix)
    if (is.matrix(totcatchmat)){
        totcatch <- rowSums(totcatchmat)
        Catch    <- rowSums(catchmat)    #Catch<-rowSums(catchmat)
        Discards <- rowSums(discardmat) #Discards<-rowSums(discardmat)
       gearCatch <- colSums(catchmat, na.rm=T)
       gearDisc  <- colSums(discardmat, na.rm=T)
    }else{
        totcatch <-totcatchmat
        Catch    <-catchmat    #Catch<-rowSums(catchmat)
        Discards <-discardmat #Discards<-rowSums(discardmat)
        gearCatch <- sum(catchmat, na.rm=T)
        gearDisc  <- sum(discardmat, na.rm=T)                     
    }   
    
    gearTot<-gearCatch+gearDisc
    model<-cbind(model, Catch, Discards)
    model<-cbind(model, totcatch)

  # flag missing pars and subset for estimation
    noB<-ifelse(is.na(model$Biomass),1,0)
    noEE<-ifelse(is.na(model$EE),1,0)
    alive<-ifelse(model$Type<2,1, 0)
    model<-cbind(model, noB, noEE, alive)

  # define detritus fate matrix
    detfate<-model[,(10+1):(10+ndead)]

  # set up and solve the system of equations for living group B or EE
    living<-subset(model, alive==1)
    Q<-living$totcatch+living$BioAcc
    A<-matrix(0, nliving, nliving)
    diag(A)<-ifelse(living$noEE==1, living$Biomass*living$PB, living$PB*living$EE)
    QBDC<-as.matrix(nodetrdiet)*living$QB[col(as.matrix(nodetrdiet))]
    dimnames(QBDC) <- list(NULL, NULL)
    QBDC[is.na(QBDC)]<-0
    QBDCa<-as.matrix(QBDC)*living$noB[col(as.matrix(QBDC))]
    A<-A-QBDCa 
    BioQB<-living$Biomass*living$QB
    cons<-as.matrix(nodetrdiet)*BioQB[col(as.matrix(nodetrdiet))]
    Q<-Q+rowSums(cons, na.rm=T)  

  # Generalized inverse does the actual solving
    pars<-ginv(A, tol=.Machine$double.eps)%*%Q
    EE<-ifelse(is.na(living$EE), pars*living$noEE, living$EE)
    B<-ifelse(is.na(living$Biomass),pars*living$noB, living$Biomass)

  # detritus EE calcs
    M0<-living$PB*(1-EE)
    QBloss<-living$QB
    QBloss[is.na(QBloss)]<-0
    loss<-c((M0*B)+(B*QBloss*living$Unassim), 
          model$DetInput[model$Type==2], 
          gearDisc)
    detinputs<-colSums(loss*detfate)
    detdiet<-diet[(nliving+1):(nliving+ndead),]
    BQB<-B*living$QB
    detcons<-as.matrix(detdiet)*BQB[col(as.matrix(detdiet))]
    detoutputs<-rowSums(detcons, na.rm=T)
    EE<-c(EE, as.vector(detoutputs/detinputs))

  # added by kya
    # if a detritus biomass is put into the spreadsheet, use that and 
    # calculate PB.  If no biomass, but a PB, use that pb with inflow to 
    # calculate biomass.  If neither, use default PB=0.5, Bio = inflow/PB  
    # This is done because Ecosim requires a detrital biomass.
    Default_Detrital_PB <- 0.5 
    inDetPB <-model$PB[(nliving+1):(nliving+ndead)] 
    inDetB  <-model$Biomass[(nliving+1):(nliving+ndead)]
    DetPB <- ifelse(is.na(inDetPB),Default_Detrital_PB, inDetPB)
    DetB  <- ifelse(is.na(inDetB), detinputs/DetPB, inDetB)
    DetPB <- detinputs/DetB

  # Trophic Level calcs
    TL<-rep(1, ngroups)
    TLcoeff<-matrix(0, ngroups, ngroups)
    diag(TLcoeff)<-rep(1, ngroups)
    gearcons<-as.matrix(totcatchmat)/gearTot[col(as.matrix(totcatchmat))]
    dimnames(gearcons) <- list(NULL, NULL)
    gearcons[is.na(gearcons)]<-0
    dietplus<-as.matrix(diet)
    dimnames(dietplus)<- list(NULL, NULL)
    dietplus<-rbind(dietplus, matrix(0,ngear, nliving))
    dietplus<-cbind(dietplus, matrix(0,ngroups,ndead), gearcons)
    TLcoeffA<-TLcoeff-dietplus
    TL<-solve(t(TLcoeffA), TL)     

  # Path outputs to csv: unbalanced flag (in name), basic ests, and morts
    wrt = ifelse(outname==FALSE, FALSE, TRUE)
    if(wrt){outname<-ifelse(max(EE, na.rm=T)>1, paste("UB_", outname, sep=""), outname)}
    #kya changed these following four lines for detritus, and removing NAs
    #to match header file format (replacing NAs with 0.0s)
    Bplus<-c(B, DetB, rep(0.0, ngear))
    PBplus<-model$PB; 
    PBplus[(nliving+1):(nliving+ndead)] <- DetPB
    PBplus[is.na(PBplus)]<-0.0
    EEplus<-c(EE, rep(0.0, ngear))
    QBplus<-model$QB; QBplus[is.na(QBplus)] = 0.0
    GE[is.na(GE)] <- 0.0
    RemPlus <- model$totcatch; RemPlus[is.na(RemPlus)] <- 0.0
    balanced<-list(Group=model$Group, TL=TL, Biomass=Bplus, PB=PBplus, QB=QBplus, EE=EEplus, GE=GE, Removals=RemPlus)

    if(wrt){write.csv(balanced, paste(outname,"_balanced.csv", sep=""))}
    M0plus<-c(M0, as.vector(detoutputs/detinputs))
    gearF<-as.matrix(totcatchmat)/B[row(as.matrix(totcatchmat))]
    newcons<-as.matrix(nodetrdiet)*BQB[col(as.matrix(nodetrdiet))]
    predM<-as.matrix(newcons)/B[row(as.matrix(newcons))]
    predM<-rbind(predM, detcons)
    morts<-list(Group=model$Group[model$Type<3], PB=model$PB[model$Type<3], M0=M0plus, F=gearF[1:(nliving+ndead),], M2=predM)
    if(wrt){write.csv(morts, paste(outname,"_morts.csv", sep=""))}
     
  # cleanup before sending to sim -- C code wants 0 as missing value, not NA
    balanced$Biomass[is.na(balanced$Biomass)]<-0
    balanced$PB[is.na(balanced$PB)]<-0
    balanced$QB[is.na(balanced$QB)]<-0
    balanced$EE[is.na(balanced$EE)]<-0
    balanced$GE[is.na(balanced$GE)]<-0
    model$BioAcc[is.na(model$BioAcc)]<-0
    model$Unassim[is.na(model$Unassim)]<-0
    dietm<-as.matrix(diet)
    dimnames(dietm) <- list(NULL, NULL)
    dietm[is.na(dietm)]<-0
    catchmatm<-as.matrix(catchmat)
    dimnames(catchmatm) <- list(NULL, NULL)
    catchmatm[is.na(catchmatm)]<-0
    discardmatm<-as.matrix(discardmat)
    dimnames(discardmatm) <- list(NULL, NULL)
    discardmatm[is.na(discardmatm)]<-0
    detfatem<-as.matrix(detfate)
    dimnames(detfatem) <- list(NULL, NULL)
    detfatem[is.na(detfatem)]<-0
    pedm<-as.matrix(ped)
    dimnames(pedm) <- list(NULL, NULL)
    pedm[is.na(pedm)]<-0

  # list structure for sim inputs
    ecosim <- list(
    NUM_GROUPS = ngroups , #define NUM_GROUPS 80  INCLUDES GEAR
    NUM_LIVING = nliving , #define NUM_LIVING 60
    NUM_DEAD = ndead   , #define NUM_DEAD 3
    NUM_GEARS= ngear   , #define NUM_GEARS 17
    spname = as.character(balanced$Group),
    type = model$Type,
    TL   = TL,
    BB = balanced$Biomass,       #float path_BB[1..NUM_GROUPS] vector
    PB = balanced$PB,       #float path_PB[1..NUM_GROUPS] vector
    QB = balanced$QB,       #float path_QB[1..NUM_GROUPS] vector
    EE = balanced$EE,       #float path_EE[1..NUM_GROUPS] vector
    BA = model$BioAcc,       #float path_BA[1..NUM_GROUPS] vector
    GS = model$Unassim,       #float path_GS[1..NUM_GROUPS] vector
    GE = balanced$GE,       #float path_GS[1..NUM_GROUPS] vector
    pedigree = (pedm), #float pedigree[B,PB,QB,Diet,1..NUM_GEARS][1..NUM_LIVING+NUM_DEAD]  matrix
    DC = (dietm),       #float path_DC[1..NUM_GROUPS][1..NUM_GROUPS]  matrix in [prey][pred] order     NUM_LIVING?
    DetFate = (detfatem),  #float path_DetFate[1..NUM_DEAD][1..NUM_GROUPS]  matrix in [det][groups] order
    Catch = (catchmatm),    #float path_Catch[1..NUM_GEARS][1..NUM_GROUPS]  matrix
    Discards = (discardmatm)  #float path_Discards[1..NUM_GEARS][1..NUM_GROUPS] matrix
    )
    
return(ecosim)
}
