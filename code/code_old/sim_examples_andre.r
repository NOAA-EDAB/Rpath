################################################################################
# EXAMPLES OF USING ECOPATH AND ECOSIM in R and C 
# written by Kerim Aydin and Sarah Gaichas
# v0.04 3/12/13
################################################################################

# Directory with all the files
  setwd("C:/Users/kaydin/d/src/sims/sim_v04")
  setwd("C:\\src\\sim_v04")
# R code files with Ecopath and Ecosim routines
  source("ecopath_r_v0_04.r")
  source("ecosim_r_v0_04.r")
  
# Load ecosim dll                                             
  ll() 

# Ecopath input parameters are kept in csv files.  See examples for format.
# Note: there is very limited checking to ensure files have the right number
# of columns, likely outcome of poor csv files is an error.

# Example is a reduced version of the Eastern Bering Sea (Aydin et al. 2007)

  # Basic Ecopath parameters (Biomass, Production, Consumption, Fishing)
    BaseFile <- "data/EBS_andre_base.csv"
  # Ecopath Diet Matrix
    DietFile <- "data/EBS_andre_diet.csv"
  # Data pedigree (if no pedigree, use file of shown format, with 1.0 for all values) 
    PedFile  <- "data/EBS_andre_ped.csv"
  # Juvenile adult split groups.  IF NO JUV/ADU IN YOUR MODEL, use file with the 
  # headers in the example, but no rows.  This is required as a placeholder 
  # (minor bug needs fixing)
    JuvFile  <- "data/EBS_andre_juvs.csv"


# EXAMPLES

# EXAMPLE 1:  BALANCE ECOPATH, SENDING RESULTS TO NEW CSV FILE and set of vectors
#             (leave filename blank or FALSE for vectors only)
  
  path <- ecopathR(BaseFile,DietFile,PedFile,"EBS_balanced")  
  summary(path) 
  
# EXAMPLE 2:  Load Ecosim parameters, declaring enough vector length for 100 years.
#             Run once in baseline (equilibrium), run with fishing change.

  # Prepare for run by loading from csv files (includes balancing Ecopath)
    base_sim <- load_ecosim( Years = 100, BaseFile, DietFile, PedFile, JuvFile)
    #long_sim <- load_ecosim( Years = 500, BaseFile, DietFile, PedFile, JuvFile)
  # Base run from year 0 to year 100
    run0 <- ecosim_run(base_sim,0,100)  

  # resulting output biomass
	summary(run0$out_BB)

  # Set up a fishing scenario, set a high F Rate on Pollock (note: this doesn't 
  # remove the fishing on gears) for 90 years.  
	fish_sim <- base_sim 	
	fish_sim$ratelist$NoIntegrate[]  <- 100
	fish_sim$ratelist$NoIntegrate[2] <- 1
	fish_sim$ratelist$NoIntegrate[3] <- 2
    fish_sim$fishlist[] <- 0
    
    fishrates <- seq(0,1.2,0.01)
    KK   <- length(fishrates)
    polF <- codF <- polB <- codB <- polC <- codC <- rep(0,KK*KK)
    
    mat_polB <- mat_polC <- mat_codB <- mat_codC <- matrix(0,KK,KK)
    nn <- 0
    for (i in 1:KK){
        fpol <- fishrates[i]
        for (j in 1:KK){	
            fcod <- fishrates[j] 
            nn <- nn + 1
            cat(nn,fpol,fcod,"\n"); flush.console()
	        fish_sim$FORCED_FRATE$'Pollock'[2:101] <- fpol
	        fish_sim$FORCED_FRATE$'Cod'[2:101]     <- fcod
	        run1 <- ecosim_run(fish_sim,0,100) 
	        polF[nn] <- fpol
	        codF[nn] <- fcod
	        polB[nn] <- run1$out_BB$'Pollock'[1200]
	        polC[nn] <- run1$out_CC$'Pollock'[1200]
	        codB[nn] <- run1$out_BB$'Cod'[1200]
	        codC[nn] <- run1$out_CC$'Cod'[1200]	     
            mat_polB[i,j] <- polB[nn]   
            mat_polC[i,j] <- polC[nn]   
            mat_codB[i,j] <- codB[nn]   
            mat_codC[i,j] <- codC[nn]   
	    }
	}

    fmin <- 1
    fmax <- 81
    
    s_polB<-mat_polB[fmin:fmax,fmin:fmax]
    s_polC<-mat_polC[fmin:fmax,fmin:fmax]
    s_codB<-mat_codB[fmin:fmax,fmin:fmax]
    s_codC<-mat_codC[fmin:fmax,fmin:fmax]
    f_ind <- fishrates[fmin:fmax]
    sk    <- length(f_ind)
    
    prec<-crec<-pOFL<-cOFL<-rep(0,sk)
    for (cc in 1:sk){
        poldif    <- abs(s_polB[,cc]/max(s_polB[,cc]) - 0.4) 
        coddif    <- abs(s_codB[cc,]/max(s_codB[cc,]) - 0.4)        
        prec[cc]  <- f_ind[which(poldif==min(poldif))[1]] 
        crec[cc]  <- f_ind[which(coddif==min(coddif))[1]] 
    }

    sscod <- sspol <-matrix(0,sk,sk)
    for (cc in 1:sk){
        poldif    <- (s_polB[,cc]/max(s_polB[,cc])) 
        coddif    <- (s_codB[cc,]/max(s_codB[cc,]))        
        sspol[,cc]  <- poldif
        sscod[cc,]  <- coddif
    }    
     

    
    par(mfrow=c(2,2))
    
    image.plot(f_ind,f_ind,s_polB,xlab="Pollock fishing rate",ylab="Cod fishing rate")
    title("(A) Pollock Biomass")
    image.plot(f_ind,f_ind,s_codB,xlab="Pollock fishing rate",ylab="Cod fishing rate")
    title("(B) Cod Biomass")        
    image.plot(f_ind,f_ind,s_polC*12,xlab="Pollock fishing rate",ylab="Cod fishing rate")
    title ("(C) Pollock Catch")
    image.plot(f_ind,f_ind,s_codC*12,xlab="Pollock fishing rate",ylab="Cod fishing rate")
    title ("(D) Cod Catch")    

findpoint <- function(x1,y1,x2,y2,f_in){      
    mat1 <- matrix(c(x1,y1),length(x1),2)
    mat2 <- matrix(c(x2,y2),length(x2),2)
    dmat <- as.matrix(nearest.dist(mat1,mat2))
    xclose <-    mean( mat1[which(dmat==min(dmat),arr.ind=T)[1],][1],
                       mat2[which(dmat==min(dmat),arr.ind=T)[2],][1])
    yclose <-    mean( mat1[which(dmat==min(dmat),arr.ind=T)[1],][2],
                       mat2[which(dmat==min(dmat),arr.ind=T)[2],][2])                       
    xref <- which(abs(f_in-xclose)==min(abs(f_in-xclose)))
    yref <- which(abs(f_in-yclose)==min(abs(f_in-yclose)))
    return(matrix(c(xref,yref),1,2))
}

  #image.plot(f_ind,f_ind,(s_polB/max(s_polB)<=0.35)|(s_codB/max(s_codB)<=0.35)) 
   plot(f_ind,f_ind,type='n',xlab="Pollock fishing rate (1/year)",ylab="Cod fishing rate (1/year)")
   
  .filled.contour(f_ind,f_ind,pmin(sspol,sscod),levels=c(-1,0.35,2),col=c("gray85","white"))
 # Multispecies F40   
    test1     <- contourLines(f_ind,f_ind,s_polB/max(s_polB),levels=0.4) 
    test2     <- contourLines(f_ind,f_ind,s_codB/max(s_codB),levels=0.4)
    outpoint  <- findpoint(test1[[1]]$x,test1[[1]]$y, test2[[1]]$x,test2[[1]]$y, f_ind)
    contour(f_ind,f_ind,s_codB/max(s_codB),levels=0.4,col='blue',drawlabels=F,add=T)
    contour(f_ind,f_ind,s_polB/max(s_polB),levels=0.4,col='blue',drawlabels=F,add=T)
    points(f_ind[outpoint[,1]],f_ind[outpoint[,2]],col="blue",pch=16,cex=2)
    outpoint1<-outpoint
    
 # Single Species F40
    ptt <- predict(loess(prec~f_ind),f_ind)
    ctt <- predict(loess(crec~f_ind),f_ind)
    outpoint  <- findpoint(ptt,f_ind,f_ind,ctt, f_ind)
    lines(f_ind,ctt,col='red')    
    lines(ptt,f_ind,col='red')    
    points(f_ind[outpoint[,1]],f_ind[outpoint[,2]],col="red",pch=16,cex=2)
    outpoint2<-outpoint
    
 # M proportional
    Mbase_pol <- 0.7435447   
    Mbase_cod <- 0.3335637
    polFrac = Mbase_pol/(Mbase_pol+Mbase_cod) 
    Bsumfrac <- (s_codB+s_polB)/(s_codB[1,1]+s_polB[1,1])
    lines(f_ind,polFrac*f_ind,col="green")
    contour(f_ind,f_ind,Bsumfrac,levels=0.4,col='green',drawlabels=F,add=T)
    test     <- contourLines(f_ind,f_ind,Bsumfrac,levels=0.4)                                                        
    outpoint <- findpoint(test[[1]]$x,test[[1]]$y, f_ind,polFrac*f_ind, f_ind)
    points(f_ind[outpoint[,1]],f_ind[outpoint[,2]],col="green",pch=16,cex=2)
    outpoint3<-outpoint    
    
 # Max Total Catch   
    outpoint<- which(s_polC+s_codC==max(s_polC+s_codC),arr.ind=T)
    points(f_ind[outpoint[,1]],f_ind[outpoint[,2]],col="black",pch=16,cex=2)
    outpoint4<-outpoint
    
    par(mfrow=c(3,2))
    
    barplot(c(s_polB[1,1],s_polB[outpoint2],s_polB[outpoint1],s_polB[outpoint3],s_polB[outpoint4]),
    names.arg=c("Unfished","A","B","C","D")); title("Pollock Biomass")    
    barplot(c(s_codB[1,1],s_codB[outpoint2],s_codB[outpoint1],s_codB[outpoint3],s_codB[outpoint4]),
    names.arg=c("Unfished","A","B","C","D")); title("Cod Biomass")    
    barplot(12*c(0,s_polC[outpoint2],s_polC[outpoint1],s_polC[outpoint3],s_polC[outpoint4]),
    names.arg=c("Unfished","A","B","C","D")); title("Pollock Catch")    
    barplot(12*c(0,s_codC[outpoint2],s_codC[outpoint1],s_codC[outpoint3],s_codC[outpoint4]),
    names.arg=c("Unfished","A","B","C","D")); title("Cod Catch")    
    barplot(c(0.0,0.48,0.39,0.48,0.49),
    names.arg=c("Unfished","A","B","C","D")); title("Pollock Fishing Rate")    
    barplot(c(0.0,0.21,0.16,0.33,0.80),
    names.arg=c("Unfished","A","B","C","D")); title("Cod Fishing Rate")    



#> f_ind[outpoint1]
#[1] 0.39 0.16
#> f_ind[outpoint2]
#[1] 0.48 0.21
#> f_ind[outpoint3]
#[1] 0.48 0.33
#> f_ind[outpoint4]
#[1] 0.49 0.80

#A - Convergence - SS option
#B - option Multispecies
#C - scalar
#D - systemwide ms
#
#D1 (D2 and D3 similar)
#
#
#3.3 sum of stocks
#3.4 system MSY


