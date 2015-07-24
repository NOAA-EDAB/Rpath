  
library(Rpath)
 
Ebase <- "data/EBS_andre_base.csv"
Ediet <- "data/EBS_andre_diet.csv"
Eped  <- "data/EBS_andre_ped.csv"
Ejuv  <- "data/EBS_andre_juvs.csv"
EBS   <- ecopath(Ebase, Ediet, Eped, eco.name = 'E. Bering')

#Ebase <- "data/ECS_eis_base_July2015.csv"
#Ediet <- "data/ECS_eis_diet_Jun2015.csv"
#Eped  <- "data/ECS_eis_ped_Jun2015.csv"
#Ejuv  <- "data/ECS_eis_juv_July2015.csv"
#EBS   <- ecopath(Ebase, Ediet, Eped, eco.name = 'Chukchi')

EBASE  <- ecosim.init(EBS,Ejuv)
EBASE$FORCED_FRATE[1:30,2]<-0.05
ERUN <- ecosim.run(EBASE,0,100)
plot(ERUN$out_BB[1:1200,2])


TBASE  <- ecotest.init(EBS)
TBASE$fishing$FRATE[1:30,2]<-0.05
TRUN <- ecotest.run(TBASE)
plot(TRUN$out_BB[1:1200,2])


etest.plot(TRUN$out_BB)

etest.plot<-function(outBB) {
  N <- length(outBB[1,])
  plot( c(1,1201), c(0,2) ,type='n') 
  for (i in 1:N) {
    lines(1:1201, outBB[,i]/outBB[1,i])    
   }  
}



test <- Adamstest(EBASE)
tderiv<-ecotest(EBASE,1,1,1)

microbenchmark()


#
as.data.frame(ecotest(ERUN,1,1,1))

require(microbenchmark)
microbenchmark(ecosim.run(EBASE,0,100),times=100L)

,Adamstest(EBASE),times=100L)

#EBS_0 <- ecosim.init(EBS,Ejuv,YEARS=100)
Estate <- ecosim.state(Epar) #optional argument for state  #return Rpath.sim.state
Eforce <- ecosim.forcing(Epar) #optional argument for forcing  #return Rpath.sim.forcing

Erun   <- ecosim.run1(Epar,Estate,Eforce,)
