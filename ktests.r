  
library(Rpath)
 
Ebase <- "data/EBS_andre_base.csv"
Ediet <- "data/EBS_andre_diet.csv"
Eped  <- "data/EBS_andre_ped.csv"
Ejuv  <- "data/EBS_andre_juvs.csv"

Ebase <- "data/ECS_eis_base_July2015.csv"
Ediet <- "data/ECS_eis_diet_Jun2015.csv"
Eped  <- "data/ECS_eis_ped_Jun2015.csv"
Ejuv  <- "data/ECS_eis_juv_July2015.csv"
EBS   <- ecopath(Ebase, Ediet, Eped, eco.name = 'Chukchi')

EBASE  <- ecosim.init(EBS,Ejuv)
TBASE  <- ecotest.init(EBS)


ERUN <- ecosim.run(EBASE,0,100)
TRUN <- ecotest.run(TBASE)

etest.plot<-function(TB) {
  N <- length(TB$out_BB[1,])
  plot( c(1,1201), c(0,2) ,type='n') 
  for (i in 1:N) {
    lines(1:1201, TB$out_BB[,i]/TB$out_BB[1,i])    
   }  
}

etest.plot(TRUN)

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
