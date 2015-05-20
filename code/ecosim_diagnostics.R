#Ecosim group plots
#'Write function for Ecopath object
#'
#'Outputs basic parameters or mortalities to a .csv file.
#'
#'@family Rpath functions
#'
#'@param Rpath.sim.obj Rpath.sim object created by the ecosim.run() function.
#'@param append object name if appending from a previous run.  
#'
#'@return Creates a list of diagnostic tables from an ecosim run.
#'@export
ecosim.group.plot <- function(Rpath.sim.obj, append = NA){
  x <- copy(Rpath.sim.obj)
  n <- x$NUM_LIVING + x$NUM_DEAD
  
  #To/From tables
  qmat <- data.table(Prey = x$PreyFrom, Pred = x$PreyTo, Q = x$QQ)
  fmat <- data.table(Target = x$FishFrom, Gear = x$FishThrough, Fate = x$FishTo, Q = x$FishQ)
  
  #Create output list
  out <- c()
  for(i in 1:n){
    group <- i
    #Biomass
    group.table <- data.table(biomass = as.numeric(as.matrix(x$out_BB[i +1])))
    
    #QB
    Q <- qmat[Pred == group, sum(Q)]
    group.table[, Q := Q]
    group.table[, QB := Q / biomass]
    group.table[, Q := NULL]
    
    #M2
    qmat.pred <- qmat[Prey == group, ]
    if(nrow(qmat.pred) > 0){
      pred <- qmat.pred[, Pred]
      for(j in 1:length(pred)){
        group.table[, M2.pred := qmat.pred[j, Q] / biomass]
        setnames(group.table, 'M2.pred', paste('M2.', pred[j], sep = ''))
      }
      group.table[, M2 := rowSums(group.table[, 3:(length(pred) + 2), with = F])]
    } else {
      pred <- NULL
      group.table[, M2 := 0]
    }
    ncol <- 3 + length(pred)
    
    #M0
    group.table[, M0 := x$MzeroMort[i +1]]
    
    #M
    group.table[, M := M0 + M2]
    ncol <- ncol + 2
    
    #Landings 
    target <- fmat[Fate == 0 & Target == group, ]
    if(nrow(target) > 0){
      gear <- target[, Gear]
      for(j in 1:length(gear)){
        group.table[, C.gear := target[j, Q] * biomass]
        setnames(group.table, 'C.gear', paste('C.', gear[j], sep = ''))
      }
      group.table[, Land := rowSums(group.table[, (ncol + 1):(length(gear) + ncol), with = F])]
    } else {
      gear <- NULL
      group.table[, Land := 0]
    }
    ncol <- ncol + length(gear) + 1
    
    #Discards
    discard <- fmat[Fate != 0 & Target == group, ]
    if(nrow(discard) > 0){
      gear <- discard[, Gear]
      for(j in 1:length(gear)){
        group.table[, D.gear := discard[j, Q] * biomass]
        setnames(group.table, 'D.gear', paste('D.', gear[j], sep = ''))
      }
      group.table[, Disc := rowSums(group.table[, (ncol + 1):(length(gear) + ncol), with = F])]
    } else {
      gear <- NULL
      group.table[, Disc := 0]
    }
    ncol <- ncol + length(gear) + 1
    
    group.table[, Catch := Land + Disc]
    
    #F
    FRate <- fmat[Target == group, ]
    FRate <- FRate[, sum(Q), by = Gear]
    setnames(FRate, 'V1', 'Q')
    if(nrow(FRate) > 0){
      gear <- FRate[, Gear]
      for(j in 1:length(gear)){
        group.table[, F.gear := FRate[Gear == gear[j], Q]]
        setnames(group.table, 'F.gear', paste('F.', gear[j], sep = ''))
      }
    }
    
    group.table[, F := Catch / biomass]
      
    group.table[, Z := M + F]
    
    out$group.table <- group.table
    names(out)[i] <- x$spname[i + 1]
  }
}







 

#plot - six panels
n <- nrow(group.table)

#Biomass
ymax <- max(group.table[, biomass]) + .1 * max(group.table[, biomass])
plot(1:n, group.table[, biomass], type = 'l', main = 'Biomass', xlab = 'Month', 
     ylab = expression('t km'^-2), ylim = c(0, ymax))

#Q/B
ymax <- max(group.table[, QB]) + .1 * max(group.table[, QB])
plot(1:n, group.table[, QB], type = 'l', main = 'QB', xlab = 'Month', 
     ylab = expression('year'^-1), ylim = c(0, ymax))

#M2
M2.table <- group.table[, grep('M2.', names(group.table)), with = F]
ymax <- max(M2.table) + .1 * max(M2.table)
plot(1:n, as.matrix(M2.table[, 1, with = F]), type = 'l', main = 'M2', xlab = 'Month',
     ylab = expression('year'^-1), ylim = c(0, ymax), col = colors[pred[1]])
if(ncol(M2.table) > 1){
  for(i in 2:ncol(M2.table)) lines(1:n, as.matrix(M2.table[, i, with = F]), col = colors[pred[i]])
}

#Catch
C.table <- group.table[, grep('C.', names(group.table)), with = F]
ymax <- max(C.table) + .1 * max(C.table)
plot(1:n, as.matrix(C.table[, 1, with = F]), type = 'l', main = 'Landings', xlab = 'Month',
     ylab = expression('year'^-1), ylim = c(0, ymax), col = colors[pred[1]])
if(ncol(C.table) > 1){
  for(i in 2:ncol(C.table)) lines(1:n, as.matrix(C.table[, i, with = F]), col = colors[pred[i]])
}

#F
F.table <- group.table[, grep('F.', names(group.table)), with = F]
ymax <- max(F.table) + .1 * max(F.table)
plot(1:n, as.matrix(F.table[, 1, with = F]), type = 'l', main = 'F rates', xlab = 'Month',
     ylab = expression('year'^-1), ylim = c(0, ymax), col = colors[pred[1]])
if(ncol(F.table) > 1){
  for(i in 2:ncol(F.table)) lines(1:n, as.matrix(F.table[, i, with = F]), col = colors[pred[i]])
}

#Mortalities
