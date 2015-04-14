#Ecosim group plots
x <- REco.s1
group <- 5
colors <- rainbow(x$NUM_LIVING + x$NUM_DEAD)

#Biomass
group.table <- data.table(biomass = as.numeric(as.matrix(x$out_BB[group + 1])))

#QB
qmat <- data.table(Prey = x$PreyFrom, Pred = x$PreyTo, Q = x$QQ)
Q <- qmat[Pred == group, sum(Q)]
group.table[, Q := Q]
group.table[, QB := Q / biomass]

#M2
qmat.pred <- qmat[Prey == group, ]
if(nrow(qmat.pred) > 0){
  pred <- qmat.pred[, Pred]
  for(i in 1:length(pred)){
    group.table[, M2.pred := qmat.pred[i, Q] / biomass]
    setnames(group.table, 'M2.pred', paste('M2.', pred[i], sep = ''))
  }
  group.table[, M2 := rowSums(group.table[, 4:(length(pred) + 3), with = F])]
} else group.table[, M2 := 0]
ncol <- 4 + length(pred)

#M0



#Catch - 
fmat <- data.table(Target = x$FishFrom, Gear = x$FishThrough, Fate = x$FishTo, Q = x$FishQ)
target <- fmat[Fate == 0 & Target == group, ]
if(nrow(target) > 0){
  gear <- target[, Gear]
  for(i in 1:length(gear)){
    group.table[, C.gear := target[i, Q] * biomass]
    setnames(group.table, 'C.gear', paste('C.', gear[i], sep = ''))
  }
  group.table[, Totland := rowSums(group.table[, (ncol + 1):(length(gear) + ncol), with = F])]
} else group.table[, Totland := 0]
ncol <- ncol + length(gear) + 1

discard <- fmat[Fate != 0 & Target == group, ]
if(nrow(discard) > 0){
  gear <- discard[, Gear]
  for(i in 1:length(gear)){
    group.table[, D.gear := discard[i, Q] * biomass]
    setnames(group.table, 'D.gear', paste('D.', gear[i], sep = ''))
  }
  group.table[, Totdisc := rowSums(group.table[, (ncol + 1):(length(gear) + ncol), with = F])]
} else group.table[, Totdisc := 0]
ncol <- ncol + length(gear) + 1

group.table[, Totcatch := Totland + Totdisc]

#F
catch <- fmat[Target == group, ]
catch <- catch[, sum(Q), by = Gear]
setnames(catch, 'V1', 'Q')
if(nrow(catch) > 0){
  gear <- catch[, Gear]
  for(i in 1:length(gear)){
    group.table[, F.gear := discard[i, Q] * biomass]
    setnames(group.table, 'D.gear', paste('D.', gear[i], sep = ''))
  }
  
group.table[, F := Totcatch / biomass]

group.table[, Z := M + F]

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
