# replicate Marchinko
# April 5, 2018

replicate_marchinko <- function(marchinko){
  # columns are
  # population
  # family
  # treatment {no = -insect, pred = +insect}
  # individual
  # standlength
  # ant.dorspine
  # sec.dorspine
  # pelspine.len
  # pelgirdle.len
  # plate.number
  # eda.genotype {AA, Aa, aa}
  # pelgirdle.presence {1=yes, 0 = no}
  # pelspine.presence {1 = yes, 0 = no}
  
  # Because spine and girdle lengths grow with body size, I corrected these traits for size using residuals from an ordinary least squares regressions of each trait on standard length.
  
  # Standardized selection differentials (i) were calculated according to equation (6.1) in Endler (1986), i = \bar{X}_a - \bar{X}_b)/\sqrt{var_b} where \bar{X}_a and \bar{X}_b were the mean trait values of fish from a single family measured at the end of the predation and control treatments, respectively, and var_b is trait variance in the control treatment. Selection differentials were calculated for each family separately and may be found in the Supporting Table S1
  
  # standlength replicates but...
  # regression ignoring treatment level does not replicate
  # regression within treatment levels does not replicate
  # regression within family replicates. Lots of noise here, a multilevel model would be better
  # *** recode any length=0.0 to NA
  
  # replace 0.0 with NA to replicate
  marchinko[ant.dorspine==0,ant.dorspine:=NA]
  marchinko[sec.dorspine==0,sec.dorspine:=NA]
  marchinko[pelspine.len==0,pelspine.len:=NA]
  marchinko[pelgirdle.len==0,pelgirdle.len:=NA]
  
  # *** ant.dorspine family 6 does not replicate
  marchinko[, ant.dorspine.s:=size_corrected(ant.dorspine, standlength), by=family]
  marchinko[, sec.dorspine.s:=size_corrected(sec.dorspine, standlength), by=family]
  marchinko[, pelspine.len.s:=size_corrected(pelspine.len, standlength), by=family]
  marchinko[, pelgirdle.len.s:=size_corrected(pelgirdle.len, standlength), by=family]
  
  working_table <- marchinko[, .(standlength=mean(standlength, na.rm=TRUE),
                                 ant.dorspine.s=mean(ant.dorspine.s, na.rm=TRUE),
                                 sec.dorspine.s=mean(sec.dorspine.s, na.rm=TRUE),
                                 pelspine.len.s=mean(pelspine.len.s, na.rm=TRUE),
                                 pelgirdle.len.s=mean(pelgirdle.len.s, na.rm=TRUE),
                                 sd.standlength=sd(standlength, na.rm=TRUE),
                                 sd.ant.dorspine.s=sd(ant.dorspine.s, na.rm=TRUE),
                                 sd.sec.dorspine.s=sd(sec.dorspine.s, na.rm=TRUE),
                                 sd.pelspine.len.s=sd(pelspine.len.s, na.rm=TRUE),
                                 sd.pelgirdle.len.s=sd(pelgirdle.len.s, na.rm=TRUE)
  ), by=.(family, treatment)]
  
  working_table.no <- working_table[treatment=='no']
  working_table.pred <- working_table[treatment=='pred']
  mucols <- c('standlength', 'ant.dorspine.s', 'sec.dorspine.s', 'pelspine.len.s', 'pelgirdle.len.s')
  sdcols <- c('sd.standlength', 'sd.ant.dorspine.s', 'sd.sec.dorspine.s', 'sd.pelspine.len.s', 'sd.pelgirdle.len.s')
  
  # checks - ant.dorspine family 6 does not replicate
  round(working_table.no[, .SD, .SDcols=mucols], 3)
  round(working_table.pred[, .SD, .SDcols=mucols], 3)
  y <- marchinko[family=='6', ant.dorspine]
  sl <- marchinko[family=='6', standlength]
  y[y==0] <- NA
  dt <- na.omit(data.table(y=y, sl=sl, treatment=marchinko[family=='6',treatment]))
  dt[, y.s := residuals(lm(y~sl))]
  dt[, .(mean=mean(y.s)), by=treatment]
  
  # Marchinko table 1
  diff_table <- (working_table.pred[, .SD, .SDcols=mucols] - working_table.no[, .SD, .SDcols=mucols])/working_table.no[, .SD, .SDcols=sdcols]
  table_1 <- apply(diff_table, 2, median)
  return(table_1)
}

size_corrected <- function(x, sl, centered=TRUE){
  # size correction via residual of regression of x on standard length sl.
  # adds the residual to mean(x) so that value is centered at original
  # centered = TRUE returns the residuals, FALSE returns residuals + mean(x)
  # to correct within groups, simply use by=.(grouping) in the data.table
  # e.g.
  # marchinko[, ant.dorspine.s:=size_corrected(ant.dorspine, standlength), by=family]
  
  x.s <- residuals(lm(x ~ sl, na.action = na.exclude))
  if(centered==FALSE){
    x.s <- x.s + mean(x)
  }
  return(x.s)
}