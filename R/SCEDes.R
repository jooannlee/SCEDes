# install.packages("devtools")
# install.packages("roxygen2")
# library("devtools")
# install_github("devtools", "hadley")
# load_all()

#'Effect sizes for single case reversal designs.
#'
#' Calculates various effect sizes for single case reversal designs.
#' @param dat A data frame with 3 columns (time, phase, y). For an example, see simulated data below:
#'
#' ### example dataset
#'
#' set.seed(135)
#'
#' ## simulation parameters
#'
#' n = c(5,5,5,5)
#'
#' phase = data.frame(phase=c(rep('A1', n[1]), rep('B1', n[2]),
#'                            rep('A2', n[3]), rep('B2', n[4])))
#'
#' dummy = model.matrix( ~0 + phase, model.frame(phase))[,-1L]
#'
#' phi = 0.2
#'
#' b = c(0, 0, 0, 0.3, 0.3, 0, 0, 0)
#'
#' x = c(0,0,0) # initialize x
#'
#' y = numeric(sum(n)) # initialize y
#'
#' MA1 = 0  # MA1 start value
#'
#'
#' ## simulate data
#'
#' set.seed(135)
#'
#' # start data generating model loop
#'
#' for(t in -1000:sum(n)){
#'
#'  if(t > 0) x = dummy[t, ]
#'
#'  u = rnorm(1, mean = 0, sd = 1)
#'
#'  epsilon = (phi*MA1) + u
#'
#'  yVal =  b[1] + (b[2]*t) + (b[3]*x[1]) + (b[4]*x[2]) + (b[5]*x[3]) +
#'  (b[6]*(t - (n[1] + n[2] + 1))*x[1]) +
#'    (b[7]*(t - (n[1] + 1))*x[2]) +
#'    (b[8]*(t - (n[1] + n[2] + n[3] + 1))*x[3]) +
#'    epsilon
#'
#'  MA1 = epsilon
#'
#'  if(t > 0) y[t] = yVal
#'
#'  if(t > 0) x = dummy[t, ]
#'
#'} # end data generating model loop
#'
#'dat <- data.frame(time = 1:sum(n), phase=phase, y=y)
#'
#'head(dat)
#'
#' @return Calculated various effect sizes for single case reversal designs:
#'
#' R2 = R-squared from multiple linear regression model,
#'
#' Beta12.mean = average of partial regression coefficients for phases A1B1 and A2B2,
#'
#' PND12.mean = average of Percentage Non Overlap statistics for phases A1B1 and A2B2,
#'
#' d12.mean = average of standardized mean difference for phases A1B1 and A2B2,
#'
#' dp12.mean = average of pooled standardized mean difference for phases A1B1 and A2B2,
#'
#' R2AG12.mean = average of Allison and Gorman's detrended R-squared for phases A1B1 and A2B2.
#'
#' @export


SCEDes <- function(dat){

  ### plot
  #windows()
  plot(dat$time,dat$y,xlab="Time", ylab="Value", pch=20); abline(v=c(length(dat$phase[phase=="A1"]),
                                                                     length(dat$phase[phase=="A1" | phase=="B1"]),
                                                                     length(dat$phase[phase=="A1" | phase=="B1" | phase=="A2"]),
                                                                     length(dat$phase[phase=="A1" | phase=="B1" | phase=="A2" | phase=="B2"])))


  ### effect sizes
  # PND
  maxA1 = max(dat$y[phase=="A1"])
  maxA2 = max(dat$y[phase=="A2"])
  PND1 = (sum(as.logical(which(dat$y[phase=="B1"] > maxA1)))/sum(phase=="B1"))*100
  PND2 = (sum(as.logical(which(dat$y[phase=="B2"] > maxA2)))/sum(phase=="B2"))*100
  PND3 = (sum(as.logical(which(dat$y[phase=="B2"] > maxA1)))/sum(phase=="B2"))*100
  PND12.mean = mean(PND1,PND2) # average of phase B1 compared to A1, B2 to A1
  PND13.mean = mean(PND1,PND3) # average of phase B1 compared to A1, B2 to A2

  # R2
  mod = lm(dat$y ~ dat$time + dat$phase + dat$time*dat$phase)
  R2 = summary(mod)$r.squared

  # partial regression coefficients
  Beta = summary(mod)$coefficients[4:5,1]
  Beta12.mean = mean(Beta)

  # SMD | prep
  mu <- tapply(dat$y, dat$phase, mean)
  SD <- tapply(dat$y[phase=="A1" | phase=="A2"], dat$phase[phase=="A1" | phase=="A2"], sd)

  SD2 <- tapply(dat$y, dat$phase, sd)
  n1<- length(dat$phase[phase=="A1"])
  n2<- length(dat$phase[phase=="A2"])
  n3<- length(dat$phase[phase=="B1"])
  n4<- length(dat$phase[phase=="B2"])
  SDpooled<- (((n1-1)*SD2[1])+((n2-1)*SD2[2])+((n3-1)*SD2[3])+((n4-1)*SD2[4]))/((n1+n2+n3+n4)-(length(unique(dat$phase))))

  # SMD
  d11 <- (mu[3] - mu[1]) / SD[1]
  d12 <- (mu[3] - mu[1]) / SD[2]

  d21 <- (mu[4] - mu[2]) / SD[1]
  d22 <- (mu[4] - mu[2]) / SD[2]

  d31 <- (mu[4] - mu[1]) / SD[1]
  d32 <- (mu[4] - mu[1]) / SD[2]

  d11p <- (mu[3] - mu[1]) / SDpooled
  d21p <- (mu[4] - mu[2]) / SDpooled
  d31p <- (mu[4] - mu[1]) / SDpooled

  d12.mean = mean(d11,d21) # average of phase B1 compared to A1, B2 to A1
  dp12.mean = mean(d11p,d21p) # average of phase B1 compared to A1, B2 to A1 with pooled std deviation

  # allison and gorman R2 | first set of AB
  fitag <- lm(dat$y[phase=="A1"] ~ dat$time[phase=="A1"])
  residfitag <- resid(fitag)
  p <- predict(fitag)

  newdata <- data.frame(time = dat$time[phase=="B1"])
  pfitag <- predict.lm(fitag, newdata)
  residpfitag <- dat$y[phase=="B1"] - pfitag

  ydetrended <- as.numeric(c(residfitag,residpfitag))

  dummy <- data.frame(model.matrix(mod))
  dummy <- cbind(dat$phase,dummy)

  xtag <- (dummy$dat.time.dat.phaseB1[phase=="A1" | phase=="B1"])
  cor1 <- cor(dummy$dat.phaseB1[phase=="A1" | phase=="B1"],ydetrended)
  cor2 <- cor(xtag,ydetrended)
  signcor1 <- sign(cor1)
  signcor2 <- sign(cor2)

  if (signcor1==1 & signcor2==1) {
    fitagsign <- lm(ydetrended ~ dummy$dat.phaseB1[phase=="A1" | phase=="B1"] + dummy$dat.time.dat.phaseB1[phase=="A1" | phase=="B1"])
  } else {
    fitagsign <- lm(ydetrended ~ dummy$dat.phaseB1[phase=="A1" | phase=="B1"])
  }

  AGrs1 <- summary(fitagsign)$r.squared
  AGrsadj1 <- summary(fitagsign)$adj.r.squared

  # allison and gorman R2 | second set of AB
  fitag2 <- lm(dat$y[phase=="A2"] ~ dat$time[phase=="A2"])
  residfitag2 <- resid(fitag2)
  p2 <- predict(fitag2)

  newdata2 <- data.frame(time = dat$time[phase=="B2"])
  pfitag2 <- predict.lm(fitag2, newdata2)
  residpfitag2 <- dat$y[phase=="B2"] - pfitag2

  ydetrended2 <- as.numeric(c(residfitag2,residpfitag2))

  xtag2 <- (dummy$dat.time.dat.phaseB2[phase=="A2" | phase=="B2"])
  cor12 <- cor(dummy$dat.phaseB2[phase=="A2" | phase=="B2"],ydetrended2)
  cor22 <- cor(xtag2,ydetrended2)
  signcor12 <- sign(cor12)
  signcor22 <- sign(cor22)

  if (signcor12==1 & signcor22==1) {
    fitagsign2 <- lm(ydetrended2 ~ dummy$dat.phaseB2[phase=="A2" | phase=="B2"] + dummy$dat.time.dat.phaseB2[phase=="A2" | phase=="B2"])
  } else {
    fitagsign2 <- lm(ydetrended2 ~ dummy$dat.phaseB2[phase=="A2" | phase=="B2"])
  }

  AGrs2 <- summary(fitagsign2)$r.squared
  AGrsadj2 <- summary(fitagsign2)$adj.r.squared

  R2AG12.mean = mean(c(AGrsadj1,AGrsadj2)) # average of phase B1 compared to A1, B2 to A1
  R2AGadj12.mean = mean(c(AGrsadj1,AGrsadj2)) # average of phase B1 compared to A1, B2 to A1

  final = rbind(R2,Beta12.mean,PND12.mean,d12.mean,dp12.mean,R2AG12.mean)
  final = data.frame(final)
  colnames(final) <- c("Effect Size Value")

  return(final)

}
