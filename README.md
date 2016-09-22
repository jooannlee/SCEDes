# SCEDes

#This package calculates various effect sizes for single case reversal designs:
R2 = R-squared from multiple linear regression model,

Beta12.mean = average of partial regression coefficients for phases A1B1 and A2B2,

PND12.mean = average of Percentage Non Overlap statistics for phases A1B1 and A2B2,

d12.mean = average of standardized mean difference for phases A1B1 and A2B2,

dp12.mean = average of pooled standardized mean difference for phases A1B1 and A2B2,

R2AG12.mean = average of Allison and Gorman's detrended R-squared for phases A1B1 and A2B2.

##Usage
SCEDes(dat)

##Arguments
###dat	
A data frame with 3 columns (time, phase, y). For an example, see simulated data below:

#### example dataset
#### simulation parameters
n = c(5,5,5,5)

phase = data.frame(phase=c(rep('A1', n[1]), rep('B1', n[2]), rep('A2', n[3]), rep('B2', n[4])))

dummy = model.matrix( ~0 + phase, model.frame(phase))[,-1L]

phi = 0.2

b = c(0, 0, 0, 0.3, 0.3, 0, 0, 0)

x = c(0,0,0) # initialize x

y = numeric(sum(n)) # initialize y

MA1 = 0 # MA1 start value

set.seed(135)

#### #start data generating model loop

for(t in -1000:sum(n)){

if(t > 0) x = dummy[t, ]

u = rnorm(1, mean = 0, sd = 1)

epsilon = (phi*MA1) + u

yVal = b[1] + (b[2]*t) + (b[3]*x[1]) + (b[4]*x[2]) + (b[5]*x[3]) + (b[6]*(t - (n[1] + n[2] + 1))*x[1]) + (b[7]*(t - (n[1] + 1))*x[2]) + (b[8]*(t - (n[1] + n[2] + n[3] + 1))*x[3]) + epsilon

MA1 = epsilon

if(t > 0) y[t] = yVal

if(t > 0) x = dummy[t, ]}

#### #end data generating model loop

dat <- data.frame(time = 1:sum(n), phase=phase, y=y)

head(dat)
