rm(list=ls())

nsamp      = 50					# length of 1st time series
nsamp.star = 55					# length of 2nd time series
nharm      = 3					# number of harmonics
phi        = c(0.9,-0.3,0.0)	# AP(p) parameters for 1st time series
phi.star   = c(0.9,-0.3,0.0)	# AP(p) parameters for 2nd time series
mu         = 0	 				# intercept for 1st time series
mu.star    = 0					# intercept for 2nd time series

set.seed(3)

harm       = rnorm(2*nharm)		# annual cycle parameters for 1st time series
harm.star  = harm + rnorm(2*nharm,sd=0)	# annual cycle parameters for 2nd time series
period     = 12					# period of the cycle

dir.Rlib = '/Users/delsole/R/delsole_tools/'	# directory of R functions
source(paste(dir.Rlib,'timeseries2ar.cycle.R',sep=''))
source(paste(dir.Rlib,'diff.regression.nested.R',sep=''))
source(paste(dir.Rlib,'diff.ar.cycle.R',sep=''))
source(paste(dir.Rlib,'aicm.R',sep=''))


########################################
####### METADATA AND SANITY CHECKS
########################################
p.order = length(phi)
if (p.order != length(phi.star)) stop('the number of AR parameters should be the same')

if (length(harm)      != 2*nharm) stop('number of annual cycle parameters should be 2*nharm')
if (length(harm.star) != 2*nharm) stop('number of annual cycle parameters should be 2*nharm')

if (2*nharm > period) stop('number of harmonics cannot exceed nyquist frequency')


########################################
####### GENERATE SYNTHETIC DATA
########################################
ts1 = numeric(nsamp)
ts2 = numeric(nsamp.star)

y.cyc       = NULL
y.cyc.star  = NULL
ncyc        = 1:nsamp
ncyc.star   = 1:nsamp.star
for (nh in 1:nharm) y.cyc      = cbind(y.cyc     ,cos(2*pi*nh*ncyc     /period),sin(2*pi*nh*ncyc     /period))
for (nh in 1:nharm) y.cyc.star = cbind(y.cyc.star,cos(2*pi*nh*ncyc.star/period),sin(2*pi*nh*ncyc.star/period))

for (n in (1+p.order):nsamp     ) ts1[n] = sum(phi      * ts1[n-1:p.order]) + sum(y.cyc     [n,] * harm     ) + rnorm(1) + mu
for (n in (1+p.order):nsamp.star) ts2[n] = sum(phi.star * ts2[n-1:p.order]) + sum(y.cyc.star[n,] * harm.star) + rnorm(1) + mu.star

yrange = range(ts1,ts2)
plot(ts1,type='l',xlab='time',ylab='',ylim=yrange,lwd=2)
lines(ts2,col='red',lwd=2)

########################################
####### TEST EQUALITY OF PROCESSES
########################################
diff.list = diff.ar.cycle(ts1,ts2,p.order,nharm,test.equal.intercept=FALSE)
print(diff.list$dev.table)

# source(paste(dir.Rlib,'timeseries2ar.cycle.old1.R',sep=''))
# source(paste(dir.Rlib,'diff.regression.nested.old1.R',sep=''))
# source(paste(dir.Rlib,'diff.ar.cycle.old1.R',sep=''))
# diff.list.old1 = diff.ar.cycle.old1(ts1,ts2,p.order,nharm)
# print(diff.list.old1$dev.table)
# print(diff.list$dev.table - diff.list.old1$dev.table)

########################################
####### COMPARE TO ANOVA FUNCTION
########################################
# for (nb in 2:length(diff.list$lm.omega)) print(anova(diff.list$lm.omega[[nb]],diff.list$lm.omega[[nb-1]]))

# if (length(diff.list$lm.omega) >= 3 ) print(anova(diff.list$lm.omega[[3]],diff.list$lm.omega[[2]],diff.list$lm.omega[[1]],test="F"))

