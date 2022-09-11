aicm = function(x,y,nsamp,x.includes.intercept=FALSE,m.fixed=1) {
########################################
### COMPUTE AICM/AICR/AICC/AIC FOR Y = XB + E
### M.FIXED: NUMBER OF FIXED/SAME PREDICTORS.  DEFAULT = 1 FOR INTERCEPT
########################################
if (length(y) %% nsamp != 0) stop('y dimensioned incorrectly')
if (length(x) %% nsamp != 0) stop('y dimensioned incorrectly')
ydim = length(y) / nsamp
xdim = length(x) / nsamp

xmat = x
if (!x.includes.intercept) {
	xmat = cbind(rep(1,nsamp),xmat)
	mdim = dim(xmat)[2]
} else {
	mdim = xdim
}

m.random = mdim - m.fixed

x.svd     = svd(xmat)
dum       = t(x.svd$u) %*% y
cov.noise = t(y) %*% y - t(dum) %*% dum
cov.noise = as.matrix(cov.noise/nsamp)

aicr      = log(det(cov.noise)) + ydim * (nsamp+1) * ( 1 + (mdim-1)/(nsamp-mdim-1))/(nsamp-mdim-ydim-1)
aicc      = log(det(cov.noise)) + ydim * (nsamp+mdim)/(nsamp-mdim-ydim-1)
# aic       = log(det(cov.noise)) + ydim * (nsamp+mdim)/nsamp
aic       = log(det(cov.noise)) + ydim + (2 * ydim * mdim + ydim*(ydim+1))/nsamp
aicm      = log(det(cov.noise)) + ydim * (nsamp + m.fixed) / (nsamp - mdim - ydim - 1) * ( 1 + m.random/(nsamp - mdim - 1))

# aicr      = aicr + ydim + ydim * log(2*pi)
# aicc      = aicc + ydim + ydim * log(2*pi)
# aic       = aic  + ydim + ydim * log(2*pi)

aicr      = aicr * nsamp 
aic       = aic  * nsamp
aicc      = aicc * nsamp
aicm      = aicm * nsamp

list(aicr=aicr,aicc=aicc,aic=aic,aicm=aicm)
	
}