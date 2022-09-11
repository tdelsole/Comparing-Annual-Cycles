diff.regression.nested = function(yvec,yvec.star,xmat,xmat.star,x.break,alpha=0.05,ntrials=10000) {
###### TEST EQUALITY OF REGRESSION MODELS THROUGH NESTED HYPOTHESES:
###### YVEC = X2  B2  + ... + XL  BL  + E
###### YVEC = X2* B2* + ... + XL* BL* + E*
###### IMPORTANT: USER SHOULD INCLUDE INTERCEPT IN XMAT AND XMAT.STAR
###### HYPOTHESES: 
######	OMEGA.0: UNRESTRICTED
######  OMEGA.1: EQUALITY OF NOISE VARIANCE (VAR[E] = VAR[E*]):  X.COMM[1] = 0
######  OMEGA.2: OMEGA.1 AND B2 = B2*; LENGTH(B2) = X.BREAK[1] = X.COMM[2]
######  OMEGA.3: OMEGA.2 AND B3 = B3*; LENGTH(B3) = X.BREAK[2] = X.COMM[3]
######  ...
###### 	OMEGA.L: OMEGA.L-1 AND BL = BL*; TESTED ONLY IF X.BREAK[L-1] IS SPECIFIED
###### INPUT:
###### 	YVEC     [NTOT]: 			Y      TIMESERIES
###### 	YVEC.STAR[NTOT.STAR]: 		Y.STAR TIMESERIES
######	XMAT     [NTOT,KTOT]:		MATRIX [X2 ,X3 ,...,XL ]
######	XMAT.STAR[NTOT,KTOT.STAR]:	MATRIX [X2*,X3*,...,XL*]
######  X.BREAK  [LTOT]:			NUMBER OF COLUMNS IN X2, X3, ... XL; SUM(X.BREAK) <= KTOT
######	ALPHA:						SIGNIFICANCE LEVEL; DEFAULT IS 5%
######	NTRIALS:					NUMBER OF MONTE CARLO TRIALS FOR ESTIMATING DEVIANCE [OMEGA.0 - OMEGA.1]


##############################
###### METADATA
##############################
x.comm    = c(0,x.break) # append 0 at the start, corresponding to no common predictors
xmat      = as.matrix(xmat)
xmat.star = as.matrix(xmat.star)

ncomm     = length(x.comm)
ntot      = length(yvec)
ntot.star = length(yvec.star)
if (sum(x.comm) > dim(xmat)     [2]) stop('x.comm inconsistent with xmat')
if (sum(x.comm) > dim(xmat.star)[2]) stop('x.comm inconsistent with xmat.star')

alpha.step  = 1 - (1-alpha)^(1/ncomm)



########################################
####### EQUALITY OF VARIANCES (OMEGA_1)
########################################
lm1 = lm(yvec      ~ xmat       - 1)
lm2 = lm(yvec.star ~ xmat.star  - 1)
sse1 = sum(residuals(lm1)^2)
sse2 = sum(residuals(lm2)^2)
dof1 = lm1$df.residual
dof2 = lm2$df.residual

########################################
####### TEST NESTED HYPOTHESES
########################################
y.all         = c(yvec,yvec.star)
sse.omega     = as.numeric(rep(NA,ncomm))
dof.omega     = as.numeric(rep(NA,ncomm))
dev.omega     = as.numeric(rep(NA,ncomm))
dev.crit      = as.numeric(rep(NA,ncomm))
fval          = as.numeric(rep(NA,ncomm))
pval          = as.numeric(rep(NA,ncomm))
param.num     = as.numeric(rep(NA,ncomm))
dev.crit.asym = as.numeric(rep(NA,ncomm))
pval.asym     = as.numeric(rep(NA,ncomm))

for (nb in 1:ncomm) {
	ncol.comm      = sum(x.comm[1:nb])
	ncol.diff      = ncol(xmat)      - ncol.comm
	ncol.diff.star = ncol(xmat.star) - ncol.comm
	
	if (ncol.comm == 0) {
		x.all = NULL
	} else {
		x.all  = rbind(xmat[,1:ncol.comm,drop=FALSE],xmat.star[,1:ncol.comm,drop=FALSE])
	}
	
	if (ncol.diff > 0) {
		x.diff      = cbind(xmat     [,1:ncol.diff      + ncol.comm])
		x.all       = cbind(x.all,rbind(x.diff,array(0,dim=c(ntot.star,ncol.diff.star))))
	}
	if (ncol.diff.star > 0) {
		x.diff.star = cbind(xmat.star[,1:ncol.diff.star + ncol.comm])		
		x.all       = cbind(x.all,rbind(array(0,dim=c(ntot,ncol.diff)),x.diff.star))
	}
	
	lm.nest = lm(y.all ~ x.all - 1)
	sse.omega[nb] = sum(residuals(lm.nest)^2)
	dof.omega[nb] = lm.nest$df.residual
	if (nb == 1) lm.omega = list(lm.nest) else lm.omega[nb] = list(lm.nest)
	param.num[nb] = ncol(x.all) + 1
}

#### DEVIANCE 0:1
var.ratio       = (sse1/dof1) / (sse2/dof2)
var.ratio.upper = qf(alpha.step/2,dof1,dof2,lower.tail=FALSE)
var.ratio.lower = qf(alpha.step/2,dof1,dof2,lower.tail=TRUE)
fval[1]         = var.ratio
if (var.ratio > 1) {
	pval[1] = pf(  var.ratio,dof1,dof2,lower.tail=FALSE)
} else {
	pval[1] = pf(1/var.ratio,dof2,dof1,lower.tail=FALSE)
}


## MAXIMUM LIKELIHOOD (BIASED)
dev.omega.check = (ntot+ntot.star)*log(sse.omega[1]/(ntot+ntot.star)) - ntot * log(sse1/ntot) - ntot.star * log(sse2/ntot.star)
dev.omega[1]    = ntot * log(1 + sse2/sse1) + ntot.star * log(1+sse1/sse2) + ntot * log(ntot/(ntot+ntot.star)) + ntot.star * log(ntot.star/(ntot+ntot.star))
dev.omega[1]    = ntot * log(1 + sse2/sse1) + ntot.star * log(1+sse1/sse2) + ntot * log(ntot) + ntot.star * log(ntot.star) - (ntot+ntot.star) * log(ntot+ntot.star)

f.random        = rf(ntrials,dof1,dof2) * dof1 / dof2
dev.mc          = ntot * log(1 + 1/f.random) + ntot.star * log(1+f.random) + ntot * log(ntot) + ntot.star * log(ntot.star) - (ntot+ntot.star) * log(ntot+ntot.star)
dev.crit[1]     = quantile(dev.mc,prob=1-alpha.step)

#### OTHER DEVIANCES
for (nb in 2:ncomm) dev.omega[nb] = (ntot + ntot.star) * log(sse.omega[nb]/sse.omega[nb-1])
for (nb in 2:ncomm) dev.crit[nb] = 
  (ntot + ntot.star) * log( 1 + qf(alpha.step,dof.omega[nb]-dof.omega[nb-1],dof.omega[nb-1],lower.tail=FALSE) * (dof.omega[nb]-dof.omega[nb-1])/dof.omega[nb-1])
  
for (nb in 2:ncomm) fval[nb] = (sse.omega[nb] - sse.omega[nb-1])/sse.omega[nb-1] * dof.omega[nb-1]/(dof.omega[nb] - dof.omega[nb-1])
for (nb in 2:ncomm) pval[nb] = pf(fval[nb],dof.omega[nb] - dof.omega[nb-1],dof.omega[nb-1],lower.tail=FALSE)

#### ASYMPTOTICS
dev.crit.asym[1] = qchisq(alpha.step,1,lower.tail=FALSE)
for (nb in 2:ncomm) dev.crit.asym[nb] = qchisq(alpha.step,param.num[nb-1] - param.num[nb],lower.tail=FALSE)

pval.asym[1] = pchisq(dev.omega[1],1,lower.tail=FALSE)
for (nb in 2:ncomm) pval.asym[nb] = pchisq(dev.omega[nb],param.num[nb-1] - param.num[nb],lower.tail=FALSE)

#### TOTAL DEVIANCE
dev.total = (ntot + ntot.star) * log(sse.omega[ncomm]/(ntot + ntot.star)) - ntot * log(sse1/ntot) - ntot.star * log(sse2/ntot.star)
dev.total.crit = qchisq(alpha,param.num[1] + 1 - param.num[ncomm],lower.tail=FALSE)
dev.total.pval = pchisq(dev.total,param.num[1] + 1 - param.num[ncomm],lower.tail=FALSE)


########################################
#### SUMMARY TABLE
########################################
dev.table = cbind(dev.omega,dev.crit.asym,pval.asym,alpha.step,dev.crit,pval,fval)
dev.table = rbind(dev.table,c(dev.total,dev.total.crit,dev.total.pval,alpha,rep(NA,ncol(dev.table)-4)))
rownames(dev.table) = c(paste('D[',1:ncomm-1,':',1:ncomm,']',sep=''),paste('D[0:',ncomm,']',sep=''))
colnames(dev.table) = c('deviance','crit(chisq)','pval(chisq)','alpha','crit(F)','pval(F)','F-value')

########################################
#### OUTPUT THE RESULTS
########################################
list(dev.table = dev.table, lm.omega = lm.omega, lm1=lm1, lm2=lm2, alpha.step=alpha.step, 
  f.lower = var.ratio.lower, f.upper = var.ratio.upper)

}