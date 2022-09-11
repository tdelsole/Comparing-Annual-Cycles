diff.ar.cycle = function(ts1,ts2,p.order,nharm,period=12,nens1=1,nens2=1,test.equal.intercept=FALSE,first.step=c(1,1),alpha=0.05) {
### TESTS EQUALITY OF AR(P) PLUS ANNUAL/DIURNAL CYCLE MODEL

if (nharm != 0 & p.order != 0 ) {
	x.break       = c(p.order,2*nharm)
	dev.row.names = c('D[0:1]; equal variance','D[1:2]; equal AR model','D[2:3]; equal cycle')
} else if (nharm == 0) {
	x.break       = p.order
	dev.row.names = c('D[0:1]; equal variance','D[1:2]; equal AR model')
} else if (p.order == 0) {
	x.break       = 2*nharm
	dev.row.names = c('D[0:1]; equal variance','D[1:2]; equal cycle')
} else stop('unexpected logic')

if (test.equal.intercept) {
	x.break       = c(x.break,1)
	dev.row.names = c(dev.row.names,paste('D[',length(dev.row.names),':',length(dev.row.names)+1,']; intercept',sep=''))
}
dev.row.names = c(dev.row.names,paste('D[0:',length(dev.row.names),']; total',sep=''))


ts1     = as.matrix(ts1)
ts1.lhs = NULL
ts1.lag = NULL
ts1.cyc = NULL
ts1.jvc = NULL
for (ne in 1:nens1) {
	ts1.list  = timeseries2ar.cycle(ts1[,ne],p.order,nharm,period,first.step=first.step[1])
	if (p.order > 0) ts1.lag = rbind(ts1.lag,as.matrix(ts1.list$y.lag))
	if (nharm   > 0) ts1.cyc = rbind(ts1.cyc,as.matrix(ts1.list$y.cyc))
	ts1.jvc  = rbind(ts1.jvc,as.matrix(ts1.list$jvec))
	ts1.lhs  = rbind(ts1.lhs,as.matrix(ts1.list$y.lhs))
}

ts2 = as.matrix(ts2)
ts2.lag   = NULL
ts2.cyc   = NULL
ts2.jvc   = NULL
ts2.lhs   = NULL
for (ne in 1:nens2) {
	ts2.list  = timeseries2ar.cycle(ts2[,ne],p.order,nharm,period,first.step=first.step[2])
	if (p.order > 0) ts2.lag  = rbind(ts2.lag,as.matrix(ts2.list$y.lag))
	if (nharm   > 0) ts2.cyc  = rbind(ts2.cyc,as.matrix(ts2.list$y.cyc))
	ts2.jvc  = rbind(ts2.jvc,as.matrix(ts2.list$jvec))
	ts2.lhs  = rbind(ts2.lhs,as.matrix(ts2.list$y.lhs))
}

# ts1.list  = timeseries2ar.cycle(ts1,p.order,nharm,period,first.step=first.step[1])
# ts2.list  = timeseries2ar.cycle(ts2,p.order,nharm,period,first.step=first.step[2])

diff.list = diff.regression.nested(ts1.lhs,ts2.lhs,
   cbind(ts1.lag,ts1.cyc,ts1.jvc),
   cbind(ts2.lag,ts2.cyc,ts2.jvc),
   x.break=x.break,alpha=alpha)
   
rownames(diff.list$dev.table) = dev.row.names

aicm1 = aicm(ts1.lhs,cbind(ts1.lag,ts1.cyc),length(ts1.lhs),m.fixed=nharm)
aicm2 = aicm(ts2.lhs,cbind(ts2.lag,ts2.cyc),length(ts2.lhs),m.fixed=nharm)


### COMPUTE FORCED RESPONSE
forced.cycle      = as.numeric(rep(0,period))
forced.cycle.star = as.numeric(rep(0,period))
if (nharm != 0) {
	phi           = coefficients(diff.list$lm1)[1:p.order ]
	phi.star      = coefficients(diff.list$lm2)[1:p.order ]
	beta          = coefficients(diff.list$lm1)[  p.order + 1:(2*nharm) ]
	beta.star     = coefficients(diff.list$lm2)[  p.order + 1:(2*nharm) ]
	np            = 1:period
	for (nh in 1:nharm) {
		forced.cycle      = forced.cycle      + Re((beta     [2*nh-1] - 1i * beta     [2*nh])/(1-sum(phi      * exp(-1i * 2 * pi * nh * (1:p.order) / period))) * exp(2i * pi * np * nh / period))
		forced.cycle.star = forced.cycle.star + Re((beta.star[2*nh-1] - 1i * beta.star[2*nh])/(1-sum(phi.star * exp(-1i * 2 * pi * nh * (1:p.order) / period))) * exp(2i * pi * np * nh / period))
	}	
}


list(dev.table = diff.list$dev.table, 
  lm.omega = diff.list$lm.omega, 
  lm1 = diff.list$lm1, lm2 = diff.list$lm2,
  ts1.cyc = ts1.list$y.cyc, ts2.cyc = ts2.list$y.cyc,
  ts1.lag = ts1.list$y.lag, ts2.lag = ts2.list$y.lag,
  ts1.lhs = ts1.list$y.lhs, ts2.lhs = ts2.list$y.lhs,
  forced.cycle = forced.cycle, forced.cycle.star = forced.cycle.star,
  aicm1=aicm1, aicm2=aicm2)

}