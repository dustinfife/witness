fitOpWit <-
function(param.values, param.form, dataset, ...){
	if (length(which(param.values>1 | param.values<0))){
		return(NA)
	} else {
		witness.able = witness.starting.params(param.form, param.values)
		fit = witness.est(witness.able, optim=TRUE, data=dataset, ...)
		return(fit)
	}
}
