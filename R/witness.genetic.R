#'Fit the WITNESS model using a Genetic Algorithm
#'
#'After generating arbitrary starting values, the genetic algorithm selects the
#'best parameter combinations to "mate." At each iteration, mutations are added
#'with a specified probability, as well as totally new parameter values.
#'
#'
#'@param parameter.form a matrix where the rows are the experiments, and the
#'columns are the parameters a, ssp, sfs, c, and wa. To fix parameters across
#'rows, one would simply input the same values. For example, if one wanted all
#'the encoding values (a) to be the same, the first row could consist of only
#'ones.
#'@param data a matrix that has the same number of rows as parameter.form.
#'Again, the rows are the experiments but the columns correspond to 1. Target
#'Identification, 2.  Foil Identification, and 3. No Identification.
#'@param N the number of iterations the WITNESS model will do to simulate the
#'eyewitness procedure.
#'@param generations how many generations the genetic algorithm will iterate
#'through before quitting.
#'@param pcross probability of cross "breeding," or the probability that a
#'parent will breed.
#'@param mutProb probability of randomly mutating a parameter (i.e., the
#'probability that noise will be added to the next generation)
#'@param prop.random to avoid local minima, the algorithm allows for "aliens"
#'to be introduced (or completely new parameter values).
#'@param gradient convergence criteria. Currently not implemented.
#'@param ... other parameters passed to the WITNESS model.
#'@seealso \code{\link{witness.optim}}, \code{\link{witness.est}}
#'@references Clark, S. E. (2003). A memory and decision model for eyewitness
#'identification. Applied Cognitive Psychology, 17, 629-654.
#'@examples
#'
#'dataMatrix = matrix(c(.471, .230, .350,
#'							.208, .137, .513,
#'							.396, .431, .242,
#'							.166, .081, .669), nrow=4, byrow=TRUE)
#'# create an parameter.form
#'params = matrix(c(rep("e", times=4),
#'		0, 0, 0, 0,
#'		"sfs1", "sfs1", "sfs2", "sfs2",
#'		"cr1", "cr2", "cr1", "cr2",
#'		1, 1, 1, 1), nrow=4)
#'
#'		####### do genetic algorithm (commented b/c it takes a while)
#'#fit = witness.genetic(params, data=dataMatrix, N=10, generations=1, sample.size=100, meth="WITC")
#'@export
witness.genetic <-function(parameter.form, data=NULL, 
				N=1000, generations=100, pcross=.8, 
				mutProb=.001, prop.random=.001, gradient=.01,...){
			
	###### generate starting parameters N times
	parent.array = array(dim=c(dim(parameter.form), N)) 
	fit.matrix = 1:N

	####### generate initial fits
	for (i in 1:N){
		parent.array[,,i] = witness.starting.params(parameter.form)
		fit.matrix[i] = witness.est(parent.array[,,i], data=data, optim=TRUE,...)
	}	

		
	######## record best ever fits
	bestEv = min(fit.matrix)
	bestEv.par = parent.array[,,which(fit.matrix==min(fit.matrix))]
	k=0
	
	######## optimize the best ever fits
	opfit = witness.optim(bestEv.par, data.set=data, ...)
	print(opfit)	
	bestEv = opfit$fit
	bestEv.par = witness.starting.params(parameter.form, values=opfit$par)
	parent.array[,,which(fit.matrix==min(fit.matrix))] = bestEv.par

	
		print(paste("optimal fit:", round(min(fit.matrix), digits=4)))
		print(paste("average fit:", round(mean(fit.matrix), digits=4)))	
		print(paste("generation:", k))	
	
	######## start generations loop
	while(k<generations){
		k = k + 1

		##### come up with PDF based on Fits
		fit.matrix = 1/fit.matrix
		pdf = fit.matrix/sum(fit.matrix)
		
		##### sample based on PDF
		samp = sample(1:N, replace=TRUE, size=N, prob=pdf)
		
		##### copy first generation
		nxtgen = parent.array[,,samp]
		
		#### randomly select parents
		randDraw = runif(N)
		parents = parentsFun(N, pcross*length(which(randDraw<pcross))) #### randomly match parents
		
		##### randomly pair parameters
		alleles = array(round(runif(nrow(parents)*5*nrow(parameter.form),0,1)), dim=c(nrow(parameter.form), 5, nrow(parents)))
		nxtgen[,,1:nrow(parents)] = nxtgen[,,parents[,1]] * alleles + nxtgen[,,parents[,1]] * (1-alleles)
		
		###### randomly mutate 
		randDraw = runif(N)
		if (length(which(randDraw<mutProb))>0){
			rand.mut.amt = runif(length(which(randDraw<mutProb)), -.2, .2)
			nxtgen[,,randDraw<mutProb] = nxtgen[,,randDraw<mutProb] + rand.mut.amt
		}
		
		###### randomly replace one and optimize
		samp = sample(1:length(N), 1)
		nxtgen[,,samp] = witness.starting.params(parameter.form)
		opfit = witness.optim(nxtgen[,,samp], data.set=data, ...)
		nxtgen[,,samp] = witness.starting.params(parameter.form, values=opfit$par)	
		
		###### now refit
		for (i in 1:N){
			parent.array[,,i] = nxtgen[,,i]
			fit.matrix[i] = witness.est(parent.array[,,i], data=data, optim=TRUE, ...)
		}	
		fit.matrix[samp]
			

		###### see if new fit surpasses previous fit
		new.Best = min(fit.matrix)
		print(new.Best)
	print(paste("3:", bestEv))
		if (new.Best<bestEv){
			bestEv = new.Best
			bestEv.par = nxtgen[,,which(fit.matrix==min(fit.matrix))]
		}
		
		####### make this generation the new parents
		parent.array = nxtgen
		
		####### report best fit
		cat(paste("\n\noptimal fit:", round(bestEv, digits=4)))
		cat(paste("\naverage fit:", round(mean(fit.matrix), digits=4)))	
		cat(paste("\ngeneration:", k))			

	}
	list(fit=bestEv, params=bestEv.par)
	

}
