witness.est <-
function(parameters, guilty=NULL, perp.removed=NULL, data = NULL, ilen = 100, lsize=6, meth="WITC", suspectMatch=FALSE, seedit=TRUE, sample.size=1000, optim = FALSE, goodFoil=FALSE){
	
		###### parameter matrix of the following form: 
		###### rows=dataset, cols=c, ssp, ssf, criteria, and wa
		###### Fixed values across rows means that the parameters are fixed to be
		###### the same between datasets.
		
		###### Example:
		###### parameters: .2  .5  .2  .06  .5
		######			   .2  .5  .2  .03  .5
		######			   .2  .5  .2  .01  .5
		######			   .2  .5  .2  .07  .5
		######			   .2  .5  .2  .07  .5						
				
		###### the datasets are fit in such a way that c, ssp, ssf, and wa are
		###### fixed across datasets, while criteria varies
		
		###### reset random seed
		if (exists(".Random.seed", .GlobalEnv)){
			seeder = TRUE
			save.seed <- get(".Random.seed", .GlobalEnv)
		} else {
			seeder = FALSE
		}
		
		
		###### check if parameters are in form of numbers
		if (mode(parameters) == "character"){
			parameters = witness.starting.params(parameters)
		}
				
		###### make sure user understands default of guilty vector
		if (is.null(guilty)){
			guilty = rep(c(T,F), length.out=nrow(parameters))
			if (!optim){
				cat("Note: You failed to specify a guilty vector. We will assume that every 
				even row is 'Target Absent' and every odd row is 'Target Present'.\n")
			}
		}
			
		p = nrow(parameters)
		
		##### don't bark if perp removed is null, just reformulate
		if (is.null(perp.removed)){
			perp.removed = rep(FALSE, times=p)
		} 

		### set up useless data vector
		if (is.null(data)) data=matrix(rep(c(0,0,1), times=p), nrow=p, byrow=T)
	
		
		###### preallocate response vectors
		responses = matrix(0, nrow=nrow(data), ncol=3)
		
		#parameters = matrix(parameters, ncol=5, byrow=TRUE)
		###### designate parameters
		if (meth=="WITC"){
			parameters[,5] = 1
		}
		a = parameters[,1]
		ssp = parameters[,2]
		ssf = parameters[,3]
		c = parameters[,4]
		wa = parameters[,5]
		

		###### begin loop of responses
		for (i in 1:sample.size){
			
			###### see if results need to be repeatable
			if (seedit){
				set.seed(i)
			}
			
			####### create feature vectors
			perp = matrix(runif(ilen*p, -1, 1), nrow=p, ncol=ilen)
			mem = matrix(runif(ilen*p, -1, 1), nrow=p, ncol=ilen)
						
			###### match according to memory (for each condition of a)
			random = matrix(runif(ilen*p, 0, 1), nrow=p, ncol=ilen)
			mem[random<a] = perp[random<a]

			#### match innocent suspect to perpetrator
			##### only do it where there's an inn suspect on dataset
			inn.rows = which(!guilty)	
			guilt.rows = which(guilty)	
			inctsp = matrix(runif(ilen*p, -1, 1), 
						nrow= p, ncol=ilen)
			random = matrix(runif(ilen*p, 0, 1), 
						nrow= p, ncol=ilen)
			ssp = ssp[-inn.rows]
			inctsp[random <ssp] = perp[random < ssp]

			##### lineup position array (rows of data, lineup position,			
			##### make foil vectors
			lineup = array(data=runif(p*ilen*lsize, -1, 1), c(p,lsize,ilen), 
			c("Dataset.Num", "Position", "Features"))
			for (k in 1:lsize){
				random = matrix(runif(ilen*p, 0, 1), nrow=p, ncol=ilen)
				##### if there's a good foil, make their ssf higher
				if (goodFoil & k==3){
					lineup[,k,][random < (ssf+.15)] = perp[random < (ssf+.15)]
				}
				
				if (!suspectMatch){
					lineup[,k,][random < ssf] = perp[random < ssf]
				} else {
					lineup[,k,][random < ssf] = inctsp[random < ssf]
				}
			}		
			
			###### put perps in position one
			lineup[guilt.rows,1,] = perp[guilt.rows,]
			
			###### put innocent suspects in position
			lineup[inn.rows,1,] = inctsp[inn.rows,]
			
			###### remove perp where applicable
			if (length(which(perp.removed))>0){
				lineup[which(perp.removed), 1,] = NA
			}
				
			####### create a vector with match scores
			matchVect = matrix(nrow=p, ncol=lsize)
			for (k in 1:lsize){
				matchVect[,k] = diag((mem) %*% t(lineup[,k,]))/ilen
			}
			
			####### find best/second best match for each dataset
			best.match = rep(NA, times=p)
			best.match.position = rep(NA, times=p)			
			second.best = rep(NA, times=p)
			second.best.position = rep(NA, times=p)
			
			best.match = apply(matchVect, 1, max, na.rm=TRUE)
			best.match.position = apply(matchVect, 1, match.position, pos=1)
			
			second.best = apply(matchVect, 1, match.position, pos=2, position=FALSE)
			second.best.position = apply(matchVect, 1, match.position, pos=2, position=TRUE)
			
			####### weight the decision
			wr= 1-wa
			evaluation = wr*(best.match-second.best) + (wa * best.match)
			
			####### check if best match exceeds criteria
			makeid = evaluation > c
	
			####### increment "fail to ID" if false
			responses[which(!makeid),3] = responses[which(!makeid),3] + 1
			
			####### incrememnt "suspect ID" if true and in position 1
			condition = which(makeid & best.match.position==1)
			responses[condition,1] = responses[condition,1]+1
			
			####### incrememnt "False ID" if true and not in position 1
			condition = which((makeid) & best.match.position!=1)
			responses[condition,2] = responses[condition,2]+1
		}	
		
		#### reset random seed
		if (seeder){
			assign(".Random.seed", save.seed, .GlobalEnv)
		}

		####### summarize results
		responses = responses/sample.size				
		
		####### compare to actual dataset (chi square)
		sse= sum((data-responses)^2)
		rmse = sqrt(sse/length(responses))

		
		
		####### report results
		parameters = data.frame(parameters)
		names(parameters) = c("a", "ssp", "sfs", "c", "wa")
		responses=data.frame(responses)
		names(responses) = c("Suspect ID", "Foil ID", "No ID")
		
		if (optim){
			return(rmse)
		} else {
			list(parameters=parameters, predicted = responses)
		}	
		
		
}
