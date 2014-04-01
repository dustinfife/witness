#'Estimate identification probabilities using WITNESS
#'
#'\code{witness} takes the parameters provided by the user then simulates the
#'WITNESS procedure and outputs probabilities
#'
#'The WITNESS model (Clark, 2003) is a direct-access matching model that has
#'been adapted for eyewitness situations. WITNESS uses numerical
#'representations of features as items in the matching process. More
#'specifically, each lineup member is represented as a vector of n arbitrary
#'features (the features do not map to literal features of a person, such as
#'eyes or a nose, however).  The value n is typically fixed to some arbitrary
#'constant such as 100.  Each position in the vector (or each feature) is
#'assigned a value between -1 and 1.
#'
#'To begin, WITNESS generates a perpetrator vector (PERP) by randomly
#'generating numbers on the interval -1, 1. This vector serves as the basis for
#'all of the subsequent lineup members. The model next "encodes" the features
#'in the PERP vector to a memory vector. The parameter a governs the degree to
#'which the Memory (MEM) vector matches the PERP vector; a is the probability
#'that each individual feature in PERP will be successfully encoded to MEM. For
#'example, if a = .4, then each feature of PERP will have a 40\% chance of
#'being represented within MEM. Conversely, there is a 1-a (60\%, in this case)
#'chance that the feature will be replaced with another arbitrary random
#'feature. Put differently, the correlation coefficient between PERP and MEM
#'will be approximately a. This results in a noisy representation in memory,
#'which accounts for phenomena such as the failure to encode details of the
#'perpetrator and post-event interference.
#'
#'WITNESS next creates the lineup members for comparison to memory. For target
#'absent lineups (i.e., an innocent suspect replaces the guilty suspect), a new
#'vector is created to represent the innocent suspect (SUSP). This vector is
#'governed by the parameter SSP, or Similarity of the Suspect to the
#'Perpetrator. SSP is the probability that each feature of SUSP will match the
#'corresponding feature of PERP. As SSP approaches 1, SUSP will be more similar
#'to PERP (and likewise, less similar as it approaches 0).
#'
#'Next, WITNESS uses one of two more parameters to create the remaining lineup
#'members (i.e., the foils), each to simulate a different method of foil
#'selection. We refer to these vectors of foils as FOILS. Foils can be selected
#'either because they match the description of the perpetrator
#'(description-matched foils) or because they match the appearance of the
#'suspect (suspect-matched foils). The parameter SFP, or Similarity of the
#'Foils to the Perpetrator, simulates a description-matched lineup. SFP is the
#'probability that a given feature of a given foil will match a corresponding
#'feature in PERP (similar to SSP, where 0 yields no shared features and 1
#'replicates PERP). Note that this results in the same foils being used in
#'target-present and target-absent lineups. SFS, Similarity of the Foils to the
#'Suspect, is a similar parameter; however, it uses the suspect in each given
#'lineup (PERP or SUSP) to generate the foils. Note that this results in
#'different foils being used in target-present and target-absent lineups. In
#'either case, the final result is a lineup constructed of vector
#'representations of a suspect (guilty or innocent) and five foils, all of
#'which share common features.
#'
#'Following the construction of the lineup, the model must simulate the actual
#'lineup procedure. WITNESS accomplishes this by comparing each lineup vector
#'(i.e., PERP, SUSP, and FOILS) to MEM to create match values, or assessments
#'of the degree to which two vectors overlap, for each lineup member. These
#'match values are the dot products of each lineup vector to MEM divided by the
#'total number of features. Larger dot products indicate a closer match to
#'memory. After computing these match values, the values are used to execute
#'the decision aspect of the lineup.
#'
#'In order to model relative/absolute contributions, WITNESS employs two
#'parameters: wa, the decision weight for the absolute contribution, and wr,
#'the decision weight for the relative contribution. These parameters are
#'proportionally complimentary in that they are constrained to sum to 1 (i.e.,
#'wa + wr = 1). WITNESS uses these weights to determine the contributions of
#'the two lineup members with the largest match values. wa governs the
#'contribution of the best match to MEM (BEST), whereas the contribution of the
#'second best match (NEXT) is governed by wr . When making its decision,
#'WITNESS will choose BEST if the evidence [EV = wa * BEST + wr * (BEST -
#'NEXT)] exceeds c, the decision criterion for simultaneous lineups. Thus, if
#'wa = 1 (and so wr = 0), the decision would be made entirely based on BEST's
#'match to memory (absolute contribution only); if wr = 1, the decision would
#'be made based only on the magnitude of the difference between BEST and NEXT
#'(relative contribution only). Both relative and absolute judgments can
#'contribute to identification decisions (meaning that wa could take any value
#'within the interval of 0 to 1). If the resulting EV value does not exceed c,
#'the lineup is rejected, meaning that no individual is selected from the
#'lineup.
#'
#'The user must specify the parameters of the model in matrix form. The columns
#'of the matrix are as follows: e (encoding parameter), ssp (similarity of
#'suspect to perp), sfs (similarity of foil to suspect), c (decision criteria),
#'and wa (relative vs. absolute criteria). The rows correspond to each data row
#'to be fit. See examples.
#'
#'@param parameters A parameters matrix. See details.
#'@param guilty A vector of logical values (i.e, TRUE or FALSE). This vector
#'indicates which rows in the dataset contain the guilty suspect.
#'@param perp.removed A vector of logical values indicating which rows have the
#'perp removed from the lineup.
#'@param sequential A vector of logical values indicating which rows should be fitted
#'with a sequential lineup. 
#'@param data Optional. Provide the dataset to be fit and this function will
#'compute the objective function
#'@param ilen The length of the features vector. Default is 100
#'@param lsize Lineup size. Default is 6.
#'@param meth Either "WITC" or "original". "WITC" will fix the wa/wr
#'parameters, while "original" will allow them to vary.
#'@param suspectMatch Logical. Should the foils be matched to the suspect?
#'@param seedit Logical. Should the random seed be fixed?  This should be set
#'at true when using a fitting algorithm so that differences between iterations
#'are due to difference in fit alone, not because of random noise.
#'@param sample.size Defines how many iterations the WITNESS function performs.
#'@param optim Logical. When TRUE, the function returns a single value (the fit
#'of the model to the data). This must be specified when using a optimization
#'routine (e.g., genetic algorithm).
#'@param goodFoil Logical. Should one of the foils be a "dead ringer"?
#'@examples
#'
#'# create data matrix
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
#'
#'est = witness.est(parameters=params, data=dataMatrix, meth="WITC", sample.size=100)
#'est
#'@export
witness.est <- function(parameters, guilty=NULL, perp.removed=NULL, sequential = NULL, data = NULL, ilen = 100, lsize=6, meth="WITC", suspectMatch=FALSE, seedit=TRUE, sample.size=1000, optim = FALSE, goodFoil=FALSE){
	
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
		# if (exists(".Random.seed", .GlobalEnv)){
			# seeder = TRUE
			# save.seed <- get(".Random.seed", .GlobalEnv)
		# } else {
			# seeder = FALSE
		# }
		
		
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
		
		###### make sure user understands default of sequential vector
		if (is.null(sequential)){
			sequential = rep(c(F), length.out=nrow(parameters))
			if (!optim){
				cat("Note: You failed to specify a sequential vector. We will assume that every 
				you want to use a simultaneous lineup.")
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
			
			####### create random feature vectors
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
			wr= 1-as.numeric(wa)
			evaluation = wr*(best.match-second.best) + (as.numeric(wa) * best.match)
			
			####### check if best match exceeds criteria
			makeid = evaluation > c
			
			####### randomize the order of positions for sequential
			randOrd = matrix(rep(1:6, times=length(a)), nrow=length(a), byrow=T)
			randOrd = t(apply(randOrd, 1, sample))
			
			######  see which one is above threshold
			firstAbove = rep(0, times=length(a))

			for (m in 1:length(a)){
				ff = (which(matchVect[m,randOrd[m,]]>c[m]))[1]
				if (!is.na(ff)){
					### if position of suspect is first
					if (which(randOrd[m,]==1)==ff){
						firstAbove[m] = 1
					} else {
						firstAbove[m] = 2
					}
				}
			}
			
			condIt = data.frame(seq=sequential, firstAbove=firstAbove, guilty=guilty, makeid=makeid,
					best.match.position=best.match.position)

	
			responseFunc = function(condition, lsize=6){
				resp = rep(0, times=3)
				if (condition["seq"]){
					if (condition['firstAbove']==1){
						resp = c(1,0,0)
					} else if (condition["firstAbove"]==0){
						resp = c(0,0,1)
					} else {
						resp = c(0,1,0)
					}
				} else {
					if (!condition["makeid"]){
						resp = c(0,0,1)
					} else if (condition["best.match.position"]==1){
						resp = c(1,0,0)
					} else {
						resp = c(0,1,0)
					}
				}
				
				return(resp)
			}	
			responses = responses + t(apply(condIt, 1, responseFunc))

		}	
		
		#### reset random seed
		# if (seeder){
			# assign(".Random.seed", save.seed, .GlobalEnv)
		# }
		if (seedit) {
			rm(.Random.seed, envir = .GlobalEnv)
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
