#'Generate random starting values
#'
#'This function takes a parameter specification matrix and returns a matrix of
#'random starting values
#'
#'
#'@param parameter.form A matrix where the rows are the experiments, and the
#'columns are the parameters a, ssp, sfs, c, and wa. To fix parameters across
#'rows, one would simply input the same values. For example, if one wanted all
#'the encoding values (a) to be the same, the first row could consist of only
#'ones. See example.
#'@param values A vector of starting values. If the user wants random starting
#'values, leave as NULL. Otherwise, specify what the starting values should be.
#'@examples
#'
#'###### prepare parameter vector
#'params = matrix(c(rep("e", times=4),
#'				0, 0, 0, 0,
#'				"sfs1", "sfs1", "sfs2", "sfs2",
#'				"cr1", "cr2", "cr1", "cr2",
#'				1, 1, 1, 1), nrow=4)
#'		# note: fixes ssp and wa encoding across experiments, but allows sfs/cr to vary
#'		###### generate initial parameters
#'st = witness.starting.params(params)
#'st
#'@export
witness.starting.params <-function(parameter.form, values=NULL){
	parameter.start.values = matrix(nrow=nrow(parameter.form), ncol=5)
	s=1
	for (i in 1:5){
		num.starting.values = length(unique(parameter.form[,i]))
		if (i ==4){ top.num=.2 } else {top.num=1}
		if (!is.null(values)){
				random.value = values[s:(s+num.starting.values-1)]
				s = s + num.starting.values
			} else {
			random.value = runif(num.starting.values, 0, top.num)
			}
		un.ones = unique(parameter.form[,i])
		new = 1:nrow(parameter.form)
		for (j in 1:length(un.ones)){
			num = which(parameter.form[,i]==un.ones[j])
			parameter.start.values[num,i] = random.value[j]
		}					
	}
	return(parameter.start.values)
}
