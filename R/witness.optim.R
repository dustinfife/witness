#'General purpose fitting for WITNESS
#'
#'\code{optim.witness} uses a steepest-descent algorithm to fit empirical data
#'using the WITNESS model
#'
#'The user must specify the parameters of the model in matrix form. The columns
#'of the matrix are as follows: e (encoding parameter), ssp (similarity of
#'suspect to perp), sfs (similarity of foil to suspect), c (decision criteria),
#'and wa (relative vs. absolute criteria). The rows correspond to each data row
#'to be fit. See examples.
#'
#'@param param.form A parameters matrix. See details.
#'@param data.set The dataset to be fit. Columns correspond to Target Chosen,
#'Foil Chosen, and Lineup Rejected. Rows are for each experiment.
#'@param generate.new If starting values are not provided, making this TRUE
#'will generate new starting values.
#'@param ... Other parameters passed to the WITNESS model
#'@seealso \code{\link{witness.genetic}}, \code{\link{witness.est}}
#'@return This function returns two objects:
#'@return \item{fitted}{A matrix in the same form as \code{param.form} that contains the fitted values.}
#'@return \item{fit}{The rmse of the fitted object.}
#'@examples
#'# create an matrix with data to be fit
#'dataMatrix = matrix(c(.471, .230, .350,
#'							.208, .137, .513,
#'							.396, .431, .242,
#'							.166, .081, .669), nrow=4, byrow=TRUE)
#'# specify parameter form
#'params = matrix(c(rep("e", times=4),
#'		0, 0, 0, 0,
#'		"sfs1", "sfs1", "sfs2", "sfs2",
#'		"cr1", "cr2", "cr1", "cr2",
#'		1, 1, 1, 1), nrow=4)
#'# find optimal parameters (currently commented to save time)
#'# fit = witness.optim(params, data.set=dataMatrix, sample.size=100, meth="WITC")
#'# fit
#'@export
witness.optim <-function(param.form, data.set, generate.new=TRUE,...){
	if (generate.new){
		random.params = witness.starting.params(param.form)
	} else { random.params = param.form }
	optimal = optim(par=unique(as.numeric(random.params)), fn=fitOpWit, param.form=random.params, dataset=data.set, ...)
	give.back = witness.starting.params(param.form, values=optimal$par)
	list(fitted = give.back, fit=optimal$value)
}
