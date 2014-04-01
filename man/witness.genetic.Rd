\name{witness.genetic}
\alias{witness.genetic}
\title{Fit the WITNESS model using a Genetic Algorithm}
\usage{
  witness.genetic(parameter.form, data = NULL, N = 1000,
    generations = 100, pcross = 0.8, mutProb = 0.001,
    prop.random = 0.001, gradient = 0.01, ...)
}
\arguments{
  \item{parameter.form}{a matrix where the rows are the
  experiments, and the columns are the parameters a, ssp,
  sfs, c, and wa. To fix parameters across rows, one would
  simply input the same values. For example, if one wanted
  all the encoding values (a) to be the same, the first row
  could consist of only ones.}

  \item{data}{a matrix that has the same number of rows as
  parameter.form. Again, the rows are the experiments but
  the columns correspond to 1. Target Identification, 2.
  Foil Identification, and 3. No Identification.}

  \item{N}{the number of iterations the WITNESS model will
  do to simulate the eyewitness procedure.}

  \item{generations}{how many generations the genetic
  algorithm will iterate through before quitting.}

  \item{pcross}{probability of cross "breeding," or the
  probability that a parent will breed.}

  \item{mutProb}{probability of randomly mutating a
  parameter (i.e., the probability that noise will be added
  to the next generation)}

  \item{prop.random}{to avoid local minima, the algorithm
  allows for "aliens" to be introduced (or completely new
  parameter values).}

  \item{gradient}{convergence criteria. Currently not
  implemented.}

  \item{...}{other parameters passed to the WITNESS model.}
}
\description{
  After generating arbitrary starting values, the genetic
  algorithm selects the best parameter combinations to
  "mate." At each iteration, mutations are added with a
  specified probability, as well as totally new parameter
  values.
}
\examples{
dataMatrix = matrix(c(.471, .230, .350,
							.208, .137, .513,
							.396, .431, .242,
							.166, .081, .669), nrow=4, byrow=TRUE)
# create an parameter.form
params = matrix(c(rep("e", times=4),
		0, 0, 0, 0,
		"sfs1", "sfs1", "sfs2", "sfs2",
		"cr1", "cr2", "cr1", "cr2",
		1, 1, 1, 1), nrow=4)

		####### do genetic algorithm (commented b/c it takes a while)
#fit = witness.genetic(params, data=dataMatrix, N=10, generations=1, sample.size=100, meth="WITC")
}
\references{
  Clark, S. E. (2003). A memory and decision model for
  eyewitness identification. Applied Cognitive Psychology,
  17, 629-654.
}
\seealso{
  \code{\link{witness.optim}}, \code{\link{witness.est}}
}

