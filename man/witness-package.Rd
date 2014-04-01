\docType{package}
\name{witness-package}
\alias{witness}
\alias{witness-package}
\title{This package fits eyewitness data using Clark's (2003) WITNESS model}
\description{
  WITNESS is a computational memory model for eyewitness
  identification. This package is able to estimate
  probabilities based on user-inputted parameters as well
  as find optimal parameters using a steepest-descent
  and/or a genetic algorithm.
}
\details{
  \tabular{ll}{ Package: \tab witness\cr Type: \tab
  Package\cr Version: \tab 1.0\cr Date: \tab 2012-10-25\cr
  License: \tab GNU 2\cr }
}
\examples{
# create data matrix
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


est = witness.est(parameters=params, data=dataMatrix, meth="WITC", sample.size=100)
est
}
\author{
  Maintainer: Dustin Fife <fife.dustin@gmail.com>
}
\references{
  Clark, S. E. (2003). A memory and decision model for
  eyewitness identification. Applied Cognitive Psychology,
  17, 629-654.
}
\seealso{
  NONE
}
\keyword{package}

