\name{witness.starting.params}
\alias{witness.starting.params}
\title{Generate random starting values}
\usage{
  witness.starting.params(parameter.form, values = NULL)
}
\arguments{
  \item{parameter.form}{A matrix where the rows are the
  experiments, and the columns are the parameters a, ssp,
  sfs, c, and wa. To fix parameters across rows, one would
  simply input the same values. For example, if one wanted
  all the encoding values (a) to be the same, the first row
  could consist of only ones. See example.}

  \item{values}{A vector of starting values. If the user
  wants random starting values, leave as NULL. Otherwise,
  specify what the starting values should be.}
}
\description{
  This function takes a parameter specification matrix and
  returns a matrix of random starting values
}
\examples{
###### prepare parameter vector
params = matrix(c(rep("e", times=4),
				0, 0, 0, 0,
				"sfs1", "sfs1", "sfs2", "sfs2",
				"cr1", "cr2", "cr1", "cr2",
				1, 1, 1, 1), nrow=4)
		# note: fixes ssp and wa encoding across experiments, but allows sfs/cr to vary
		###### generate initial parameters
st = witness.starting.params(params)
st
}

