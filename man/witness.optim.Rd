\name{witness.optim}
\alias{witness.optim}
\title{General purpose fitting for WITNESS}
\usage{
  witness.optim(param.form, data.set, generate.new = TRUE,
    ...)
}
\arguments{
  \item{param.form}{A parameters matrix. See details.}

  \item{data.set}{The dataset to be fit. Columns correspond
  to Target Chosen, Foil Chosen, and Lineup Rejected. Rows
  are for each experiment.}

  \item{generate.new}{If starting values are not provided,
  making this TRUE will generate new starting values.}

  \item{...}{Other parameters passed to the WITNESS model}
}
\value{
  This function returns two objects:

  \item{fitted}{A matrix in the same form as
  \code{param.form} that contains the fitted values.}

  \item{fit}{The rmse of the fitted object.}
}
\description{
  \code{optim.witness} uses a steepest-descent algorithm to
  fit empirical data using the WITNESS model
}
\details{
  The user must specify the parameters of the model in
  matrix form. The columns of the matrix are as follows: e
  (encoding parameter), ssp (similarity of suspect to
  perp), sfs (similarity of foil to suspect), c (decision
  criteria), and wa (relative vs. absolute criteria). The
  rows correspond to each data row to be fit. See examples.
}
\examples{
# create an matrix with data to be fit
dataMatrix = matrix(c(.471, .230, .350,
							.208, .137, .513,
							.396, .431, .242,
							.166, .081, .669), nrow=4, byrow=TRUE)
# specify parameter form
params = matrix(c(rep("e", times=4),
		0, 0, 0, 0,
		"sfs1", "sfs1", "sfs2", "sfs2",
		"cr1", "cr2", "cr1", "cr2",
		1, 1, 1, 1), nrow=4)
# find optimal parameters (currently commented to save time)
# fit = witness.optim(params, data.set=dataMatrix, sample.size=100, meth="WITC")
# fit
}
\seealso{
  \code{\link{witness.genetic}}, \code{\link{witness.est}}
}

