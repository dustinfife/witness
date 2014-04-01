\name{witness.est}
\alias{witness.est}
\title{Estimate identification probabilities using WITNESS}
\usage{
  witness.est(parameters, guilty = NULL,
    perp.removed = NULL, sequential = NULL, data = NULL,
    ilen = 100, lsize = 6, meth = "WITC",
    suspectMatch = FALSE, seedit = TRUE,
    sample.size = 1000, optim = FALSE, goodFoil = FALSE)
}
\arguments{
  \item{parameters}{A parameters matrix. See details.}

  \item{guilty}{A vector of logical values (i.e, TRUE or
  FALSE). This vector indicates which rows in the dataset
  contain the guilty suspect.}

  \item{perp.removed}{A vector of logical values indicating
  which rows have the perp removed from the lineup.}

  \item{sequential}{A vector of logical values indicating
  which rows should be fitted with a sequential lineup.}

  \item{data}{Optional. Provide the dataset to be fit and
  this function will compute the objective function}

  \item{ilen}{The length of the features vector. Default is
  100}

  \item{lsize}{Lineup size. Default is 6.}

  \item{meth}{Either "WITC" or "original". "WITC" will fix
  the wa/wr parameters, while "original" will allow them to
  vary.}

  \item{suspectMatch}{Logical. Should the foils be matched
  to the suspect?}

  \item{seedit}{Logical. Should the random seed be fixed?
  This should be set at true when using a fitting algorithm
  so that differences between iterations are due to
  difference in fit alone, not because of random noise.}

  \item{sample.size}{Defines how many iterations the
  WITNESS function performs.}

  \item{optim}{Logical. When TRUE, the function returns a
  single value (the fit of the model to the data). This
  must be specified when using a optimization routine
  (e.g., genetic algorithm).}

  \item{goodFoil}{Logical. Should one of the foils be a
  "dead ringer"?}
}
\description{
  \code{witness} takes the parameters provided by the user
  then simulates the WITNESS procedure and outputs
  probabilities
}
\details{
  The WITNESS model (Clark, 2003) is a direct-access
  matching model that has been adapted for eyewitness
  situations. WITNESS uses numerical representations of
  features as items in the matching process. More
  specifically, each lineup member is represented as a
  vector of n arbitrary features (the features do not map
  to literal features of a person, such as eyes or a nose,
  however).  The value n is typically fixed to some
  arbitrary constant such as 100.  Each position in the
  vector (or each feature) is assigned a value between -1
  and 1.

  To begin, WITNESS generates a perpetrator vector (PERP)
  by randomly generating numbers on the interval -1, 1.
  This vector serves as the basis for all of the subsequent
  lineup members. The model next "encodes" the features in
  the PERP vector to a memory vector. The parameter a
  governs the degree to which the Memory (MEM) vector
  matches the PERP vector; a is the probability that each
  individual feature in PERP will be successfully encoded
  to MEM. For example, if a = .4, then each feature of PERP
  will have a 40\% chance of being represented within MEM.
  Conversely, there is a 1-a (60\%, in this case) chance
  that the feature will be replaced with another arbitrary
  random feature. Put differently, the correlation
  coefficient between PERP and MEM will be approximately a.
  This results in a noisy representation in memory, which
  accounts for phenomena such as the failure to encode
  details of the perpetrator and post-event interference.

  WITNESS next creates the lineup members for comparison to
  memory. For target absent lineups (i.e., an innocent
  suspect replaces the guilty suspect), a new vector is
  created to represent the innocent suspect (SUSP). This
  vector is governed by the parameter SSP, or Similarity of
  the Suspect to the Perpetrator. SSP is the probability
  that each feature of SUSP will match the corresponding
  feature of PERP. As SSP approaches 1, SUSP will be more
  similar to PERP (and likewise, less similar as it
  approaches 0).

  Next, WITNESS uses one of two more parameters to create
  the remaining lineup members (i.e., the foils), each to
  simulate a different method of foil selection. We refer
  to these vectors of foils as FOILS. Foils can be selected
  either because they match the description of the
  perpetrator (description-matched foils) or because they
  match the appearance of the suspect (suspect-matched
  foils). The parameter SFP, or Similarity of the Foils to
  the Perpetrator, simulates a description-matched lineup.
  SFP is the probability that a given feature of a given
  foil will match a corresponding feature in PERP (similar
  to SSP, where 0 yields no shared features and 1
  replicates PERP). Note that this results in the same
  foils being used in target-present and target-absent
  lineups. SFS, Similarity of the Foils to the Suspect, is
  a similar parameter; however, it uses the suspect in each
  given lineup (PERP or SUSP) to generate the foils. Note
  that this results in different foils being used in
  target-present and target-absent lineups. In either case,
  the final result is a lineup constructed of vector
  representations of a suspect (guilty or innocent) and
  five foils, all of which share common features.

  Following the construction of the lineup, the model must
  simulate the actual lineup procedure. WITNESS
  accomplishes this by comparing each lineup vector (i.e.,
  PERP, SUSP, and FOILS) to MEM to create match values, or
  assessments of the degree to which two vectors overlap,
  for each lineup member. These match values are the dot
  products of each lineup vector to MEM divided by the
  total number of features. Larger dot products indicate a
  closer match to memory. After computing these match
  values, the values are used to execute the decision
  aspect of the lineup.

  In order to model relative/absolute contributions,
  WITNESS employs two parameters: wa, the decision weight
  for the absolute contribution, and wr, the decision
  weight for the relative contribution. These parameters
  are proportionally complimentary in that they are
  constrained to sum to 1 (i.e., wa + wr = 1). WITNESS uses
  these weights to determine the contributions of the two
  lineup members with the largest match values. wa governs
  the contribution of the best match to MEM (BEST), whereas
  the contribution of the second best match (NEXT) is
  governed by wr . When making its decision, WITNESS will
  choose BEST if the evidence [EV = wa * BEST + wr * (BEST
  - NEXT)] exceeds c, the decision criterion for
  simultaneous lineups. Thus, if wa = 1 (and so wr = 0),
  the decision would be made entirely based on BEST's match
  to memory (absolute contribution only); if wr = 1, the
  decision would be made based only on the magnitude of the
  difference between BEST and NEXT (relative contribution
  only). Both relative and absolute judgments can
  contribute to identification decisions (meaning that wa
  could take any value within the interval of 0 to 1). If
  the resulting EV value does not exceed c, the lineup is
  rejected, meaning that no individual is selected from the
  lineup.

  The user must specify the parameters of the model in
  matrix form. The columns of the matrix are as follows: e
  (encoding parameter), ssp (similarity of suspect to
  perp), sfs (similarity of foil to suspect), c (decision
  criteria), and wa (relative vs. absolute criteria). The
  rows correspond to each data row to be fit. See examples.
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

