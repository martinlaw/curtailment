% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find2stageDesigns.R
\name{find2stageDesigns}
\alias{find2stageDesigns}
\title{Find two-stage designs}
\usage{
find2stageDesigns(
  nmin,
  nmax,
  p0,
  p1,
  alpha,
  power,
  maxthetaF = NA,
  benefit = FALSE
)
}
\arguments{
\item{nmin}{Minimum permitted sample size. Should be a multiple of block size or number of stages.}

\item{nmax}{Maximum permitted sample size. Should be a multiple of block size or number of stages.}

\item{p0}{Probability for which to control the type-I error-rate}

\item{p1}{Probability for which to control the power}

\item{alpha}{Significance level}

\item{power}{Required power (1-beta)}

\item{maxthetaF}{Maximum permitted conditional power for futility stopping. Optional.}

\item{benefit}{Allow the trial to end for a go decision and reject the null hypothesis at the interim analysis (i.e., the design of Mander and Thompson)}
}
\value{
A list of class "curtailment_simon" containing two data frames. The first data frame, $input,
has a single row and contains all the inputted values. The second data frame, $all.des, contains one
row for each design realisation, and contains the details of each design, including sample size,
stopping boundaries and operating characteristics. To see a diagram of any obtained design realisation
and its corresponding stopping boundaries, simply call the function drawDiagram with this output as the only argument.
}
\description{
This function finds two-stage designs for a given set of design parameters, allowing
stopping for benefit at the interim (Mander and Thompson's design) or no stopping
for benefit at the interim (Simon's design). It returns not
only the optimal and minimax design realisations, but all design realisations that could
be considered "best" in terms of expected sample size under p=p0 (EssH0), expected
sample size under p=p1 (Ess), maximum sample size (n) or any weighted combination of these
three optimality criteria.
}
\examples{
\donttest{
find2stageDesigns(nmin=23,
 nmax=27,
 p0=0.75,
 p1=0.92,
 alpha=0.22,
 power=0.95,
 benefit=TRUE)
 }
}
\references{
\doi{10.1016/j.cct.2010.07.008}A.P. Mander, S.G. Thompson,
Two-stage designs optimal under the alternative hypothesis for phase II cancer clinical trials,
Contemporary Clinical Trials,
Volume 31, Issue 6,
2010,
Pages 572-578

\doi{10.1016/0197-2456(89)90015-9}Richard Simon,
Optimal two-stage designs for phase II clinical trials,
Controlled Clinical Trials,
Volume 10, Issue 1,
1989,
Pages 1-10
}
\author{
Martin Law, \email{martin.law@mrc-bsu.cam.ac.uk}
}
