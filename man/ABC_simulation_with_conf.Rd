% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABC_simulation_with_conf.R
\name{ABC_simulation_with_conf}
\alias{ABC_simulation_with_conf}
\title{Do ABC simulation with a configuration table}
\usage{
ABC_simulation_with_conf(
  npod,
  conf,
  time.range,
  buffer = 0,
  concentrationscale,
  concentrationShape,
  prefix = "temp",
  BayeSSCallocation,
  do.parallel = 1,
  write.reference.file = F
)
}
\arguments{
\item{npod}{An integer, the number of simulation replicates}

\item{conf}{A data frame with configurations}

\item{time.range}{A vector containing two numeric elements specifying the minimal and maximum event time.}

\item{buffer}{The minimal amount of time separating two events. The default is 0, when events' time are sampled randomly}

\item{concentrationscale}{A numeric. The \code{scale} for the gamma distribution. See 'Details'}

\item{concentrationShape}{A numeric. The \code{shape} for the gamma distribution. See 'Details'}

\item{prefix}{BayeSSC will generate many intermediate files, this function will create temporary folders for storing these files. One folder for one row in the configuration table, and the folders are named as \code{prefix}\code{simulation replicates}_row number.If folders with the same name already exist, user will be prompt for providing another \code{prefix} or overwriting the folders.}

\item{BayeSSCallocation}{A character string providing the location of the BayeSSC executable (including the file name)}

\item{do.parallel}{A integer for the number of parallel threads to use. if equals \code{1},no parallel execution. The parallel R package required}

\item{write.reference.file}{\code{TRUE} or \code{FALSE}(Default). If \code{TRUE}, the simulated hyperstat will be append to a file \code{prefix}_reference_table. Useful for running large number of replicates.If \code{False}, the simulated hyperstat will be collected and returned as a data frame}
}
\value{
If \code{write.reference.file} is \code{FALSE}, this function returns a data frame containing the simulated hyperstats. The first column is 'uid' (\code{prefix}+the replicates number), the second column is the number of events, and the third column is the total number of species. Following are columns of the "real" expansion time for each species in the simulation, and then, the columns of hyperstats.If \code{write.reference.file} is \code{TRUE}, this function returns a character string providing the path to the reference file.
}
\description{
A function for doing abc simulation with a configuration table. see Readme on github for required columns in the table
}
\details{
For the prior on the number of co-expansion events, hBayeSSC used a flat distribution. In the paper, we adopt the [PyMsbayes-style](http://joaks1.github.io/PyMsBayes/) prior: a gamma distribution for the alpha parameter of the dirichelete process. User need to select the two parameters for gamma distribution: concentrationShape and consentrationScale. Check [prior selection in PyMsbayes](http://joaks1.github.io/PyMsBayes/tutorials/selecting-priors.html) for how to select these two parameters.
}
