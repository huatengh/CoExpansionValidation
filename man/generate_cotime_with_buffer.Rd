% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_cotime_with_buffer.R
\name{generate_cotime_with_buffer}
\alias{generate_cotime_with_buffer}
\title{Generate co-expansion time}
\usage{
generate_cotime_with_buffer(time.range, nco.events, buffer = 0)
}
\arguments{
\item{time.range}{A numeric vector of length 2 giving the minimal and maximum time.}

\item{nco.events}{An integer. The number of co-expansion events. See 'Details'.}

\item{buffer}{A numeric. The minimal amount of time separating two events. The default is 0, when events' time are sampled randomly}
}
\value{
A numeric vector of length \code{nco.events} containing the sampled expansion time
}
\description{
This function randomly draw events' time from a time range.
User can specify the minimal time length (buffer) between events
}
\details{
Note that the \code{buffer} is the time length on each side of a expansion event. If the \code{buffer} is non-zero,the maximum number of events (max \code{nco.event})=\code{time.range}/\code{buffer} +2. However, for a large \code{nco.events}, even if it does not exceed the maximum, this function might take a while to find a set of time that are at least \code{buffer} apart from each other.
}
