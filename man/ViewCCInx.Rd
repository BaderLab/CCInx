% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShinyInx.R
\name{ViewCCInx}
\alias{ViewCCInx}
\title{Run the Shiny app to interactively explore cell-cell interactions.}
\usage{
ViewCCInx(INX, ...)
}
\arguments{
\item{INX}{The list containing edgelist and node metadata as generated by
\code{\link{BuildCCInx}}.}

\item{...}{Named options that should be passed to the
\code{\link[shiny]{runApp}} call (these can be any of the following:
"port", "launch.browser", "host", "quiet", "display.mode" and "test.mode").}
}
\description{
Run the Shiny app to interactively explore cell-cell interactions.
}
