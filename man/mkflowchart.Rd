% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohorttools.R
\name{mkflowchart}
\alias{mkflowchart}
\title{Function makes flowchart in graphviz}
\usage{
mkflowchart(N, text.M, text.P, type = 1)
}
\arguments{
\item{N}{Population sizes}

\item{text.M}{Text for exclusions, length one less than N}

\item{text.P}{Text for main boxes, must be same length with N}

\item{type}{flowchart type (1 or 2)}
}
\value{
Character string, graphviz language
}
\description{
Function makes flowchart in graphviz
}
\examples{
DiagrammeR::grViz(mkflowchart(N=c(743,32,20),
text.M=c("Excluded","Excluded \n other with reasons"),
text.P=c("Studies","Relevant studies","Included in final review"),type=1))
}
\author{
Jari Haukka \email{jari.haukka@helsinki.fi}
}
