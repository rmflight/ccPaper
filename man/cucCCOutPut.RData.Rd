\docType{data}
\name{cucCCOutPut.RData}
\alias{cucCCOutPut.RData}
\title{cucCCOutPut.RData}
\source{
  interacting with GO nodes in Cytoscape to partition into
  groups
}
\description{
  cucCCOutPut.RData
}
\details{
  After performing \code{ccCompare} on \code{cucRomer}, any
  edges =< 0.8 similarity were removed, and the graph
  loaded into Cytoscape using \code{ccOutCyt}. Individual
  groups of GO nodes were determined by Robert M Flight
  visualy and exploratively, and labels assigned based on
  biological knowledge of Robert M Flight. Each GO group is
  defined in \code{cucNodeGroups}, with the list of GO
  terms that are part of that group, along with their
  attributes.
}
\keyword{datasets}

