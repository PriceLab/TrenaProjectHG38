% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrenaProjectHG38.R
\docType{class}
\name{TrenaProjectHG38-class}
\alias{TrenaProjectHG38-class}
\alias{.TrenaProjectHG38}
\alias{TrenaProjectHG38}
\title{TrenaProjectHG38}
\usage{
TrenaProjectHG38(
  projectName,
  supportedGenes,
  footprintDatabaseHost,
  footprintDatabaseNames,
  footprintDatabasePort = 5432,
  packageDataDirectory,
  quiet
)
}
\arguments{
\item{supportedGenes}{a vector of character strings}

\item{footprintDatabaseHost}{Character string (e.g., "khaleesi.systemsbiology.net")}

\item{footprintDatabaseNames}{Character string (e.g., "hint_brain_20")}

\item{packageDataDirectory}{A string pointing to the parent of a more-or-less standard set of data subdirectories}

\item{quiet}{A logical indicating whether or not the Trena object should print output}
}
\value{
An object of the TrenaProjectHG38 class
}
\description{
TrenaProjectHG38 and its (projected) subclasses provide convenient containers in which to collect
 trena-related aggregation of a gene's (a hierarchy of classes) including expression data,
transcript and variant info, genomic and epigenomic context, trena models and/or the means to create them
}
