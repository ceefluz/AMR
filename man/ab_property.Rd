% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ab_property.R
\name{ab_property}
\alias{ab_property}
\alias{ab_atc}
\alias{ab_official}
\alias{ab_name}
\alias{ab_trivial_nl}
\alias{ab_certe}
\alias{ab_umcg}
\alias{ab_tradenames}
\title{Property of an antibiotic}
\usage{
ab_property(x, property = "official")

ab_atc(x)

ab_official(x, language = NULL)

ab_name(x, language = NULL)

ab_trivial_nl(x)

ab_certe(x)

ab_umcg(x)

ab_tradenames(x)
}
\arguments{
\item{x}{a (vector of a) valid \code{\link{atc}} code or any text that can be coerced to a valid atc with \code{\link{as.atc}}}

\item{property}{one of the column names of one of the \code{\link{antibiotics}} data set, like \code{"atc"} and \code{"official"}}

\item{language}{language of the returned text, defaults to English (\code{"en"}) and can be set with \code{\link{getOption}("AMR_locale")}. Either one of \code{"en"} (English) or \code{"nl"} (Dutch).}
}
\value{
A vector of values. In case of \code{ab_tradenames}, if \code{x} is of length one, a vector will be returned. Otherwise a \code{\link{list}}, with \code{x} as names.
}
\description{
Use these functions to return a specific property of an antibiotic from the \code{\link{antibiotics}} data set, based on their ATC code. Get such a code with \code{\link{as.atc}}.
}
\examples{
ab_atc("amcl")         # J01CR02
ab_name("amcl")        # Amoxicillin and beta-lactamase inhibitor
ab_name("amcl", "nl")  # Amoxicilline met enzymremmer
ab_trivial_nl("amcl")  # Amoxicilline/clavulaanzuur
ab_certe("amcl")       # amcl
ab_umcg("amcl")        # AMCL
}
\seealso{
\code{\link{antibiotics}}
}
