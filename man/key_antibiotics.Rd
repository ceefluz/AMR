% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/key_antibiotics.R
\name{key_antibiotics}
\alias{key_antibiotics}
\alias{key_antibiotics_equal}
\title{Key antibiotics for first \emph{weighted} isolates}
\usage{
key_antibiotics(tbl, col_mo = "mo", universal_1 = "amox",
  universal_2 = "amcl", universal_3 = "cfur", universal_4 = "pita",
  universal_5 = "cipr", universal_6 = "trsu", GramPos_1 = "vanc",
  GramPos_2 = "teic", GramPos_3 = "tetr", GramPos_4 = "eryt",
  GramPos_5 = "oxac", GramPos_6 = "rifa", GramNeg_1 = "gent",
  GramNeg_2 = "tobr", GramNeg_3 = "coli", GramNeg_4 = "cfot",
  GramNeg_5 = "cfta", GramNeg_6 = "mero", warnings = TRUE,
  col_bactid = "bactid")

key_antibiotics_equal(x, y, type = c("keyantibiotics", "points"),
  ignore_I = TRUE, points_threshold = 2, info = FALSE)
}
\arguments{
\item{tbl}{table with antibiotics coloms, like \code{amox} and \code{amcl}.}

\item{col_mo}{column name of the unique IDs of the microorganisms, see \code{\link{mo}}. If this column has another class than \code{"mo"}, values will be coerced using \code{\link{as.mo}}.}

\item{universal_1, universal_2, universal_3, universal_4, universal_5, universal_6}{column names of \strong{broad-spectrum} antibiotics, case-insensitive}

\item{GramPos_1, GramPos_2, GramPos_3, GramPos_4, GramPos_5, GramPos_6}{column names of antibiotics for \strong{Gram positives}, case-insensitive}

\item{GramNeg_1, GramNeg_2, GramNeg_3, GramNeg_4, GramNeg_5, GramNeg_6}{column names of antibiotics for \strong{Gram negatives}, case-insensitive}

\item{warnings}{give warning about missing antibiotic columns, they will anyway be ignored}

\item{col_bactid}{(deprecated, use \code{col_mo} instead)}

\item{x, y}{characters to compare}

\item{type}{type to determine weighed isolates; can be \code{"keyantibiotics"} or \code{"points"}, see Details}

\item{ignore_I}{logical to determine whether antibiotic interpretations with \code{"I"} will be ignored when \code{type = "keyantibiotics"}, see Details}

\item{points_threshold}{points until the comparison of key antibiotics will lead to inclusion of an isolate when \code{type = "points"}, see Details}

\item{info}{print progress}
}
\description{
These function can be used to determine first isolates (see \code{\link{first_isolate}}). Using key antibiotics to determine first isolates is more reliable than without key antibiotics. These selected isolates will then be called first \emph{weighted} isolates.
}
\details{
The function \code{key_antibiotics} returns a character vector with 12 antibiotic results for every isolate. These isolates can then be compared using \code{key_antibiotics_equal}, to check if two isolates have generally the same antibiogram. Missing and invalid values are replaced with a dot (\code{"."}). The \code{\link{first_isolate}} function only uses this function on the same microbial species from the same patient. Using this, an MRSA will be included after a susceptible \emph{S. aureus} (MSSA) found within the same episode (see \code{episode} parameter of \code{\link{first_isolate}}). Without key antibiotic comparison it wouldn't.

  At default, the antibiotics that are used for \strong{Gram positive bacteria} are (colum names): \cr
  \code{"amox"}, \code{"amcl"}, \code{"cfur"}, \code{"pita"}, \code{"cipr"}, \code{"trsu"} (until here is universal), \code{"vanc"}, \code{"teic"}, \code{"tetr"}, \code{"eryt"}, \code{"oxac"}, \code{"rifa"}.

  At default, the antibiotics that are used for \strong{Gram negative bacteria} are (colum names): \cr
  \code{"amox"}, \code{"amcl"}, \code{"cfur"}, \code{"pita"}, \code{"cipr"}, \code{"trsu"} (until here is universal), \code{"gent"}, \code{"tobr"}, \code{"coli"}, \code{"cfot"}, \code{"cfta"}, \code{"mero"}.


  The function \code{key_antibiotics_equal} checks the characters returned by \code{key_antibiotics} for equality, and returns a logical vector.
}
\section{Key antibiotics}{

    There are two ways to determine whether isolates can be included as first \emph{weighted} isolates which will give generally the same results: \cr

    \strong{1. Using} \code{type = "keyantibiotics"} \strong{and parameter} \code{ignore_I} \cr
    Any difference from S to R (or vice versa) will (re)select an isolate as a first weighted isolate. With \code{ignore_I = FALSE}, also differences from I to S|R (or vice versa) will lead to this. This is a reliable method and 30-35 times faster than method 2. \cr

    \strong{2. Using} \code{type = "points"} \strong{and parameter} \code{points_threshold} \cr
    A difference from I to S|R (or vice versa) means 0.5 points, a difference from S to R (or vice versa) means 1 point. When the sum of points exceeds \code{points_threshold}, an isolate will be (re)selected as a first weighted isolate.
}

\examples{
# septic_patients is a dataset available in the AMR package
?septic_patients
my_patients <- septic_patients

library(dplyr)
# set key antibiotics to a new variable
my_patients <- my_patients \%>\%
  mutate(keyab = key_antibiotics(.)) \%>\%
  mutate(
    # now calculate first isolates
    first_regular = first_isolate(., "date", "patient_id", "mo"),
    # and first WEIGHTED isolates
    first_weighted = first_isolate(., "date", "patient_id", "mo",
                                   col_keyantibiotics = "keyab")
  )

# Check the difference, in this data set it results in 7\% more isolates:
sum(my_patients$first_regular, na.rm = TRUE)
sum(my_patients$first_weighted, na.rm = TRUE)


# output of the `key_antibiotics` function could be like this:
strainA <- "SSSRR.S.R..S"
strainB <- "SSSIRSSSRSSS"

key_antibiotics_equal(strainA, strainB)
# TRUE, because I is ignored (as are missing values)

key_antibiotics_equal(strainA, strainB, ignore_I = FALSE)
# FALSE, because I is not ignored and so the 4th value differs
}
\seealso{
\code{\link{first_isolate}}
}
