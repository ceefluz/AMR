% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{septic_patients}
\alias{septic_patients}
\title{Data set with 2000 blood culture isolates of septic patients}
\format{A \code{\link{tibble}} with 2,000 observations and 49 variables:
\describe{
  \item{\code{date}}{date of receipt at the laboratory}
  \item{\code{hospital_id}}{ID of the hospital, from A to D}
  \item{\code{ward_icu}}{logical to determine if ward is an intensive care unit}
  \item{\code{ward_clinical}}{logical to determine if ward is a regular clinical ward}
  \item{\code{ward_outpatient}}{logical to determine if ward is an outpatient clinic}
  \item{\code{age}}{age of the patient}
  \item{\code{sex}}{sex of the patient}
  \item{\code{patient_id}}{ID of the patient, first 10 characters of an SHA hash containing irretrievable information}
  \item{\code{mo}}{ID of microorganism, see \code{\link{microorganisms}}}
  \item{\code{peni:rifa}}{40 different antibiotics with class \code{rsi} (see \code{\link{as.rsi}}); these column names occur in \code{\link{antibiotics}} data set and can be translated with \code{\link{abname}}}
}}
\usage{
septic_patients
}
\description{
An anonymised data set containing 2,000 microbial blood culture isolates with their full antibiograms found in septic patients in 4 different hospitals in the Netherlands, between 2001 and 2017. It is true, genuine data. This \code{data.frame} can be used to practice AMR analysis. For examples, press F1.
}
\examples{
# ----------- #
# PREPARATION #
# ----------- #

# Save this example data set to an object, so we can edit it:
my_data <- septic_patients

# load the dplyr package to make data science A LOT easier
library(dplyr)

# Add first isolates to our data set:
my_data <- my_data \%>\%
  mutate(first_isolates = first_isolate(my_data, "date", "patient_id", "mo"))

# -------- #
# ANALYSIS #
# -------- #

# 1. Get the amoxicillin resistance percentages (p)
#     and numbers (n) of E. coli, divided by hospital:

my_data \%>\%
  filter(mo == guess_mo("E. coli"),
         first_isolates == TRUE) \%>\%
  group_by(hospital_id) \%>\%
  summarise(n = n_rsi(amox),
            p = portion_IR(amox))


# 2. Get the amoxicillin/clavulanic acid resistance
#    percentages of E. coli, trend over the years:

my_data \%>\%
  filter(mo == guess_mo("E. coli"),
         first_isolates == TRUE) \%>\%
  group_by(year = format(date, "\%Y")) \%>\%
  summarise(n = n_rsi(amcl),
            p = portion_IR(amcl, minimum = 20))
}
\keyword{datasets}
