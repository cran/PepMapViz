% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_modification.R
\name{obtain_mod_Comet}
\alias{obtain_mod_Comet}
\title{Obtain modification information from Peptide data generated by Comet}
\usage{
obtain_mod_Comet(
  data,
  mod_column,
  PTM_table = NULL,
  PTM_annotation = FALSE,
  PTM_mass_column
)
}
\arguments{
\item{data}{A data.table with a column containing PTM information.}

\item{mod_column}{The name of the column containing the modified peptide sequences.}

\item{PTM_table}{A data.table with columns 'PTM_mass' and 'PTM_type' containing PTM annotation information.}

\item{PTM_annotation}{A logical value indicating whether to include PTM annotation information in the result.}

\item{PTM_mass_column}{The name of the column containing the PTM mass information}
}
\value{
A data.table with 'PTM_position', 'PTM_type', 'reps', and other columns.
}
\description{
This function takes Comet output containing a column with modified peptide
sequences including PTM information and converts it into a new dataframe with the
desired format of peptide sequences and associated PTM information.
}
\examples{
library(data.table)
data <- data.table(
  modified_peptide = c(
    "AAM[15.9949]Q[-0.98]RGSLYQCDYSTGSC[57.02]EPIR",
    "K.AAQQTGKLVHANFGT.K",
    "K.[-0.98]AATVTGKLVHANFGT.K"
  ),
  plain_peptide = c(
    "AAMQRGSLYQCDYSTGSCEPIR",
    "AAQQTGKLVHANFGT",
    "AATVTGKLVHANFGT"
  ),
  Condition = c("A", "B", "B")
)
PTM_table <- data.table(
  PTM_mass = c(57.02, -0.98, 15.9949),
  PTM_type = c("Cam", "Amid", "Ox")
)
mod_column <- 'modified_peptide'
PTM_mass_column <- "PTM_mass"
converted_data <- obtain_mod_Comet(data, mod_column, PTM_table,
PTM_annotation = TRUE, PTM_mass_column)

}
