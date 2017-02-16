\name{hyperLOPITU2OS2017}
\alias{hyperLOPITU2OS2017}
\alias{hyperLOPITU2OS2017b}

\docType{data}

\title{ 2017 hyperLOPIT on U2OS cells }

\description{

  This data contains 2 replicated hyperLOPIT experiments/gradients (see
  Christoforou et al. 2016 for details on the technology) on U2OS
  cells. The quantitation was obtain using two sets of TMT 10-plex per
  replicate, producing a total of 40 quanitation channels.

  \code{hyperLOPITU2OS2017} contains all 40 fractions while in
  \code{hyperLOPITU2OS2017b}, 3 fractions with low protein yields have
  been remove (see example below).
  
}


\usage{data("hyperLOPITU2OS2017")}

\format{ An object of class \code{MSnSet}, defined in the \code{MSnbase}
  package.}


\details{

  The data (expression and feature variable) contain:

  \itemize{
    \item UniProt Accession for Protein Group (no isoform information):
    Unique UniProt accession for quantified protein group reported by
    Proteome Discoverer (1\% FDR) - isoform information not retained.
    \item Normalized TMT 10-plex Reporter Ion Distribution: ReplicateX
    TMT SetX-126 Normalized TMT 10-plex reporter ion values,
    representing the distribution of each protein across the
    fractionation scheme for each experiment. Protein-level reporter ion
    values were calculated by taking the median of all quantifiable PSMs
    for the protein group, then normalized so that the sum of all 10
    channels was equal to 1. The numeric value in the tag name
    corresponds to the nominal mass of each TMT reporter ion. The N and
    C suffixes differentiates between the 15N or 13C isotopologue
    variants of TMT 10-plex reporter ions of the same nominal mass.
    \item UniProt Accession for Protein Group (with isoform
    information): Unique UniProt accession for quantified protein group
    reported by Proteome Discoverer (1\% FDR) - isoform information
    retained.
    \item UniProt Protein Description: UniProt description for protein
    accession.
    \item Coverage: Percentage of protein sequence covered by identified
    peptides.
    \item Quantified Proteins: Number of quantified protein groups.
    \item Quantified Unique Peptides: Number of unique quantified
    peptides. Only these peptides were used for quantification.
    \item Quantified Peptides: Number of quantified peptides. Only peptides
    that were unique to a single protein group were used for
    quantification.
    \item Quantified PSMs: Number of quantified peptide-spectrum matches.
    \item Score - ReplicateX TMT SetX: Total score of identified protein group
    for each experiment. This score is equal to the sum of the
    individual peptide scores.
    \item Coverage - ReplicateX TMT SetX: Percentage of protein sequence
    covered by identified peptides for each experiment.
    \item Quantified Peptides - ReplicateX TMT SetX: Number of quantified
    peptides for each experiment. Only peptides that were unique to a
    single protein group were used for quantification.
    \item Quantified PSMs - ReplicateX TMT SetX: Number of quantified
    peptide-spectrum matches for each experiment.
    \item SVM Marker Set: Final marker set used for SVM classification of
    protein subcellular localization to 14 subcellular compartments.
    \item SVM Classification: Subcellular class to which the protein group
    was assigned by SVM classification. All proteins are assigned to a
    single class by SVM.
    \item SVM Score: Confidence score for localization assignment, ranging
    from 0 to 1. A score close to 0 represents a very low confidence
    assignment, whereas a score of 1 indicates a very high confidence
    assignment.
    \item Final SVM Classification (5\% FDR) (assignment): Predicted
    localization, with SVM score thresholds determined empirically by
    comparison to GO annotation and protein database annotation. The
    SVM score thresholds were set individually for each class so that
    the false discovery rate of the SVM classification was equal or
    lower than 5\%.
  }  
}



\examples{
data(hyperLOPITU2OS2017)
hyperLOPITU2OS2017

library("pRoloc")
plot2D(hyperLOPITU2OS2017, addLegend = "bottomleft")

## removing low intensity fractions
sort(colSums(exprs(hyperLOPITU2OS2017)))
i <- order(colSums(exprs(hyperLOPITU2OS2017)))[1:3]
x <- hyperLOPITU2OS2017[, -i]
plot2D(x, mirrorY = TRUE)

data(hyperLOPITU2OS2017b)
## only difference if subsetting date
all.equal(hyperLOPITU2OS2017b, x)
processingData(hyperLOPITU2OS2017b)
processingData(x)
}
\keyword{datasets}
