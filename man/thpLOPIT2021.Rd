\name{thpLOPIT_lps_mulvey2021}
\alias{thpLOPIT_unstimulated_mulvey2021}

\alias{thpLOPIT_unstimulated_rep1_mulvey2021}
\alias{thpLOPIT_unstimulated_rep2_mulvey2021}
\alias{thpLOPIT_unstimulated_rep3_mulvey2021}
\alias{thpLOPIT_lps_rep1_mulvey2021}
\alias{thpLOPIT_lps_rep2_mulvey2021}
\alias{thpLOPIT_lps_rep3_mulvey2021}

\alias{psms_thpLOPIT_lps_rep1_set1}
\alias{psms_thpLOPIT_lps_rep1_set2}
\alias{psms_thpLOPIT_lps_rep2_set1}
\alias{psms_thpLOPIT_lps_rep2_set2}
\alias{psms_thpLOPIT_lps_rep3_set1}
\alias{psms_thpLOPIT_lps_rep3_set2}
\alias{psms_thpLOPIT_unstim_rep1_set1}
\alias{psms_thpLOPIT_unstim_rep1_set2}
\alias{psms_thpLOPIT_unstim_rep2_set1}
\alias{psms_thpLOPIT_unstim_rep2_set2}
\alias{psms_thpLOPIT_unstim_rep3_set1}
\alias{psms_thpLOPIT_unstim_rep3_set2}

\docType{data}

\title{

  Protein and PMS-level hyperLOPIT datasets from THP-1 human leukaema
  cells

}

\description{
  These are spatial proteomics datasets from a hyperLOPIT experimental
  design on the THP-1 human leukaema cell line. 
}

\usage{
data(thpLOPIT_lps_mulvey2021)
data(thpLOPIT_unstimulated_mulvey2021)
data(thpLOPIT_unstimulated_rep1_mulvey2021)
data(thpLOPIT_unstimulated_rep2_mulvey2021)
data(thpLOPIT_unstimulated_rep3_mulvey2021)
data(thpLOPIT_lps_rep1_mulvey2021)
data(thpLOPIT_lps_rep2_mulvey2021)
data(thpLOPIT_lps_rep3_mulvey2021)
data(psms_thpLOPIT_lps_rep1_set1)
data(psms_thpLOPIT_lps_rep1_set2)
data(psms_thpLOPIT_lps_rep2_set1)
data(psms_thpLOPIT_lps_rep2_set2)
data(psms_thpLOPIT_lps_rep3_set1)
data(psms_thpLOPIT_lps_rep3_set2)
data(psms_thpLOPIT_unstim_rep1_set1)
data(psms_thpLOPIT_unstim_rep1_set2)
data(psms_thpLOPIT_unstim_rep2_set1)
data(psms_thpLOPIT_unstim_rep2_set2)
data(psms_thpLOPIT_unstim_rep3_set1)
data(psms_thpLOPIT_unstim_rep3_set2)
}

\format{
  These data are instances of class \code{MSnSet} from package
  \code{MSnbase}.
}

\details{

  These datasets are from a spatiotemporal proteomic profiling experiment
  (Mulvey et al 2021) in which the dynamic pro-inflammatory response to 
  lipopolysaccharide (LPS) in the THP-1 human leukaemia cell line is mapped.
  
  Triplicate hyperLOPIT experiments using subcellular fractionation  
  were conducted and fractions were digested at either 0h-LPS or 
  12h-LPS post-stimulation. For each replicate 20 fractions were selected 
  and labelled with 2 x TMT 10plex then acquired with SPS-MS3 
  acquisition on the Orbitrap Fusion Lumos Tribrid instrument (Thermo). 
  The resulting datasets were processed by Proteome Discoverer 
  software (v2.1, Thermo) and analysed using the Bioconductor 
  packages \code{MSnbase} and \code{pRoloc}. 

  PSM-level csv files are available in the code{extdata} directory
  and have been imported as \code{MSnSet} instances using 
  \code{readMSnSet2}. In total there are 12 PSM-level datasets, 
  6 for each coniditon and 2 sets per replicate (2 x TMT-10plex).
  A maximum to 2 missing values per PSM was allowed and are
  designated as NA.
  
  There are 8 protein-level datasets; 3 replicates for each condition,
  e.g. \code{thpLOPIT_lps_rep1_mulvey2021}, 
  \code{thpLOPIT_lps_rep2_mulvey2021},  etc. then a final 2 datasets 
  \code{thpLOPIT_lps_mulvey2021} and \code{thpLOPIT_unstim_mulvey2021}
  in which the 3 replicates per condition have been concatenated. 
  This last 2 datasets form part of the main analysis in the 
  manuscript from Mulvey et al 2021. 
  
  The protein-level data was generated from the PSM-data in the 
  following steps: any PSMs with missing values were assessed 
  by examining if there was any trend in missing values. 
  Barplots of the data suggest the few missing values that
  appear accumulate in the first few fractions and we deduce
  they are not missing at random and on the whole reflect the 
  gradient distributions. A left-cenored method was used for
  imputation using \code{MinDet} in \code{MSnbase}. PSMs
  were normalised by sum across the fractions and then 
  combined to protein according to the median PSM per protein 
  group. 
  
  Marker proteins were annotated based on the combined 
  protein-level data in \code{thpLOPIT_lps_mulvey2021} and
  \code{thpLOPIT_unstimulated_mulvey2021} and can be found
  in the \code{fData} slot. A list of well 
  annotated, unambiguous resident organelle marker proteins 
  from 11 subcellular niches: mitochondria, ER, Golgi apparatus, 
  lysosome, peroxisome, PM, nucleus, nucleolus, chromatin, 
  ribosome and cytosol, were curated from the Uniprot database 
  (The UniProt, 2017), Gene Ontology (Ashburner et al., 2000) 
  and from mining the literature. Only proteins known to localise 
  to a single location were included as markers. The
  processing script is \code{scripts/thpLOPIT2021.R}.
  
  The \code{fData} slot of the two combined datasets also contains
  the results from the data analysis as described in Mulvey et al
  (2021). 
  
  
   The data (expression and feature variable) contain:

  \itemize{
    \item UniProt Accession: found in the \code{featureNames} e.g.
    \code{featureNames(thpLOPIT_lps_mulvey2021)}. The Protein Group 
    (no isoform information):Unique UniProt accession for quantified 
    protein group reported by Proteome Discoverer (1\% FDR).
    \item Expression data slot: Normalized TMT 10-plex reporter 
    distributions representing the normalised abundance of each protein 
    across the fractionation scheme for each experiment. 
    Protein-level reporter ion values were calculated by taking the 
    median of all quantifiable PSMs for the protein group, then normalized 
    so that the sum of all 10 channels was equal to 1 before being
    concatenated across replicates.
    \item GN: Gene name for protein accession.
    \item Description: UniProt description for protein
    accession.
    \item Confidence_x: The confidence level of protein identification 
    FDR determined in hyperLOPIT experiment. Only proteins with Medium 
    (Q ² 5 \%) and High (Q ² 1 \%) FDR confidence levels were retained; 
    Percolator v. 2.05 was used to determine FDR; for details, see 
    L. KSll et al., Nat. Methods 2007, 4, 923-925, L. KSll et al., 
    J. Proteome Res. 2008, 7, 29-34, and L. Kall et al., Bioinformatics 
    2008, 24, i42-i48.
    \item Coverage: Percentage of protein sequence covered by 
    identified peptides for each experiment.
    \item Quantified Peptides:  Number of quantified peptides for each 
    experiment. Only peptides that were unique to a single protein
    group were used for quantification.
    \item Quantified PSMs: Number of quantified peptide-spectrum 
    matches for each experiment.
    \item tagm.xxx.xx: \code{TAGM} allocation results from the 
    Baysian T-augmented Gaussian Mixture modelling approach as 
    described in Crook et al. (2018). See \code{??tagmMcmcTrain}.
    \item L2_distance: the natural L2 distance  between the 
    TAGM joint posterior probabilities
    \item non_movers: proteins predicted to not change location
    \item type1_translocation: proteins predicted from one organelle 
    class in the unstimulated condition to a different organelle 
    class in the LPS-stimulated dataset i.e. organelle to organelle
    \item type2_translocation: proteins precicted to move from 
    an unknown localisation in the unstimulated dataset to a 
    predicted organelle class in the LPS-stimulated dataset 
    i.e. unknown to organelle
    \item type3_translocation: proteins predicted to move from a 
    organelle localisation in the unstimulated dataset to an 
    unknown location in the LPS-stimulated dataset i.e. 
    organelle to unknown
    \item type4_translocation: a translocation event within the 
    unknown class i.e. a protein that exhibits a large change 
    between posterior probabilities in both conditions and is 
    classified to an unknown location. For more information
    please see Mulvey et al 2021.
    \item markers: annotated protein location based on the combined 
    protein-level data used for training the \code{TAGM} MCMC
    classifier.
    \item localisation.prob: assignment score, product of the
    \code{tagm.mcmc.probability} * 1 - \code{tagm.mcmc.outlier}.
    \item localisation.pred: predicted localisation as filtered 
    by the \code{localisation.prob}. Proteins were assigned 
    the localisation predicted from \code{tagm.mcmc.allocation}
    if their \code{localisation.prob} was lower than .99 (1\% FDR)
 }
}

\source{
  The data was generated by C. Mulvey and L. Breckels in the
  Cambridge Centre for Proteomics.
  \url{http://www.bio.cam.ac.uk/proteomics/}.
}

\references{

  \emph{Using hyperLOPIT to perform high-resolution mapping of the spatial proteome}
  Mulvey, Claire M and Breckels, Lisa M and Geladaki, Aikaterini and Britovšek, 
  Nina Kočevar and Nightingale, Daniel J H and Christoforou, Andy and Elzek, 
  Mohamed and Deery, Michael J and Gatto, Laurent and Lilley, 
  Kathryn S Nature Protocols 12, 1110–1135 (2017). https://doi.org/10.1038/nprot.2017.026

  \emph{A Bayesian Mixture Modelling Approach For Spatial Proteomics}
  Oliver M Crook, Claire M Mulvey, Paul D. W. Kirk, Kathryn S Lilley,
  Laurent Gatto bioRxiv 282269; doi: https://doi.org/10.1101/282269
}

\examples{
## load a THP data (unstimulated combined tripcate)
data(thpLOPIT_unstimulated_mulvey2021)
thpLOPIT_unstimulated_mulvey2021
pData(thpLOPIT_unstimulated_mulvey2021)
head(exprs(thpLOPIT_unstimulated_mulvey2021))

## simple protocol for combining psm to protein 
## load LPS stimulated data for replicate 1

library("MSnbase")
data(psms_thpLOPIT_lps_rep1_set1)
data(psms_thpLOPIT_lps_rep1_set2)

## impute missing values with "MinDet"
set1 <- impute(psms_thpLOPIT_lps_rep1_set1, method = "MinDet")
set2 <- impute(psms_thpLOPIT_lps_rep1_set2, method = "MinDet")

## normalise to sum
set1 <- normalise(set1, "sum")
set2 <- normalise(set2, "sum")

## combine to protein
set1 <- combineFeatures(set1,
		     groupBy = fData(set1)$Master.Protein.Accessions,
		     method = median)
set2 <- combineFeatures(set2,
		     groupBy = fData(set2)$Master.Protein.Accessions,
		     method = median)
		     
## update fvarLabels for set 2 to differentiate them from set 1   
set2 <- updateFvarLabels(set2)

## combine sets to form one replicate
xx <- combine(set1, set2)

## keep on proteins common in both sets
xx <- filterNA(xx)

library("pRoloc")
plot2D(thpLOPIT_unstimulated_mulvey2021, main = "Protein-level unstimulated data")
}

\keyword{datasets}
