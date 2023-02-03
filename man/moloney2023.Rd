\name{moloneyTbBSF}
\alias{moloneyTbBSF}
\alias{moloneyTbPCF}
\alias{moloneyTcBSF}
\alias{moloneyTcPCF}


\docType{data}

\title{

  Spatial proteomics datasets from two African trypanosome species

}

\description{
  Spatial proteomics datasets from a hyperLOPIT experimental design on two on
  two African trypanosome species, \emph{Trypanosoma brucei} and 
  \emph{Trypanosoma congolense}, which have been mapped across two life-stages.
}

\usage{
data(moloneyTbBSF)
data(moloneyTbBSF)
data(moloneyTcBSF)
data(moloneyTcPCF)
}

\format{
  These data are instances of class \code{MSnSet} from package
  \code{MSnbase}.
}

\details{
  Protein function is intimately linked with localisation and as a consequence the
  subcellular distribution of a protein provides information on its role in the
  cell. We have optimised a method for resolving subcellular compartments in
  \emph{Trypanosoma brucei} and \emph{Trypanosoma congolense} and implemented it
  in the spatial proteomics strategy of hyperLOPIT (hyperplexed localisation of
  organelle proteins by isotope tagging) (Christoforou et al. 2016; Mulvey et al.
  2017). Between the vertebrate and insect stages of these parasites, represented
  by bloodstream and procyclic forms respectively, we have detected over 7000
  proteins in each species across three biological iterations. Of these, 6182
  \emph{T. brucei} proteins and 6324 \emph{T. congolense} proteins are included in
  a spatial proteome characterisation (Trotter et al., 2010). Classification to
  19-23 subcellular compartments was performed using a machine learning approach
  based on a T-augmented Gaussian mixture model (Crook et al. 2019; Crook et al.
  2018). With 713-852 compartment marker proteins, this has yielded localisation
  information for 2504-2795 proteins in each organism.
  
  The data (expression and feature variable) contain:

  \itemize{
    \item MW:	TriTrypDB based molecular mass (Aslett et al. 2010)
    \item Signal_peptide:	TriTrypDB based signal peptide prediction.
    \item TM:	TriTrypDB based predicted number of transmembrane domains.
    \item Curated_GO_Processes:	TriTrypDB based curated gene ontology biological 
    process term.
    \item Computed_GO_Processes:	TriTrypDB based computed gene ontology biological 
    process term.
    \item Curated_GO_Components:	TriTrypDB based curated gene ontology cellular 
    component term.
    \item Computed_GO_Components:	TriTrypDB based computed gene ontology cellular 
    component term.
    \item PFam_Description:	TriTrypDB based Pfam description.
    \item HDBSCAN_cluster:	Cluster number according to HDBSCAN unsupervised 
    clustering (Campello et al. 2013)
    \item HDBSCAN_cluster_probability:	Probability of membership associated with 
    HDBSCAN clustering.
    \item NetGPI:	Binary prediction of GPI anchor according to NetGPI (Gíslason 
    et al. 2021)
    \item DeepLoc_location:	Subcellular localisation predicted by DeepLoc (Almagro
    Armenteros et al. 2017)
    \item computed_pI:	pI computed by pIR (Perez-Riverol et al. 2012, Audain et al. 
    2016)
    \item markers:	The marker set used for TAGM-MAP classification of protein 
    subcellular localisation.
    \item tagm.map.allocation:	The TAGM-MAP prediction for the most probable 
    subcellular allocation.
    \item tagm.map.probability:	The posterior probability for the master protein
    subcellular allocations computed by TAGM-MAP.
    \item tagm.map.outlier:	The posterior probability for the master protein to 
    belong to the outlier component rather than any of the annotated components.
    \item tagm.map.localisation.probability:	The localisation probability for the 
    master protein to belong a subcellular class; defined as the product of the
    \item tagm.map.probability: and 1 - tagm.map.outlier
    \item tagm.map.localisation.prediction:	The final prediction of the master 
    protein subcellular localisation based on its localisation probability; only 
    proteins with a localisation probability of greater than 99.9 percent and outlier
    probability of less than 5E-5 were retained.
    \item tagm.mcmc.allocation:	The TAGM-MCMC prediction for the most probable 
    subcellular allocation.
    \item tagm.mcmc.probability:	The mean posterior probability for the master 
    protein subcellular allocations computed by TAGM-MCMC.
    \item tagm.mcmc.probability.lowerquantile:	The lower boundary to the equitailed 
    95-credible interval of tagm.mcmc.probability.
    \item tagm.mcmc.probability.upperquantile:	The upper boundary to the equitailed 
    95-credible interval of tagm.mcmc.probability.
    \item tagm.mcmc.mean.shannon:	A Monte-Carlo averaged Shannon entropy, which is a
    measure of uncertainty in the allocations.
    \item tagm.mcmc.outlier:	The posterior probability for the master protein to 
    belong to the outlier component rather than any of the annotated components. 
    \item tagm.mcmc.joint:	The posterior probability for the master protein 
    allocation to each of the subcellular classes determined by TAGM-MCMC.
 }
}

\source{
  The data was generated by N. Moloney at the Cambridge Centre for Proteomics.
  \url{http://www.bio.cam.ac.uk/proteomics/}.
}

\references{

  Christoforou, A., Mulvey, C.M., Breckels, L.M., Geladaki, A., Hurrell, T.,
  Hayward, P.C., Naake, T., Gatto, L., Viner, R., Martinez Arias, A., and Lilley,
  K.S. (2016). A draft map of the mouse pluripotent stem cell spatial proteome.
  Nat Commun 7, 8992. 10.1038/ncomms9992.

  Crook, O.M., Breckels, L.M., Lilley, K.S., Kirk, P.D.W., and Gatto, L. (2019). A
  Bioconductor workflow for the Bayesian analysis of spatial proteomics.
  F1000Research 8, 446. 10.12688/f1000research.18636.1.

  Crook, O.M., Mulvey, C.M., Kirk, P.D.W., Lilley, K.S., and Gatto, L. (2018). A
  Bayesian mixture modelling approach for spatial proteomics. PLOS Computational
  Biology 14, e1006516. 10.1371/journal.pcbi.1006516.

  Mulvey, C.M., Breckels, L.M., Geladaki, A., Britovsek, N.K., Nightingale,
  D.J.H., Christoforou, A., Elzek, M., Deery, M.J., Gatto, L., and Lilley, K.S.
  (2017). Using hyperLOPIT to perform high-resolution mapping of the spatial
  proteome. Nat Protoc 12, 1110-1135. 10.1038/nprot.2017.026.

  Trotter, M.W., Sadowski, P.G., Dunkley, T.P., Groen, A.J., and Lilley, K.S.
  (2010). Improved sub‐cellular resolution via simultaneous analysis of organelle
  proteomics data across varied experimental conditions. Proteomics 10, 4213-4219.

}

\examples{
  ## load the data
  data(moloneyTbBSF)
  
  ## View a summary of the data MSnSet data container
  moloneyTbBSF
  
  ## PCA plot of the data
  plot2D(moloneyTbBSF, fcol = "tagm.map.localisation.prediction",
         main = "PCA plot map of T. Brucei bloodstream form")
         
  ## Pass pre-computed t-SNE coords and plot the data
  tsne_coords <- as.matrix(fData(moloneyTbBSF)[, c("Dim.1", "Dim.2")])
  plot2D(tsne_coords, fcol = "tagm.map.localisation.prediction",
         method = "none", methargs = list(moloneyTbBSF),
         main = "t-SNE map of T. Brucei bloodstream form")
}

\keyword{datasets}
