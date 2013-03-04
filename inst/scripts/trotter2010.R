library("MSnbase")

shallow <- read.csv("../extdata/pmic_201000359_sm_Shallow.csv", sep = ";", row.names = 1)
steep <- read.csv("../extdata/pmic_201000359_sm_Steep.csv", row.names = 1)

shallow <- new("MSnSet", exprs = as.matrix(shallow[, 1:8]),
               featureData = new("AnnotatedDataFrame",
                 data = shallow[, 9, drop = FALSE]))

steep <- new("MSnSet", exprs = as.matrix(steep[, 1:8]),
             featureData = new("AnnotatedDataFrame",
               data = steep[, 9, drop = FALSE]))

shallow@experimentData <-
  new("MIAPE",
      title = c(
        "Improved sub-cellular resolution via simultaneous analysis of organelle proteomics data across varied experimental conditions.",
        "Mapping the Arabidopsis organelle proteome."),
      abstract = c(
        "Spatial organisation of proteins according to their function plays an important role in the specificity of their molecular interactions. Emerging proteomics methods seek to assign proteins to sub-cellular locations by partial separation of organelles and computational analysis of protein abundance distributions among partially separated fractions. Such methods permit simultaneous analysis of unpurified organelles and promise proteome-wide localisation in scenarios wherein perturbation may prompt dynamic re-distribution. Resolving organelles that display similar behavior during a protocol designed to provide partial enrichment represents a possible shortcoming. We employ the Localisation of Organelle Proteins by Isotope Tagging (LOPIT) organelle proteomics platform to demonstrate that combining information from distinct separations of the same material can improve organelle resolution and assignment of proteins to sub-cellular locations. Two previously published experiments, whose distinct gradients are alone unable to fully resolve six known protein-organelle groupings, are subjected to a rigorous analysis to assess protein-organelle association via a contemporary pattern recognition algorithm. Upon straightforward combination of single-gradient data, we observe significant improvement in protein-organelle association via both a non-linear support vector machine algorithm and partial least-squares discriminant analysis. The outcome yields suggestions for further improvements to present organelle proteomics platforms, and a robust analytical methodology via which to associate proteins with sub-cellular organelles.",
        "A challenging task in the study of the secretory pathway is the identification and localization of new proteins to increase our understanding of the functions of different organelles. Previous proteomic studies of the endomembrane system have been hindered by contaminating proteins, making it impossible to assign proteins to organelles. Here we have used the localization of organelle proteins by the isotope tagging technique in conjunction with isotope tags for relative and absolute quantitation and 2D liquid chromatography for the simultaneous assignment of proteins to multiple subcellular compartments. With this approach, the density gradient distributions of 689 proteins from Arabidopsis thaliana were determined, enabling confident and simultaneous localization of 527 proteins to the endoplasmic reticulum, Golgi apparatus, vacuolar membrane, plasma membrane, or mitochondria and plastids. This parallel analysis of endomembrane components has enabled protein steady-state distributions to be determined. Consequently, genuine organelle residents have been distinguished from contaminating proteins and proteins in transit through the secretory pathway."),
      pubMedIds = c("21058340", "16618929"),
      samples = list(species = "Arabidopsis thaliana",
        tissue = "Callus",
        gradient = "shallow"),
      lab = "Cambridge Centre for Proteomics (CCP)",
      name = "Kathryn S. Lilley",
      email = "k.s.lilley@bioc.cam.ac.uk",
      url = "http://www.bio.cam.ac.uk/proteomics/")

steep@experimentData <-
  new("MIAPE",
      title = c(
        "Improved sub-cellular resolution via simultaneous analysis of organelle proteomics data across varied experimental conditions.",
        "Sub-cellular localization of membrane proteins."),
      abstract = c(
        "Spatial organisation of proteins according to their function plays an important role in the specificity of their molecular interactions. Emerging proteomics methods seek to assign proteins to sub-cellular locations by partial separation of organelles and computational analysis of protein abundance distributions among partially separated fractions. Such methods permit simultaneous analysis of unpurified organelles and promise proteome-wide localisation in scenarios wherein perturbation may prompt dynamic re-distribution. Resolving organelles that display similar behavior during a protocol designed to provide partial enrichment represents a possible shortcoming. We employ the Localisation of Organelle Proteins by Isotope Tagging (LOPIT) organelle proteomics platform to demonstrate that combining information from distinct separations of the same material can improve organelle resolution and assignment of proteins to sub-cellular locations. Two previously published experiments, whose distinct gradients are alone unable to fully resolve six known protein-organelle groupings, are subjected to a rigorous analysis to assess protein-organelle association via a contemporary pattern recognition algorithm. Upon straightforward combination of single-gradient data, we observe significant improvement in protein-organelle association via both a non-linear support vector machine algorithm and partial least-squares discriminant analysis. The outcome yields suggestions for further improvements to present organelle proteomics platforms, and a robust analytical methodology via which to associate proteins with sub-cellular organelles.",
        "In eukaryotes, numerous complex sub-cellular structures exist. The majority of these are delineated by membranes. Many proteins are trafficked to these in order to be able to carry out their correct physiological function. Assigning the sub-cellular location of a protein is of paramount importance to biologists in the elucidation of its role and in the refinement of knowledge of cellular processes by tracing certain activities to specific organelles. Membrane proteins are a key set of proteins as these form part of the boundary of the organelles and represent many important functions such as transporters, receptors, and trafficking. They are, however, some of the most challenging proteins to work with due to poor solubility, a wide concentration range within the cell and inaccessibility to many of the tools employed in proteomics studies. This review focuses on membrane proteins with particular emphasis on sub-cellular localization in terms of methodologies that can be used to determine the accurate location of membrane proteins to organelles. We also discuss what is known about the membrane protein cohorts of major organelles."),
      pubMedIds = c("21058340", "18780351"),
      samples = list(species = "Arabidopsis thaliana",
        tissue = "Callus",
        gradient = "Steep"),
      lab = "Cambridge Centre for Proteomics (CCP)",
      name = "Kathryn S. Lilley",
      email = "k.s.lilley@bioc.cam.ac.uk",
      url = "http://www.bio.cam.ac.uk/proteomics/")

shallow <- updateSampleNames(shallow)
steep <- updateSampleNames(steep)

if (validObject(shallow))
  save(shallow, file = "../../data/trotter2010shallow.rda",
       compress = "xz", compression_level = 9)

if (validObject(steep))
    save(steep, file = "../../data/trotter2010steep.rda",
         compress = "xz", compression_level = 9)




