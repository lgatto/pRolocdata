Spatial proteomics datasets
============================

## The `pRolocdata` package

`pRolocdata` is a [Bioconductor](http://bioconductor.org/)
[experiment package](http://bioconductor.org/packages/release/BiocViews.html#___ExperimentData)
([releast](http://bioconductor.org/packages/release/data/experiment/html/pRolocdata.html)
and
[devel](http://bioconductor.org/packages/devel/data/experiment/html/pRolocdata.html)
pages) that collects published (mainly, although some unpublished
datasets are also available) mass spectrometry-based spatial/organelle
dataset. The data are distributed as `MSnSet` instances (see the
[`MSnbase`](http://www.bioconductor.org/packages/release/bioc/html/MSnbase.html)
for details) and are used throughout the
[`pRoloc`](http://bioconductor.org/packages/release/data/experiment/html/pRolocdata.html)
and
[`pRolocGUI`](http://bioconductor.org/packages/devel/bioc/html/pRolocGUI.html)
software for spatial proteomics data analysis and visualisation.

### Installation




```r
source("http://bioconductor.org/biocLite.R")
biocLite("pRolocdata")
```

Once installed, the package needs to be loaded


```r
library("pRolocdata")
```

### Available datasets


```r
pRolocdata()
```


|Data               |Description                                                                                                   |
|:------------------|:-------------------------------------------------------------------------------------------------------------|
|andy2011           |Data from LOPIT experiment on Human Embryonic Kidney fibroblast cells                                         |
|at_chloro          |The AT_CHLORO data base                                                                                       |
|dunkley2006        |LOPIT data from Dunkley et al. 2006                                                                           |
|foster2006         |PCP data from Foster et al, 2006                                                                              |
|groen2014cmb       |Data from LOPIT experiments on Arabidopsis thaliana callus roots, taken from Groen et al (Accepted, Dec 2013) |
|groen2014r1        |Data from LOPIT experiments on Arabidopsis thaliana callus roots, taken from Groen et al (Accepted, Dec 2013) |
|groen2014r2        |Data from LOPIT experiments on Arabidopsis thaliana callus roots, taken from Groen et al (Accepted, Dec 2013) |
|groen2014r3        |Data from LOPIT experiments on Arabidopsis thaliana callus roots, taken from Groen et al (Accepted, Dec 2013) |
|hall2009           |LOPIT data from Hall et al. 2009                                                                              |
|nikolovski2012     |Data from Nikolovski et al. 2012                                                                              |
|nikolovski2012imp  |Data from Nikolovski et al. 2012                                                                              |
|nikolovski2014     |Data from Nikolovski et al. 2014                                                                              |
|tan2009r1          |LOPIT data from Tan et al. 2009                                                                               |
|tan2009r2          |LOPIT data from Tan et al. 2009                                                                               |
|tan2009r3          |LOPIT data from Tan et al. 2009                                                                               |
|trotter2010shallow |LOPIT data sets used in Trotter et al. 2010.                                                                  |
|trotter2010steep   |LOPIT data sets used in Trotter et al. 2010.                                                                  |
### Loading data

Data is loaded into the `R` session using the `load` function; for
instance, to get the data from
[Dunkley et al (2006)](http://www.pnas.org/content/103/17/6518.abstract),
one would type


```r
print(data(dunkley2006))
```

```
## [1] "dunkley2006"
```

To get more information about a given dataset, see its manual page


```r
?dunkley2006
```

## The datasets

Each data object in `pRolocdata` is available as an `MSnSet`
instance. The instances contain the actual quantitative data, sample
and feautres annotations (see `pData` and `fData`
respectively). Additional MIAPE data
[[1](https://en.wikipedia.org/wiki/Minimum_Information_About_a_Proteomics_Experiment),
[2](http://www.nature.com/nbt/journal/v25/n8/abs/nbt1329.html)]
experimental data is available in the `experimentData` slot, as
described below.

### Required metadata

#### Species
Documented in `experimentData(.)@samples$species`

#### Tissue

Documented in `experimentData(.)@samples$tissue`. If tissue is `Cell
line`, then details about the cell line are available in
`experimentData(.)@samples$cellLine`.

#### PMID
Documented in `pubMedIds(.)`.

#### Spatial proteomics experiment annotation

Documented in `experimentData(.)@other$pRolocMetadata`:
  - `$MS`: iTRAQ8, iTRAQ4, TMT6, LF, SC, ...
  - `$spatexp`: LOPIT, LOPIMS, substractive, PCP, other, PCP-SILAC, ...
  - `$type`: new, meta, both
  - `$markers.fcol`: default is `markers`
  - `$prediction.fcol`

### Adding new data
