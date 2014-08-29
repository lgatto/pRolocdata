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

### Loading data

Data is loaded into the `R` session using the `load` function; for
instance, to get the data from
[Dunkley et al (2006)](http://www.pnas.org/content/103/17/6518.abstract),
one would type

```r
data(dunkley2006)
```

To get more information about a given dataset, see its manual page

```r
?dunkley2006
```

## The datasets

Each data object in `pRolocdata` is available as an `MSnSet` instance, that 

* pRolocmetadata
** species 
   see experimentData(.)@samples
** tissue
   see experimentData(.)@samples
   If tissue is Cell, then 

- PMID: see pubMedIds(.)

- experimentData(.)@other$pRolocMetadata
  - $MS: iTRAQ8, iTRAQ4, TMT6, LF, SC, ...
  - $spatexp: LOPIT, LOPIMS, substractive, PCP, other, PCP-SILAC, ...
  - $type: new, meta, both
  - $fcol: default is markers

### Adding new data
