# ccPaper

This is the work for all the analysis in the **categoryCompare** publication and make it easier to work with. To allow testing and incorporation of new features for the reviewers and readers, a specific branch of the **categoryCompare** Bioconductor package now exists in the public repository based on the most recent Bioconductor release. This can be installed using 

## Vignettes

If you are just looking for the supplementary data that is available as `vignettes`, they are available as web-pages for this site.

* [Hypothetical Example](http://rmflight.github.io/ccPaper/inst/doc/hypotheticalExample_vignette.html)
* [CROHNS vs UC](http://rmflight.github.io/ccPaper/inst/doc/uc_crohns_analysis.html)
* [Skin vs Muscle Denervation](http://rmflight.github.io/ccPaper/inst/doc/skin_muscle_preprocessing.html)

## Installation 

```r
source("http://bioconductor.org/biocLite.R")
library(BiocInstaller)
biocLite(c("KEGG.db", "GO.db", "categoryCompare")) # install current bioconductor version
install_github("categoryCompare", username="rmflight", ref="paperRev", quick=TRUE) # install the version for ccPaper
```

Note that the `quick=TRUE` part is required if running an up-to-date version of `devtools`.

This will install the package **categoryComparePaper**, thus named so that it does not change any other version of **categoryCompare** previously installed. Note that the **categoryComparePaper** package should be installed prior to installing this package, as this one depends on it.

## Installation

```r
library(devtools)
library(BiocInstaller)
biocLite(c("ALL", "hgu95av.db", "hgu133plus2.db", "limma")) # other required packages
library(devtools)
install_github("ccPaper", "rmflight", quick=TRUE) # install this package
```

