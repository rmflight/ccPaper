# ccPaperRev

This is work to address reviewer comments on the **categoryCompare** publication and make it easier to work with. To allow testing and incorporation of new features for the reviewers, a specific branch of the **categoryCompare** Bioconductor package now exists in the public repository based on the most recent Bioconductor release. This can be installed using 

```r
library(devtools)
install_github("categoryCompare", username="rmflight", ref="paperRev")
```

This will install the package **categoryComparePaperRev**, thus named so that it does not change any other version of **categoryCompare** previously installed. Note that the **categoryComparePaperRev** package should be installed prior to installing this package, as this one depends on it.

## Installation

```r
source("http://bioconductor.org/biocLite.R")
biocLite("categoryCompare") # makes sure all normal CC prereqs are there

library(devtools)
install_github("categoryCompare", username="rmflight", ref="paperRev") # install the dependency package
install_github("ccPaperRev", "rmflight") # install this package
```

