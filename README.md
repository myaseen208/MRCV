
## `MRCV`: Methods for Analyzing Multiple Response Categorical Variables (MRCVs)

###### Version : [0.5-0](https://myaseen208.com/MRCV/); Copyright (C) 2013-2024: License: [GPL-3](https://www.r-project.org/Licenses/)

##### *Natalie Koziol<sup>1</sup>, and Chris Bilder<sup>2</sup>*

1.  [Nebraska Center for Research on Children, Youth, Families and
    Schools (CYFS), University of Nebraska Lincoln, NE,
    USA](https://cyfs.unl.edu/personnel/bios/koziol-natalie.php)
2.  [Department of Statistics, University of Nebraska-Lincoln, Lincoln
    NE, USA](https://www.chrisbilder.com/)

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/MRCV)](https://cran.r-project.org/package=MRCV)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/MRCV?color=green)](https://CRAN.R-project.org/package=MRCV)
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/myaseen208/MRCV) -->

[![develVersion](https://img.shields.io/badge/devel%20version-0.5-0-orange.svg)](https://github.com/myaseen208/MRCV)

<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/myaseen208/MRCV/total.svg)] -->

[![Project Status:
WIP](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--10--24-yellowgreen.svg)](https://github.com/myaseen208/MRCV)
\*\*\*

## Description

Provides functions for analyzing the association between one single
response categorical variable (SRCV) and one multiple response
categorical variable (MRCV), or between two or three MRCVs. A modified
Pearson chi-square statistic can be used to test for marginal
independence for the one or two MRCV case, or a more general loglinear
modeling approach can be used to examine various other structures of
association for the two or three MRCV case. Bootstrap- and
asymptotic-based standardized residuals and model-predicted odds ratios
are available, in addition to other descriptive information. Statisical
methods implemented are described in Bilder et al. (2000)
[doi:10.1080/03610910008813665\>, Bilder and Loughin (2004)
\<doi:10.1111/j.0006-341X.2004.00147.x\>, Bilder and Loughin (2007)
\<doi:10.1080/03610920600974419\>, and Koziol and Bilder (2014)
\<https://journal.r-project.org/articles/RJ-2014-014/](https://doi.org/10.1080/03610910008813665%3E,%20Bilder%20and%20Loughin%20(2004)%20%3Cdoi:10.1111/j.0006-341X.2004.00147.x%3E,%20Bilder%20and%20Loughin%20(2007)%20%3Cdoi:10.1080/03610920600974419%3E,%20and%20Koziol%20and%20Bilder%20(2014)%20%3Chttps://journal.r-project.org/articles/RJ-2014-014/).

## Installation

The package can be installed from CRAN as follows:

``` r
install.packages("MRCV", dependencies = TRUE)
```

## What’s new

To know whats new in this version type:

``` r
news(package = "MRCV")
```

## Links

[CRAN page](https://cran.r-project.org/package=MRCV)

## Citing `MRCV`

To cite the package use:

``` r
citation("MRCV")
Please, support this project by citing it in your publications!

  Koziol, A. N, Bilder, R. C (2014). _MRCV: Methods for
  Analyzing Multiple Response Categorical Variables (MRCVs)_,
  volume 6. https://rjournal.github.io/.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {MRCV: Methods for Analyzing Multiple Response Categorical Variables (MRCVs)},
    author = {{Koziol} and Natalie A. and {Bilder} and Christopher R.},
    journal = {The R Journal},
    year = {2014},
    note = {https://rjournal.github.io/},
    volume = {6},
    issue = {1},
    issn = {2073-4859},
    pages = {144-150},
  }
```
