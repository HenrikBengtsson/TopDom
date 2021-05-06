

<div id="badges"><!-- pkgdown markup -->

<a href="https://github.com/HenrikBengtsson/TopDom/actions?query=workflow%3AR-CMD-check"><img border="0" src="https://github.com/HenrikBengtsson/TopDom/workflows/R-CMD-check/badge.svg?branch=develop" alt="Build status"></a></a>
<a href="https://travis-ci.org/HenrikBengtsson/TopDom"><img border="0" src="https://travis-ci.org/HenrikBengtsson/TopDom.svg" alt="Build status"></a></a>
<a href="https://ci.appveyor.com/project/HenrikBengtsson/topdom"><img border="0" src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/TopDom?svg=true" alt="Build status"></a></a>
<a href="https://codecov.io/gh/HenrikBengtsson/TopDom"><img border="0" src="https://codecov.io/gh/HenrikBengtsson/TopDom/branch/develop/graph/badge.svg" alt="Coverage Status"></a></a>

</div>

# TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes 

The TopDom method identifies topological domains in genomes from Hi-C sequence data (Shin et al., 2016).  The authors published an implementation of their method as an R script (two different versions; also available in this package).  This package originates from those original TopDom R scripts and provides help pages adopted from the original TopDom PDF documentation.  It also provides a small number of bug fixes to the original code.


## Citing TopDom - the method and the package

Whenever using the TopDom method, please cite:

* Hanjun Shin, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, Xianghong Jasmine Zhou; TopDom: an efficient and deterministic method for identifying topological domains in genomes, Nucleic Acids Research, Volume 44, Issue 7, 20 April 2016, Pages e70, [10.1093/nar/gkv1505](https://doi.org/10.1093/nar/gkv1505). PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/), PMID: 26704975.

Whenever using the **TopDom** package, please cite:

* Henrik Bengtsson, Hanjun Shin, Harris Lazaris, Gangqing Hu and Xianghong Zhou (2020). R Package TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes. R package version 0.10.0. https://github.com/HenrikBengtsson/TopDom

The above information is also available as plain text as well as BibTeX entries via `citation("TopDom")`.


## Background

Please note that I, Henrik Bengtsson, is _not_ the developer of the TopDom method (Shin et al. 2016), or the (original) TopDom code available from the TopDom authors.   As I needed the TopDom method in a project, I ended up putting their `TopDom_v0.0.*.R` scripts into a proper R package (this **TopDom** package) so I could validate the[ir] code via `R CMD check` and so on. Then I started to add a bit of help documentation to make my own life easier.  I also found a few bugs that I fixed and did some improvements but I really tried not to diverge from the original functionality.

Now, since I find it a waste of resources if someone else has to go through the same efforts that I had to, I decided to make this package public.  The goal is also to submit the **TopDom** package to CRAN so that the TopDom method is properly archived for reproducible purposes.  In order to do this, the original authors agreed on releasing their original TopDom scripts under the GPL license, which is also the license of the **TopDom** package.

Having said this, please note that I don't intend to do user support for the TopDom method, how to use it, add new features, and so on.  This is simply because I don't have the resources to do that.  More importantly, I am not the best person to answer questions on the TopDom method and its implementation.  Instead, I recommend that you reach out to the original TopDom authors to get your questions answered.


## Archaeology

The original R script (versions 0.0.1 and 0.0.2) and the TopDom manual were made available on a website of the TopDom author (http<span/>://zhoulab.usc.edu/TopDom/).  As of 2019, that website is no longer available.  Instead, TopDom is now listed as a paragraph on <https://labs.dgsom.ucla.edu/zhou/pages/software> with a link to one of the author GitHub repository [jasminezhoulab/TopDom](https://github.com/jasminezhoulab/TopDom), which holds the TopDom manual and version 0.0.2 of the script but not version 0.0.1.



[R]: https://www.r-project.org/
[TopDom]: https://github.com/HenrikBengtsson/TopDom/
[TopDomData]: https://github.com/HenrikBengtsson/TopDomData/

## Installation
R package TopDom is only available via [GitHub](https://github.com/HenrikBengtsson/TopDom) and can be installed in R as:
```r
remotes::install_github("HenrikBengtsson/TopDom", ref="master")
```


### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/TopDom", ref="develop")
```
This will install the package from source.  

<!-- pkgdown-drop-below -->

