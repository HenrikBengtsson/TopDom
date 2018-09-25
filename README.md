# TopDom - An efficient and Deterministic Method for identifying Topological Domains in Genomes

This is an R package that originates from the original [TopDom](http://zhoulab.usc.edu/TopDom/) R script by Shin et al.  This package was put together by Henrik Bengtsson independently of the TopDom authors - please see the [TopDom website] for the original code.  The TopDom article is:

Shin et al., TopDom: an efficient and deterministic method for identifying topological domains in genomes, Nucleic Acids Res. 2016 Apr 20; 44(7): e70., 2015.
doi: 10.1093/nar/gkv1505, PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/), PMID: 26704975.

_Note_, the license for the original TopDom code is ~~unknown~~ to be decided by Shin et al. (private conversation on 2017-11-21).


## Installation

This package is only available from GitHub - it is neither on CRAN nor on Bioconductor. To install this package, use:
```r
remotes::install_github("HenrikBengtsson/TopDom")
```

To get the example data, install the [TopDomData] package.


## Contributions

This Git repository uses the [Git Flow](http://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/TopDom/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/TopDom) branch contains the code of the latest release.

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [TopDom repository](https://github.com/HenrikBengtsson/TopDom).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/TopDom">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/TopDom">AppVeyor CI</a> when the PR is submitted.


## Software status

| Resource:     | CRAN                | Travis CI       | Appveyor         |
| ------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | | <a href="https://travis-ci.org/HenrikBengtsson/TopDom"><img src="https://travis-ci.org/HenrikBengtsson/TopDom.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/TopDom"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/TopDom?svg=true" alt="Build status"></a> |
| Test coverage | | <a href="https://codecov.io/gh/HenrikBengtsson/TopDom"><img src="https://codecov.io/gh/HenrikBengtsson/TopDom/branch/develop/graph/badge.svg" alt="Coverage Status"/></a> | |


[R]: https://www.r-project.org/
[TopDom]: https://github.com/HenrikBengtsson/TopDom/
[TopDomData]: https://github.com/HenrikBengtsson/TopDomData/
[TopDom website]: http://zhoulab.usc.edu/TopDom/
