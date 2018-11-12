# TopDom - An efficient and Deterministic Method for identifying Topological Domains in Genomes

This is an R package that originates from the original [TopDom](http://zhoulab.usc.edu/TopDom/) R script by Shin et al.  This package was put together by Henrik Bengtsson independently of the TopDom authors - please see the [TopDom website] for the original code.

When using the TopDom method, or this TopDom package, please cite:

* Hanjun Shin, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, Xianghong Jasmine Zhou; TopDom: an efficient and deterministic method for identifying topological domains in genomes, Nucleic Acids Research, Volume 44, Issue 7, 20 April 2016, Pages e70, [10.1093/nar/gkv1505](https://doi.org/10.1093/nar/gkv1505). PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/), PMID: 26704975.

This package is licensed under the [GNU Public License (GPL)](https://www.gnu.org/licenses/gpl.html) (per email from XJ Zhou on 2018-11-12).


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
