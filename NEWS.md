# Version (development version)

## Documentation

 * Fix equation format issue and update one URL.
 

# Version 0.10.1 [2020-05-04]

## CRAN Policies

 * Update URL that otherwise redirects HTTP to HTTPS.
 

# Version 0.10.0 [2020-05-03]

## New Features

 * Package now includes mouse Chr19 data from the TopDom study. They
   can be found in the `system.file("exdata", package = "TopDom")`
   folder.

 * The orignal TopDom scripts `TopDom_v0.0.1.R` and `TopDom_v0.0.2.R`
   are now distributed part of the package as-is.  They can be found
   in the `system.file("original-scripts", package = "TopDom")`
   folder.

 * `readHiC()` gained arguments `...` which is passed as-is to
   `read.table()`.
 

# Version 0.9.1 [2020-04-02]

## Bug Fixes

 * `ggCountHeatmap()` for `TopDomData` could produce a warning on a
   partial argument name.


# Version 0.9.0 [2020-02-19]

## Significant Changes

 * The list data.frame elements returned by `overlapScores()` now has
   column `chromosome` as the first position.  The data.frame:s are of
   kind `tibble`.
 
## New Features

 * Add `as_tibble()` for `TopDomOverlapScores`.
 
## Documentation

 * Add further documentation on the `window.size` parameter.

 * Add reference to Hanjun Shin's PhD thesis.
 

# Version 0.8.2 [2019-12-17]

## Documentation

 * Improved help on `overlapScores()` and `TopDom()`.

 * Provide a reference for the default value for `window.size` of
   `TopDom()`.
 

# Version 0.8.1 [2019-06-19]

## Software Quality

 * Fix two cases of partial argument name.


# Version 0.8.0 [2019-04-11]

## New Features

 * The TopDom object returned by `TopDom()` now has an attribute
   `parameters` which records the value of arguments `window.size` and
   `statFilter`.

 * Made `TopDom()` faster and more memory efficient by lower the
   number of replicated computations.
 

# Version 0.7.1 [2019-04-09]

## Robustness

 * Add internal sanity checks to `TopDom()` asserting that the
   intermediate and final results are of proper length and does not
   contain missing values.
 
## Bug Fixes

 * Internal `Convert.Bin.To.Domain.TMP()` used by `TopDom()` could
   produce `Error in `[<-.data.frame`(`*tmp*`, , "to.coord", value =
   c(NA, 2500, 247500 : replacement has 3 rows, data has 1`, because
   it assumed at least one domain was identified.


# Version 0.7.0 [2019-03-15]

## Significant Changes

 * Renamed fields returned by `overlapScores()` to be in singular
   form, e.g.  `best_score` instead of `best_scores`.

## New Features

 * Now `overlapScores()` returns also the lengths of the reference
   domains.
 

# Version 0.6.0 [2019-01-08]

## Significant Changes

 * Renamed and swapped the order of the first two arguments of
   `overlapScores()` and renamed the second argument to `reference`.
   This was done in order to make it clear which set of topological
   domains the overlap scores are calculated relative to.
 

# Version 0.5.0 [2018-11-12]

## Significant Changes

 * Lead TopDom author Xianghong Jasmine Zhou confirms by email that
   the original TopDom scripts, and thereby this package, may be
   released under the GNU Public License (GPL).
 

# Version 0.4.0 [2018-10-23]

## New Features

 * Add `countsPerRegion()` for calculating the total contact-frequency
   counts per region specified, e.g. per domain.

 * Add `print()`, `dim()`, `[()`, and `subsetByRegion()` for TopDom
   objects where the number of rows in the dimension reflect the
   number of TopDom domains.
 

# Version 0.3.0 [2018-10-08]

## New Features

 * The legacy `TopDom()` functions, available via `legacy()`, also
   accept `TopDomData` objects as returned by `readHiC()`.  This is
   supported mostly to make it possible to efficiently compare the
   different implementations.

 * Added `[()` for `TopDomData` objects, e.g. `tdd[1:100]`.
 
 * Added `subsetByRegion()` for `TopDomData` objects.

 * Added `ggCountHeatmap()`, `ggDomain()`, and `ggDomainLabel()` for
   `TopDomData` objects.


# Version 0.2.0 [2018-07-26]

## Significant Changes

 * See 'BUG FIXES' below.
 
## New Features

 * Add function `legacy()` for access to the original TopDom v0.0.1
   and TopDom v0.0.2 implementations,
   e.g. `TopDom::legacy("0.0.1")$TopDom()`.

## Bug Fixes

 * All TopDom functions except `TopDom()` itself were the ones from
   TopDom v0.0.2.


# Version 0.1.2 [2018-06-28]

## Significant Changes

 * Released on GitHub.
 
## Miscellaneous

 * The license of the underlying TopDom R code/scripts is still unknown,
   i.e. to be decided by the original authors of TopDom.  Any mentioning
   of code licenses in this package / git repository history is invalid.
   
## New Features

 * Add `overlapScores()`.

 * Add `image()` for `TopDomData`.
 
 * List returned by `TopDom()` gained class `TopDom`.

 * Added logical option `TopDom.debug`, which controls whether
   functions produce debug output or not.  The default is FALSE.

## Documentation

 * Updated the 'Value' section of `help("TopDom")` with details from
   the TopDom Manual (an online PDF) provided by Shin et al.
 

# Version 0.1.1 [2018-05-04]

## New Features

 * Add `print()` method for `TopDomData` object.

 * Reference the TopDom paper (Shin et al., 2016) in the help and the
   README.
 

# Version 0.1.0 [2018-04-23]

## Significant Changes

 * Turned the original TopDom R script into a package.

 * All progress messages are outputted done to standard error.

 * Add `readHiC()`.

## New Features

 * `TopDom()` can now read, via `readHiC()`, a pure count matrix file
   without bin information.  To read such files, specify what
   chromosome is being read (argument `chr`) and the bin size of the
   count matrix (argument `binSize`).

 * If the matrix file is not of a known format, then `TopDom()`
   produces an informative error. Previously it gave a message on
   stdout and returned 0.

## Code Style

 * Tidied up code.

## Software Quality

 * Added package tests.
 

# Version 0.0.2 [2016-07-08]

 * TopDom v0.0.2 script from http://zhoulab.usc.edu/TopDom/ with the below
   entries from the official release note:
   
 * Gap Identification module is changed.
   
 * Minor bug related to Change Points identification in very small
   regions is fixed.
     
 * bed format support.
   

# Version 0.0.1 [2015-11-09]

 * TopDom v0.0.1 script from http://zhoulab.usc.edu/TopDom/.
