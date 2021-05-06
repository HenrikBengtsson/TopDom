#' Data for the TopDom Package
#' 
#' The \file{exdata/} folder of this package provides a example data set
#' used in examples.  The data are also used to validate the \pkg{TopDom}
#' implementation toward the original TopDom scripts.
#' 
#' @section Origin:
#' The data herein contain a tiny subset of the HiC and TopDom data used
#' in the TopDom study (Shin et al., 2016).
#' More precisely, it contains:
#' 
#' 1. A TopDom file \file{mESC_5w_chr19.nij.HindIII.comb.40kb.domain}, which
#' is part of the \file{mESC_5w_domain.zip} file
#' (5,504 bytes; md5 ffb19996f681a4d35d5c9944f2c44343) from the
#' Supplementary Materials of Shin et al. (2016).
#' These data were downloaded from the
#' TopDom website (http://zhoulab.usc.edu/TopDom/ - now defunct).
#' 
#' 2. A normalized HiC-count matrix file \file{nij.chr19.gz}, where the
#' non-compressed version is part of the \file{mESC.norm.tar.gz} file
#' (1,305,763,679 bytes; md5 2e79d0f57463b5b7c4bf86b187086d3c) originally
#' downloaded from the
#' [UCSD Ren Lab](http://chromosome.sdsc.edu/mouse/hi-c/download.html).
#' It is a tab-delimited file containing a 3250-by-3250 numeric matrix
#' non-negative decimal values. The underlying HiC sequence data is
#' available from
#' [GSE35156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35156)
#' on GEO and was published part of Dixon, et al. (2012).
#'
#' @references
#' 1. Dixon JR, Selvaraj S, Yue F, Kim A, et al. Topological domains in
#'    mammalian genomes identified by analysis of chromatin interactions.
#'    Nature 2012 Apr 11; 485(7398):376-80,
#'    doi: 10.1038/nature11082,
#'    PMCID: [PMC3356448](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3356448/),
#'    PMID: 22495300.
#' 
#' 2. Shin, et al., TopDom: an efficient and deterministic method for
#'    identifying topological domains in genomes,
#'    Nucleic Acids Res. 2016 Apr 20; 44(7): e70., 2016.
#'    doi: 10.1093/nar/gkv1505,
#'    PMCID: [PMC4838359](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4838359/),
#'    PMID: 26704975.
#'
#' @name TopDom-data
#' @aliases TopDom-data
#' @docType data
#' @keywords data
NULL
