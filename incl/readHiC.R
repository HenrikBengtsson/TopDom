path <- system.file("exdata", package = "TopDom", mustWork = TRUE)

## Original count data
chr <- "chr19"
pathname <- file.path(path, sprintf("nij.%s.gz", chr))
data <- readHiC(pathname, chr = chr, binSize = 40e3)
print(data)
str(data)
