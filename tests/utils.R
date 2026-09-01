library("TopDom")

stop_if_not <- TopDom:::stop_if_not
mcat <- TopDom:::mcat

## stop_if_not() - passing cases return invisibly (NULL)
res <- withVisible(stop_if_not(TRUE))
stopifnot(is.null(res$value), isFALSE(res$visible))
stop_if_not(TRUE, 1 == 1, !FALSE)

## stop_if_not() - failing cases signal an error, and the message contains
## a deparsed version of the failing expression
res <- tryCatch(stop_if_not(1 == 2), error = identity)
stopifnot(inherits(res, "error"))
stopifnot(grepl("1 == 2", conditionMessage(res)))

## Non-scalar / NA / non-logical arguments all fail
res <- tryCatch(stop_if_not(c(TRUE, TRUE)), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(stop_if_not(NA), error = identity)
stopifnot(inherits(res, "error"))

res <- tryCatch(stop_if_not(TRUE, FALSE), error = identity)
stopifnot(inherits(res, "error"))

## Long expressions are truncated with a trailing '....'
res <- tryCatch(
  stop_if_not(is.numeric("a really really really really really long expression here")),
  error = identity
)
stopifnot(inherits(res, "error"))

## mcat() writes to stderr and returns invisibly
err <- utils::capture.output(out <- withVisible(mcat("hello", " ", "world")), type = "message")
stopifnot(isFALSE(out$visible))
stopifnot(any(grepl("hello world", err)))
