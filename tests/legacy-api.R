library("TopDom")

## Known versions return an environment with the legacy API ------------
api1 <- legacy("0.0.1")
api2 <- legacy("0.0.2")
stopifnot(is.environment(api1), is.environment(api2))
stopifnot(is.function(api1$TopDom), is.function(api2$TopDom))
stopifnot(is.function(api1$Detect.Local.Extreme))

## Unknown version -> error -------------------------------------------
res <- tryCatch(legacy("9.9.9"), error = identity)
stopifnot(inherits(res, "error"),
          grepl("Unknown TopDom legacy version", conditionMessage(res)))
