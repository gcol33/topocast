pkg <- "C:/Users/GillesC/Documents/dev/topocast"

Rcpp::compileAttributes(pkg)
devtools::document(pkg, quiet = TRUE)

res <- devtools::test(pkg, stop_on_failure = FALSE)
df <- as.data.frame(res)
print(df[, intersect(c("test", "nb", "failed", "error", "warning"), names(df))])
