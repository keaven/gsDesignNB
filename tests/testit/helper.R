# TRUE if x and y differ (negation of all.equal(), which returns a message when unequal)
not_equal <- function(x, y, ...) !isTRUE(all.equal(x, y, ...))

# R CMD check sets a relative R_TESTS that child R processes (e.g. future
# workers) can't resolve once test_pkg() changes the working directory
Sys.unsetenv('R_TESTS')
