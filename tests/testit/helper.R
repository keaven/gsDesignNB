# TRUE if x and y differ (negation of all.equal(), which returns a message when unequal)
not_equal <- function(x, y, ...) !isTRUE(all.equal(x, y, ...))
