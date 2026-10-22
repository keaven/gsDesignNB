# `expect_equal()` analog: TRUE/FALSE for use inside `assert()`
is_equal <- function(x, y, ...) isTRUE(all.equal(x, y, ...))
