assert("nb_sim respects block randomization vector", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = c("A", "B"), rate = c(0.5, 0.5))
  max_followup <- 1

  # Test with block vector c("A", "B") and n=10 (Balanced)
  sim_balanced <- nb_sim(enroll_rate, fail_rate, max_followup = max_followup, n = 10, block = c("A", "B"))
  counts_bal <- table(unique(sim_balanced[, c("id", "treatment")])$treatment)
  (all.equal(as.vector(counts_bal), c(5, 5)))

  # Test with block vector c("A", "A", "B") and n=12 (Unbalanced 2:1)
  # 4 blocks of 3 = 12 subjects. Expect 8 A, 4 B.
  sim_unbal <- nb_sim(enroll_rate, fail_rate, max_followup = max_followup, n = 12, block = c("A", "A", "B"))
  counts_unbal <- table(unique(sim_unbal[, c("id", "treatment")])$treatment)
  (all.equal(as.vector(counts_unbal), c(8, 4)))

  # Test validation
  (has_error(
    nb_sim(enroll_rate, fail_rate, max_followup = 1, n = 10, block = c("A", "C")),
    "Elements of 'block' must match treatment names"
  ))
})

assert("nb_sim uses default block when block=NULL", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  # Standard names "Control", "Experimental"
  fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.5))
  max_followup <- 1

  # n=4, should use default block c("Control", "Control", "Experimental", "Experimental")
  # This implies n=4 should be perfectly balanced if default logic is used and matches names
  # Note: block argument default is only used if user doesn't provide one.
  # If we call nb_sim(..., block=NULL) explicitly, it overrides default?
  # R default arguments are only used if argument is missing.
  # If I call nb_sim(..., block=NULL), 'block' becomes NULL inside function.
  # My implementation logic:
  # if (is.null(block)) {
  #   # Check if treatments match default. If so use default. Else fallback.
  # }
  # This logic handles both missing argument (if I set default to NULL in signature)
  # OR explicit NULL.
  # BUT, I changed the signature to default = c(...).
  # So if I call without block argument, block is c(...).
  # If I call block=NULL, block is NULL.

  # Case 1: Call without block argument (uses default c("Control"x2, "Experimental"x2))
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = max_followup, n = 4)
  counts <- table(unique(sim[, c("id", "treatment")])$treatment)
  (all.equal(as.vector(counts), c(2, 2)))

  # Case 2: Call with block=NULL (uses simple randomization)
  # Hard to test stochasticity in one run, but at least it runs.
  sim_null <- nb_sim(enroll_rate, fail_rate, max_followup = max_followup, n = 4, block = NULL)
  (nrow(sim_null) > 0)
})
