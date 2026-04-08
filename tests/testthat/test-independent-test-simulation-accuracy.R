# Monte Carlo vs analytical expected events (nSurv $d) — regression on simulation accuracy

test_that("Monte Carlo mean total events matches nSurv expected d (scalar hazards)", {
  simulate_lf_total_events <- function(x, n_subj = NULL, seed = 1L) {
    if (is.null(n_subj)) {
      n_subj <- round(x$n)
    }
    ratio <- x$ratio
    R <- x$R
    T <- x$T
    lamC <- as.numeric(x$lambdaC)[1L]
    hr <- x$hr
    etaC <- as.numeric(x$etaC)[1L]
    etaE <- as.numeric(x$etaE)[1L]
    accr <- sum(R)
    set.seed(seed)
    enroll <- stats::runif(n_subj, 0, accr)
    is_exp <- stats::rbinom(n_subj, 1L, prob = ratio / (1 + ratio)) == 1L
    lam_fail <- ifelse(is_exp, lamC * hr, lamC)
    eta <- ifelse(is_exp, etaE, etaC)
    Tf <- stats::rexp(n_subj, rate = lam_fail)
    Td <- stats::rexp(n_subj, rate = eta)
    te <- pmin(Tf, Td)
    fail_first <- Tf <= Td
    cal <- enroll + te
    as.integer(sum(fail_first & (cal <= T)))
  }

  x <- gsDesign::nSurv(
    lambdaC = 0.002,
    eta = 5e-4,
    etaE = 5e-4,
    gamma = 10,
    R = 18,
    S = NULL,
    T = 24,
    minfup = 8,
    hr = 0.7,
    hr0 = 1,
    ratio = 1,
    alpha = 0.025,
    beta = 0.1
  )
  nsim <- 1500L
  vals <- vapply(
    seq_len(nsim),
    function(i) simulate_lf_total_events(x, seed = i),
    integer(1)
  )
  m <- mean(vals)
  mc_se <- stats::sd(vals) / sqrt(nsim)
  d <- unname(x$d)
  # Within ~4 MC standard errors of analytical expectation (very likely under correct model)
  expect_true(abs(m - d) < 4 * max(mc_se, 1e-6))
  # Relative error small (design-dependent scale)
  expect_true(abs(m - d) < 0.02 * d)
})

test_that("thin_interevent_gap enforces minimum spacing between counted events", {
  thin <- function(times, delta) {
    if (length(times) == 0L) {
      return(numeric(0))
    }
    out <- times[1L]
    for (i in seq_len(length(times) - 1L) + 1L) {
      if (times[i] >= max(out) + delta) {
        out <- c(out, times[i])
      }
    }
    out
  }
  t <- c(0, 5, 10, 12, 35)
  expect_equal(thin(t, delta = 20), c(0, 35))
  expect_equal(thin(c(100), delta = 20), 100)
})
