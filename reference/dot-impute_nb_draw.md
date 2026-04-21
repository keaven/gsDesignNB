# Gamma–Poisson compound draw (internal)

Draws count observations from the Gamma–Poisson mixture that corresponds
to the negative binomial distribution used throughout **gsDesignNB**:
\$\$\lambda_i \sim \text{Gamma}(1/k,\\ \mu_i k), \quad Y_i^{(m)} \mid
\lambda_i \sim \text{Poisson}(\lambda_i).\$\$

## Usage

``` r
.impute_nb_draw(mu, k, n_imp = 1L)
```

## Arguments

- mu:

  Numeric vector of predicted means (\> 0).

- k:

  Numeric scalar. NB dispersion: Var(Y) = mu + k \* mu^2.

- n_imp:

  Integer. Number of independent draws per element of `mu`.

## Value

Integer matrix, `length(mu)` rows × `n_imp` columns.
