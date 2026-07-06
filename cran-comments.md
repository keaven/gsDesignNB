## Resubmission

This is a resubmission. In this version I have reworded the DESCRIPTION to
address the CRAN incoming spell-check NOTE for "unblinded".
I also reduced the number of illustrative simulation replicates in one vignette
to address the Windows incoming checktime NOTE.

## Test environments

* Local macOS Tahoe 26.5.1, R 4.6.1

## R CMD check results

0 errors | 0 warnings | 1 note

* checking HTML version of manual ... NOTE
  Skipping checking HTML validation: local `tidy` does not look like a recent
  enough HTML Tidy installation. This is a local system-tooling note.

## Package size

The source tarball is approximately 4.8 MB. The previous package-size NOTE was
resolved by replacing CRAN-bundled raw simulation caches with compact
precomputed summaries used by the simulation vignettes.

Installed size reported by local `R CMD check --as-cran`:

* total: 7.5 MB
* doc: 5.2 MB
* extdata: 1.9 MB

## Changes

This release updates score-test sizing guidance, group sequential simulation
documentation, and sample-size examples. It also adds support for harm-bound
group sequential designs available in `gsDesign` 3.10.0 and reduces the CRAN
package size by using compact summary caches for simulation articles.
