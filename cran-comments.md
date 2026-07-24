## Submission

This is a new submission of topocast to CRAN.

## Test environments

* local Windows 11, R 4.6.0
* win-builder (r-devel)
* GitHub Actions: ubuntu-latest (r-devel, r-release, r-oldrel), macOS-latest (r-release), windows-latest (r-release)

## R CMD check results

0 errors | 0 warnings | 0 notes locally.

win-builder r-devel reports one NOTE, the expected new-submission one. The
incoming checks also list possibly misspelled words in the Description field, all
of which are correct as written:

* "Downscales", "Downscaling" -- standard terms in the climate and remote-sensing
  literature for deriving a fine-resolution field from a coarse one. They are
  ordinary English words here, not software names, so they are left unquoted.
* "Karger", "et", "al" -- the author surname and the standard abbreviation from
  the method reference in the Description field.

## Notes for the reviewer

The package compiles OpenMP-parallel C++ (Rcpp, RcppArmadillo). The fitting loop
requests two threads when `_R_CHECK_LIMIT_CORES_` is set, and the number of
available cores otherwise.

The method reference in the Description field is given as a DOI
(<doi:10.1038/sdata.2017.122>).
