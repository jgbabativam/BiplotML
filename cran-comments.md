## Resubmission notes (v1.1.1)

This is a resubmission after the package was archived on 2023-10-29.

### Root cause of archiving
The `optimr` package was removed from CRAN on 2023-10-29 (merged into `optimx`).
All calls to `optimr()` now use `optimx::optimr()`. The `optimr` package has
been removed from `Imports` entirely.

### Dependency changes
* `optimr` removed from `Imports` (package no longer on CRAN; function absorbed
  into `optimx`).
* `optimr`, `optimx`, and `shapes` moved from `Depends` to `Imports`.
* `RSpectra` moved from `Suggests` to `Imports` (required for large matrices).

### Documentation and code fixes
* `DESCRIPTION` Title in title case; Description includes `<doi:>` reference.
* `utils::globalVariables()` moved to `R/zzz.R` (previously caused a roxygen
  block error in `plotBLB.R`).
* `importFrom(stats,svd)` and `importFrom(stats,eigen)` removed — these
  functions belong to base R, not `stats`.
* All ggplot2 `aes()` calls use `.data[[]]` tidy evaluation to avoid undefined
  global variable NOTEs.
* Non-ASCII characters (em-dash, en-dash, e-grave) replaced with ASCII in all
  R source files.
* `bootBLB()`: fixed column-name bug where `Ahat` received `Dimb1`/`Dimb2`
  instead of `Dim1`/`Dim2`, causing `plotBLB()` to fail with
  "Column 'Dim1' not found".
* `bootBLB()`: suppressed verbose optimr CG control messages printed to console.
* Fixed `Sensitivy` typo -> `Sensitivity` in confusion matrix output.
* `inst/CITATION` updated from deprecated `citEntry()` to `bibentry()`.

### Timestamp note
The "future timestamps" warning on `LICENSE` and `README.md` is caused by
git not preserving file modification times. Files were re-touched before
submission. This is a known artefact of git-based workflows and does not
affect package functionality.

## R CMD check results

0 errors | 0 warnings | 1 note

* NOTE: New submission (package was previously archived on 2023-10-29).
