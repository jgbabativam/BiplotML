## Resubmission (v1.1.1)

This is a resubmission after the package was archived on 2023-10-29.
All changes are described in NEWS.md.

## Test environments

* Windows 11 x64, R 4.5.1 (local), devtools::check(cran = TRUE)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Notes on previous check rounds

**`badhess.txt` (resolved):** `optimx` writes this file to the working
directory when it encounters an indefinite Hessian during CG optimization.
It has been suppressed by routing all `optimr()` and `optimx()` calls through
an internal wrapper (`.run_optimr()`) that changes to `tempdir()` before the
call and explicitly deletes `badhess.txt` from both `tempdir()` and the
original working directory in an `on.exit()` handler.

**`future file timestamps: NEWS.md` (infrastructure note):** The check machine
reported "unable to verify current time", meaning it could not reach a time
server to confirm whether the timestamp is truly in the future. This is a
transient infrastructure issue on the check server and is not caused by the
package. `README.md` has been added to `.Rbuildignore` so it is not included
in the built package.

**`optimr` dependency (resolved):** The `optimr` package was removed from CRAN
on 2023-10-29. All calls now use `optimx::optimr()`.
