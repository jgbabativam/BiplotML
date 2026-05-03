## Resubmission (v1.1.1)

This is a resubmission after the package was archived on 2023-10-29.
All changes are described in NEWS.md.

## Test environments

* Windows 11 x64, R 4.5.1 (local), devtools::check(cran = TRUE)
* win-builder (R-devel)

## R CMD check results

0 errors | 0 warnings | 1 note

**NOTE: New submission / Package was archived on CRAN**
This note is expected for any resubmission of a previously archived package.
The package was archived on 2023-10-29 because it depended on the `optimr`
package, which was removed from CRAN at the maintainer's request. All calls
now use `optimx::optimr()`, which is the current home of that function.
The CRAN repository comment (`X-CRAN-Comment`) is automatically appended by
CRAN infrastructure and cannot be suppressed.

## Notes on previous check rounds

**`badhess.txt` (resolved):** `optimx` writes this file to the working
directory when it encounters an indefinite Hessian during CG optimization.
It has been suppressed by routing all `optimr()` and `optimx()` calls through
an internal wrapper (`.run\\\_optimr()`) that changes to `tempdir()` before the
call and explicitly deletes `badhess.txt` from both `tempdir()` and the
original working directory in an `on.exit()` handler.

**`future file timestamps` (resolved):** The `Date` field in `DESCRIPTION`
has been updated to the current date. In a previous round the check machine
also reported "unable to verify current time"; that is a transient
infrastructure issue on the check server and is not caused by the package.

