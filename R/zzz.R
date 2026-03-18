# Package-level declarations
# This file is loaded last (alphabetically) and contains
# internal declarations that must not be processed as
# documentation blocks by roxygen2.

# Suppress R CMD check NOTEs for ggplot2 tidy-evaluation column
# names used inside aes() in plotBLB().
utils::globalVariables(c(
  ".data", "label",
  "x.50",  "y.50",
  "x.end", "y.end",
  "ind",   "resample",
  "Dim1",  "Dim2",
  "Dimb1", "Dimb2",
  "b0",    "b1",    "b2",
  "bb0",   "bb1",   "bb2",
  "Dim1c", "Dim2c"
))
