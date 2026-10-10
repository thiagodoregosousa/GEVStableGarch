## Resubmission

This is a resubmission of GEVStableGarch, which was archived on CRAN on
2020-10-22. Version 2.0.0 is a complete rewrite with a new interface; the
check problems that led to the archival no longer apply.

## R CMD check results

0 errors | 0 warnings | 1 note. On win-builder (release and R-devel) the only
note is the CRAN incoming feasibility note: "New submission", "Package was
archived on CRAN", and possibly misspelled words in DESCRIPTION, which are the
model acronyms APARCH, GAt and GEV and the word "pluggable" — all intentional.

A local `R CMD check --as-cran` additionally reported one WARNING (the LaTeX
package `inconsolata.sty` was not installed locally) and two NOTEs (HTML Tidy
too old; `V8` unavailable). These are local tooling gaps, not package issues,
and do not occur on the CRAN check machines.

## Test environments

* local: macOS, R release
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (R-devel, release, oldrel-1)

## Reverse dependencies

None (the package is currently archived).
