## Resubmission

This is a resubmission of GEVStableGarch, which was archived on CRAN on
2020-10-22. Version 2.0.0 is a complete rewrite with a new interface; the
check problems that led to the archival no longer apply.

I am the same maintainer and author as the archived version (Thiago do Rego
Sousa). My email has changed from the old address used for version 1.x
(thiagoestatistico@gmail.com) to my current address
(thiagodoregosousa@gmail.com), which is the maintainer address in DESCRIPTION.

## Changes after the manual review (Konstanze Lauseker)

* DESCRIPTION: all acronyms (ARMA, GARCH, APARCH, ARCH, GEV, GAt) are now
  spelled out.
* Added `\value` to `gs_fit-methods.Rd`, describing the class and meaning of
  each method's result.
* The `GAT.Rd` example now restores `par()` after changing it.
* The internal distribution check no longer seeds the RNG nor writes to
  `.GlobalEnv`; it uses the caller's RNG stream as is.

## R CMD check results

0 errors | 0 warnings | 1 note. On win-builder (release and R-devel) the only
note is the CRAN incoming feasibility note: "New submission", "Package was
archived on CRAN", and possibly misspelled words in DESCRIPTION, which are the
model acronyms APARCH, GAt and GEV and the words "heteroskedasticity" and
"pluggable" — all intentional.

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
