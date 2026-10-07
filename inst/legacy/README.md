# Legacy analysis scripts

The original flat R scripts that this package was built from, kept for
provenance:

- `Analysis.R` — real-data (SCR/SNCR ozone) analysis.
- `Real_DAPS_and_Recover_U.R` — simulation study (data-generating process).
- `Analysis_plots.R` — figures for the report.
- `analysis_dat.RData` — the ozone analysis data in R format.

The same ozone data in spreadsheet form is bundled at
`inst/extdata/analysis_dat.xlsx` and is what the package's application vignette
loads via `system.file("extdata", "analysis_dat.xlsx", package = "spaci")`.

These scripts are superseded by the package's documented, tested functions
(`idaps()`, `recoverUplus()`, etc.); they are not run as part of the package.
