# RegionScan

## Introduction
`RegionScan` is a R package for comprehensive and scalable region-level genome-wide association testing of alternative region-level multiple-variant and single-variant statistics and visualization of the results. It implements various state-of-the-art region-level tests to improve signal detection under heterogeneous genetic architectures and comparison of multiple-variant region-level and single-variant test results. It leverages LD-based genomic partitioning for LD-adaptive region definition. `RegionScan` is compatible with VCF input file format, and accommodates parallel region-level processing and analysis to improve computational time and memory efficiency. It also provide options for analysis of multi-allelic variants, unbalanced binary phenotypes, with detailed outputs for results interpretation, and utility functions for visualization, comparison, and interpretation.

The detailed information about `RegionScan` can be found in our paper '?submitted to biorxiv'.

## Download and install
Use the `devtools` package (available from
[CRAN](http://cran-r.c3sl.ufpr.br/web/packages/devtools/index.html)) to
install automatically from this GitHub repository:

```{r, eval=TRUE}
library(devtools)
install_github("brossardMyriam/RegionScan")
```

## Illustration of `RegionScan` capabilities
[!C:\Users\bross\OneDrive\Haplotypes_CLQD_MLC\RegionScan_paper\sent_to_coauthors_May30/Fig1_version1_sept13.pdf]

## Authors
- Myriam Brossard (brossard@lunenfeld.ca), author and main developer
- Yun J Yoo, contributor

## Documentation 
The reference and manuel PDF can be found here: ...

## Citation
RegionScan....
MLC
gpart

## License
This package is released under the [GNU General Public License (GPL) v3.0.](https://www.gnu.org/licenses/gpl-3.0.html)
