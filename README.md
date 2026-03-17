# IndivSTATIS

[![R-CMD-check](https://github.com/juchiyu/IndivSTATIS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/juchiyu/IndivSTATIS/actions/workflows/R-CMD-check.yaml)\
[![GitHub release](https://img.shields.io/github/v/release/juchiyu/IndivSTATIS)](https://github.com/juchiyu/IndivSTATIS/releases) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**IndivSTATIS** introduces a multivariate framework to analyze brain networks with individualized parcellation that allows the numbers of parcels/networks to vary across individuals.

IndivSTATIS extracts a shared component space, called the *compromise* to quantify functional brain network organization of individually-defined brain nodes and networks.

### Preprint:

[https://www.biorxiv.org/content/10.64898/2025.12.19.695601v1](#0){.uri}

------------------------------------------------------------------------

## IndivSTATIS

![](workflow_IndivSTATIS.png)

## Installation

You can install the development version directly from GitHub:

``` r
install.packages("remotes")
remotes::install_github("juchiyu/IndivSTATIS")

library(IndivSTATIS)
```

## Quick Example

``` r
library(IndivSTATIS)

# Load example data
data(wines2012List)

# Perform Individual STATIS
result <- statis(
  wines2012List,
  Norm = 'MFA',
  center = TRUE,
  scale = "SS1",
  RV = TRUE,
  BigShareDim = FALSE)

# Inspect results
print(result$Cmat) # component space of the RV matrix, which measures similarity between data tables
print(result$Splus) # component space of the compromise
```

## Main Function

`statis()`

Performs STATIS analysis on a list of matrices

**Arguments**

-   `X` — list of correlation matrices between the time series derived from a common reference (e.g., Gordon atlas, the Schaefer 400 atlas, mean time series of common networks, vertices) and from the individualized parcellation. These matrices should be in a size of $I$ (number of the reference nodes) x $J_k$ (number of individualized parcels of the $k$th individual).

-   `Norm` — normalization methods for each table, e.g., "MFA" or "SUMPCA"

-   `center` — logical; whether to center (subtract the means) from each column

-   `scale` — whether to normalize the columns. `TRUE` will normalize the columns by their standard deviations; `"SS1"` will normalize the columns so that each of them will have a sum of squares = 1; and `FALSE` means no normalization.

-   `RV` — logical; whether to use RV coefficients to measure matrix similarity.

-   `BigShareDim` — logical; whether the common reference has high dimensionality, e.g., vertices

**Returns**

-   `Cmat` — the RV matrix and its eigenvalue decomposition. The weights for each table is derived from the first eigenvalue of the RV matrix

-   `Splus` — the compromise and its SVD output, including singular values, proportion of explained variance (tau), row factor scores (Fi) (and partial factor scores), column factor scores (Fj), and column loadings (Q; generalized singular vectors)

### Example Dataset

`wines2012List`

A small example sensory dataset used to demonstrate the method.

-   **Format:** A list of matrices

-   **Rows:** wines

-   **Columns:** sensory attributes

### Dependencies

-   `ExPosition` – core multivariate functions

-   `DistatisR` (github package: `HerveAbdi/DistatisR`) – base computation functions for STATIS-related computation

-   `rsvd` – faster SVD for big dataset

### Notes for Users

-   Ensure your data is formatted as a **list of matrices** with the same number of rows (observations or a common reference).

-   The `statis()` function returns a comprehensive list of matrices and scores for further visualization or analysis.

## License

This package is licensed under **MIT License**.

## Development

This package is under active development. Contributions, suggestions, and bug reports are welcome via GitHub issues.

**GitHub Actions** are set up to automatically run R CMD check on each push.

### 
