
<!-- README.md is generated from README.Rmd. Please edit that file -->

# abtrac

<!-- badges: start -->
<!-- badges: end -->

`abtrac` is a small companion package to our paper Prevalent and
persistent new-onset autoantibodies in mild to severe COVID-19[^1]. The
package demonstrates the clustering approach taken for classification of
new-onset autoantibody trajectories. For details on the rationale and
use case, see the paper.

The clustering method of `abtrac` is based on the Paritioning Around
Medioids (PAM) algorithm implemented in
[`cluster::pam()`](https://cran.r-project.org/web/packages/cluster/index.html)[^2],
combined with the `cosine*euclidean` distance metric. Prior to using
`abtrac`, antibody trajectories are measured using bead arrays and data
are acquired as median fluorescent intensity, in arbitrary units (MFI
\[AU\]). Using a set point of reference, e.g., seroconversion, MFI are
converted to fold changes (FC). FC are used to cluster trajectories in
`abtrac`.

Please see the vignette (accessible using `browseVignettes("abtrac")`)
for an example of running `abtrac` on the included dataset for antibody
trajectory clustering.

## Installation

### `abtrac`

To install the package:

``` r
if(!requireNamespace("remotes")) install.packages("remotes")

remotes::install_github("jernbom/abtrac")
```

To load the package:

``` r
library(abtrac)
```

[^1]: August Jernbom Falk, Lovisa Skoglund, Elisa Pin, Ronald Sjöberg,
    Hanna Tegel, Sophia Hober, Elham Rostami, Annica Rasmusson, Janet L.
    Cunningham, Sebastian Havervall, Charlotte Thålin, Anna Månberg,
    Peter Nilsson. Prevalent and persistent new-onset autoantibodies in
    mild to severe COVID-19. medRxiv 2024.02.15.24302857; doi:
    <https://doi.org/10.1101/2024.02.15.24302857>

[^2]: Maechler M, Rousseeuw P, Struyf A, Hubert M, Hornik K (2023).
    cluster: Cluster Analysis Basics and Extensions. R package version
    2.1.6
