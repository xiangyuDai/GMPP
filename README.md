# GMPP

Genomic Mating to predict expected progeny genetic values.

## Installation

``` r
install.packages("devtools")

# Install from GitHub
devtools::install_github("xiangyuDai/GMPP")

# Or install from Git
devtools::install_git("https://github.com/xiangyuDai/GMPP")
```

If the above methods fail, you can manually download the latest source package from the [Releases page](https://github.com/xiangyuDai/GMPP/releases) and install it locally.

## Quick Start

``` r
library(GMPP)

# Load example data
data(simulated_mate)

# Predict expected progeny genetic values for all sire × dam crosses
epv <- GMPP(mating_sire = sire_data, mating_dam = dam_data, method = "additive")
head(epv)
```

## Core Functions

-   **GRM()** — Compute additive and dominance genomic relationship matrices (VanRaden method)
-   **Marker_effect()** — Estimate marker effects via GBLUP (calls `gaston::lmm.aireml`)
-   **GMPP()** — Predict expected progeny genetic values for candidate parental mating combinations, supporting purebred, crossbred, and multi-breed designs

## Input Format

Input matrices have phenotypes in column 1 and SNP genotypes (coded 0/1/2) in the remaining columns. Rows represent individuals; row names are optional.

## References

Aliloo, H. et al. (2017). Including nonadditive genetic effects in mating programs to maximize dairy farm profitability. *Journal of Dairy Science*, 100(2), 1203–1222.

## License

GPL (\>= 3)
