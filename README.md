# PermafrostTools
R package for permafrost

## Installation

### From source
```bash
git pull https://github.com/geocryology/PermafrostTools
cd PermafrostTools
Rscript -e "devtools::install()"
```

### From github with `devtools`

Generate a personal authentication token at [https://github.com/settings/tokens](https://github.com/settings/tokens)
```r
library(devtools)
install_github("geocryology/PermafrostTools", ref="main")
```
