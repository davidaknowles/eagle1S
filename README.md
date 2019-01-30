# eagle1S

Re-implementation of EAGLE using Stan. 

## Original EAGLE 

The original EAGLE codebase is here: https://github.com/davidaknowles/eagle
and the accompanying paper is here: http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.4298.html

## Getting started

The easiest way to install will be using `devtools`
```
if (!require("devtools")) install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_github("davidaknowles/eagle1S/eagle1S")
```

You will need the packages `tidyverse`, `magrittr`, `doMC` and `intervals` to run the example code. 

Then clone this repo somewhere and look in `example/`. `simulate.Rmd` simulates a synthetic dataset on chromosome 22 and `fit_per_snp.Rmd` runs EAGLE on the simulated data. This will give you appropriate input files to try to match for your data. 

