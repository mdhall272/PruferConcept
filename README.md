# Lighting Up Random Trees

Code and preprint for "Lighting up random trees: the combinatorics of partially observed transmission forests" (Hall, 2026).

## Overview

In infectious disease genomic epidemiology, sampled individuals form clusters and singletons based on transmission linkage. This project develops the combinatorial framework for counting labelled rooted forests with an independence constraint -- where a designated set of *n* "sampled" nodes cannot be adjacent to each other -- and uses it to perform Bayesian inference on:

- **N**: the total size of the samplable infected population
- **k**: the number of independent lineage introductions

given observed cluster size data.

The key mathematical contribution is a modified Prufer sequence representation for forests with independence constraints, yielding a closed-form forest count and enabling uniform sampling over forest structures. This is wrapped in an importance sampling framework for posterior inference over N and k.

## Interactive dashboard

A Shiny dashboard is hosted at:

**https://mdhall-lshtm.shinyapps.io/prufer-concept/**

Enter comma-separated cluster sizes (e.g. `3,3,2,1,1,1` for two clusters of size 3, one of size 2, and three singletons), configure priors on N and k, and run importance sampling inference.

To run locally:

```r
# From the code/ directory
shiny::runApp("app.R")
```

### Requirements

- R (>= 4.0)
- Packages: `shiny`, `tidyverse`, `Rcpp`
- A C++ compiler (for Rcpp)

## Repository structure

```
code/
  app.R              Shiny dashboard
  forest_v2.5.R      Core R functions (sampling, likelihood, posterior)
  forest_rcpp.cpp    C++ backend (Rcpp) for performance
preprint/
  LURT_preprint.pdf  Preprint
  LURT_preprint.tex  LaTeX source
  prufer_refs.bib    Bibliography
```

## Citation

If you use this code, please cite:

> Hall, M. (2026). Lighting up random trees: the combinatorics of partially observed transmission forests. *Preprint.*

## License

See [LICENSE](LICENSE).
