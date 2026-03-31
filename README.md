# Lighting Up Random Trees

Code and preprint for "Lighting up random trees: the combinatorics of partially observed transmission forests" (Hall, 2026).

## Overview

In infectious disease genomic epidemiology, sampled individuals form clusters and singletons based on transmission linkage. This project develops the combinatorial framework for counting labelled rooted forests with an independence constraint -- where a designated set of *n* "sampled" nodes cannot be adjacent to each other -- and uses it to perform Bayesian inference on:

- **N**: the total size of the samplable infected population
- **k**: the number of independent lineage introductions

given observed cluster size data.

The key mathematical contribution is a modified Prufer sequence representation for forests with independence constraints, yielding a closed-form forest count and enabling uniform sampling over forest structures. This is wrapped in an importance sampling framework for posterior inference over N and k.

## Installation

```r
# install.packages("remotes")
remotes::install_github("mdhall272/PruferConcept")
```

## Usage

```r
library(PruferConcept)

data <- c(3, 3, 2, 1, 1, 1)

# Define a proposal distribution
samp_dist <- sampling_distribution(nb.mean = 18, nb.size = 0.5, components.mean = 3)

# Generate importance samples
samples <- generate_samples(10000, n.sampled = length(data), sampling.dist = samp_dist)

# Compute posterior weights (Poisson priors on N-n and k-1)
result <- posterior(data, samples, prior.N.mean = 18, prior.k.mean = 2)

# Normalise weights
lw <- result$log.weight - max(result$log.weight)
w <- exp(lw)
w <- w / sum(w)

# Posterior summaries
sum(w * result$actual.size)  # E[N]
sum(w * result$components)   # E[k]
1 / sum(w^2)                 # Effective sample size
```

For faster inference on larger datasets, use the C++ batch functions directly:

```r
samples <- generate_samples_batch(50000L, n_sampled = 6L,
    nb_mean = 18, nb_size = 0.5, components_mean = 3)

result <- compute_posterior_batch(as.integer(data), samples, list(
    prior_N_dist = "poisson", prior_N_mean = 18, prior_N_size = 1,
    prior_k_dist = "poisson", prior_k_mean = 2, prior_k_size = 1,
    nb_mean = 18, nb_size = 0.5, components_mean = 3,
    fix_N = FALSE, fix_k = FALSE, cluster = FALSE
))
```

## Interactive dashboard

A Shiny dashboard is hosted at:

**https://mdhall-lshtm.shinyapps.io/prufer-concept/**

The dashboard source is in `code/` and can be run locally with `shiny::runApp("code/app.R")`.

## Repository structure

```
R/                   Package R source
src/                 C++ backend (Rcpp)
man/                 Documentation
code/
  app.R              Shiny dashboard (standalone, not part of the package)
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
