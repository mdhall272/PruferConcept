library(tidyverse)

get.type <- function(i, N, n, k1, k2){
  case_when(
    i <= k1 ~ "IR",
    i > k1 & i <= n ~ "IN",
    i > n & i <= n + k2 ~ "DR",
    TRUE ~ "DN"
  )
}

prufer.sequence.generator <- function(N, n, k1, k2) {

  k <- k1 + k2 # k is the number of components

  if (N < 2) {
    stop("N must be at least 2")
  }
  if (n < 0) {
    stop("n must be at least 0")
  }
  if (n > N) {
    stop("n must be at most N")
  }

  # Special case: no IN (n = k1)
  # Standard rooted forest on N nodes with k roots
  # Count: k * N^(N-k-1)
  if (n == k1) {
    if (N - k - 1 < 0) {
      # A forest of roots
      return(list(prufer = list(c()), type = "standard"))
    } else {
      if (N - k - 1 == 0) {
        seq <- c()
      } else {
        seq <- sample(1:N, N - k - 1, replace = TRUE)
      }
      # Last element must be a root
      ir_roots <- if(k1 > 0) 1:k1 else integer(0)
      dr_roots <- if(k2 > 0) (n+1):(n+k2) else integer(0)
      roots <- c(ir_roots, dr_roots)
      if (length(roots) == 1) {
        last <- roots
      } else {
        last <- sample(roots, 1)
      }
      return(list(prufer = list(c(seq, last)), type = "standard"))
    }
  }

  # Special case: no DN (k2 = N - n)
  # Each IN picks a DR to attach to
  # Count: k2^(n-k1)
  if (k2 == N - n) {
    if (n - k1 == 0) {
      seq <- c()
    } else if (k2 == 1) {
      seq <- rep(n + 1, n - k1)
    } else {
      seq <- sample((n+1):(n+k2), n - k1, replace = TRUE)
    }
    return(list(prufer = list(seq), type = "star"))
  }

  # General case below

  full.seq3.table <- tibble(s3.1 = 1:N)

  ir_roots <- if(k1 > 0) 1:k1 else integer(0)
  dr_roots <- if(k2 > 0) (n+1):(n+k2) else integer(0)
  full.seq3.table$s3.2 <- list(c(ir_roots, dr_roots))

  seq3.choices <- full.seq3.table %>%
    unnest(s3.2) %>%
    mutate(type.1 = get.type(s3.1, N, n, k1, k2)) %>%
    mutate(type.2 = get.type(s3.2, N, n, k1, k2)) %>%
    filter(!(type.1 == "IN" & type.2 == "IR"))


  # to get uniform sampling, I think we need to just choose which category we are in ahead of time

  weights <- c(k1*k2*(n-k1), (n-k1)*(k*N-k1*(n+k2)), (N-n)*(k*N-k1*(n-k1)))

  type <- sample(c("small", "medium", "large"), 1, prob = weights)

  # Empty seq2 is always medium (condition 2 holds vacuously)
  if(N - n - k2 - 1 == 0) type <- "medium"

  if(n-k1-1  == 0){
    seq1 <- c()
  } else {
    seq1 <- sample((n + 1):N, n-k1-1, replace = TRUE)
  }

if(type == "small"){
    acceptable <- FALSE
    while(!acceptable){
      seq2 <- sample(1:N, N-n-k2-1, replace = TRUE)

      types <- get.type(seq2, N, n, k1, k2)

      if(!any(types == "IN")){
        if(types[1] == "IR"){
          acceptable <- TRUE
        }
        next
      }

      final.IN <- which(types == "IN") %>% max()

      if(final.IN == length(types)){
        next
      }

      if(types[final.IN + 1] == "IR"){
        acceptable <- TRUE
        next
      }
    }

    seq3.choices <- seq3.choices %>%
      filter(type.1 == "IN" & type.2 == "DR")

    chosen.row <- sample(1:nrow(seq3.choices), 1)

    seq3 <- c(seq3.choices$s3.1[chosen.row], seq3.choices$s3.2[chosen.row])

  } else if(type == "medium"){

    if(N-n-k2-1 == 0){
      seq2 <- c()
    } else {
      if(N-n-k2-2 == 0){
        seq2 <- c()
      } else {
      seq2 <- sample(1:N, N-n-k2-2, replace = TRUE)
      }
      seq2 <- c(seq2, sample((k1+1):n, 1))
    }

    seq3.choices <- seq3.choices %>%
      filter(type.1 == "IN" & type.2 == "DR" |
               type.1 == "DN" & type.2 == "IR" |
               type.1 == "DN" & type.2 == "DR" |
               type.1 == "DR" & type.2 == "IR" |
               type.1 == "DR" & type.2 == "DR")

    chosen.row <- sample(1:nrow(seq3.choices), 1)

    seq3 <- c(seq3.choices$s3.1[chosen.row], seq3.choices$s3.2[chosen.row])


    } else if(type == "large"){

    acceptable <- FALSE

    while(!acceptable){
      if(N-n-k2-2 == 0){
        seq2 <- c()
      } else {
        seq2 <- sample(1:N, N-n-k2-2, replace = TRUE)
      }
      if(k1==0){
        seq2 <- c(seq2, sample((n+1):N, 1))
      } else {
        seq2 <- c(seq2, sample(c(1:k1, (n+1):N), 1))
      }

      types <- get.type(seq2, N, n, k1, k2)

      if(types[1] == "IR" & !any(types == "IN")){
        next
      }

      IN.class <- which(types == "IN")

      if(length(IN.class) == 0){
        acceptable <- TRUE
        next
      }

      final.IN <- IN.class %>% max()

      if(types[final.IN + 1] == "IR"){
        next
      }

      acceptable <- TRUE

    }

    chosen.row <- sample(1:nrow(seq3.choices), 1)

    seq3 <- c(seq3.choices$s3.1[chosen.row], seq3.choices$s3.2[chosen.row])

  }

  return(list(prufer = list(seq1, seq2, seq3), type = type))
}

#' Create a sampling distribution specification
#'
#' @param nb.mean Mean of negative binomial for forest size (N - n). Ignored if fix.N is set.
#' @param nb.size Size parameter of negative binomial (larger = closer to Poisson). Ignored if fix.N is set.
#' @param components.mean Mean number of components. Ignored if fix.k is set.
#' @param fix.N If not NULL, fix the collapsed forest size to this value.
#' @param fix.k If not NULL, fix the number of components to this value.
#' @return A list with class "sampling_dist" containing the parameters
sampling_distribution <- function(nb.mean = NULL, nb.size = NULL, components.mean = NULL,
                                  fix.N = NULL, fix.k = NULL) {
  if (is.null(fix.N) && (is.null(nb.mean) || is.null(nb.size))) {
    stop("nb.mean and nb.size are required when N is not fixed")
  }
  if (is.null(fix.k) && is.null(components.mean)) {
    stop("components.mean is required when k is not fixed")
  }
  structure(
    list(
      nb.mean = nb.mean,
      nb.size = nb.size,
      components.mean = components.mean,
      fix.N = fix.N,
      fix.k = fix.k
    ),
    class = "sampling_dist"
  )
}

#' Generate Prüfer sequences with roots using importance sampling
#'
#' @param Nreps Number of samples to generate
#' @param n.sampled Number of independent (sampled) nodes
#' @param sampling.dist A sampling_dist object from sampling_distribution()
#' @return A tibble with samples and attached sampling distribution
generate.prufer.sequences.with.roots <- function(
    Nreps,
    n.sampled,
    sampling.dist
) {

  if (!inherits(sampling.dist, "sampling_dist")) {
    stop("sampling.dist must be created with sampling_distribution()")
  }

  if (n.sampled < 0) {
    stop("n must be at least 0")
  }
  if (n.sampled == 0) {
    n.sampled <- 1 # these are equivalent and this is easier
  }

  if (is.null(sampling.dist$fix.N)) {
    if (sampling.dist$nb.mean <= 0) {
      stop("sampling.dist$nb.mean must be positive")
    }
    sizes <- rnbinom(Nreps, mu = sampling.dist$nb.mean, size = sampling.dist$nb.size)
    forest.sizes <- sizes + n.sampled
  } else {
    forest.sizes <- rep(as.integer(sampling.dist$fix.N), Nreps)
  }

  if (is.null(sampling.dist$fix.k)) {
    components <- map_int(forest.sizes, function(sz){
      if(sz == n.sampled){
        return(as.integer(n.sampled))
      } else {
        rbinom(1, sz-1, (sampling.dist$components.mean - 1) / (sz-1)) + 1L
      }
    })
  } else {
    components <- rep(as.integer(sampling.dist$fix.k), Nreps)
  }

  samples <- tibble(size = forest.sizes, components = components)

  # We now need to determine which of the n.sampled nodes are roots.
  # This is just a hypergeometric draw

  samples <- samples %>% mutate(k1 = map2_int(size, components, function(sz, comp){
    rhyper(1, n.sampled, sz - n.sampled, comp)
  })) %>%
    mutate(k2 = components - k1)

  samples <- samples %>% mutate(sequence = pmap(
    list(size, k1, k2),
    function(sz, k1, k2) {

      prufer.sequence.generator(sz, n.sampled, k1, k2)
    }
  )) %>% unnest_wider(sequence)

  # Attach sampling distribution as attribute
  attr(samples, "sampling_dist") <- sampling.dist
  attr(samples, "n.sampled") <- n.sampled

  return(samples)
}

classify_seq2 <- function(seq2, N, n, k1, k2) {
  if(length(seq2) == 0){
    return("medium")
  }

  pat <- get.type(seq2, N, n, k1, k2)
  L <- length(pat)

  # Medium: last element is IN
  if (pat[L] == "IN") {
    return("medium")
  }

  # Small: (no IN AND first is IR) OR (last IN immediately followed by IR)
  in_positions <- which(pat == "IN")

  if (length(in_positions) == 0) {
    if (pat[1] == "IR") {
      return("small")
    } else {
      return("large")
    }
  }

  last_in_pos <- max(in_positions)
  if (pat[last_in_pos + 1] == "IR") {
    return("small")
  }

  return("large")
}


log.forest.count <- function(N, n, k1, k2, roots.known = TRUE) {
  k <- k1 + k2

  if(!roots.known){
    root.choices <- log(choose(n, k1)) + log(choose(N - n, k2))
  } else {
    root.choices <- 0
  }

  # Special case: no IN (n = k1)
  # Count: k * N^(N-k-1)
  if (n == k1) {
    return(root.choices + log(k) + (N - k - 1) * log(N))  # Works for N=k too: log(k) - log(k) = 0
  }

  # Special case: no DN (k2 = N - n)
  # Count: k2^(n-k1)
  if (k2 == N - n) {
    return(root.choices + (n - k1) * log(k2))
  }

  # General case
  # Count: (N-n)^(n-k1-1) * N^(N-n-k2-1) * (kN - k1(n+k2))
  return(
    root.choices +
    (n - k1 - 1) * log(N - n) +
    (N - n - k2 - 1) * log(N) +
    log(k * N - k1 * (n + k2))
  )


}


log.likelihood <- function(data, sequence, components, seq.type, cluster = FALSE) {

  independent.nodes <- length(data)
  n <- independent.nodes

  if (seq.type == "standard") {
    
    # n = k1 case: single sequence of length N-k

    full.sequence <- sequence[[1]]

    collapsed.nodes <- length(full.sequence) + components


  } else if (seq.type == "star") {
    # k2 = N - n case: sequence of length n-k1

    full.sequence <- sequence[[1]]

    collapsed.nodes <- length(full.sequence) + components  # N = n + k2 = seq_len + k

  } else {
    # General case: seq1, seq2, seq3

    full.sequence <- c(sequence[[1]], sequence[[2]], sequence[[3]])

    collapsed.nodes <- length(full.sequence) + components

  }

  k <- components
  N <- collapsed.nodes

  expanded.nodes <- collapsed.nodes + sum(data - 1)

  log.denominator <- log(k) +
    (expanded.nodes - k - 1) * log(expanded.nodes) +
    log(choose(expanded.nodes, k))

  log.numerator <- 0

  for(component in 1:length(data)) {

    nodes <- data[component]

    # number of ways of opening up this collapsed node. If cluster == FALSE this cancels

    if(cluster == TRUE){
      log.numerator <- log.numerator + (nodes - 1) * log(nodes)
    }
    
    
    # number of ways of connecting the edges incident to this node

    occurrences <- length(which(full.sequence == component))

    log.numerator <- log.numerator + occurrences * log(nodes)
  }

  log.likelihood <- log.numerator - log.denominator

  return(log.likelihood)
}


#' Compute posterior weights for importance sampling
#'
#' @param data Vector of subtree sizes at independent nodes
#' @param samples Samples from generate.prufer.sequences.with.roots()
#' @param prior.N.dist Prior distribution for N-n: "poisson" or "nbinom"
#' @param prior.N.mean Prior mean for N - n (number of unsampled individuals)
#' @param prior.N.size Dispersion parameter for NB prior on N-n (ignored if Poisson)
#' @param prior.k.dist Prior distribution for k-1: "poisson" or "nbinom"
#' @param prior.k.mean Prior mean for k - 1 (extra introductions)
#' @param prior.k.size Dispersion parameter for NB prior on k-1 (ignored if Poisson)
#' @param fix.N If not NULL, fix N to this value (prior on N is omitted)
#' @param fix.k If not NULL, fix k to this value (prior on k is omitted)
#' @param cluster If TRUE, marginalize over tree structures within subtrees
#' @return The samples tibble with added log.weight column
posterior <- function(data, samples,
  prior.N.dist = "poisson",
  prior.N.mean,
  prior.N.size = 2,
  prior.k.dist = "poisson",
  prior.k.mean,
  prior.k.size = 2,
  cluster = FALSE) {

  # Extract sampling distribution from samples attribute
  sampling.dist <- attr(samples, "sampling_dist")
  if (is.null(sampling.dist)) {
    stop("samples must be generated by generate.prufer.sequences.with.roots() with sampling distribution attached")
  }

  n <- length(data)
  min.size <- sum(data)

  # N and k fixing is now handled in sampling_distribution()
  # Read fix.N and fix.k from the sampling distribution
  fix.N <- sampling.dist$fix.N
  fix.k <- sampling.dist$fix.k

  samples$data.order <- map(1:nrow(samples), function(i){
    data[sample.int(length(data))] # FIX: sample(x) with scalar x samples 1:x
  })

  is.N.fixed <- !is.null(fix.N)
  is.k.fixed <- !is.null(fix.k)

  samples <- samples %>%
    mutate(actual.size = size + sum(data - 1)) %>%
    mutate(
      likelihood = pmap_dbl(
        list(data.order, prufer, components, type),
        function(do, prf, cmp, tp) {
          log.likelihood(do, prf, cmp, tp, cluster)
        }
      )
    ) %>%
    mutate(
      prior = map2_dbl(size, components, function(sz, cmp) {
        expanded <- sz + sum(data - 1)
        prior <- 0

        # Prior on N-n (omitted if N is fixed)
        if (!is.N.fixed) {
          if (prior.N.dist == "poisson") {
            prior <- prior + dpois(expanded - min.size, prior.N.mean, log = TRUE)
          } else {
            prior <- prior + dnbinom(expanded - min.size, mu = prior.N.mean,
                                     size = prior.N.size, log = TRUE)
          }
        }

        # Prior on k-1 (omitted if k is fixed)
        if (!is.k.fixed) {
          if (prior.k.dist == "poisson") {
            prior <- prior + dpois(cmp - 1, prior.k.mean, log = TRUE)
          } else {
            prior <- prior + dnbinom(cmp - 1, mu = prior.k.mean,
                                     size = prior.k.size, log = TRUE)
          }
        }

        # Uniform prior over forest structures given N and k
        # prior <- prior - log(cmp) - (expanded - cmp - 1) * log(expanded)
        prior
      })
    ) %>%
    mutate(
      samp.prob = pmap_dbl(list(size, components, k1, k2), function(sz, cmp, k1, k2) {
        log_prob <- 0

        # Size sampling probability (omitted if N is fixed)
        if (!is.N.fixed) {
          log_prob <- log_prob + dnbinom(sz - n, mu = sampling.dist$nb.mean,
                                         size = sampling.dist$nb.size, log = TRUE)
        }

        # Components sampling probability (omitted if k is fixed)
        if (!is.k.fixed) {
          if (sz > n) {
            log_prob <- log_prob + dbinom(cmp - 1, sz - 1,
                                          (sampling.dist$components.mean - 1) / (sz - 1),
                                          log = TRUE)
          }
        }

        sampled.root.assignment.log.probability <- log(factorial(n))
        unsampled.root.assignment.log.probability <- log(choose(sz - n, k2))

        # k1 and sequence probabilities (always included)
        log_prob <- log_prob + dhyper(k1, n, sz - n, cmp, log = TRUE)
        log_prob <- log_prob - log.forest.count(sz, n, k1, k2, roots.known = TRUE)
        log_prob <- log_prob - sampled.root.assignment.log.probability
        log_prob <- log_prob - unsampled.root.assignment.log.probability

        log_prob
      })
    ) %>%
    mutate(log.weight = likelihood + prior - samp.prob)

  # Preserve the sampling_dist attribute
  attr(samples, "sampling_dist") <- sampling.dist

  return(samples)
}

#' Format a Prüfer sequence list as a string
#'
#' @param prufer A list of sub-sequences (from prufer.sequence.generator)
#' @param type The sequence type: "standard", "star", or "small"/"medium"/"large"
#' @return A string representation: sub-sequences comma-separated, joined by ";"
format.prufer <- function(prufer, type) {
  parts <- map_chr(prufer, function(s) {
    if (is.null(s) || length(s) == 0) return("")
    paste(s, collapse = ",")
  })
  paste(parts, collapse = ";")
}

#' Identify which sampled nodes are roots given the data order
#'
#' The data order assigns cluster sizes to node positions 1..n.
#' Nodes 1..k1 are independent roots (IR). This function returns
#' the original data indices that were assigned to root positions.
#'
#' @param data.order The permuted data vector
#' @param data The original data vector
#' @param k1 Number of independent roots
#' @return An integer vector of 1-based indices into the original data
#'         that are assigned to root positions
root.indices <- function(data.order, data, k1) {
  if (k1 == 0) return(integer(0))
  # data.order is data[perm], so data.order[i] = data[perm[i]]
  # Root positions are 1..k1. We want to know which original data
  # entries ended up there. Match by position in the permutation.
  # Since data.order = data[perm], we need perm. But we don't store perm.
  # Instead, report the cluster sizes at root positions.
  # Actually: return the root cluster sizes directly, since that's
  # what matters for the likelihood.
  data.order[1:k1]
}

#' Build a sample log tibble from posterior output
#'
#' @param result Output from posterior()
#' @param data The original data vector
#' @return A tibble with one row per sample, including Prüfer strings,
#'         root cluster sizes, and all weight components
sample.log <- function(result, data) {
  result %>%
    mutate(
      prufer.string = map2_chr(prufer, type, format.prufer),
      root.cluster.sizes = map2_chr(data.order, k1, function(do, k1) {
        if (k1 == 0) return("")
        paste(do[1:k1], collapse = ",")
      }),
      nonroot.cluster.sizes = map2_chr(data.order, k1, function(do, k1) {
        n <- length(do)
        if (k1 >= n) return("")
        paste(do[(k1+1):n], collapse = ",")
      })
    ) %>%
    select(size, actual.size, components, k1, k2, type,
           prufer.string, root.cluster.sizes, nonroot.cluster.sizes,
           likelihood, prior, samp.prob, log.weight)
}
