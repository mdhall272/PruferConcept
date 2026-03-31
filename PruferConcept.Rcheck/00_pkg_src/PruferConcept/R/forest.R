#' @importFrom dplyr mutate filter group_by summarise select case_when %>%
#' @importFrom purrr map map_int map2_int map_chr map2_chr pmap pmap_dbl map2_dbl
#' @importFrom tibble tibble
#' @importFrom tidyr unnest unnest_wider
#' @importFrom stats rnbinom rbinom rhyper dpois dnbinom dbinom dhyper
#' @importFrom Rcpp sourceCpp
#' @useDynLib PruferConcept, .registration = TRUE
NULL

# Suppress R CMD check notes for NSE variables used in dplyr/tidyr
utils::globalVariables(c(
  "s3.1", "s3.2", "type.1", "type.2",
  "size", "components", "k1", "k2", "prufer", "type",
  "data.order", "likelihood", "prior", "samp.prob"
))

get_type <- function(i, N, n, k1, k2) {
  case_when(
    i <= k1 ~ "IR",
    i > k1 & i <= n ~ "IN",
    i > n & i <= n + k2 ~ "DR",
    TRUE ~ "DN"
  )
}

prufer_sequence_generator <- function(N, n, k1, k2) {

  k <- k1 + k2

  if (N < 2) stop("N must be at least 2")
  if (n < 0) stop("n must be at least 0")
  if (n > N) stop("n must be at most N")

  # Special case: no IN (n = k1)
  if (n == k1) {
    if (N - k - 1 < 0) {
      return(list(prufer = list(c()), type = "standard"))
    } else {
      if (N - k - 1 == 0) {
        seq <- c()
      } else {
        seq <- sample(1:N, N - k - 1, replace = TRUE)
      }
      ir_roots <- if (k1 > 0) 1:k1 else integer(0)
      dr_roots <- if (k2 > 0) (n + 1):(n + k2) else integer(0)
      roots <- c(ir_roots, dr_roots)
      last <- if (length(roots) == 1) roots else sample(roots, 1)
      return(list(prufer = list(c(seq, last)), type = "standard"))
    }
  }

  # Special case: no DN (k2 = N - n)
  if (k2 == N - n) {
    if (n - k1 == 0) {
      seq <- c()
    } else if (k2 == 1) {
      seq <- rep(n + 1, n - k1)
    } else {
      seq <- sample((n + 1):(n + k2), n - k1, replace = TRUE)
    }
    return(list(prufer = list(seq), type = "star"))
  }

  # General case
  full_seq3_table <- tibble(s3.1 = 1:N)
  ir_roots <- if (k1 > 0) 1:k1 else integer(0)
  dr_roots <- if (k2 > 0) (n + 1):(n + k2) else integer(0)
  full_seq3_table$s3.2 <- list(c(ir_roots, dr_roots))

  seq3_choices <- full_seq3_table %>%
    unnest(s3.2) %>%
    mutate(type.1 = get_type(s3.1, N, n, k1, k2)) %>%
    mutate(type.2 = get_type(s3.2, N, n, k1, k2)) %>%
    filter(!(type.1 == "IN" & type.2 == "IR"))

  weights <- c(k1 * k2 * (n - k1), (n - k1) * (k * N - k1 * (n + k2)), (N - n) * (k * N - k1 * (n - k1)))
  type <- sample(c("small", "medium", "large"), 1, prob = weights)

  if (N - n - k2 - 1 == 0) type <- "medium"

  if (n - k1 - 1 == 0) {
    seq1 <- c()
  } else {
    seq1 <- sample((n + 1):N, n - k1 - 1, replace = TRUE)
  }

  if (type == "small") {
    acceptable <- FALSE
    while (!acceptable) {
      seq2 <- sample(1:N, N - n - k2 - 1, replace = TRUE)
      types <- get_type(seq2, N, n, k1, k2)
      if (!any(types == "IN")) {
        if (types[1] == "IR") acceptable <- TRUE
        next
      }
      final_IN <- max(which(types == "IN"))
      if (final_IN == length(types)) next
      if (types[final_IN + 1] == "IR") {
        acceptable <- TRUE
        next
      }
    }
    seq3_choices <- seq3_choices %>% filter(type.1 == "IN" & type.2 == "DR")
    chosen_row <- sample(1:nrow(seq3_choices), 1)
    seq3 <- c(seq3_choices$s3.1[chosen_row], seq3_choices$s3.2[chosen_row])

  } else if (type == "medium") {
    if (N - n - k2 - 1 == 0) {
      seq2 <- c()
    } else {
      if (N - n - k2 - 2 == 0) {
        seq2 <- c()
      } else {
        seq2 <- sample(1:N, N - n - k2 - 2, replace = TRUE)
      }
      seq2 <- c(seq2, sample((k1 + 1):n, 1))
    }
    seq3_choices <- seq3_choices %>%
      filter(type.1 == "IN" & type.2 == "DR" |
               type.1 == "DN" & type.2 == "IR" |
               type.1 == "DN" & type.2 == "DR" |
               type.1 == "DR" & type.2 == "IR" |
               type.1 == "DR" & type.2 == "DR")
    chosen_row <- sample(1:nrow(seq3_choices), 1)
    seq3 <- c(seq3_choices$s3.1[chosen_row], seq3_choices$s3.2[chosen_row])

  } else if (type == "large") {
    acceptable <- FALSE
    while (!acceptable) {
      if (N - n - k2 - 2 == 0) {
        seq2 <- c()
      } else {
        seq2 <- sample(1:N, N - n - k2 - 2, replace = TRUE)
      }
      if (k1 == 0) {
        seq2 <- c(seq2, sample((n + 1):N, 1))
      } else {
        seq2 <- c(seq2, sample(c(1:k1, (n + 1):N), 1))
      }
      types <- get_type(seq2, N, n, k1, k2)
      if (types[1] == "IR" & !any(types == "IN")) next
      IN_class <- which(types == "IN")
      if (length(IN_class) == 0) {
        acceptable <- TRUE
        next
      }
      final_IN <- max(IN_class)
      if (types[final_IN + 1] == "IR") next
      acceptable <- TRUE
    }
    chosen_row <- sample(1:nrow(seq3_choices), 1)
    seq3 <- c(seq3_choices$s3.1[chosen_row], seq3_choices$s3.2[chosen_row])
  }

  list(prufer = list(seq1, seq2, seq3), type = type)
}


#' Create a sampling distribution specification
#'
#' @param nb.mean Mean of negative binomial for forest size (N - n).
#'   Ignored if fix.N is set.
#' @param nb.size Size (dispersion) parameter of negative binomial.
#'   Ignored if fix.N is set.
#' @param components.mean Mean number of components.
#'   Ignored if fix.k is set.
#' @param fix.N If not NULL, fix the collapsed forest size to this value.
#' @param fix.k If not NULL, fix the number of components to this value.
#' @return A list with class \code{sampling_dist}
#' @export
sampling_distribution <- function(nb.mean = NULL, nb.size = NULL,
                                  components.mean = NULL,
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


#' Generate importance samples from the forest proposal distribution
#'
#' @param Nreps Number of samples to generate
#' @param n.sampled Number of independent (sampled) nodes
#' @param sampling.dist A \code{sampling_dist} object from
#'   \code{\link{sampling_distribution}}
#' @return A tibble with columns \code{size}, \code{components}, \code{k1},
#'   \code{k2}, \code{prufer}, \code{type}, with the sampling distribution
#'   attached as an attribute
#' @export
generate_samples <- function(Nreps, n.sampled, sampling.dist) {

  if (!inherits(sampling.dist, "sampling_dist")) {
    stop("sampling.dist must be created with sampling_distribution()")
  }
  if (n.sampled < 0) stop("n must be at least 0")
  if (n.sampled == 0) n.sampled <- 1

  if (is.null(sampling.dist$fix.N)) {
    if (sampling.dist$nb.mean <= 0) stop("nb.mean must be positive")
    sizes <- rnbinom(Nreps, mu = sampling.dist$nb.mean, size = sampling.dist$nb.size)
    forest_sizes <- sizes + n.sampled
  } else {
    forest_sizes <- rep(as.integer(sampling.dist$fix.N), Nreps)
  }

  if (is.null(sampling.dist$fix.k)) {
    components <- map_int(forest_sizes, function(sz) {
      if (sz == n.sampled) {
        return(as.integer(n.sampled))
      } else {
        rbinom(1, sz - 1, (sampling.dist$components.mean - 1) / (sz - 1)) + 1L
      }
    })
  } else {
    components <- rep(as.integer(sampling.dist$fix.k), Nreps)
  }

  samples <- tibble(size = forest_sizes, components = components)

  samples <- samples %>%
    mutate(k1 = map2_int(size, components, function(sz, comp) {
      rhyper(1, n.sampled, sz - n.sampled, comp)
    })) %>%
    mutate(k2 = components - k1)

  samples <- samples %>%
    mutate(sequence = pmap(
      list(size, k1, k2),
      function(sz, k1, k2) {
        prufer_sequence_generator(sz, n.sampled, k1, k2)
      }
    )) %>%
    unnest_wider(sequence)

  attr(samples, "sampling_dist") <- sampling.dist
  attr(samples, "n.sampled") <- n.sampled
  samples
}


#' Log forest count
#'
#' @param N Number of nodes in the collapsed forest
#' @param n Number of independent nodes
#' @param k1 Number of independent roots
#' @param k2 Number of dependent roots
#' @param roots.known If TRUE, root identities are known (no combinatorial factor)
#' @return Log of the number of valid forests
#' @export
log_forest_count <- function(N, n, k1, k2, roots.known = TRUE) {
  k <- k1 + k2

  if (!roots.known) {
    root_choices <- log(choose(n, k1)) + log(choose(N - n, k2))
  } else {
    root_choices <- 0
  }

  if (n == k1) {
    return(root_choices + log(k) + (N - k - 1) * log(N))
  }
  if (k2 == N - n) {
    return(root_choices + (n - k1) * log(k2))
  }

  root_choices +
    (n - k1 - 1) * log(N - n) +
    (N - n - k2 - 1) * log(N) +
    log(k * N - k1 * (n + k2))
}


#' Log-likelihood of cluster data given a forest structure
#'
#' @param data Integer vector of cluster sizes
#' @param sequence List of Prufer sub-sequences
#' @param components Number of components (k)
#' @param seq.type Sequence type: "standard", "star", "small", "medium", or "large"
#' @param cluster If TRUE, marginalize over subtree structures within clusters
#' @return Log-likelihood value
#' @export
log_likelihood <- function(data, sequence, components, seq.type, cluster = FALSE) {

  n <- length(data)

  if (seq.type == "standard") {
    full_sequence <- sequence[[1]]
    collapsed_nodes <- length(full_sequence) + components
  } else if (seq.type == "star") {
    full_sequence <- sequence[[1]]
    collapsed_nodes <- length(full_sequence) + components
  } else {
    full_sequence <- c(sequence[[1]], sequence[[2]], sequence[[3]])
    collapsed_nodes <- length(full_sequence) + components
  }

  k <- components
  N <- collapsed_nodes
  expanded_nodes <- collapsed_nodes + sum(data - 1)

  log_denom <- log(k) +
    (expanded_nodes - k - 1) * log(expanded_nodes) +
    log(choose(expanded_nodes, k))

  log_num <- 0
  for (i in seq_along(data)) {
    nodes <- data[i]
    if (cluster) {
      log_num <- log_num + (nodes - 1) * log(nodes)
    }
    occurrences <- length(which(full_sequence == i))
    log_num <- log_num + occurrences * log(nodes)
  }

  log_num - log_denom
}


#' Compute posterior importance weights
#'
#' @param data Integer vector of cluster sizes
#' @param samples Output from \code{\link{generate_samples}}
#' @param prior.N.dist Prior on N-n: "poisson" or "nbinom"
#' @param prior.N.mean Prior mean for N - n
#' @param prior.N.size Dispersion for NB prior on N-n (ignored if Poisson)
#' @param prior.k.dist Prior on k-1: "poisson" or "nbinom"
#' @param prior.k.mean Prior mean for k - 1
#' @param prior.k.size Dispersion for NB prior on k-1 (ignored if Poisson)
#' @param cluster If TRUE, marginalize over subtree structures
#' @return The samples tibble with added columns: \code{actual.size},
#'   \code{data.order}, \code{likelihood}, \code{prior}, \code{samp.prob},
#'   \code{log.weight}
#' @export
posterior <- function(data, samples,
                      prior.N.dist = "poisson",
                      prior.N.mean,
                      prior.N.size = 2,
                      prior.k.dist = "poisson",
                      prior.k.mean,
                      prior.k.size = 2,
                      cluster = FALSE) {

  sampling.dist <- attr(samples, "sampling_dist")
  if (is.null(sampling.dist)) {
    stop("samples must be generated by generate_samples()")
  }

  n <- length(data)
  min_size <- sum(data)

  fix.N <- sampling.dist$fix.N
  fix.k <- sampling.dist$fix.k

  samples$data.order <- map(1:nrow(samples), function(i) {
    data[sample.int(length(data))]
  })

  is_N_fixed <- !is.null(fix.N)
  is_k_fixed <- !is.null(fix.k)

  samples <- samples %>%
    mutate(actual.size = size + sum(data - 1)) %>%
    mutate(
      likelihood = pmap_dbl(
        list(data.order, prufer, components, type),
        function(do, prf, cmp, tp) {
          log_likelihood(do, prf, cmp, tp, cluster)
        }
      )
    ) %>%
    mutate(
      prior = map2_dbl(size, components, function(sz, cmp) {
        expanded <- sz + sum(data - 1)
        prior <- 0
        if (!is_N_fixed) {
          if (prior.N.dist == "poisson") {
            prior <- prior + dpois(expanded - min_size, prior.N.mean, log = TRUE)
          } else {
            prior <- prior + dnbinom(expanded - min_size, mu = prior.N.mean,
                                     size = prior.N.size, log = TRUE)
          }
        }
        if (!is_k_fixed) {
          if (prior.k.dist == "poisson") {
            prior <- prior + dpois(cmp - 1, prior.k.mean, log = TRUE)
          } else {
            prior <- prior + dnbinom(cmp - 1, mu = prior.k.mean,
                                     size = prior.k.size, log = TRUE)
          }
        }
        prior
      })
    ) %>%
    mutate(
      samp.prob = pmap_dbl(list(size, components, k1, k2), function(sz, cmp, k1, k2) {
        log_prob <- 0
        if (!is_N_fixed) {
          log_prob <- log_prob + dnbinom(sz - n, mu = sampling.dist$nb.mean,
                                         size = sampling.dist$nb.size, log = TRUE)
        }
        if (!is_k_fixed) {
          if (sz > n) {
            log_prob <- log_prob + dbinom(cmp - 1, sz - 1,
                                          (sampling.dist$components.mean - 1) / (sz - 1),
                                          log = TRUE)
          }
        }
        log_prob <- log_prob + dhyper(k1, n, sz - n, cmp, log = TRUE)
        log_prob <- log_prob - log_forest_count(sz, n, k1, k2, roots.known = TRUE)
        log_prob <- log_prob - log(factorial(n))
        log_prob <- log_prob - log(choose(sz - n, k2))
        log_prob
      })
    ) %>%
    mutate(log.weight = likelihood + prior - samp.prob)

  attr(samples, "sampling_dist") <- sampling.dist
  samples
}


#' Format a Prufer sequence list as a string
#'
#' @param prufer A list of sub-sequences
#' @param type The sequence type
#' @return A string with sub-sequences comma-separated, joined by ";"
#' @export
format_prufer <- function(prufer, type) {
  parts <- map_chr(prufer, function(s) {
    if (is.null(s) || length(s) == 0) return("")
    paste(s, collapse = ",")
  })
  paste(parts, collapse = ";")
}
