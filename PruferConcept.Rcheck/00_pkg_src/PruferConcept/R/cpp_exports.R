#' Generate samples in batch using C++ backend
#'
#' @param Nreps Number of samples
#' @param n_sampled Number of sampled nodes
#' @param nb_mean Mean of NegBinom proposal for N-n
#' @param nb_size Size of NegBinom proposal
#' @param components_mean Mean number of components in proposal
#' @param fix_N If > 0, fix collapsed forest size to this value
#' @param fix_k If > 0, fix number of components to this value
#' @return A list with elements size, components, k1, k2, prufer, type
#' @export
generate_samples_batch <- function(Nreps, n_sampled, nb_mean = 0,
                                   nb_size = 0, components_mean = 0,
                                   fix_N = 0L, fix_k = 0L) {
  generate_samples_batch_cpp(
    as.integer(Nreps), as.integer(n_sampled),
    as.double(nb_mean), as.double(nb_size), as.double(components_mean),
    as.integer(fix_N), as.integer(fix_k)
  )
}

#' Compute posterior weights in batch using C++ backend
#'
#' @param data Integer vector of cluster sizes
#' @param samples_list Output from \code{\link{generate_samples_batch}}
#' @param params List with elements: prior_N_dist, prior_N_mean, prior_N_size,
#'   prior_k_dist, prior_k_mean, prior_k_size, nb_mean, nb_size,
#'   components_mean, fix_N, fix_k, cluster
#' @return A list with elements: actual.size, likelihood, prior, samp.prob,
#'   log.weight, data.order, prufer.string, root.cluster.sizes,
#'   nonroot.cluster.sizes
#' @export
compute_posterior_batch <- function(data, samples_list, params) {
  compute_posterior_batch_cpp(as.integer(data), samples_list, params)
}

#' Compute log forest count using C++ backend
#'
#' @param N Number of nodes
#' @param n Number of independent nodes
#' @param k1 Number of independent roots
#' @param k2 Number of dependent roots
#' @param roots_known Whether root identities are known
#' @return Log forest count
#' @export
log_forest_count_cpp_wrapper <- function(N, n, k1, k2, roots_known = TRUE) {
  log_forest_count_cpp(
    as.integer(N), as.integer(n), as.integer(k1), as.integer(k2), roots_known
  )
}

#' Compute log-likelihood using C++ backend
#'
#' @param data Integer vector of cluster sizes
#' @param sequence List of Prufer sub-sequences
#' @param components Number of components
#' @param seq_type Sequence type string
#' @param cluster Whether to marginalize over subtree structures
#' @return Log-likelihood
#' @export
log_likelihood_cpp_wrapper <- function(data, sequence, components, seq_type,
                                       cluster = FALSE) {
  log_likelihood_cpp(
    as.integer(data), sequence, as.integer(components), seq_type, cluster
  )
}
