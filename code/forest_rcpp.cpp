#include <Rcpp.h>
using namespace Rcpp;

// Node type: 0=IR, 1=IN, 2=DR, 3=DN
inline int get_type(int i, int N, int n, int k1, int k2) {
  if (i <= k1) return 0;        // IR
  if (i <= n) return 1;         // IN
  if (i <= n + k2) return 2;    // DR
  return 3;                      // DN
}

// Sample one integer from [a, b] inclusive
inline int sample_one(int a, int b) {
  if (a == b) return a;
  return a + (int)(R::runif(0, 1) * (b - a + 1));
}

// Sample n integers from [a, b] with replacement
IntegerVector sample_replace(int a, int b, int n) {
  IntegerVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = sample_one(a, b);
  }
  return result;
}

// Sample one element from a vector
inline int sample_from_vec(IntegerVector v) {
  if (v.size() == 1) return v[0];
  int idx = (int)(R::runif(0, 1) * v.size());
  if (idx >= v.size()) idx = v.size() - 1;
  return v[idx];
}

// Fisher-Yates shuffle of an IntegerVector (in place)
inline void shuffle_inplace(IntegerVector &v) {
  int n = v.size();
  for (int i = n - 1; i > 0; i--) {
    int j = (int)(R::runif(0, 1) * (i + 1));
    if (j > i) j = i;
    int tmp = v[i];
    v[i] = v[j];
    v[j] = tmp;
  }
}

// [[Rcpp::export]]
List prufer_sequence_generator_cpp(int N, int n, int k1, int k2) {
  int k = k1 + k2;

  if (N < 2) stop("N must be at least 2");
  if (n < 0) stop("n must be at least 0");
  if (n > N) stop("n must be at most N");

  // Special case: no IN (n = k1)
  if (n == k1) {
    if (N - k - 1 < 0) {
      // Forest of roots
      return List::create(
        Named("prufer") = List::create(IntegerVector(0)),
        Named("type") = "standard"
      );
    }

    IntegerVector seq;
    if (N - k - 1 == 0) {
      seq = IntegerVector(0);
    } else {
      seq = sample_replace(1, N, N - k - 1);
    }

    // Build roots vector
    IntegerVector roots(k);
    int idx = 0;
    for (int i = 1; i <= k1; i++) roots[idx++] = i;
    for (int i = n + 1; i <= n + k2; i++) roots[idx++] = i;

    int last = sample_from_vec(roots);

    // Combine seq and last
    IntegerVector full_seq(seq.size() + 1);
    for (int i = 0; i < seq.size(); i++) full_seq[i] = seq[i];
    full_seq[seq.size()] = last;

    return List::create(
      Named("prufer") = List::create(full_seq),
      Named("type") = "standard"
    );
  }

  // Special case: no DN (k2 = N - n)
  if (k2 == N - n) {
    IntegerVector seq;
    if (n - k1 == 0) {
      seq = IntegerVector(0);
    } else if (k2 == 1) {
      seq = IntegerVector(n - k1, n + 1);
    } else {
      seq = sample_replace(n + 1, n + k2, n - k1);
    }

    return List::create(
      Named("prufer") = List::create(seq),
      Named("type") = "star"
    );
  }

  // General case

  // Build seq3 choices: pairs (s3_1, s3_2) where s3_2 is a root and not (IN, IR)
  std::vector<std::pair<int, int>> seq3_small;   // (IN, DR) pairs
  std::vector<std::pair<int, int>> seq3_medium;  // medium category pairs
  std::vector<std::pair<int, int>> seq3_large;   // large category pairs (all except IN,IR)

  for (int s1 = 1; s1 <= N; s1++) {
    int t1 = get_type(s1, N, n, k1, k2);

    // IR roots
    for (int s2 = 1; s2 <= k1; s2++) {
      int t2 = 0; // IR
      if (!(t1 == 1 && t2 == 0)) { // not (IN, IR)
        seq3_large.push_back({s1, s2});
        // Medium: (DN,IR), (DR,IR)
        if ((t1 == 3 && t2 == 0) || (t1 == 2 && t2 == 0)) {
          seq3_medium.push_back({s1, s2});
        }
      }
    }

    // DR roots
    for (int s2 = n + 1; s2 <= n + k2; s2++) {
      int t2 = 2; // DR
      // (IN, DR) is always allowed
      seq3_large.push_back({s1, s2});

      if (t1 == 1 && t2 == 2) { // (IN, DR)
        seq3_small.push_back({s1, s2});
        seq3_medium.push_back({s1, s2});
      } else if ((t1 == 3 && t2 == 2) || (t1 == 2 && t2 == 2)) {
        // (DN, DR) or (DR, DR)
        seq3_medium.push_back({s1, s2});
      }
    }
  }

  // Compute weights for category selection
  double w_small = (double)k1 * k2 * (n - k1);
  double w_medium = (double)(n - k1) * (k * N - k1 * (n + k2));
  double w_large = (double)(N - n) * (k * N - k1 * (n - k1));
  double total_w = w_small + w_medium + w_large;

  // Choose category
  double r = R::runif(0, 1) * total_w;
  std::string type_str;
  if (r < w_small) {
    type_str = "small";
  } else if (r < w_small + w_medium) {
    type_str = "medium";
  } else {
    type_str = "large";
  }

  // Generate seq1
  IntegerVector seq1;
  if (n - k1 - 1 == 0) {
    seq1 = IntegerVector(0);
  } else {
    seq1 = sample_replace(n + 1, N, n - k1 - 1);
  }

  // Generate seq2 and seq3
  IntegerVector seq2;
  IntegerVector seq3(2);

  if (N - n - k2 - 1 == 0) {
    seq2 = IntegerVector(0);
    // Empty seq2 treated as medium category (condition 2 holds vacuously)
    int idx = (int)(R::runif(0, 1) * seq3_medium.size());
    if (idx >= (int)seq3_medium.size()) idx = seq3_medium.size() - 1;
    seq3[0] = seq3_medium[idx].first;
    seq3[1] = seq3_medium[idx].second;

  } else if (type_str == "small") {
    // Rejection sampling for small category
    bool acceptable = false;
    while (!acceptable) {
      seq2 = sample_replace(1, N, N - n - k2 - 1);

      // Check types
      bool has_IN = false;
      int final_IN = -1;
      for (int i = 0; i < seq2.size(); i++) {
        if (get_type(seq2[i], N, n, k1, k2) == 1) {
          has_IN = true;
          final_IN = i;
        }
      }

      if (!has_IN) {
        if (get_type(seq2[0], N, n, k1, k2) == 0) { // first is IR
          acceptable = true;
        }
        continue;
      }

      if (final_IN == seq2.size() - 1) continue;

      if (get_type(seq2[final_IN + 1], N, n, k1, k2) == 0) { // followed by IR
        acceptable = true;
      }
    }

    int idx = (int)(R::runif(0, 1) * seq3_small.size());
    if (idx >= (int)seq3_small.size()) idx = seq3_small.size() - 1;
    seq3[0] = seq3_small[idx].first;
    seq3[1] = seq3_small[idx].second;

  } else if (type_str == "medium") {
    // Last element must be IN
    if (N - n - k2 - 2 == 0) {
      seq2 = IntegerVector(0);
    } else {
      seq2 = sample_replace(1, N, N - n - k2 - 2);
    }
    // Append random IN
    int last_in = sample_one(k1 + 1, n);
    IntegerVector new_seq2(seq2.size() + 1);
    for (int i = 0; i < seq2.size(); i++) new_seq2[i] = seq2[i];
    new_seq2[seq2.size()] = last_in;
    seq2 = new_seq2;

    int idx = (int)(R::runif(0, 1) * seq3_medium.size());
    if (idx >= (int)seq3_medium.size()) idx = seq3_medium.size() - 1;
    seq3[0] = seq3_medium[idx].first;
    seq3[1] = seq3_medium[idx].second;

  } else { // large
    bool acceptable = false;
    while (!acceptable) {
      if (N - n - k2 - 2 == 0) {
        seq2 = IntegerVector(0);
      } else {
        seq2 = sample_replace(1, N, N - n - k2 - 2);
      }

      // Last element from IR or dependent nodes
      int last_elem;
      if (k1 == 0) {
        last_elem = sample_one(n + 1, N);
      } else {
        // Sample from 1:k1 or (n+1):N
        int total = k1 + (N - n);
        int r = sample_one(1, total);
        if (r <= k1) {
          last_elem = r;
        } else {
          last_elem = n + (r - k1);
        }
      }

      IntegerVector new_seq2(seq2.size() + 1);
      for (int i = 0; i < seq2.size(); i++) new_seq2[i] = seq2[i];
      new_seq2[seq2.size()] = last_elem;
      seq2 = new_seq2;

      // Check validity
      bool has_IN = false;
      int final_IN = -1;
      for (int i = 0; i < seq2.size(); i++) {
        if (get_type(seq2[i], N, n, k1, k2) == 1) {
          has_IN = true;
          final_IN = i;
        }
      }

      // Reject if first is IR and no IN
      if (get_type(seq2[0], N, n, k1, k2) == 0 && !has_IN) continue;

      if (!has_IN) {
        acceptable = true;
        continue;
      }

      // Reject if last IN followed by IR
      if (final_IN < seq2.size() - 1 && get_type(seq2[final_IN + 1], N, n, k1, k2) == 0) {
        continue;
      }

      acceptable = true;
    }

    int idx = (int)(R::runif(0, 1) * seq3_large.size());
    if (idx >= (int)seq3_large.size()) idx = seq3_large.size() - 1;
    seq3[0] = seq3_large[idx].first;
    seq3[1] = seq3_large[idx].second;
  }

  return List::create(
    Named("prufer") = List::create(seq1, seq2, seq3),
    Named("type") = type_str
  );
}

// [[Rcpp::export]]
double log_forest_count_cpp(int N, int n, int k1, int k2, bool roots_known = true) {
  int k = k1 + k2;

  double root_choices = 0.0;
  if (!roots_known) {
    root_choices = R::lchoose(n, k1) + R::lchoose(N - n, k2);
  }

  // Special case: n = k1
  if (n == k1) {
    return root_choices + log((double)k) + (N - k - 1) * log((double)N);
  }

  // Special case: k2 = N - n
  if (k2 == N - n) {
    return root_choices + (n - k1) * log((double)k2);
  }

  // General case
  return root_choices +
    (n - k1 - 1) * log((double)(N - n)) +
    (N - n - k2 - 1) * log((double)N) +
    log((double)(k * N - k1 * (n + k2)));
}

// [[Rcpp::export]]
double log_likelihood_cpp(IntegerVector data, List sequence, int components,
                          std::string seq_type, bool cluster = false) {
  int n = data.size();

  IntegerVector full_sequence;
  int collapsed_nodes;

  if (seq_type == "standard") {
    IntegerVector seq0 = sequence[0];
    full_sequence = seq0;
    collapsed_nodes = seq0.size() + components;

  } else if (seq_type == "star") {
    IntegerVector seq0 = sequence[0];
    full_sequence = seq0;
    collapsed_nodes = seq0.size() + components;

  } else {
    IntegerVector seq0 = sequence[0];
    IntegerVector seq1 = sequence[1];
    IntegerVector seq2 = sequence[2];

    int total_len = seq0.size() + seq1.size() + seq2.size();
    full_sequence = IntegerVector(total_len);
    int idx = 0;
    for (int i = 0; i < seq0.size(); i++) full_sequence[idx++] = seq0[i];
    for (int i = 0; i < seq1.size(); i++) full_sequence[idx++] = seq1[i];
    for (int i = 0; i < seq2.size(); i++) full_sequence[idx++] = seq2[i];

    collapsed_nodes = total_len + components;
  }

  int k = components;
  int N = collapsed_nodes;

  int data_sum = 0;
  for (int i = 0; i < n; i++) data_sum += data[i];
  int expanded_nodes = collapsed_nodes + data_sum - n;

  double log_denom = log((double)k) +
    (expanded_nodes - k - 1) * log((double)expanded_nodes) +
    R::lchoose(expanded_nodes, k);

  double log_num = 0.0;

  for (int comp = 0; comp < n; comp++) {
    int nodes = data[comp];

    if (cluster) {
      log_num += (nodes - 1) * log((double)nodes);
    }

    // Count occurrences of (comp+1) in full_sequence (1-indexed)
    int occurrences = 0;
    for (int i = 0; i < full_sequence.size(); i++) {
      if (full_sequence[i] == comp + 1) occurrences++;
    }

    log_num += occurrences * log((double)nodes);
  }

  return log_num - log_denom;
}

// [[Rcpp::export]]
double log_sampling_probability_cpp(int N, int n, int k1, int k2,
                                    double nb_mean, double nb_size,
                                    double components_mean) {
  int drawn_size = N - n;
  int k = k1 + k2;

  double size_log_prob = R::dnbinom_mu(drawn_size, nb_size, nb_mean, true);

  double comp_log_prob;
  if (N == n) {
    comp_log_prob = 0.0;
  } else {
    comp_log_prob = R::dbinom(k - 1, N - 1, (components_mean - 1) / (N - 1), true);
  }

  double k1_log_prob = R::dhyper(k1, n, N - n, k, true);

  double seq_log_prob = -log_forest_count_cpp(N, n, k1, k2, true);

  double root_assign_log_prob = R::lgammafn(n + 1);  // log(n!)

  return size_log_prob + comp_log_prob + k1_log_prob + seq_log_prob - root_assign_log_prob;
}

// Classify a seq2 sequence into small/medium/large category
// Returns 0=small, 1=medium, 2=large
// [[Rcpp::export]]
int classify_seq2_cpp(IntegerVector seq2, int N, int n, int k1, int k2) {
  if (seq2.size() == 0) {
    return 1; // medium
  }

  int L = seq2.size();

  // Medium: last element is IN
  if (get_type(seq2[L - 1], N, n, k1, k2) == 1) {
    return 1; // medium
  }

  // Find IN positions
  int last_in_pos = -1;
  bool has_IN = false;
  for (int i = 0; i < L; i++) {
    if (get_type(seq2[i], N, n, k1, k2) == 1) {
      has_IN = true;
      last_in_pos = i;
    }
  }

  if (!has_IN) {
    if (get_type(seq2[0], N, n, k1, k2) == 0) { // first is IR
      return 0; // small
    } else {
      return 2; // large
    }
  }

  // last IN followed by IR?
  if (get_type(seq2[last_in_pos + 1], N, n, k1, k2) == 0) {
    return 0; // small
  }

  return 2; // large
}

// Format a Prufer sequence list as a string
// Sub-sequences are comma-separated, joined by ";"
std::string format_prufer(List prufer, std::string type) {
  std::string result;
  for (int i = 0; i < prufer.size(); i++) {
    if (i > 0) result += ";";
    if (prufer[i] == R_NilValue) continue;
    IntegerVector seq = as<IntegerVector>(prufer[i]);
    for (int j = 0; j < seq.size(); j++) {
      if (j > 0) result += ",";
      result += std::to_string(seq[j]);
    }
  }
  return result;
}

// Format cluster sizes for root or non-root positions
std::string format_cluster_sizes(IntegerVector perm_data, int k1, bool roots) {
  int n = perm_data.size();
  std::string result;
  if (roots) {
    for (int i = 0; i < k1; i++) {
      if (i > 0) result += ",";
      result += std::to_string(perm_data[i]);
    }
  } else {
    for (int i = k1; i < n; i++) {
      if (i > k1) result += ",";
      result += std::to_string(perm_data[i]);
    }
  }
  return result;
}

// Generate multiple samples at once (batch version)
// fix_N: if > 0, fix collapsed forest size to this value
// fix_k: if > 0, fix number of components to this value
// nb_mean, nb_size, components_mean: ignored when corresponding parameter is fixed
// [[Rcpp::export]]
List generate_samples_batch_cpp(int Nreps, int n_sampled,
                                 double nb_mean = 0, double nb_size = 0,
                                 double components_mean = 0,
                                 int fix_N = 0, int fix_k = 0) {
  IntegerVector sizes(Nreps);
  IntegerVector components(Nreps);
  IntegerVector k1_vec(Nreps);
  IntegerVector k2_vec(Nreps);
  List prufer_list(Nreps);
  CharacterVector type_vec(Nreps);

  // Convert mu/size parametrization to prob parametrization for rnbinom
  double prob = 0;
  if (fix_N <= 0 && nb_size > 0) {
    prob = nb_size / (nb_size + nb_mean);
  }

  for (int i = 0; i < Nreps; i++) {
    int sz;
    if (fix_N > 0) {
      sz = fix_N;
    } else {
      sz = (int)R::rnbinom(nb_size, prob) + n_sampled;
    }
    sizes[i] = sz;

    int comp;
    if (fix_k > 0) {
      comp = fix_k;
    } else if (sz == n_sampled) {
      comp = n_sampled;
    } else {
      comp = R::rbinom(sz - 1, (components_mean - 1) / (sz - 1)) + 1;
    }
    components[i] = comp;

    int k1 = R::rhyper(n_sampled, sz - n_sampled, comp);
    int k2 = comp - k1;
    k1_vec[i] = k1;
    k2_vec[i] = k2;

    List result = prufer_sequence_generator_cpp(sz, n_sampled, k1, k2);
    prufer_list[i] = result["prufer"];
    type_vec[i] = as<std::string>(result["type"]);
  }

  return List::create(
    Named("size") = sizes,
    Named("components") = components,
    Named("k1") = k1_vec,
    Named("k2") = k2_vec,
    Named("prufer") = prufer_list,
    Named("type") = type_vec
  );
}

// Batch computation of posterior weights
// This replaces the slow row-by-row pmap_dbl calls in the R posterior() function
//
// samples_list: list with elements size, components, k1, k2, prufer, type
// params: list with elements:
//   prior_N_dist, prior_N_mean, prior_N_size,
//   prior_k_dist, prior_k_mean, prior_k_size,
//   nb_mean, nb_size, components_mean,
//   fix_N, fix_k, cluster
// [[Rcpp::export]]
List compute_posterior_batch_cpp(IntegerVector data, List samples_list, List params) {

  IntegerVector sizes = samples_list["size"];
  IntegerVector components_vec = samples_list["components"];
  IntegerVector k1_vec = samples_list["k1"];
  IntegerVector k2_vec = samples_list["k2"];
  List prufer_list = samples_list["prufer"];
  CharacterVector type_vec = samples_list["type"];

  std::string prior_N_dist = as<std::string>(params["prior_N_dist"]);
  double prior_N_mean = as<double>(params["prior_N_mean"]);
  double prior_N_size = as<double>(params["prior_N_size"]);
  std::string prior_k_dist = as<std::string>(params["prior_k_dist"]);
  double prior_k_mean = as<double>(params["prior_k_mean"]);
  double prior_k_size = as<double>(params["prior_k_size"]);
  double nb_mean = as<double>(params["nb_mean"]);
  double nb_size = as<double>(params["nb_size"]);
  double components_mean = as<double>(params["components_mean"]);
  bool fix_N = as<bool>(params["fix_N"]);
  bool fix_k = as<bool>(params["fix_k"]);
  bool cluster = as<bool>(params["cluster"]);

  int Nreps = sizes.size();
  int n = data.size();
  int data_sum = 0;
  for (int i = 0; i < n; i++) data_sum += data[i];
  int min_size = data_sum;

  NumericVector actual_size_vec(Nreps);
  NumericVector likelihood_vec(Nreps);
  NumericVector prior_vec(Nreps);
  NumericVector samp_prob_vec(Nreps);
  NumericVector log_weight_vec(Nreps);
  List data_order_list(Nreps);
  CharacterVector prufer_str_vec(Nreps);
  CharacterVector root_sizes_vec(Nreps);
  CharacterVector nonroot_sizes_vec(Nreps);

  for (int rep = 0; rep < Nreps; rep++) {
    int sz = sizes[rep];
    int cmp = components_vec[rep];
    int k1 = k1_vec[rep];
    int k2 = k2_vec[rep];

    // Permute data (Fisher-Yates shuffle)
    IntegerVector perm_data = clone(data);
    shuffle_inplace(perm_data);
    data_order_list[rep] = perm_data;

    // Actual size (expanded)
    actual_size_vec[rep] = sz + data_sum - n;

    // Likelihood
    // Convert NULL sub-sequences to empty IntegerVectors (R stores c() as NULL)
    List prufer = as<List>(prufer_list[rep]);
    for (int j = 0; j < prufer.size(); j++) {
      if (prufer[j] == R_NilValue) {
        prufer[j] = IntegerVector(0);
      }
    }
    std::string tp = as<std::string>(type_vec[rep]);
    likelihood_vec[rep] = log_likelihood_cpp(perm_data, prufer, cmp, tp, cluster);

    // Log fields
    prufer_str_vec[rep] = format_prufer(prufer, tp);
    root_sizes_vec[rep] = format_cluster_sizes(perm_data, k1, true);
    nonroot_sizes_vec[rep] = format_cluster_sizes(perm_data, k1, false);

    // Prior
    double prior = 0.0;
    int expanded = sz + data_sum - n;
    if (!fix_N) {
      if (prior_N_dist == "poisson") {
        prior += R::dpois(expanded - min_size, prior_N_mean, true);
      } else {
        prior += R::dnbinom_mu(expanded - min_size, prior_N_size, prior_N_mean, true);
      }
    }
    if (!fix_k) {
      if (prior_k_dist == "poisson") {
        prior += R::dpois(cmp - 1, prior_k_mean, true);
      } else {
        prior += R::dnbinom_mu(cmp - 1, prior_k_size, prior_k_mean, true);
      }
    }
    prior_vec[rep] = prior;

    // Sampling probability
    double log_prob = 0.0;
    if (!fix_N) {
      log_prob += R::dnbinom_mu(sz - n, nb_size, nb_mean, true);
    }
    if (!fix_k) {
      if (sz > n) {
        log_prob += R::dbinom(cmp - 1, sz - 1, (components_mean - 1.0) / (sz - 1), true);
      }
    }
    log_prob += R::dhyper(k1, n, sz - n, cmp, true);
    log_prob -= log_forest_count_cpp(sz, n, k1, k2, true);
    log_prob -= R::lgammafn(n + 1);  // log(n!) for root assignment
    samp_prob_vec[rep] = log_prob;

    // Weight
    log_weight_vec[rep] = likelihood_vec[rep] + prior_vec[rep] - samp_prob_vec[rep];
  }

  return List::create(
    Named("actual.size") = actual_size_vec,
    Named("likelihood") = likelihood_vec,
    Named("prior") = prior_vec,
    Named("samp.prob") = samp_prob_vec,
    Named("log.weight") = log_weight_vec,
    Named("data.order") = data_order_list,
    Named("prufer.string") = prufer_str_vec,
    Named("root.cluster.sizes") = root_sizes_vec,
    Named("nonroot.cluster.sizes") = nonroot_sizes_vec
  );
}

