// The data-driven model of cases and deaths, fitted together (Section "A
// data-driven model" and the appendix on it). Each series has its own delay
// distribution (for cases, with a tail beyond d_c common to all onset days),
// report effects, overdispersion and removals. The new count of an onset day
// between two reports is NB additions less Poisson removals. The series
// are joined through the expected eventual counts: the expected deaths of
// onset day t are its expected cases times a case fatality ratio,
//   lambda_d(t) = cfr(t) * lambda_c(t),
// with logit cfr(t) a random walk over onset days. The deaths of recent onset
// days, most of which are still to come, then follow the cases of those days
// (which are more complete) and a ratio which changes slowly, rather than the
// deaths' own random walk.
//
// Suffix _c is for cases and _d for deaths. The likelihood of each series is
// summed in parallel over its new counts (reduce_sum), and the delay
// probabilities are computed where they are used rather than kept as
// transformed parameters, so that the saved draws stay small. day_target is
// data, so the draws can be projected afterwards to any day by standalone
// generated quantities.
functions {
  // log-probability of a new count n which is NB additions less Poisson
  // removals, the removals summed out up to w_max
  real nb_minus_poisson_lpmf(int n, real log_mu_a, real phi, real mu_w, int w_max) {
    int w0 = n < 0 ? -n : 0;
    if (w0 > w_max)
      return negative_infinity();
    vector[w_max - w0 + 1] terms;
    for (w in w0:w_max)
      terms[w - w0 + 1] = neg_binomial_2_log_lpmf(n + w | log_mu_a, phi) +
                          poisson_lpmf(w | mu_w);
    return log_sum_exp(terms);
  }
  // the log-likelihood of the new counts start..end of a series, given the
  // log of their expected additions, for reduce_sum
  real partial_sum_lpmf(array[] int n_slice, int start, int end, vector log_mu_a,
                        array[] int w_max, real phi, real rho) {
    real lp = 0;
    for (i in start:end)
      lp += nb_minus_poisson_lpmf(n_slice[i - start + 1] | log_mu_a[i], phi,
                                  rho * exp(log_mu_a[i]), w_max[i]);
    return lp;
  }
  // the log delay distribution of each onset day: a softmax over delays of
  // the reference, the average shape and the departure of the onset day,
  // split at d_c into a part which changes with onset time and a common tail
  // with probability q if d_c < D
  matrix delay_log_p(int T, int D, int Kt, int Kd, matrix Bt, matrix Bd,
                     vector log_p_ref, vector beta_g_free, matrix z_h,
                     real tau_t, int d_c, array[] real q) {
    matrix[Kt, Kd - 1] steps = append_row(z_h[1], tau_t * z_h[2:Kt]);
    matrix[Kt, Kd - 1] cum = steps;
    for (k in 2:Kt)
      cum[k] = cum[k - 1] + steps[k];
    matrix[Kt, Kd] beta_h = append_col(rep_vector(0, Kt), cum);
    vector[Kd] beta_g = append_row(0, beta_g_free);
    vector[D + 1] eta_common = log_p_ref + Bd * beta_g;
    matrix[T, D + 1] eta = rep_matrix(eta_common', T) + Bt * beta_h * Bd';
    matrix[T, D + 1] log_p;
    if (d_c >= D) {
      for (i in 1:T)
        log_p[i] = log_softmax(eta[i]')';
    } else {
      row_vector[D - d_c] log_tail = (log(q[1]) + log_softmax(eta_common[(d_c + 2):(D + 1)]))';
      for (i in 1:T)
        log_p[i] = append_col(log1m(q[1]) + log_softmax(eta[i, 1:(d_c + 1)]')', log_tail);
    }
    return log_p;
  }
  // the count each onset day will show by day_target: its count in the last
  // report fitted plus the additions less removals of the delays still to
  // come, each with the effect of a future report
  array[] int project_rng(int T, int D, array[] int count_last, array[] int delay_last,
                          int day_target, vector log_lambda, matrix log_p,
                          vector u_future, real phi, real rho) {
    array[T] int projected;
    for (i in 1:T) {
      projected[i] = count_last[i];
      if (day_target > 0) {
        int d_end = min(D, day_target - i);
        for (d in (delay_last[i] + 1):d_end) {
          real log_mu_a = log_lambda[i] + log_p[i, d + 1] + u_future[i + d];
          if (log_mu_a > -20) {
            real mu_a = exp(log_mu_a);
            projected[i] += neg_binomial_2_rng(mu_a, phi) - poisson_rng(rho * mu_a);
          }
        }
      }
    }
    return projected;
  }
}
data {
  int<lower=1> T;                      // onset days, shared
  int<lower=1> D;                      // longest delay, shared
  int<lower=2> Kt;
  matrix[T, Kt] Bt;
  int<lower=1> R;                      // reports, shared
  int<lower=0> day_target;
  // cases
  int<lower=3> Kd_c;
  matrix[D + 1, Kd_c] Bd_c;
  int<lower=0> N_c;
  array[N_c] int<lower=1, upper=T> t_c;
  array[N_c] int<lower=0, upper=D> lo_c;
  array[N_c] int<lower=0, upper=D> hi_c;
  array[N_c] int<lower=1, upper=R> r_c;
  array[N_c] int n_c;
  array[N_c] int<lower=0> w_max_c;
  array[T] int count_last_c;
  array[T] int<lower=0, upper=D> delay_last_c;
  vector[D + 1] log_p_ref_c;
  int<lower=0, upper=D> d_c_c;         // D or more for no common tail
  real<lower=0> q_alpha_c;
  real<lower=0> q_beta_c;
  real<lower=0> sigma_lambda_scale_c;
  real<lower=0> beta_g_scale_c;
  real log_phi_mu_c;
  real<lower=0> log_phi_sd_c;
  real<lower=0> tau_t_scale_c;
  real<lower=0> sigma_u_scale_c;
  real<lower=0> rho_scale_c;
  // deaths
  int<lower=3> Kd_d;
  matrix[D + 1, Kd_d] Bd_d;
  int<lower=0> N_d;
  array[N_d] int<lower=1, upper=T> t_d;
  array[N_d] int<lower=0, upper=D> lo_d;
  array[N_d] int<lower=0, upper=D> hi_d;
  array[N_d] int<lower=1, upper=R> r_d;
  array[N_d] int n_d;
  array[N_d] int<lower=0> w_max_d;
  array[T] int count_last_d;
  array[T] int<lower=0, upper=D> delay_last_d;
  vector[D + 1] log_p_ref_d;
  int<lower=0, upper=D> d_c_d;
  real<lower=0> q_alpha_d;
  real<lower=0> q_beta_d;
  real<lower=0> beta_g_scale_d;
  real log_phi_mu_d;
  real<lower=0> log_phi_sd_d;
  real<lower=0> tau_t_scale_d;
  real<lower=0> sigma_u_scale_d;
  real<lower=0> rho_scale_d;
  // the case fatality ratio
  real logit_cfr_mu;
  real<lower=0> logit_cfr_sd;
  real<lower=0> sigma_cfr_scale;
  int<lower=1> grainsize;              // of reduce_sum
}
transformed data {
  int has_tail_c = d_c_c < D;
  int has_tail_d = d_c_d < D;
}
parameters {
  real log_lambda_1;
  vector[T - 1] z_lambda;
  real<lower=0> sigma_lambda;
  real logit_cfr_1;
  vector[T - 1] z_cfr;
  real<lower=0> sigma_cfr;
  // cases
  vector[Kd_c - 1] beta_g_free_c;
  matrix[Kt, Kd_c - 1] z_h_c;
  real<lower=0> tau_t_c;
  vector[R] z_u_c;
  real<lower=0> sigma_u_c;
  real log_phi_c;
  real<lower=0> rho_c;
  array[has_tail_c] real<lower=0, upper=1> q_c;
  // deaths
  vector[Kd_d - 1] beta_g_free_d;
  matrix[Kt, Kd_d - 1] z_h_d;
  real<lower=0> tau_t_d;
  vector[R] z_u_d;
  real<lower=0> sigma_u_d;
  real log_phi_d;
  real<lower=0> rho_d;
  array[has_tail_d] real<lower=0, upper=1> q_d;
}
transformed parameters {
  real phi_c = exp(log_phi_c);
  real phi_d = exp(log_phi_d);
  vector[T] log_lambda_c = log_lambda_1 +
    append_row(0, cumulative_sum(sigma_lambda * z_lambda));
  vector[T] logit_cfr = logit_cfr_1 + append_row(0, cumulative_sum(sigma_cfr * z_cfr));
  vector[T] log_lambda_d = log_lambda_c + log_inv_logit(logit_cfr);
  vector[R] u_c = sigma_u_c * z_u_c;
  vector[R] u_d = sigma_u_d * z_u_d;
}
model {
  log_lambda_1 ~ normal(log(5), 3);
  z_lambda ~ std_normal();
  sigma_lambda ~ normal(0, sigma_lambda_scale_c);
  logit_cfr_1 ~ normal(logit_cfr_mu, logit_cfr_sd);
  z_cfr ~ std_normal();
  sigma_cfr ~ normal(0, sigma_cfr_scale);
  beta_g_free_c ~ normal(0, beta_g_scale_c);
  to_vector(z_h_c) ~ std_normal();
  tau_t_c ~ normal(0, tau_t_scale_c);
  z_u_c ~ std_normal();
  sigma_u_c ~ normal(0, sigma_u_scale_c);
  log_phi_c ~ normal(log_phi_mu_c, log_phi_sd_c);
  rho_c ~ normal(0, rho_scale_c);
  q_c ~ beta(q_alpha_c, q_beta_c);
  beta_g_free_d ~ normal(0, beta_g_scale_d);
  to_vector(z_h_d) ~ std_normal();
  tau_t_d ~ normal(0, tau_t_scale_d);
  z_u_d ~ std_normal();
  sigma_u_d ~ normal(0, sigma_u_scale_d);
  log_phi_d ~ normal(log_phi_mu_d, log_phi_sd_d);
  rho_d ~ normal(0, rho_scale_d);
  q_d ~ beta(q_alpha_d, q_beta_d);
  {
    matrix[T, D + 1] log_p_c = delay_log_p(T, D, Kt, Kd_c, Bt, Bd_c, log_p_ref_c,
                                           beta_g_free_c, z_h_c, tau_t_c, d_c_c, q_c);
    matrix[T, D + 1] log_p_d = delay_log_p(T, D, Kt, Kd_d, Bt, Bd_d, log_p_ref_d,
                                           beta_g_free_d, z_h_d, tau_t_d, d_c_d, q_d);
    vector[N_c] log_mu_c;
    vector[N_d] log_mu_d;
    for (i in 1:N_c)
      log_mu_c[i] = log_lambda_c[t_c[i]] +
        log_sum_exp(log_p_c[t_c[i], (lo_c[i] + 1):(hi_c[i] + 1)]) + u_c[r_c[i]];
    for (i in 1:N_d)
      log_mu_d[i] = log_lambda_d[t_d[i]] +
        log_sum_exp(log_p_d[t_d[i], (lo_d[i] + 1):(hi_d[i] + 1)]) + u_d[r_d[i]];
    target += reduce_sum(partial_sum_lpmf, n_c, grainsize, log_mu_c, w_max_c, phi_c, rho_c);
    target += reduce_sum(partial_sum_lpmf, n_d, grainsize, log_mu_d, w_max_d, phi_d, rho_d);
  }
}
generated quantities {
  // a future report's effect is shared by every onset day it counts, and
  // drawn separately for cases and deaths
  array[T] int projected_c;
  array[T] int projected_d;
  vector[T] mean_delay_c;
  vector[T] mean_delay_d;
  vector[T] cfr = inv_logit(logit_cfr);
  {
    matrix[T, D + 1] log_p_c = delay_log_p(T, D, Kt, Kd_c, Bt, Bd_c, log_p_ref_c,
                                           beta_g_free_c, z_h_c, tau_t_c, d_c_c, q_c);
    matrix[T, D + 1] log_p_d = delay_log_p(T, D, Kt, Kd_d, Bt, Bd_d, log_p_ref_d,
                                           beta_g_free_d, z_h_d, tau_t_d, d_c_d, q_d);
    vector[max(day_target, 1)] u_future_c;
    vector[max(day_target, 1)] u_future_d;
    for (s in 1:max(day_target, 1)) {
      u_future_c[s] = normal_rng(0, sigma_u_c);
      u_future_d[s] = normal_rng(0, sigma_u_d);
    }
    projected_c = project_rng(T, D, count_last_c, delay_last_c, day_target,
                              log_lambda_c, log_p_c, u_future_c, phi_c, rho_c);
    projected_d = project_rng(T, D, count_last_d, delay_last_d, day_target,
                              log_lambda_d, log_p_d, u_future_d, phi_d, rho_d);
    for (i in 1:T) {
      mean_delay_c[i] = dot_product(exp(log_p_c[i]), linspaced_vector(D + 1, 0, D));
      mean_delay_d[i] = dot_product(exp(log_p_d[i]), linspaced_vector(D + 1, 0, D));
    }
  }
}
