setwd("D:/randomslopes/simulations_randomslopes")
library(dplyr)

kl_div_cov_matrix=function(Sigma_hat,Sigma_true){
  q = ncol(Sigma_true)
  eigen_values_Sigma_hat=eigen(Sigma_hat)$values
  eigen_values_Sigma_true=(eigen(Sigma_true))$values
  return((sum(diag(solve(Sigma_hat)%*%(Sigma_true))) + sum(log(eigen_values_Sigma_hat)) - sum(log(eigen_values_Sigma_true)) - q)/2)    
}


p = 3
Sigmaa = diag(0.3, p + 1)
Sigmab = diag(0.1, p + 1)

diagonal = 1
Gam = seq(3.2, 5, 0.2)
conv_thres_exp = 6
lmer = 1
final_summary = NULL
for(gam in Gam){
  N = round(10^gam)
  R = round(N^(0.6))
  niter = 100
  file = sprintf("randomslopes_after_review_N=%d_R=%d_p=%d_seed=1_niter=%d_conv_thres_exp=%d_diagonal_%d_lmer_%d.RDS", N,R,p,niter,conv_thres_exp, diagonal, lmer)
  results = readRDS(file)
  niter = length(results)
  df = NULL
  for(i in 1:niter){
    a = results[[i]][1:43]
    df = rbind(df, as.numeric(a))
  }
  
  df = data.frame(df)
  final_summary = rbind(final_summary, colMeans(df))
  
}
required_values = c("N","mom_status", "moment_estimation_user_time", "backfitting_user_time", "var_backfitting_user_time", "clubbed_backfitting_user_time", "clubbed_var_backfitting_user_time",
                    "moment_estimation_elapsed_time", "backfitting_elapsed_time", "var_backfitting_elapsed_time", "clubbed_backfitting_elapsed_time", "clubbed_var_backfitting_elapsed_time",
                    "niter_backfit", "niter_var_backfit", "niter_clubbed_backfit", "niter_clubbed_var_backfit", 
                    "err_beta_backfit", "err_beta_var_backfit", "err_beta_clubbed_backfit", "err_beta_clubbed_var_backfit",
                    "err_sigma2e_mom", "err_sigma2e_var_backfit", "err_sigma2e_clubbed_var_backfit", 
                    "kl_err_Sigma_a_mom", "kl_err_Sigma_a_var_backfit", "kl_err_Sigma_a_clubbed_var_backfit", 
                    "kl_err_Sigma_b_mom", "kl_err_Sigma_b_var_backfit", "kl_err_Sigma_b_clubbed_var_backfit", 
                    "frob_err_Sigma_a_mom", "frob_err_Sigma_a_var_backfit", "frob_err_Sigma_a_clubbed_var_backfit", 
                    "frob_err_Sigma_b_mom", "frob_err_Sigma_b_var_backfit", "frob_err_Sigma_b_clubbed_var_backfit", 
                    "max_lik_user_time", "max_lik_elapsed_time", "niter_mle", "err_beta_mle", "err_sigma2e_mle", "kl_err_Sigma_a_mle", "kl_err_Sigma_b_mle", "frob_err_Sigma_a_mle", "frob_err_Sigma_b_mle")

final_summary = cbind(round(10^Gam), final_summary)
colnames(final_summary) = required_values
write.csv(final_summary, "../results/final_summary_diagonal_with_lmer.csv", row.names = FALSE)


Gam = seq(5.2, 6.2, 0.2)
lmer = 0
final_summary = NULL
for(gam in Gam){
  N = round(10^gam)
  R = round(N^(0.6))
  niter = 100
  file = sprintf("randomslopes_after_review_N=%d_R=%d_p=%d_seed=1_niter=%d_conv_thres_exp=%d_diagonal_%d_lmer_%d.RDS", N,R,p,niter,conv_thres_exp, diagonal, lmer)
  results = readRDS(file)
  niter = length(results)
  df = NULL
  for(i in 1:niter){
    a = results[[i]][1:34]
    df = rbind(df, as.numeric(a))
  }
  
  df = data.frame(df)
  final_summary = rbind(final_summary, colMeans(df))
  
}
required_values = c("N","mom_status", "moment_estimation_user_time", "backfitting_user_time", "var_backfitting_user_time", "clubbed_backfitting_user_time", "clubbed_var_backfitting_user_time",
                    "moment_estimation_elapsed_time", "backfitting_elapsed_time", "var_backfitting_elapsed_time", "clubbed_backfitting_elapsed_time", "clubbed_var_backfitting_elapsed_time",
                    "niter_backfit", "niter_var_backfit", "niter_clubbed_backfit", "niter_clubbed_var_backfit", 
                    "err_beta_backfit", "err_beta_var_backfit", "err_beta_clubbed_backfit", "err_beta_clubbed_var_backfit",
                    "err_sigma2e_mom", "err_sigma2e_var_backfit", "err_sigma2e_clubbed_var_backfit", 
                    "kl_err_Sigma_a_mom", "kl_err_Sigma_a_var_backfit", "kl_err_Sigma_a_clubbed_var_backfit", 
                    "kl_err_Sigma_b_mom", "kl_err_Sigma_b_var_backfit", "kl_err_Sigma_b_clubbed_var_backfit", 
                    "frob_err_Sigma_a_mom", "frob_err_Sigma_a_var_backfit", "frob_err_Sigma_a_clubbed_var_backfit", 
                    "frob_err_Sigma_b_mom", "frob_err_Sigma_b_var_backfit", "frob_err_Sigma_b_clubbed_var_backfit")
                    
final_summary = cbind(round(10^Gam), final_summary)
colnames(final_summary) = required_values
write.csv(final_summary, "../results/final_summary_diagonal_without_lmer.csv", row.names = FALSE)


Sigmaa=diag(0.8,(p+1)) + matrix(0.2,(p+1),(p+1))
Sigmab=diag(0.8,(p+1)) + matrix(0.2,(p+1),(p+1))

diagonal = 0
Gam = seq(3.2, 5.0, 0.2)
conv_thres_exp = 6
lmer = 1
final_summary = NULL
for(gam in Gam){
  N = round(10^gam)
  R = round(N^(0.6))
  p = 3
  niter = 100
  file = sprintf("randomslopes_after_review_N=%d_R=%d_p=%d_seed=1_niter=%d_conv_thres_exp=%d_diagonal_%d_lmer_%d.RDS", N,R,p,niter,conv_thres_exp, diagonal, lmer)
  results = readRDS(file)
  niter = length(results)
  df = NULL
  for(i in 1:niter){
    a = c(results[[i]][1:43])
    df = rbind(df, as.numeric(a))
  }
  
  df = data.frame(df)
  final_summary = rbind(final_summary, colMeans(df))
  
}
required_values = c("N","mom_status", "moment_estimation_user_time", "backfitting_user_time", "var_backfitting_user_time", "clubbed_backfitting_user_time", "clubbed_var_backfitting_user_time",
                    "moment_estimation_elapsed_time", "backfitting_elapsed_time", "var_backfitting_elapsed_time", "clubbed_backfitting_elapsed_time", "clubbed_var_backfitting_elapsed_time",
                    "niter_backfit", "niter_var_backfit", "niter_clubbed_backfit", "niter_clubbed_var_backfit", 
                    "err_beta_backfit", "err_beta_var_backfit", "err_beta_clubbed_backfit", "err_beta_clubbed_var_backfit",
                    "err_sigma2e_mom", "err_sigma2e_var_backfit", "err_sigma2e_clubbed_var_backfit", 
                    "kl_err_Sigma_a_mom", "kl_err_Sigma_a_var_backfit", "kl_err_Sigma_a_clubbed_var_backfit", 
                    "kl_err_Sigma_b_mom", "kl_err_Sigma_b_var_backfit", "kl_err_Sigma_b_clubbed_var_backfit", 
                    "frob_err_Sigma_a_mom", "frob_err_Sigma_a_var_backfit", "frob_err_Sigma_a_clubbed_var_backfit", 
                    "frob_err_Sigma_b_mom", "frob_err_Sigma_b_var_backfit", "frob_err_Sigma_b_clubbed_var_backfit", 
                    "max_lik_user_time","max_lik_elapsed_time", "niter_mle", "err_beta_mle", "err_sigma2e_mle", "kl_err_Sigma_a_mle", "kl_err_Sigma_b_mle", "frob_err_Sigma_a_mle", "frob_err_Sigma_b_mle")

final_summary = cbind(round(10^Gam), final_summary)
colnames(final_summary) = required_values
write.csv(final_summary, "../results/final_summary_non_diagonal_with_lmer.csv", row.names = FALSE)


Gam = seq(5.2, 6.0, 0.2)
lmer = 0
final_summary = NULL
for(gam in Gam){
  N = round(10^gam)
  R = round(N^(0.6))
  p = 3
  niter = 100
  file = sprintf("randomslopes_after_review_N=%d_R=%d_p=%d_seed=1_niter=%d_conv_thres_exp=%d_diagonal_%d_lmer_%d.RDS", N,R,p,niter,conv_thres_exp, diagonal, lmer)
  results = readRDS(file)
  niter = length(results)
  df = NULL
  for(i in 1:niter){
    a = c(results[[i]][1:34])
    df = rbind(df, as.numeric(a))
  }
  
  df = data.frame(df)
  final_summary = rbind(final_summary, colMeans(df))
  
}
required_values = c("N","mom_status", "moment_estimation_user_time", "backfitting_user_time", "var_backfitting_user_time", "clubbed_backfitting_user_time", "clubbed_var_backfitting_user_time",
                    "moment_estimation_elapsed_time", "backfitting_elapsed_time", "var_backfitting_elapsed_time", "clubbed_backfitting_elapsed_time", "clubbed_var_backfitting_elapsed_time",
                    "niter_backfit", "niter_var_backfit", "niter_clubbed_backfit", "niter_clubbed_var_backfit", 
                    "err_beta_backfit", "err_beta_var_backfit", "err_beta_clubbed_backfit", "err_beta_clubbed_var_backfit",
                    "err_sigma2e_mom", "err_sigma2e_var_backfit", "err_sigma2e_clubbed_var_backfit", 
                    "kl_err_Sigma_a_mom", "kl_err_Sigma_a_var_backfit", "kl_err_Sigma_a_clubbed_var_backfit", 
                    "kl_err_Sigma_b_mom", "kl_err_Sigma_b_var_backfit", "kl_err_Sigma_b_clubbed_var_backfit", 
                    "frob_err_Sigma_a_mom", "frob_err_Sigma_a_var_backfit", "frob_err_Sigma_a_clubbed_var_backfit", 
                    "frob_err_Sigma_b_mom", "frob_err_Sigma_b_var_backfit", "frob_err_Sigma_b_clubbed_var_backfit")

final_summary = cbind(round(10^Gam), final_summary)
colnames(final_summary) = required_values
write.csv(final_summary, "../results/final_summary_non_diagonal_without_lmer.csv", row.names = FALSE)

non_diagonal_with_lmer = read.csv("../results/final_summary_non_diagonal_with_lmer.csv")
non_diagonal_without_lmer = read.csv("../results/final_summary_non_diagonal_without_lmer.csv")

summaries_binded = bind_rows(non_diagonal_with_lmer, non_diagonal_without_lmer)
write.csv(summaries_binded, "../results/non_diagonal_binded.csv", row.names = FALSE)


diagonal_with_lmer = read.csv("../results/final_summary_diagonal_with_lmer.csv")
diagonal_without_lmer = read.csv("../results/final_summary_diagonal_without_lmer.csv")

summaries_binded = bind_rows(diagonal_with_lmer, diagonal_without_lmer)
write.csv(summaries_binded, "../results/diagonal_binded.csv", row.names = FALSE)
