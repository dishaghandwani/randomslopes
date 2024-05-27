library(dplyr)
library(tibble)
library(mvtnorm)
library(CVXR)
library(lme4)
library(parallel)

source("randomslopes.R")

kl_div_cov_matrix=function(Sigma_hat,Sigma_true){
    q = ncol(Sigma_true)
    eigen_values_Sigma_hat=eigen(Sigma_hat)$values
    eigen_values_Sigma_true=(eigen(Sigma_true))$values
    return((sum(diag(solve(Sigma_hat)%*%(Sigma_true))) + sum(log(eigen_values_Sigma_hat)) - sum(log(eigen_values_Sigma_true)) - q)/2)    
}

results = function(N, R, p, conv_thres, diagonal = FALSE, lmer = TRUE){


    C = R
    f1 <- sample(1:R, N, replace=TRUE)
    f2 <- sample(1:C, N, replace=TRUE)
    beta = 0.1*(1:(p+1))
    sigmae = 1
    
    if(diagonal){
        Sigmaa = diag(0.3, p + 1)
        Sigmab = diag(0.1, p + 1)
    }else{
        Sigmaa=diag(0.8,(p+1)) + matrix(0.2,(p+1),(p+1))
        Sigmab=diag(0.8,(p+1)) + matrix(0.2,(p+1),(p+1))
    }
    
    X <- cbind(1,matrix(rnorm(N*p),ncol=p))
    A = rmvnorm(R,sigma=Sigmaa)
    B = rmvnorm(C,sigma=Sigmab)
    
    y = X%*%beta + rowSums(X*A[f1,]) + rowSums(X*B[f2,]) + rnorm(N,sd=sigmae)
    p_a = p + 1
    p_b = p + 1
    
    if(diagonal){
        moment_estimation_time <- system.time(mom <- method_of_moments_diagonal(y, X, X, X, f1, f2))
        moment_estimation_user_time <- as.numeric(moment_estimation_time[1])
        moment_estimation_elapsed_time <- as.numeric(moment_estimation_time[3])
        if(all(diag(mom$Sigma_a) > 0) &  all(diag(mom$Sigma_b) > 0) & mom$sigma2e > 0){
        
            print("Method of moments approach worked")
            mom_status = 1
            mm_Sigma_a <- mom$Sigma_a
            mm_Sigma_b <- mom$Sigma_b
            mm_sigma2e <- mom$sigma2e
            
            kl_err_Sigma_a_mom <- kl_div_cov_matrix(mm_Sigma_a, Sigmaa)
            frob_err_Sigma_a_mom <- norm(mm_Sigma_a - Sigmaa,type="2")
            kl_err_Sigma_b_mom <- kl_div_cov_matrix(mm_Sigma_b, Sigmab)
            frob_err_Sigma_b_mom <- norm(mm_Sigma_b - Sigmab,type="2")
            err_sigma2e_mom <- abs(mm_sigma2e - 1)
          
            Za = mom$Za
            Zb = mom$Zb
            
          
          
        }else{
            print("Method of moments approach failed")
            mom_status = 0
            mm_Sigma_a = lam = diag(p+1)
            mm_Sigma_b = alpha = diag(p+1)
            mm_sigma2e = 1
            
            
            Za = mom$Za
            Zb = mom$Zb
            
            kl_err_Sigma_a_mom = NA
            frob_err_Sigma_a_mom = NA
            kl_err_Sigma_b_mom = NA
            frob_err_Sigma_b_mom = NA
            err_sigma2e_mom <- NA
            
            
        }
    
    }else{
        moment_estimation_time <- system.time(mom <- method_of_moments_non_diagonal(y, X, X, X, f1, f2))
        moment_estimation_user_time <- as.numeric(moment_estimation_time[1])
        moment_estimation_elapsed_time <- as.numeric(moment_estimation_time[3])
        if(mom$result_status == 'optimal'){
            print("Method of moments approach worked")
            mom_status <- 1
            mm_Sigma_a <- mom$Sigma_a
            mm_Sigma_b <- mom$Sigma_b
            mm_sigma2e <- mom$sigma2e
            
            kl_err_Sigma_a_mom <- kl_div_cov_matrix(mm_Sigma_a, Sigmaa)
            frob_err_Sigma_a_mom <- norm(mm_Sigma_a - Sigmaa,type="2")
            kl_err_Sigma_b_mom <- kl_div_cov_matrix(mm_Sigma_b, Sigmab)
            frob_err_Sigma_b_mom <- norm(mm_Sigma_b - Sigmab,type="2")
            err_sigma2e_mom <- abs(mm_sigma2e - 1)

            
            Za = mom$Za
            Zb = mom$Zb                
        }else{
            print("Method of moments approach failed")
            mom_status = 0
            mm_Sigma_a = lam = diag(p_a)
            mm_Sigma_b = alpha = diag(p_b)
            mm_sigma2e = 1
            Za = mom$Za
            Zb = mom$Zb
            
            kl_err_Sigma_a_mom = NA
            frob_err_Sigma_a_mom = NA
            kl_err_Sigma_b_mom = NA
            frob_err_Sigma_b_mom = NA
            err_sigma2e_mom <- NA
            
 
        }
    
    }
    backfitting_time <- system.time(back_fit <- vanilla_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 1000, conv_thres = conv_thres))
    backfitting_user_time <- as.numeric(backfitting_time[1])
    backfitting_elapsed_time <- as.numeric(backfitting_time[3])
    niter_backfit <- back_fit$niter
    beta_backfit <- back_fit$beta
    err_beta_backfit <- sum((back_fit$beta-beta)^2)
    mean_ai_backfit <- rowMeans(back_fit$A)
    mean_bj_backfit <- rowMeans(back_fit$B)
    
    if(diagonal){
        var_backfitting_time <- system.time(var_backfit <- variational_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 1000, conv_thres = conv_thres, covariances_diag_Sigma_a = TRUE, covariances_diag_Sigma_b = TRUE))
    }else{
        var_backfitting_time <- system.time(var_backfit <- variational_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 1000, conv_thres = conv_thres))
    }
    var_backfitting_user_time <- as.numeric(var_backfitting_time[1])
    var_backfitting_elapsed_time <- as.numeric(var_backfitting_time[3])
    
    niter_var_backfit <- var_backfit$niter
    beta_var_backfit <- var_backfit$beta
    err_beta_var_backfit <- sum((var_backfit$beta - beta)^2)
    
    Sigma_a_var_backfit <- var_backfit$Sigma_a
    Sigma_b_var_backfit <- var_backfit$Sigma_b
    sigma2e_var_backfit <- var_backfit$sigma2e
    err_sigma2e_var_backfit <- abs(var_backfit$sigma2e/(sigmae^2) - 1)
    
    kl_err_Sigma_a_var_backfit <- kl_div_cov_matrix(var_backfit$Sigma_a, Sigmaa)
    kl_err_Sigma_b_var_backfit <- kl_div_cov_matrix(var_backfit$Sigma_b, Sigmab)
    
    frob_err_Sigma_a_var_backfit <- norm(var_backfit$Sigma_a - Sigmaa, type="2")
    frob_err_Sigma_b_var_backfit <- norm(var_backfit$Sigma_b - Sigmab, type="2")
    
    clubbed_backfitting_time <- system.time(club_backfit <- clubbed_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 800, conv_thres = conv_thres))
    clubbed_backfitting_user_time <- as.numeric(clubbed_backfitting_time[1])
    clubbed_backfitting_elapsed_time <- as.numeric(clubbed_backfitting_time[3])
    niter_clubbed_backfit <- club_backfit$niter
    beta_clubbed_backfit <- club_backfit$beta
    err_beta_clubbed_backfit <- sum((club_backfit$beta-beta)^2)
    
    if(diagonal){
        clubbed_var_backfitting_time <- system.time(clubbed_var_backfit <- clubbed_variational_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 800, conv_thres = conv_thres, covariances_diag_Sigma_a = TRUE, covariances_diag_Sigma_b = TRUE))
    }else{
        clubbed_var_backfitting_time <- system.time(clubbed_var_backfit <- clubbed_variational_backfitting(y, X, X, X, f1, f2, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, max_iter = 800, conv_thres = conv_thres))
    }
    clubbed_var_backfitting_user_time <- as.numeric(clubbed_var_backfitting_time[1])
    clubbed_var_backfitting_elapsed_time <- as.numeric(clubbed_var_backfitting_time[3])
    niter_clubbed_var_backfit <- clubbed_var_backfit$niter
    beta_clubbed_var_backfit <- clubbed_var_backfit$beta
    err_beta_clubbed_var_backfit <- sum((clubbed_var_backfit$beta - beta)^2)
    
    Sigma_a_clubbed_var_backfit <- clubbed_var_backfit$Sigma_a
    Sigma_b_clubbed_var_backfit <- clubbed_var_backfit$Sigma_b
    sigma2e_clubbed_var_backfit <- clubbed_var_backfit$sigma2e
    err_sigma2e_clubbed_var_backfit <- abs(clubbed_var_backfit$sigma2e/(sigmae^2) - 1)
    
    kl_err_Sigma_a_clubbed_var_backfit <- kl_div_cov_matrix(clubbed_var_backfit$Sigma_a, Sigmaa)
    kl_err_Sigma_b_clubbed_var_backfit <- kl_div_cov_matrix(clubbed_var_backfit$Sigma_b, Sigmab)
    
    frob_err_Sigma_a_clubbed_var_backfit <- norm(clubbed_var_backfit$Sigma_a - Sigmaa, type="2")
    frob_err_Sigma_b_clubbed_var_backfit <- norm(clubbed_var_backfit$Sigma_b - Sigmab, type="2")

  
    if(lmer){
        if(diagonal){
            data = data.frame(x1 = X[,1], x2 = X[, 2], x3 = X[, 3], x4 = X[,4], y = y, f1 = f1, f2 = f2)
            max_lik_time <- system.time(max_lik <- lmer(y ~ 0 + x1 + x2 + x3 + x4  + (1|f1) + (1|f2)
                + (0 + x2 | f1) + (0 + x2|f2) + (0 + x3|f1) + (0 + x3|f2) + (0 + x4|f1) + (0 + x4|f2), data, REML = FALSE))
            max_lik_user_time <- as.numeric(max_lik_time[1])
            max_lik_elapsed_time <- as.numeric(max_lik_time[3])
            
            fit_mle = summary(max_lik)
            niter_mle <- fit_mle$optinfo$feval
            beta_mle <- fit_mle$coefficients[,1]
            err_beta_mle <- sum((fit_mle$coefficients[,1] - beta)^2)
            Sigma_a_mle <- diag(as.numeric(summary(max_lik)$varcor[c(1,3,5,7)]))
            Sigma_b_mle <- diag(as.numeric(summary(max_lik)$varcor[c(2,4,6,8)]))
            kl_err_Sigma_a_mle <- kl_div_cov_matrix(Sigma_a_mle, Sigmaa)
            kl_err_Sigma_b_mle <- kl_div_cov_matrix(Sigma_b_mle, Sigmab)
            
            frob_err_Sigma_a_mle <- norm(Sigma_a_mle - Sigmaa, type = "2")
            frob_err_Sigma_b_mle <- norm(Sigma_b_mle - Sigmab, type = "2")
            sigma2e_mle <- (fit_mle$sigma)^2
            err_sigma2e_mle <- abs(sigma2e_mle/(sigmae^2) - 1)
            
        }else{
            data = data.frame(X = X, y = y)
            max_lik_time <- system.time(max_lik <- lmer(y ~ 0 + X + (0+X|f1) +(0+X|f2),data, REML = FALSE))
            max_lik_user_time <- as.numeric(max_lik_time[1])
            max_lik_elapsed_time <- as.numeric(max_lik_time[3])
            
            fit_mle = summary(max_lik)
            niter_mle <- fit_mle$optinfo$feval
            beta_mle <- fit_mle$coefficients[,1]
            err_beta_mle <- sum((fit_mle$coefficients[,1] - beta)^2)
            Sigma_a_mle <- fit_mle$varcor$f1
            Sigma_b_mle <- fit_mle$varcor$f2
            
            kl_err_Sigma_a_mle <- kl_div_cov_matrix(matrix(Sigma_a_mle[1:(p_a^2)], p_a, p_a), Sigmaa)
            kl_err_Sigma_b_mle <- kl_div_cov_matrix(matrix(Sigma_b_mle[1:(p_b^2)], p_b, p_b), Sigmab)
            
            frob_err_Sigma_a_mle <- norm(matrix(Sigma_a_mle[1:(p_a^2)], p_a, p_a)-Sigmaa, type = "2")
            frob_err_Sigma_b_mle <- norm(matrix(Sigma_b_mle[1:(p_a^2)], p_b, p_b)-Sigmab, type = "2")
            sigma2e_mle <- (fit_mle$sigma)^2
            err_sigma2e_mle <- abs(sigma2e_mle/(sigmae^2) - 1)
            
        
        }
        
    
        
        
        return(enlist(mom_status, moment_estimation_user_time, backfitting_user_time, var_backfitting_user_time, clubbed_backfitting_user_time, clubbed_var_backfitting_user_time,
            moment_estimation_elapsed_time, backfitting_elapsed_time, var_backfitting_elapsed_time, clubbed_backfitting_elapsed_time, clubbed_var_backfitting_elapsed_time,
            niter_backfit, niter_var_backfit, niter_clubbed_backfit, niter_clubbed_var_backfit,
            err_beta_backfit, err_beta_var_backfit, err_beta_clubbed_backfit, err_beta_clubbed_var_backfit, 
            err_sigma2e_mom, err_sigma2e_var_backfit, err_sigma2e_clubbed_var_backfit,
            kl_err_Sigma_a_mom, kl_err_Sigma_a_var_backfit, kl_err_Sigma_a_clubbed_var_backfit, 
            kl_err_Sigma_b_mom, kl_err_Sigma_b_var_backfit, kl_err_Sigma_b_clubbed_var_backfit, 
            frob_err_Sigma_a_mom, frob_err_Sigma_a_var_backfit, frob_err_Sigma_a_clubbed_var_backfit, 
            frob_err_Sigma_b_mom, frob_err_Sigma_b_var_backfit, frob_err_Sigma_b_clubbed_var_backfit, 
            max_lik_user_time, max_lik_elapsed_time, niter_mle, err_beta_mle, err_sigma2e_mle, kl_err_Sigma_a_mle, kl_err_Sigma_b_mle, frob_err_Sigma_a_mle, frob_err_Sigma_b_mle,
            moment_estimation_time, backfitting_time, var_backfitting_time, clubbed_backfitting_time, clubbed_var_backfitting_time, max_lik_time,
            mm_Sigma_a, mm_Sigma_b, mm_sigma2e,
            beta_backfit,
            Sigma_a_var_backfit, Sigma_b_var_backfit, sigma2e_var_backfit, beta_var_backfit,
            beta_clubbed_backfit, 
            Sigma_a_clubbed_var_backfit, Sigma_b_clubbed_var_backfit, sigma2e_clubbed_var_backfit, beta_clubbed_var_backfit,
            beta_mle, Sigma_a_mle, Sigma_b_mle, sigma2e_mle 
            # current_seed
        ))
    }else{
        return(enlist(mom_status, moment_estimation_user_time, backfitting_user_time, var_backfitting_user_time, clubbed_backfitting_user_time, clubbed_var_backfitting_user_time,
            moment_estimation_elapsed_time, backfitting_elapsed_time, var_backfitting_elapsed_time, clubbed_backfitting_elapsed_time, clubbed_var_backfitting_elapsed_time,
            niter_backfit, niter_var_backfit, niter_clubbed_backfit, niter_clubbed_var_backfit,
            err_beta_backfit, err_beta_var_backfit, err_beta_clubbed_backfit, err_beta_clubbed_var_backfit,
            err_sigma2e_mom, err_sigma2e_var_backfit, err_sigma2e_clubbed_var_backfit,
            kl_err_Sigma_a_mom, kl_err_Sigma_a_var_backfit, kl_err_Sigma_a_clubbed_var_backfit,
            kl_err_Sigma_b_mom, kl_err_Sigma_b_var_backfit, kl_err_Sigma_b_clubbed_var_backfit,
            frob_err_Sigma_a_mom, frob_err_Sigma_a_var_backfit, frob_err_Sigma_a_clubbed_var_backfit,
            frob_err_Sigma_b_mom, frob_err_Sigma_b_var_backfit, frob_err_Sigma_b_clubbed_var_backfit,
            moment_estimation_time, backfitting_time, var_backfitting_time, clubbed_backfitting_time, clubbed_var_backfitting_time, 
            mm_Sigma_a, mm_Sigma_b, mm_sigma2e,
            beta_backfit,
            Sigma_a_var_backfit, Sigma_b_var_backfit, sigma2e_var_backfit, beta_var_backfit,
            beta_clubbed_backfit, 
            Sigma_a_clubbed_var_backfit, Sigma_b_clubbed_var_backfit, sigma2e_clubbed_var_backfit, beta_clubbed_var_backfit
            # current_seed
                      
        ))    
    
    
    }
  
  
  
}


total = function(gam, diagonal = 0, lmer = 1, conv_thres_exp = 6, niter = 100){
  
    N = round(10^gam)
    R = round(N^(0.6))
    p = 3
    
    cat("N = ",N)
    cat("R = ",R)
    cat("p = ",p)
    
    
    set.seed(1)
    l = function(N){results(N, R, p, conv_thres = 10^(-conv_thres_exp), diagonal, lmer)}
    final = mclapply(rep(N,niter),l)
    saveRDS(final, file = sprintf("randomslopes_after_review_N=%d_R=%d_p=%d_seed=1_niter=%d_conv_thres_exp=%d_diagonal_%d_lmer_%d.RDS", N,R,p,niter,conv_thres_exp, diagonal, lmer))
    return(final)
}



