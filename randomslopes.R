library(dplyr)
library(tibble)
library(AHM)
library(CVXR)
library(expm)
library(parallel)

## This is a function to ignore 0/0 cases
div = function(u,v){
  ifelse(v==0,0,u/v)
}

## This is a summarise function, to summarise by sum across clients and items 

tmat = function(X, f, k = max(f)){
  ## X is nxp matrix, including intercept if appropriate
  ## f is a factor variable of length n with k classes, stored as an integer
  ## we reshape X to be a matrix if X is a vector
  if(length(X)==length(f)){
    X=matrix(X,ncol=1)
  }
  d=as_tibble(cbind(f,X))
  psum = d %>%
    group_by(f) %>%
    summarize_all(sum)
  return(data.matrix(psum[,-1]))
}




## This is method of moments function if covariance matrices are diagonal.
## X_a is covariates for which random clients slopes are used.
## X_b is covariates for which random items slopes are used.
## fa is the vector representing the client, data point corresponds to
## fb is the vector representing the client, data point corresponds to

method_of_moments_diagonal = function(y, X, X_a, X_b, fa, fb){
    n = length(y)
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    p = ncol(X) - 1
    C = length(unique(fb))
    R = length(unique(fa))

    W_a = NULL
    for(i in 1:p_a){
        for(j in 1:p_a){
            W_a = cbind(W_a, X_a[,i]*X_a[,j])
        }
    }
  
    Za = tmat(W_a,fa)
    remove(W_a)
  
    W_b = NULL
    for(i in 1:p_b){
        for(j in 1:p_b){
            W_b = cbind(W_b,X_b[,i]*X_b[,j])
        }
    }
  
    Zb = tmat(W_b,fb)
    remove(W_b)
    
    t_a = seq(1,p_a^2,(p_a + 1))
    t_b = seq(1,p_b^2,(p_b + 1))
    
    sum_xa_2 = colSums(Za[, t_a, drop = FALSE])
    sum_xb_2 = colSums(Zb[, t_b, drop = FALSE])

  
    left = NULL
    for(i in 1:p_a){
        a = c()
        for(i1 in 1:p_a){
            a = c(a, sum_xa_2[i1] - sum(div((Za[,(i-1)*p_a + i1]^2),Za[,(i-1)*p_a + i])))
        }
        left = rbind(left, a) 
    }
    right = NULL
    for(i in 1:p_a){
        a = c()
        for(i1 in 1:p_b){
            a = c(a, sum_xb_2[i1] - sum(div(tmat((X_a[, i] * X_b[, i1])^2,fa),Za[,(i-1)*p_a + i])))
        }
        right = rbind(right, a) 
    }
    
    mat = cbind(left, right, n - R)

    left = NULL
    for(i in 1:p_b){
        a = c()
        for(i1 in 1:p_a){
            a = c(a, sum_xa_2[i1] - sum(div(tmat((X_b[, i] * X_a[, i1])^2,fb),Zb[,(i-1)*p_b + i])))
        }
        left = rbind(left, a) 
    }
  
    right = NULL
    for(i in 1:p_b){
        a = c()
        for(i1 in 1:p_b){
            a = c(a, sum_xb_2[i1] - sum(div((Zb[,(i-1)*p_b + i1]^2),Zb[,(i-1)*p_b + i])))
        }
        right = rbind(right, a) 
    }
    mat = rbind(mat, cbind(left, right, n - C))
    
    last = c()
    for(i1 in 1:p_a){
        last = c(last, sum_xa_2[i1] - sum(Za[, i1]^2) / n)
    }
    for(i1 in 1:p_b){
        last = c(last, sum_xb_2[i1] - sum(Zb[, i1]^2) / n )
    }
    last = c(last, n - 1)
    mat = rbind(mat, last)
    
    beta_ols = lm(y ~ X - 1)$coef
    r_ols = matrix(y - X%*% beta_ols, ncol = 1)
    
    Za_y = tmat(X_a * matrix(rep(r_ols, p_a), ncol = p_a),fa)
    Zb_y = tmat(X_b * matrix(rep(r_ols, p_b), ncol = p_b),fb)
    right_side = sum(r_ols^2) - c(colSums(div(((Za_y)^2),Za[,t_a, drop = FALSE])),colSums(div(((Zb_y)^2), Zb[,t_b, drop = FALSE])),((sum(r_ols))^2)/n)
    cov_par = solve(mat) %*% right_side
    if(p_a == 1){
        Sigma_a = cov_par[1:p_a]
    }else{
        Sigma_a = diag(cov_par[1:p_a])
    }
    
    if(p_b == 1){
        Sigma_b = cov_par[(p_a + 1): (p_a + p_b)]
    }else{
        Sigma_b = diag(cov_par[(p_a + 1): (p_a + p_b)])
    }
    sigma2e = cov_par[p_a + p_b + 1]
    return(enlist(Sigma_a, Sigma_b, sigma2e, Za, Zb))
    
    

}

## This is method of moments function if covariance matrices are non diagonal, we although provide an option to constraint any of them to be diagonal.
## The function involves solving convex optimization problem

method_of_moments_non_diagonal = function(y, X, X_a, X_b, fa, fb, diagonal_Sigma_a = FALSE, diagonal_Sigma_b = FALSE){
    n = length(y)
    p = ncol(X) - 1
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    C = length(unique(fb))
    R = length(unique(fa))
    
    W_a = NULL
    for(i in 1:p_a){
        for(j in 1:p_a){
          W_a = cbind(W_a, X_a[,i]*X_a[,j])
        }
    }
      
    Za = tmat(W_a,fa)
    remove(W_a)
    
    W_b = NULL
    for(i in 1:p_b){
        for(j in 1:p_b){
          W_b = cbind(W_b,X_b[,i]*X_b[,j])
        }
    }
    
    Zb = tmat(W_b,fb)
    remove(W_b)
  
    t_a = seq(1,p_a^2,(p_a + 1))
    t_b = seq(1,p_b^2,(p_b + 1))
    
    sum_xa_2 = colSums(Za[, t_a, drop = FALSE])
    sum_xb_2 = colSums(Zb[, t_b, drop = FALSE])
  
  
    left = NULL
    for(i in 1:p_a){
        a = c()
        for(i1 in 1:p_a){
            for(j1 in i1:p_a){
                if(i1 == j1){
                    a = c(a, sum_xa_2[i1] - sum(div((Za[,(i-1)*p_a + i1]^2),Za[,(i-1)*p_a + i])))
                }else{
                    a = c(a, 2*(sum(Za[,(i1-1)*p_a + j1])-sum(div((Za[,(i-1)*p_a + j1]*Za[,(i-1)*p_a + i1]),(Za[,(i-1)*p_a + i])))))
                }
            }
        }
        left = rbind(left, a) 
    }
    right = NULL
    for(i in 1:p_a){
        a = c()
        for(i1 in 1:p_b){
            for(j1 in i1:p_b){
                if(i1 == j1){
                    a = c(a, sum_xb_2[i1] - sum(div(tmat((X_a[, i] * X_b[, i1])^2,fa),Za[,(i-1)*p_a + i])))
                }else{
                    a = c(a,2*(sum(Zb[,(i1-1)*p_b + j1])-sum(div(tmat(((X_a[,i]^2) * X_b[, i1] * X_b[, j1]),fa),(Za[,(i-1)*p_a + i])))))
                }
            }
        }
        right = rbind(right, a) 
    }
    
    mat = cbind(left, right, n - R)
  
    left = NULL
    for(i in 1:p_b){
        a = c()
        for(i1 in 1:p_a){
            for(j1 in i1:p_a){
                if(i1 == j1){
                    a = c(a, sum_xa_2[i1] - sum(div(tmat((X_b[, i] * X_a[, i1])^2,fb),Zb[,(i-1)*p_b + i])))
                }else{
                    a = c(a,2*(sum(Za[,(i1-1)*p_a + j1])-sum(div(tmat(((X_b[,i]^2) * X_a[, i1] * X_a[, j1]),fb),(Zb[,(i-1)*p_b + i])))))
                }
            }
        }
        left = rbind(left, a) 
    }

    right = NULL
    for(i in 1:p_b){
        a = c()
        for(i1 in 1:p_b){
            for(j1 in i1:p_b){
                if(i1 == j1){
                    a = c(a, sum_xb_2[i1] - sum(div((Zb[,(i-1)*p_b + i1]^2),Zb[,(i-1)*p_b + i])))
                }else{
                    a = c(a,2*(sum(Zb[,(i1-1)*p_b + j1])-sum(div((Zb[,(i-1)*p_b + j1]*Zb[,(i-1)*p_b + i1]),(Zb[,(i-1)*p_b + i])))))
                }
            }        
        }
        right = rbind(right, a) 
    }
    mat = rbind(mat, cbind(left, right, n - C))
    
    last = c()
    for(i1 in 1:p_a){
        for(j1 in i1:p_a){
            if(i1 == j1){
                last = c(last, sum_xa_2[i1] - sum(Za[, i1]^2) / n)
            }else{
                last = c(last, 2*sum(Za[,(i1-1)*p_a + j1]) - 2*sum((Za[,(i1)]*Za[,(j1)])) / n)
            }
        }
    }
    for(i1 in 1:p_b){
        for(j1 in i1:p_b){
            if(i1 == j1){
                last = c(last, sum_xb_2[i1] - sum(Zb[, i1]^2) / n )
            }else{
                last = c(last, 2*sum(Zb[,(i1-1)*p_b + j1]) - 2*sum((Zb[,(i1)]*Zb[,(j1)])) / n)
            }
        }
    }
    last = c(last, n - 1)
    mat = rbind(mat, last)
  
    beta_ols = lm(y ~ X - 1)$coef
    r_ols = matrix(y - X%*% beta_ols, ncol = 1)
    
    Za_y = tmat(X_a * matrix(rep(r_ols, p_a), ncol = p_a),fa)
    Zb_y = tmat(X_b * matrix(rep(r_ols, p_b), ncol = p_b),fb)

    right_side = sum(r_ols^2) - c(colSums(div(((Za_y)^2),Za[,t_a, drop = FALSE])),colSums(div(((Zb_y)^2), Zb[,t_b, drop = FALSE])),((sum(r_ols))^2)/n)
    m = mean(mat)
    mat = mat / m
    right_side = right_side / m
    
    ## removing non - required variables
    rm(list = c('Za_y', 'Zb_y', 'left', 'right', 'a', 'last'))
  
    sequence_a = c()
    sequence_a_diag = c()
    for(i in 1:p_a){
        q = ((i-1)*p_a)+(i:p_a)
        sequence_a = c(sequence_a, q)
        if(diagonal_Sigma_a & i < p_a){
            sequence_a_diag = c(sequence_a_diag, ((i-1)*p_a)+((i+1):p_a))
        }
    }
    
    sequence_b = c()
    sequence_b_diag = c()
    for(i in 1:p_b){
        q = ((i-1)*p_b)+(i:p_b)
        sequence_b = c(sequence_b,q)
        if(diagonal_Sigma_b & i < p_b){
            sequence_b_diag = c(sequence_b_diag, ((i-1)*p_b)+((i+1):p_b))
        }
    }
    
    Sigma_a=Variable(p_a, p_a, PSD=TRUE)
    Sigma_b=Variable(p_b, p_b, PSD=TRUE)
    sigma2e=Variable(1)
  
    constraints = list(Sigma_a%>>%0, Sigma_a == t(Sigma_a), Sigma_b%>>%0, sigma2e>=0, Sigma_b == t(Sigma_b), 
                   mat[,(1:(p_a*(p_a + 1)/2))]%*%(vec(Sigma_a))[sequence_a]
                   + mat[,(1+((p_a)*(p_a + 1)/2)):(((p_a)*(p_a + 1))/ 2 + ((p_b)*(p_b + 1))/ 2)]%*%(vec(Sigma_b))[sequence_b]
                   + sigma2e*mat[,1 + (((p_a)*(p_a + 1))/ 2 + ((p_b)*(p_b + 1))/ 2)] == right_side
    )
    if(diagonal_Sigma_a){
        constraints = append(constraints, vec(Sigma_a)[sequence_a_diag] == 0 )
    }
    if(diagonal_Sigma_b){
        constraints = append(constraints, vec(Sigma_b)[sequence_b_diag] == 0 )
    }
    objective = 0
    prob = Problem(Minimize(objective),constraints)
    result = solve(prob, solver = 'MOSEK')
    print(result$status)
    
    if(result$status == 'optimal'){
        return(list(result_status = result$status, Sigma_a = result$getValue(Sigma_a), Sigma_b = result$getValue(Sigma_b), sigma2e = result$getValue(sigma2e), Za = Za, Zb = Zb))
    }else{
    print("Method of moments approach failed")
        return(list(result_status = 'Failed', Za = Za, Zb = Zb))
    
    }
  
}



## This is vanilla backfitting function to estimate fixed effects. It requires an initial estimate of covariance parameters.

vanilla_backfitting = function(y, X, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, Za, Zb, max_iter = 600, conv_thres = 10^(-6)){
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    C = length(unique(fb))
    R = length(unique(fa))
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    
    xtx_lam_inverse = function(i){
        return(solve(matrix(Za[i,],p_a, p_a)+(lam)))
    }
    
    xtx_alpha_inverse = function(i){
        return(solve(matrix(Zb[i,],p_b, p_b)+(alpha)))
    }
    
    projection_a = mclapply(1:R, xtx_lam_inverse)
    
    projection_b = mclapply(1:C, xtx_alpha_inverse)
    
    rslope_a = function(r){
        Rvec = X_a * matrix(rep(r,p_a),ncol = p_a)
        
        X_R = tmat(Rvec,fa)
        
        pr = function(i){
            projection_a[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:R),nrow = p_a)
        fit = rowSums(X_a*(t(beta[,fa, drop = FALSE])))
        enlist(beta,fit)
    }
    
    rslope_b = function(r){
        Rvec = X_b * matrix(rep(r,p_b),ncol = p_b)
        
        X_R = tmat(Rvec,fb)
        
        pr = function(i){
            projection_b[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:C),nrow = p_b)
        fit = rowSums(X_b*(t(beta[,fb, drop = FALSE])))
        enlist(beta,fit)
    }
    
    
    fit0 = lm(y~X-1)
    beta = coef(fit0)
    fvo = fitted(fit0)
    faxo = fbxo = fvo*0
    diff = 1
    iter = 1
    while(iter <= max_iter & diff >= conv_thres){
        r <- y- fvo - fbxo
        fita <- rslope_a(r)
        fax <- fita$fit
        a = fita$beta
        r <- y- fvo - fax
        fitb <- rslope_b(r)
        fbx <- fitb$fit
        b = fitb$beta
        r <- y - fax - fbx
        fit0 <- lm(r~X-1)
        fv <- fitted(fit0)
        n0 <- norm(fv-fvo,type="2")/norm(fv,type="2")
        na <- norm(fax-faxo,type="2")/norm(fax,type="2")
        nb <- norm(fbx-fbxo,type="2")/norm(fbx,type="2")
        diff = max(n0,na,nb)
        cat(iter,"max change:",diff,"\n")
        fvo = fv; faxo=fax; fbxo=fbx
        iter = iter + 1
        if(iter > 50 & diff > 0.1){
            print(R, C, n, p, p_a, p_b)
        }
    }
    return(list(beta = coef(fit0), niter = iter, A = a, B = b))
   
}

## This is variational backfitting function (based on variational EM) to estimate fixed effects as well as update covariance parameters.
## It requires an initial estimate of covariance parameters. If covariance matrices are constrained to be diagonal, those constraints could be fed into the function.

variational_backfitting = function(y, X, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e,Za, Zb, max_iter = 600, conv_thres = 10^(-6), covariances_diag_Sigma_a = FALSE, covariances_diag_Sigma_b = FALSE){
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    C = length(unique(fb))
    R = length(unique(fa))
    n = length(y)
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    
    xtx_lam_inverse = function(i, lam){
        return(solve(matrix(Za[i,],p_a, p_a)+(lam)))
    }
    
    xtx_alpha_inverse = function(i, alpha){
        return(solve(matrix(Zb[i,],p_b, p_b)+(alpha)))
    }
    

    rslope_a = function(r, xtx_lam_inv){
        Rvec = X_a * matrix(rep(r,p_a),ncol = p_a)
        
        X_R = tmat(Rvec,fa)
        
        pr = function(i){
            xtx_lam_inv[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:R),nrow = p_a)
        fit = rowSums(X_a*(t(beta[,fa, drop = FALSE])))
        enlist(beta,fit)
    }
    
    rslope_b = function(r, xtx_alpha_inv){
        Rvec = X_b * matrix(rep(r,p_b),ncol = p_b)
        
        X_R = tmat(Rvec,fb)
        
        pr = function(i){
            xtx_alpha_inv[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:C),nrow = p_b)
        fit = rowSums(X_b*(t(beta[,fb, drop = FALSE])))
        enlist(beta,fit)
    }
    
    
    fit0 = lm(y~X-1)
    beta = coef(fit0)
    fvo = fitted(fit0)
    faxo = fbxo = fvo*0
    diff = 1
    iter = 1
    while(iter <= max_iter & diff >= conv_thres){
        r <- y - fvo - fbxo
        xtx_lam_inv = mclapply(1:R,function(i){xtx_lam_inverse(i,lam)})
        xtx_alpha_inv = mclapply(1:C,function(i){xtx_alpha_inverse(i,alpha)})
        fita <- rslope_a(r, xtx_lam_inv)
        fax <- fita$fit
        a = fita$beta
        Sigma_a = mclapply(1:R,function(i){(matrix(a[,i],ncol=1))%*%(matrix(a[,i],nrow=1))})
        Sigma_a = (Reduce('+', Sigma_a)+sigma2e*Reduce('+',xtx_lam_inv))/R
        if(p_a > 1 & covariances_diag_Sigma_a){
            Sigma_a = diag(diag(Sigma_a))
        }

        r <- y- fvo - fax
        fitb <- rslope_b(r, xtx_alpha_inv)
        fbx <- fitb$fit
        b = fitb$beta
        Sigma_b = mclapply(1:C,function(i){(matrix(b[,i],ncol=1))%*%(matrix(b[,i],nrow=1))})
        Sigma_b = (Reduce('+', Sigma_b)+sigma2e*Reduce('+',xtx_alpha_inv))/C
        if(p_b > 1 & covariances_diag_Sigma_b){
            Sigma_b = diag(diag(Sigma_b))
        }

        r = y - fax - fbx 
        fit0 <- lm(r~X-1)
        fv <- fitted(fit0)
        r = y - fax - fbx - fv
        
        a_comp = function(i){
            return(sum(diag(matrix(Za[i, ], ncol = p_a, nrow = p_a) %*% xtx_lam_inv[[i]])))
        }
        b_comp = function(i){
            return(sum(diag(matrix(Zb[i, ], ncol = p_b, nrow = p_b) %*% xtx_alpha_inv[[i]])))
        }
        sigma2e = (sigma2e * (sum(mcmapply(a_comp, 1:R)) + sum(mcmapply(b_comp, 1:C))) + sum(r^2))/ n
        
        
        n0 <- norm(fv-fvo,type="2")/norm(fv,type="2")
        na <- norm(fax-faxo,type="2")/norm(fax,type="2")
        nb <- norm(fbx-fbxo,type="2")/norm(fbx,type="2")
        diff = max(n0,na,nb)
        cat(iter,"max change:",diff,"\n")
        fvo = fv; faxo = fax; fbxo = fbx
        iter = iter + 1
        lam = sigma2e*solve(Sigma_a)
        alpha = sigma2e*solve(Sigma_b)
    }
    return(list(beta = coef(fit0), niter = iter, A = a, B = b, Sigma_a = Sigma_a, Sigma_b = Sigma_b, sigma2e = sigma2e))
   
}


clubbed_backfitting = function(y, X, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, Za, Zb, max_iter = 600, conv_thres = 10^(-6)){
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    n = length(y)
    p = ncol(X) - 1
    C = length(unique(fb))
    R = length(unique(fa))
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    xtx_lam_inverse = function(i){
        return(solve(matrix(Za[i,],p_a, p_a)+(lam)))
    }
    
    xtx_alpha_inverse = function(i){
        return(solve(matrix(Zb[i,],p_b, p_b)+(alpha)))
    }
    
    projection_a = mclapply(1:R, xtx_lam_inverse)
    
    projection_b = mclapply(1:C, xtx_alpha_inverse)
    
    rslope_a = function(r){
        Rvec = X_a * matrix(rep(r,p_a),ncol = p_a)
        
        X_R = tmat(Rvec,fa)
        
        pr = function(i){
        projection_a[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:R),nrow = p_a)
        fit = rowSums(X_a*(t(beta[,fa, drop = FALSE])))
        enlist(beta,fit)
    }
    
    rslope_b = function(r){
        Rvec = X_b * matrix(rep(r,p_b),ncol = p_b)
        
        X_R = tmat(Rvec,fb)
        
        pr = function(i){
        projection_b[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:C),nrow = p_b)
        fit = rowSums(X_b*(t(beta[,fb, drop = FALSE])))
        enlist(beta,fit)
    }
    
    club_a = function(r){
        Xresid = X
        for(i in 1:ncol(X)){
            Xresid[,i] = X[,i]-(rslope_a(X[,i]))$fit
        }
        beta = solve(t(X)%*%Xresid,t(Xresid)%*%r)
        fv = X%*%beta
        fita = rslope_a(r-fv)
        betaa = fita$beta
        fax = fita$fit
        enlist(beta = beta, betaa = betaa, fax = fax)
    }
    
    club_b = function(r){
        Xresid = X
        for(i in 1:ncol(X)){
            Xresid[,i] = X[,i]-(rslope_b(X[,i]))$fit
        }
        beta = solve(t(X)%*%Xresid,t(Xresid)%*%r)
        fv = X%*%beta
        fitb = rslope_b(r-fv)
        betab = fitb$beta
        fbx = fitb$fit
        enlist(beta = beta, betab = betab, fbx = fbx)
    }
    
    fvo = faxo =  fbxo = y*0
    diff = 1
    iter = 1
    while(iter <= max_iter & diff>= conv_thres){
        r = y - fbxo
        fitxa = club_a(r)
        fax = fitxa$fax
        r = y - fax
        fitxb = club_b(r)
        fbx = fitxb$fbx
        fv = X%*%fitxb$beta
        n0 = norm(fv-fvo,type="2")/norm(fv,type="2")
        na = norm(fax-faxo,type="2")/norm(fax,type="2")
        nb = norm(fbx-fbxo,type="2")/norm(fbx,type="2")
        diff = max(n0,na,nb)
        cat(iter, "max change:",diff,"\n")
        fvo = fv; faxo = fax; fbxo = fbx
        iter = iter + 1
    }
    return(list(beta = fitxb$beta, niter = iter, A = fitxa$betaa, B = fitxb$betab))
   
}
clubbed_variational_backfitting = function(y, X, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, Za, Zb, max_iter = 600, conv_thres = 10^(-6), covariances_diag_Sigma_a = FALSE, covariances_diag_Sigma_b = FALSE){
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    p = ncol(X) - 1
    n = length(y)
    C = length(unique(fb))
    R = length(unique(fa))
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    
    xtx_lam_inverse = function(i, lam){
        return(solve(matrix(Za[i,],p_a, p_a)+(lam)))
    }
    
    xtx_alpha_inverse = function(i, alpha){
        return(solve(matrix(Zb[i,],p_b, p_b)+(alpha)))
    }
    

    rslope_a = function(r, xtx_lam_inv){
        Rvec = X_a * matrix(rep(r,p_a),ncol = p_a)
        
        X_R = tmat(Rvec,fa)
        
        pr = function(i){
            xtx_lam_inv[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:R),nrow = p_a)
        fit = rowSums(X_a*(t(beta[,fa, drop = FALSE])))
        enlist(beta,fit)
    }
    
    rslope_b = function(r, xtx_alpha_inv){
        Rvec = X_b * matrix(rep(r,p_b),ncol = p_b)
        
        X_R = tmat(Rvec,fb)
        
        pr = function(i){
            xtx_alpha_inv[[i]]%*%X_R[i,] 
        }
        beta = matrix(mcmapply(pr, 1:C),nrow = p_b)
        fit = rowSums(X_b*(t(beta[,fb, drop = FALSE])))
        enlist(beta,fit)
    }

    club_a = function(r,xtx_lam_inv){
        Xresid = X
        for(i in 1:ncol(X)){
            Xresid[,i] = X[,i]-(rslope_a(X[,i], xtx_lam_inv))$fit
        }
        beta = solve(t(X)%*%Xresid,t(Xresid)%*%r)
        fv = X%*%beta
        fita = rslope_a(r - fv,xtx_lam_inv)
        betaa = fita$beta
        fax = fita$fit
        list(beta = beta,betaa = betaa, fax = fax)
    }
    
    club_b = function(r,xtx_alpha_inv){
        Xresid = X
        for(i in 1:ncol(X)){
            Xresid[,i]=X[,i]-(rslope_b(X[,i], xtx_alpha_inv))$fit
        }
        beta = solve(t(X)%*%Xresid,t(Xresid)%*%r)
        fv = X%*%beta
        fitb = rslope_b(r - fv,xtx_alpha_inv)
        betab = fitb$beta
        fbx = fitb$fit
        list(beta = beta, betab = betab, fbx = fbx)
    }

  
    fvo<-faxo <- fbxo <- y*0
    iter = 1
    diff = 1
    Changes_sigma2e = c()
    Changes_Sigma_a = c()
    Changes_Sigma_b = c()
    
    if(is.null(sigma2e)){
        sigma2e = 1
    }
    if(is.null(Sigma_a)){
        Sigma_a = diag(ncol(X_a))
    }
    
    if(is.null(Sigma_b)){
      Sigma_b = diag(ncol(X_b))
    }
    old_sigma2e = sigma2e
    old_Sigma_a = Sigma_a
    old_Sigma_b = Sigma_b
    
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    
    training_error = c()
    
    while(iter <= max_iter && diff >= conv_thres){
        xtx_lam_inv = lapply(1:R,function(i){xtx_lam_inverse(i,lam)})
        xtx_alpha_inv = lapply(1:C,function(i){xtx_alpha_inverse(i,alpha)})  
        r = y - fbxo
        fitxa = club_a(r, xtx_lam_inv)
        fax = fitxa$fax
        a = fitxa$betaa
        Sigma_a = mclapply(1 : R,function(i){(matrix(a[,i],ncol=1))%*%(matrix(a[,i],nrow=1))})
        new_Sigma_a = (Reduce('+', Sigma_a) + old_sigma2e*Reduce('+',xtx_lam_inv))/R
        if(p_a > 1 & covariances_diag_Sigma_a){
            new_Sigma_a = diag(diag(new_Sigma_a))
        }
        r = y - fax
        fitxb = club_b(r, xtx_alpha_inv)
        fbx = fitxb$fbx
        b = fitxb$betab
        Sigma_b = mclapply(1 : C,function(i){(matrix(b[,i],ncol=1))%*%(matrix(b[,i],nrow=1))})
        new_Sigma_b = (Reduce('+', Sigma_b) + old_sigma2e*Reduce('+',xtx_alpha_inv))/C
        if(p_b > 1 & covariances_diag_Sigma_b){
            new_Sigma_b = diag(diag(new_Sigma_b))
        }
        
        fv = X %*% fitxb$beta
        n0 =  norm(fv - fvo,type="2")/norm(fv,type="2")
        na = norm(fax - faxo,type="2")/norm(fax,type="2")
        nb = norm(fbx - fbxo,type="2")/norm(fbx,type="2")
        diff = max(n0,na,nb)
        fvo = fv; faxo = fax; fbxo = fbx
        r = y - fv - fax - fbx
        a_comp = function(i){
            return(sum(diag(matrix(Za[i, ], ncol = p_a, nrow = p_a) %*% xtx_lam_inv[[i]])))
        }
        b_comp = function(i){
            return(sum(diag(matrix(Zb[i, ], ncol = p_b, nrow = p_b) %*% xtx_alpha_inv[[i]])))
        }
        new_sigma2e = (old_sigma2e * (sum(mcmapply(a_comp, 1:R)) + sum(mcmapply(b_comp, 1:C))) + sum(r^2))/ n
        
        lam = new_sigma2e*solve(new_Sigma_a)
        alpha = new_sigma2e*solve(new_Sigma_b)
        
        change_sigma2e = abs(new_sigma2e/old_sigma2e - 1)
        change_Sigma_a = norm(new_Sigma_a - old_Sigma_a, type = "2")
        change_Sigma_b = norm(new_Sigma_b - old_Sigma_b, type = "2")
        Changes_sigma2e = c(Changes_sigma2e, change_sigma2e)
        Changes_Sigma_a = c(Changes_Sigma_a, change_Sigma_a)
        Changes_Sigma_b = c(Changes_Sigma_b, change_Sigma_b)
        old_sigma2e = new_sigma2e
        old_Sigma_a = new_Sigma_a
        old_Sigma_b = new_Sigma_b
        training_error = c(training_error, mean(r^2))
        cat(iter,"max change:",diff,"\n")
        iter = iter + 1
    }
    return(list(beta = fitxb$beta, Sigma_a = new_Sigma_a, Sigma_b = new_Sigma_b, sigma2e = new_sigma2e, niter = iter, A = a, B = b, training_error = training_error, Changes_sigma2e = Changes_sigma2e, Changes_Sigma_a = Changes_Sigma_a, Changes_Sigma_b = Changes_Sigma_b))
  
  
}



fit_on_random_effects = function(y, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, max_iter = 600, conv_thres = 10^(-6)){
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    C = length(unique(fb))
    R = length(unique(fa))
    
    lam = sigma2e * solve(Sigma_a)
    alpha = sigma2e * solve(Sigma_b)
    
    W_a = NULL
    for(i in 1:p_a){
        for(j in 1:p_a){
            W_a = cbind(W_a,X_a[,i]*X_a[,j])
        }
    }
    W_b = NULL
    for(i in 1:p_b){
        for(j in 1:p_b){
            W_b = cbind(W_b,X_b[,i]*X_b[,j])
        }
    }
    Za = tmat(W_a, fa)
    Zb = tmat(W_b, fb)
    xtx_lam_inverse = function(i){
        return(solve(matrix(Za[i,],p_a, p_a) + lam))
    }
    
    xtx_alpha_inverse = function(i){
        return(solve(matrix(Zb[i,],p_b, p_b) + alpha))
    }
    
    Ma = mclapply(1:R, xtx_lam_inverse)
    
    Mb = mclapply(1:C, xtx_alpha_inverse)

    rslope_a = function(r){
        Rvec = matrix(rep(r,p_a),ncol=p_a)
        Rvec = X_a*Rvec
        X_R = tmat(Rvec,fa)
        l = function(i){
            Ma[[i]]%*%X_R[i,] 
        }
      
        beta=matrix(mcmapply(l, 1:R),nrow = p_a)
        fit=rowSums(X_a*(t(beta[,fa, drop = FALSE])))
        enlist(beta,fit)
    }
    
    rslope_b = function(r){
        Rvec = matrix(rep(r,p_b),ncol=p_b)
        Rvec = X_b*Rvec
        X_R = tmat(Rvec,fb)
        l = function(i){
            Mb[[i]]%*%X_R[i,] 
        }
      
        beta=matrix(mcmapply(l, 1:C),nrow = p_b)
        fit=rowSums(X_b*(t(beta[,fb, drop = FALSE])))
        enlist(beta,fit)
    }
    
    faxo <- fbxo <- 0
    diff = 1
    iter = 1
    while(iter <= max_iter && diff >= conv_thres){
        r <- y - fbxo
        fita <- rslope_a(r)
        fax <- fita$fit
        r <- y - fax
        fitb <- rslope_b(r)
        fbx <- fitb$fit
        na <- norm(fax-faxo,type="2")/norm(fax,type="2")
        nb <- norm(fbx-fbxo,type="2")/norm(fbx,type="2")
        diff=max(na,nb)
        cat(iter,"max change:",diff,"\n")
        faxo=fax;fbxo=fbx
        iter = iter + 1
    }
    return(list(betaa=fita$beta,betab=fitb$beta,fax=fax,fbx=fbx))
    

}


covariances_beta_hat = function(y, X, X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, X_residual = NULL){
    
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    q = ncol(X)
    
    if(is.null(X_residual)){
        X_residual = X
        for(k in 1:q){
            fit_random_effects = fit_on_random_effects(X[, k], X_a, X_b, fa, fb, Sigma_a, Sigma_b, sigma2e, max_iter = 100, conv_thres = 10^(-6))
            X_residual[, k] = X[, k] - fit_random_effects$fax - fit_random_effects$fbx
    
        }
    
    }

    
    X_a_Sigma_a_half = X_a %*% sqrtm(Sigma_a)
    X_b_Sigma_b_half = X_b %*% sqrtm(Sigma_b)

    X_t_X_res_inverse = solve(t(X)%*%X_residual)
    
    X_res_X_a_Sigma_a_half = NULL
    for(i in 1:q){
        for(j in 1:p_a){
            X_res_X_a_Sigma_a_half = cbind(X_res_X_a_Sigma_a_half, X_residual[,i] * X_a_Sigma_a_half[,j])  
        }
    }
    X_res_X_a_Sigma_a_half_sum = tmat(X_res_X_a_Sigma_a_half,fa)
    
    final_a = matrix(0,q,q)
    for(i in 1:q){
        for(j in i:q){
            final_a[i,j] = sum(X_res_X_a_Sigma_a_half_sum[,((i-1)*p_a+1):(i*p_a)]*X_res_X_a_Sigma_a_half_sum[,((j-1)*p_a+1):(j*p_a)])
            final_a[j,i] = final_a[i,j]
        }
    }
    
    X_res_X_b_Sigma_b_half=NULL
    
    for(i in 1:q){
        for(j in 1:p_b){
            X_res_X_b_Sigma_b_half = cbind(X_res_X_b_Sigma_b_half, X_residual[,i]*X_b_Sigma_b_half[,j])  
        }
    }
    X_res_X_b_Sigma_b_half_sum = tmat(X_res_X_b_Sigma_b_half,fb)
    
    final_b = matrix(0,q,q)
    for(i in 1:q){
        for(j in i:q){
            final_b[i,j] = sum(X_res_X_b_Sigma_b_half_sum[,((i-1)*p_b+1):(i*p_b)]*X_res_X_b_Sigma_b_half_sum[,((j-1)*p_b+1):(j*p_b)])
            final_b[j,i] = final_b[i,j]
        }
    }
    
    ## Crosscheck if X_t_X_res_inverse is symmetric
    covariance_gls_beta_hat_gls = X_t_X_res_inverse%*%(final_a+final_b+sigma2e*t(X_residual)%*%X_residual)%*%X_t_X_res_inverse
    return(list(covariance_gls_beta_hat_gls = covariance_gls_beta_hat_gls, X_residual = X_residual))

}





scalable_crossed_random = function(y, X, X_a, X_b, fa, fb, Sigma_a_diagonal = FALSE, Sigma_b_diagonal = FALSE, variational = TRUE, clubbing = TRUE, method_of_mom = TRUE, beta_covariance = FALSE, max_iter = 600, conv_thres = 10^(-6)){

    n = length(y)
    
    range_f1_train=sort(unique(fa))
    range_f2_train=sort(unique(fb))
    
    R = length(range_f1_train)
    C = length(range_f2_train)
    
    
    bij_f1=list()
    bij_f2=list()
    
    for(i in 1:R){
        bij_f1[[range_f1_train[i]]] = i
    
    }
    for(i in 1:C){
        bij_f2[[range_f2_train[i]]] = i
    
    }
    
    
    
    for(i in 1:n){
        fa[i]=bij_f1[[fa[i]]]
    }
    
    for(i in 1:n){
        fb[i]=bij_f2[[fb[i]]]
    }
    p = ncol(X) - 1
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    
    if((Sigma_a_diagonal & Sigma_b_diagonal) | (p_a == 1 & p_b == 1)){
        mom = method_of_moments_diagonal(y, X, X_a, X_b, fa, fb)
        if(all(diag(mom$Sigma_a) > 0) &  all(diag(mom$Sigma_b) > 0) & mom$sigma2e > 0 & method_of_mom){
            print("Method of moments approach worked")
            Za = mom$Za
            Zb = mom$Zb
            mm_Sigma_a <- mom$Sigma_a
            mm_Sigma_b <- mom$Sigma_b
            mm_sigma2e <- mom$sigma2e
            mom_status <- 1
            
        }else{
            print("Method of moments approach failed")
            mm_Sigma_a = lam = diag(p_a)
            mm_Sigma_b = alpha = diag(p_b)
            mm_sigma2e = 1
            Za = mom$Za
            Zb = mom$Zb
            mom_status <- 0
        }
        
    }else{
        mom <- method_of_moments_non_diagonal(y, X, X_a, X_b, fa, fb, diagonal_Sigma_a = Sigma_a_diagonal, diagonal_Sigma_b = Sigma_b_diagonal)
        saveRDS(mom, "mom_results_stitch_fix.RDS")
        if(mom$result_status == 'optimal' & method_of_mom){
            print("Method of moments approach worked")
            mom_status <- 1
            mm_Sigma_a <- mom$Sigma_a
            mm_Sigma_b <- mom$Sigma_b
            mm_sigma2e <- mom$sigma2e
            
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
            

        }
        
    }
    
    if(variational | mom_status == 0){
        if(clubbing){
            fit <- clubbed_variational_backfitting(y, X, X_a, X_b, fa, fb, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, covariances_diag_Sigma_a = Sigma_a_diagonal, covariances_diag_Sigma_b = Sigma_b_diagonal, conv_thres = conv_thres, max_iter = max_iter)
        }
        else{
            fit <- variational_backfitting(y, X, X_a, X_b, fa, fb, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, covariances_diag_Sigma_a = Sigma_a_diagonal, covariances_diag_Sigma_b = Sigma_b_diagonal, conv_thres = conv_thres, max_iter = max_iter)
        }
        results = list(fixed_effects = fit$beta, cov_fa = fit$Sigma_a, cov_fb = fit$Sigma_b, sigma = sqrt(fit$sigma2e), random_effects_fa = fit$A, random_effects_fb = fit$B)
    }else{
        if(clubbing){
            fit <- clubbed_backfitting(y, X, X_a, X_b, fa, fb, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, conv_thres = conv_thres, max_iter = max_iter)
            
        }else{
            fit <- vanilla_backfitting(y, X, X_a, X_b, fa, fb, mm_Sigma_a, mm_Sigma_b, mm_sigma2e, Za, Zb, conv_thres = conv_thres, max_iter = max_iter)
        }
        
        results = list(fixed_effects = fit$beta, cov_fa = mm_Sigma_a, cov_fb = mm_Sigma_b, sigma = sqrt(mm_sigma2e), random_effects_fa = fit$A, random_effects_fb = fit$B)
        
    }
    
    if(beta_covariance){
        covariance_beta_hat = covariances_beta_hat(y, X, X_a, X_b, fa, fb, results$cov_fa, results$cov_fb, results$sigma^2, X_residual = NULL)
        results$covariance_beta_hat = covariance_beta_hat$covariance_gls_beta_hat_gls
        
    }
    return(results)
}




naivety_inefficiency_wrt_ols = function(y, X, X_a, X_b, fa, fb, fit, X_residual = NULL){

    q = ncol(X)
    p_a = ncol(X_a)
    p_b = ncol(X_b)
    fit_ols = summary(lm(y~X-1))
    ols_sigma2e = (fit_ols$sigma)^2
    
    X_t_X_inverse = solve(t(X) %*% X)

    if(is.null(X_residual)){
        X_residual = X
        for(k in 1:q){
            fit_random_effects = fit_on_random_effects(X[, k], X_a, X_b, fa, fb, fit$cov_fa, fit$cov_fb, fit$sigma^2, max_iter = 800, conv_thres = 10^(-6))
            X_residual[, k] = X[, k] - fit_random_effects$fax - fit_random_effects$fbx
    
        }
        saveRDS(X_residual, sprintf("/home/groups/hastie/disha/x_residual_p_a_%d_p_b_%d.RDS", p_a, p_b))
    }
    
    if(is.null(fit$covariance_beta_hat)){
        cov = covariances_beta_hat(y, X, X_a, X_b, fa, fb, fit$cov_fa, fit$cov_fb, fit$sigma^2, X_residual = X_residual)
        covariance_gls_beta_hat_gls = cov$covariance_gls_beta_hat_gls
    }else{
        covariance_gls_beta_hat_gls = fit$covariance_beta_hat
    }
    
    X_t_X_res_inverse = solve(t(X)%*%X_residual)
    covariance_ols_beta_hat_gls = (ols_sigma2e*X_t_X_res_inverse)%*%t(X_residual)%*%X_residual%*%X_t_X_res_inverse
    
    X_a_Sigma_a_half = X_a%*%sqrtm(fit$cov_fa)
    X_b_Sigma_b_half = X_b%*%sqrtm(fit$cov_fb)
    
    X_X_a_Sigma_a_half = NULL
    for(i in 1:q){
        for(j in 1:p_a){
            X_X_a_Sigma_a_half = cbind(X_X_a_Sigma_a_half,X[,i]*X_a_Sigma_a_half[,j])  
        }
    }
    X_X_a_Sigma_a_half_sum = tmat(X_X_a_Sigma_a_half, fa)
    
    final_a = matrix(0,q,q)
    for(i in 1:q){
        for(j in i:q){
            final_a[i,j] = sum(X_X_a_Sigma_a_half_sum[,(((i-1)*p_a)+1):(i*p_a)]*(X_X_a_Sigma_a_half_sum[,(((j-1)*p_a)+1):(j*p_a)]))
            final_a[j,i] = final_a[i,j]
        }
    }
    
    X_X_b_Sigma_b_half = NULL
    for(i in 1:q){
        for(j in 1:p_b){
        X_X_b_Sigma_b_half=cbind(X_X_b_Sigma_b_half,X[,i]*X_b_Sigma_b_half[,j])  
        }
    }
    X_X_b_Sigma_b_half_sum = tmat(X_X_b_Sigma_b_half, fb)
    final_b = matrix(0,q,q)
    for(i in 1:q){
        for(j in i:q){
            final_b[i,j] = sum(X_X_b_Sigma_b_half_sum[,(((i-1)*p_b)+1):(i*p_b)]*(X_X_b_Sigma_b_half_sum[,(((j-1)*p_b)+1):(j*p_b)]))
            final_b[j,i] = final_b[i,j]
        }
    }
    covariance_gls_beta_hat_ols = X_t_X_inverse %*% (final_a+final_b+(fit$sigma^2)*t(X)%*%X)%*%X_t_X_inverse
    covariance_ols_beta_hat_ols = ols_sigma2e*X_t_X_inverse
    
    naivety = diag(covariance_gls_beta_hat_ols)/diag(covariance_ols_beta_hat_ols)
    inefficiency = diag(covariance_gls_beta_hat_ols)/diag(covariance_gls_beta_hat_gls)
    
    return(list(X_residual = X_residual, covariance_ols_beta_hat_ols = covariance_ols_beta_hat_ols, covariance_gls_beta_hat_ols = covariance_gls_beta_hat_ols, 
        covariance_ols_beta_hat_gls = covariance_ols_beta_hat_gls, covariance_gls_beta_hat_gls = covariance_gls_beta_hat_gls,
            naivety = naivety, inefficiency = inefficiency))

    
}

