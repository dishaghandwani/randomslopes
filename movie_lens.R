library(softImpute)
source("randomslopes.R")
## Read movie lens 100k dataset from kaggle
user_info = read.table("/home/users/disha123/factor_analysis/movie_lens/u.user", header = FALSE, sep = "|", quote = "", stringsAsFactors = FALSE)
user_info = data.frame(user_info)
colnames(user_info) = c("user_id","age","gender","occupation","zip code")
user_info$gender = ifelse(user_info$gender == 'M', 1, 0)

movie_info = read.table("/home/users/disha123/factor_analysis/movie_lens/u.item", header = FALSE, sep = "|", quote = "", stringsAsFactors = FALSE)
movie_info = movie_info[,-4]
movie_info = data.frame(movie_info)
colnames(movie_info) = c("movie_id","movie_title","video_release_date","IMDb_URL","unknown","Action","Adventure","Animation","Children's","Comedy","Crime","Documentary","Drama","Fantasy","Film-Noir","Horror","Musical","Mystery","Romance","Sci-Fi","Thriller","War","Western")




results = function(t,xa, xb, max_iter = 150,c = 6, k = 1){
    
    train_data = read.table(sprintf("/home/users/disha123/factor_analysis/movie_lens/u%d.base",t), header = FALSE)
    test_data = read.table(sprintf("/home/users/disha123/factor_analysis/movie_lens/u%d.test",t), header = FALSE)
    train_data = data.frame(train_data)
    colnames(train_data) = c("user_id","item_id","rating","timestamp")
    test_data = data.frame(test_data)
    colnames(test_data) = c("user_id","item_id","rating","timestamp")
    
    y = train_data$rating
    train_movie = movie_info[as.numeric(train_data$item_id), ]
    train_movie = train_movie[,6:23]
    train_user = user_info[train_data$user_id,2:3]
    test_user = user_info[test_data$user_id,2:3]
    
    X = cbind(train_movie, train_user)
    X = as.matrix(X)
    N = nrow(X)
    X = cbind(1,X)
    
    
    y_test = test_data$rating
    test_movie = movie_info[as.numeric(test_data$item_id), ]
    test_movie = test_movie[,6:23]
    
    X_test = cbind(test_movie, test_user)
    X_test = as.matrix(X_test)
    X_test = cbind(1,X_test)
    n_test = nrow(X_test)
    f1 = train_data$user_id
    f2 = train_data$item_id
    
    f1_test = test_data$user_id
    f2_test = test_data$item_id
    
    
    range_f1_train=sort(unique(f1))
    range_f2_train=sort(unique(f2))
    
    train_R=length(range_f1_train)
    train_C=length(range_f2_train)
    
    
    bij_f1=list()
    bij_f2=list()
    
    for(i in 1:(train_R)){
        bij_f1[[range_f1_train[i]]]=i
    
    }
    for(i in 1:(train_C)){
        bij_f2[[range_f2_train[i]]]=i
    
    }
    
    
    
    for(i in 1:N){
    f1[i]=bij_f1[[f1[i]]]
    }
    
    for(i in 1:N){
    f2[i]=bij_f2[[f2[i]]]
    }
    
    R = max(f1)
    C = max(f2)
    
    
    for(i in 1:n_test){
    if(f1_test[i] <= length(bij_f1)){
      if(!(is.null(bij_f1[[f1_test[i]]]))){
        f1_test[i]=bij_f1[[f1_test[i]]]
      }else{
        f1_test[i] = 0
      }
    }else{
      f1_test[i] = 0}
    }
    for(i in 1:n_test){
    if(f2_test[i] <= length(bij_f2)){
      if(!(is.null(bij_f2[[f2_test[i]]]))){
        f2_test[i] = bij_f2[[f2_test[i]]]
      }else{
        f2_test[i] = 0
      }
    }else{
      f2_test[i] = 0}
    }
    
    sel = which(f1_test>0 & f2_test>0)
    X_test = X_test[sel,]
    y_test = y_test[sel]
    f1_test = f1_test[sel]
    f2_test = f2_test[sel]
    n_test = length(y_test)
    
    X_0 = matrix(rep(1,N), ncol = 1)
    X_test_0 = matrix(rep(1,n_test), ncol = 1)
    
    if(xa == 1){
        X_a = matrix(rep(1, N), ncol = 1)
        X_a_test = matrix(rep(1, n_test), ncol = 1)
    }else{
        X_a = train_movie
        X_a_test = test_movie
        X_a_test = X_a_test[sel, ]
     } 
     
    if(xb == 1){
        X_b = matrix(rep(1, N), ncol = 1)
        X_b_test = matrix(rep(1, n_test), ncol = 1)
    }else{
        X_b = cbind(1, train_user)
        X_b_test = cbind(1, test_user)
        X_b_test = X_b_test[sel, ]
    }
    
    
    if(xa == 1){
        fit_ols = lm(y ~ X - 1)$coef
        residuals_square_ols = list(sum_res_sq = sum((y_test - X_test %*% fit_ols)^2), n_test = n_test)
        saveRDS(residuals_square_ols, sprintf("movie_lens_pred_error_randomslopes_ols_t_%d.RDS", t))
        fit = scalable_crossed_random(y, X, X_a, X_b, f1, f2)
    }else{
        fit = scalable_crossed_random(y, X, X_a, X_b, f1, f2, Sigma_a_diagonal = TRUE)
    }    
    residuals = y_test - X_test %*% fit$fixed_effects
    for(i in 1:n_test){
        if(f1_test[i]==0){
            a = 0
        }else{
            a = sum(X_a_test[i, ] * fit$random_effects_fa[, f1_test[i] ])
        }
        if(f2_test[i]==0){
            b = 0
        }else{
            b = sum(X_b_test[i, ] * fit$random_effects_fb[, f2_test[i] ])
        }
        residuals[i] = residuals[i] - a - b
    }
    pred_error = list(sum_res_sq = sum((residuals)^2), n_test = n_test)
    saveRDS(pred_error, sprintf("movie_lens_pred_error_randomslopes_t_%d_xa_%d_xb_%d.RDS", t, xa, xb))
    
    return(list(pred_error = pred_error, n_test = n_test))
    
}
