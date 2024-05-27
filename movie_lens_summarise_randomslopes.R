setwd("movie_lens_results")
dir()
results = NULL
for(t in 1:5){
  x = as.numeric(readRDS(sprintf("movie_lens_pred_error_randomslopes_ols_t_%d.RDS", t)))
  x = c(x, as.numeric(readRDS(sprintf("movie_lens_pred_error_randomslopes_t_%d_xa_1_xb_1.RDS", t))[1]))
  x = c(x, as.numeric(readRDS(sprintf("movie_lens_pred_error_randomslopes_t_%d_xa_1_xb_2.RDS", t))[1]))
  x = c(x, as.numeric(readRDS(sprintf("movie_lens_pred_error_randomslopes_t_%d_xa_2_xb_1.RDS", t))[1]))
  x = c(x, as.numeric(readRDS(sprintf("movie_lens_pred_error_randomslopes_t_%d_xa_2_xb_2.RDS", t))[1]))
  results = rbind(results, x)
}
View(results)
results = data.frame(results)
colnames(results) = c('ols','n_test', 'xa_1_xb_1', 'xa_1_xb_2', 'xa_2_xb_1', 'xa_2_xb_2')
write.csv(results, "movie_lens_randomslopes.csv")
r = colSums(results)
r/r[2]
