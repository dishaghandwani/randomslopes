## randomslopes
This package provides scalable solutions to crossed random effects models with random slopes. 

Accompanying paper, (https://arxiv.org/abs/2307.12378)


``Scalable solutions to crossed random effects models with random slopes" 
D. Ghandwani, S. Ghosh, T. Hastie, A. Owen


##Contents

 - `randomslopes.R` R file containing functions to perform method of moments, backfitting (with or without clubbing ) and variational methods. We also provide code for estimating standard errors for fixed effects. 
 -  `randomslopes.py` python file for performing the same functions in Python language.
 -  `simulation.R` R file for perform the simulations reflected in the paper.
 -  `processing.R` R file to summarise the simulations.
 -   `plots.Rmd` Markdown file to generate the plots using the summarised results.
 -   `movie_lens.R` R file to perform analysis on Movie Lens 100K dataset (https://www.kaggle.com/datasets/prajitdatta/movielens-100k-dataset)
 -   `movie_lens_summarise_randomslopes.R` R file to summarise the prediction errors on the Movie Lens dataset.
 -   `simulations_randomslopes.zip` Zip folder containing the simulation results produced through simulations.R
 -   `simualations_results.zip` Zip folder containing the summarised results compiled through processing.R
 -   `movie_lens_results.zip` Zip folder containing the results on fit on the movie lens dataset. 


## Examples

```py
from randomslopes import *
import numpy as np
import cvxpy as cp
import mosek # you need a license for this
from numpy.linalg import inv, norm
from scipy.linalg import sqrtm
import pandas as pd
from sklearn.linear_model import LinearRegression
from plotnine import *
import statsmodels.api as sm

linear_model = LinearRegression()

gam = 3.0
p = 7
p_a = 3
p_b = 4

N = round(10**gam)
R = round(N**0.6)

C = R
f1 = np.random.choice(np.arange(0, R), N, replace=True)
f2 = np.random.choice(np.arange(0, C), N, replace=True)
beta = 0.1 * np.arange(1, p + 2)
sigmae = 1

Sigmaa = np.diag(np.repeat(0.8, p_a)) + 0.2 * np.ones((p_a, p_a))
Sigmab = np.diag(np.repeat(0.1, p_b))

X = np.column_stack((np.ones(N), np.random.normal(size=(N, p))))
X_a = np.column_stack((np.ones(N), np.random.normal(size=(N, p_a))))
X_b = np.column_stack((np.ones(N), np.random.normal(size=(N, p_b))))
A = np.random.multivariate_normal(np.zeros(p+1), cov=Sigmaa, size=R)
B = np.random.multivariate_normal(np.zeros(p+1), cov=Sigmab, size=C)

y = np.dot(X, beta) + np.sum(X_a * A[f1, :], axis=1) + np.sum(X_b * B[f2, :], axis=1) + np.random.normal(scale=sigmae, size=N)

fit = scalable_crossed_random(y, X, X_a, X_b, f1, f2, Sigma_a_diagonal = False, Sigma_b_diagonal = TRUE, 
                            variational = True, clubbing = True, method_of_mom = True, beta_covariance = TRUE)

```
