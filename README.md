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

import pyreadr
import numpy as np
import cvxpy as cp
import mosek
from numpy.linalg import inv, norm
from scipy.linalg import sqrtm
import pandas as pd
from sklearn.linear_model import LinearRegression
from plotnine import *
import statsmodels.api as sm

linear_model = LinearRegression()

gam = 3.0
p = 1
diagonal = False

N = round(10**gam)
R = round(N**0.6)

C = R
f1 = np.random.choice(np.arange(0, R), N, replace=True)
f2 = np.random.choice(np.arange(0, C), N, replace=True)
beta = 0.1 * np.arange(1, p + 2)
sigmae = 1

if diagonal:
    Sigmaa = np.diag(np.repeat(0.3, p + 1))
    Sigmab = np.diag(np.repeat(0.1, p + 1))
else:
    Sigmaa = np.diag(np.repeat(0.8, p + 1)) + 0.2 * np.ones((p + 1, p + 1))
    Sigmab = np.diag(np.repeat(0.8, p + 1)) + 0.2 * np.ones((p + 1, p + 1))

X = np.column_stack((np.ones(N), np.random.normal(size=(N, p))))
A = np.random.multivariate_normal(np.zeros(p+1), cov=Sigmaa, size=R)
B = np.random.multivariate_normal(np.zeros(p+1), cov=Sigmab, size=C)

y = np.dot(X, beta) + np.sum(X * A[f1, :], axis=1) + np.sum(X * B[f2, :], axis=1) + np.random.normal(scale=sigmae, size=N)
p_a = p + 1
p_b = p + 1

```
