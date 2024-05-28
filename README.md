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
