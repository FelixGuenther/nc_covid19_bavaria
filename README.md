# nc_covid19_bavaria

## Nowcast
This is the public code repository for the manuscript:

"Felix Günther, Andreas Bender, Katharina Katz, Helmut Küchenhoff, Michael Höhle: Nowcasting the COVID-19 Pandemic in Bavaria"

A preprint of the manuscript is available here:
https://www.medrxiv.org/content/10.1101/2020.06.26.20140210v2

Code was mainly written by Günther, Bender, Höhle. We thank Titus Laska for initial code on the R(t) estimation.

The folder 'data_public' contains the datafile "data_public.csv" (tab-separated) with individual-specific information on the
reporting date and (artificial) age, and disease onset date (i.e. symptom onset date) for all Bavarian COVID-19 cases up
until 2020-04-09 (registered after 2020-03-01, excluding the first known 16 cases).

The artificial data was created the following way:
We utilize for each of the 29246 cases the official reporting date at LGL and sampled an artifical age based on a truncated normal distribution with mean=50, sd=15 and truncation 5 and 100 years, corresponding roughly to the observed age distribution in the original data. Based on an individuals' age, the reporting week and weekday, we sampled the expected reporting delay from the original Weibull GAMLSS imputation model and derived the artificial disease onset date for each case.
We removed the onset date for 16175 randomly selected cases, yielding an available disease onset date for 13087 cases (approx. 45%, as in the original data).
The dataset contains reporting date at regional health authority, artificial onset date (if available, otherwise NA), and artificial age for 29246 cases.

The folder 'code_public' contains code to estimate nowcasts and the time-dependent case reproduction number R(t)
based on a hierarchical bayesian model (nowcast) and the method of Wallinga and Teunis (2004) for R(t) as described
in the manuscript and performed on the original data. The code is structured as follows:

- File 1_analysis.R contains all necessary steps for performing the nowcast and estimating R(t) by calling functions from file analysis_fun.R

- File analysis_fun.R contains the documented functions (1) summarize_data, (2) perform_imputation, (3) estimate_nowcasts, and (4) estimate_Rt for the main analysis steps of the nowcasting. (1) summarizes the observed delay distribution, (2) estimates the imputation GAMLSS and performs the imputation, (3) performs nowcasting based on a Bayesian hierarchical model implemented in rstan, and (4) performs the estimation of R(t) based on the posterior of the nowcast.

- analysis_fun.R implements two more helper functions that are called during nowcasting and further helper functions can be found in the subfolder 'general'

- 'general' contains functions for the estimation of the nowcast and R(t), most notably the stan model for the Bayesian hierarcical nowcast and the customized R-function est.R0.TD (customized R0::est.R0.TD) for estimation of R(t).

Results of the nowcast are written to the 'results_public' folder
In case of questions, feel free to contact felix.guenther(at)stat.uni-muenchen.de.

## Segmented regression
Additionally to the nowcast code we provide code for analysing the epidemic curve (estimated based on the nowcast) using segmented regression. A manuscript for the analysis of the Bavarian data (in German) can be found here:
https://www.stablab.stat.uni-muenchen.de/_assets/docs/analyse_covid19_bayern.pdf

In 'code public' we provide the file '2_bp_analysis.R' that implements the performed analysis based on the results of the nowcast on artificial data from the 'results_public' folder. The file '3_bp_analysis_bavaria_20-05-05.R' shows the code for the segmented regression analysis based on the nowcast for the officially reported Bavarian COVID-19 case data from May 5th, 2020. The nowcast results for this date can be found in results_public/nowcasting_results_bavaria_2020-05-05.csv.

Functions to estimate the segmented regression models are implemtend in 'breakpoint_fun.R'. One of the functions is based on discrete optimization over all possible combinations of breakpoints using quasi-poisson regression based on the 'glm' function. The second function estimates the segmented regression model using the 'segmented::segmented' function.

