# nc_covid19_bavaria

## Nowcast
This is the public code repository for the manuscript:

"Felix Günther, Andreas Bender, Katharina Katz, Helmut Küchenhoff, Michael Höhle: Nowcasting the COVID-19 Pandemic in Bavaria"

A preprint of the manuscript is available here:
https://www.stablab.stat.uni-muenchen.de/_assets/docs/nowcasting_covid19_bavaria.pdf

Code was mainly written by Günther, Bender, Höhle. We thank Titus Laska for initial code on the R(t) estimation.

The folder 'data_public' contains the datafile "data_public.csv" (tab-separated) with individual-specific information on the
reporting date and (artificial) age, and disease onset date (i.e. symptom onset date) for all Bavarian COVID-19 cases up
until 2020-04-09 (registered after 2020-03-01, excluding the first known 16 cases).

The artificial data was created the following way:
We utilize for each of the 29246 cases the official reporting date at LGL and sampled an artifical age based on a truncated normal distribution with mean=50, sd=15 and truncation 5 and 100 years, corresponding roughly to the observed age distribution in the original data. Based on an individuals' age, the reporting week and weekday, we sampled the expected reporting delay from the original Weibull GAMLSS imputation model and derived the artificial disease onset date for each case.
We removed the onset date for 16175 randomly selected cases, yielding an available disease onset date for 13087 cases (approx. 45%, as in the original data).
The dataset consets reporting date, artificial onset date (if available, otherwise NA), and artificial age for 29246 cases.

The folder 'code_public' contains code to estimate nowcasts and the time-dependent case reproduction number R(t)
based on a hierarchical bayesian model (nowcast) and the method of Wallinga and Teunis (2004) for R(t) as described
in the manuscript and performed on the original data. The code is structured as follows:

- File 1_analysis.R contains all necessary steps for performing the nowcast and estimating R(t) by calling functions from file analysis_fun.R

- File analysis_fun.R contains the documented functions (1) summarize_data, (2) perform_imputation and (3) estimate_nowcasts for the main analysis steps of the nowcasting. (1) summarizes the observed delay distribution, (2) estimates the imputation GAMLSS and performs the imputation, and (3) performs 3 different nowcasts.

- analysis_fun.R implements two more helper functions that are called during nowcasting and further helper functions can be found in the subfolder 'general'

- 'general' contains functions for the estimation of the nowcast and R(t), most notably the function nowcast_w.R as an extension of the nowcast function of the surveillance package (surveillance::nowcast), that implements the estimation of weekday effects (or general extensions of the 'W' matrix) in the estimation of the delay distribution for the nowcast and the customized est.R0.TD (customizin R0::est.R0.TD) for estimation of R(t).

Results of the nowcast are written to the 'results_public' folder
In case of questions, feel free to contact felix.guenther(at)stat.uni-muenchen.de.

## Segmented regression
Additionally to the nowcast code we provide code for analysing the epidemic curve (estimated based on the nowcast) using segmented regression. A manuscript for the analysis of the Bavarian data (in German) can be found here:
https://www.stablab.stat.uni-muenchen.de/_assets/docs/analyse_covid19_bayern.pdf

In 'code public' we provide the file '2_bp_analysis.R' that implements the performed analysis based on the results of the nowcast from the 'results_public' folder. 
Functions to estimate the segmented regression models are implemtend in 'breakpoint_fun.R'. One of the functions is based on discrete optimization over all possible combinations of breakpoints using quasi-poisson regression based on the 'glm' function. The second function estimates the segmented regression model using the 'segmented::segmented' function.

