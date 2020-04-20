# nc_covid19_bavaria

This is the public code repository for the manuscript:

"Felix Günther, Andreas Bender, Katharina Katz, Helmut Küchenhoff, Michael Höhle: Nowcasting the COVID-19 Pandemic in Bavaria"

The folder 'data_public' contains the file "data.public.csv" (tab-separated) with individual-specific information on the
reporting date and (artificial) age, and disease onset date (i.e. symptom onset date) for all Bavarian COVID-19 cases up
until 2020-04-09 (and after 2020-03-01, excluding the first known 16 cases).

The artificial data was created the following way:
We utilize for each of the 29246 cases the official reporting date at LGL and sample an artifical age based on a truncated normal distribution with mean=50, sd=15 and truncation 5 and 100 years, corresponding roughly to the observed age distribution in the original data. Based on an individuals' age, the reporting week and weekday, we sample the expected reporting delay from the original Weibull GAMLSS imputation model and derive the artificial disease onset date for each case.
We remove the onset date for 16175 randomly selected cases, yielding an available disease onset date for 13087 cases (approx. 45% as in the original data).
The dataset consets reporting date, onset date (if available, otherwise NA), and artificial age for 29246 cases.

The folder 'code_public' contains code to estimate nowcasts and the time-dependent case reproduction number R(t)
based on a hierarchical bayesian model (nowcast) and the method of Wallinga and Teunis (2004) for R(t) as described
in the manuscript. The code is structured as follows:

- File 1_analysis.R contains all necessary steps for performing the nowcast and estimating R(t) by calling functions from file analysis_fun.R

- File analysis_fun.R contains the documented functions (1) summarize_data, (2) perform_imputation and (3) estimate_nowcasts for the main analysis steps of the nowcasting. (1) summarizes the observed delay distribution, (2) estimates the imputation GAMLSS and performs the imputation, and (3) performs 3 different nowcasts.

- analysis_fun.R implements two more helper functions that are called during nowcasting and further helper functions can be found in the subfolder 'general'

- 'general' contains functions for the estimation of the nowcast and R(t), most notably the function nowcast_w.R as an extension of the nowcast function of the surveillance package (surveillance::nowcast), that implements the estimation of weekday effects (or general extensions of the 'W' matrix) in the estimation of the delay distribution for the nowcast and the customized est.R0.TD (customizin R0::est.R0.TD) for estimation of R(t).

In case of questions, feel free to contact felix.guenther(at)stat.uni-muenchen.de.
