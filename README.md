# A generalized phase I/II dose optimization trial design with multi-categorical and multi-graded outcomes
R codes to implement the simulation study in "A generalized phase I/II dose optimization trial design with multi-categorical and multi-graded outcomes"

# Description
We propose a phase I/II design that simultaneously considers multi-categorical toxicity and efficacy with multi-graded outcomes, measured as quasi-continuous probability based on prespecified weight matrices of clinical significance. Following keyboard design, our approach aims to screen out overly toxic doses by the toxicity probability intervals and adaptively makes dose escalation or de-escalation decisions by comparing the posterior distributions of dose desirability (utility) among the adjacent levels of the current dose. It helps to more accurately identify the optimal biological dose (OBD) in a non-monotonically increasing dose-efficacy relationship.

Here is a brief introduction about the function of each R file. Detailed instruction is provided in each file by comprehensive comments. It is necessary to set the working path accordingly for storing results and load data for tables and visualizations. Some simulation codes also provide HPC options to run some time-consuming scenarios.

* get_oc.R: to conduct simulations and obtain the operating characteristics of the proposed design under certain scenarios.
* scenario_generator_local.R: to conduct various scenario data for simulation, including true MTD and OBD, true toxicity and efficacy probabilities, etc.
* Ours_sensitivity_analysis.R: sensitivity analysis section regarding sample size cutoff and sub-interval width of efficacy unit.
* SA_WE.R: sensitivity analysis section regarding various random values of toxicity severity matrix $W$ and efficacy weighting matrix $E$.
* gBOIN-ET_poor_allication.R: the manually updated gboinet() function to obtain the poor allocation results in the simulation study.
* Random_scenario_benchmark.R: major methods regarding the random scenario generator and simulation results.
* OBD_selection_visualization_local.R: major data wrangling process for visualization of simulation results.

We will consistently provide the necessary maintenance for these major R files.
