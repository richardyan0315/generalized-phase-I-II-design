# A generalized phase I/II dose optimization trial design with multi-categorical and multi-graded outcomes
R codes to implement the simulation study in "A generalized phase I/II dose optimization trial design with multi-categorical and multi-graded outcomes" and corresponding output data for visualization reproduction.

# Description
We propose a phase I/II design that simultaneously considers multi-categorical toxicity and efficacy with multi-graded outcomes, measured as quasi-continuous probability based on prespecified weight matrices of clinical significance. Following keyboard design, our approach aims to screen out overly toxic doses by the toxicity probability intervals and adaptively makes dose escalation or de-escalation decisions by comparing the posterior distributions of dose desirability (utility) among the adjacent levels of the current dose. It helps to more accurately identify the optimal biological dose (OBD) in a non-monotonically increasing dose-efficacy relationship.

Here is a brief introduction about the function of each R file. Comprehensive comments provide detailed instructions in each file. Setting the relative working path for storing results and loading data for tables and visualizations is necessary. All the necessary output files for quick reproduction of visualization can be accessed via **`Outputs`**

* getoc.R: to conduct simulations and obtain the operating characteristics of the proposed design under certain scenarios.
* scenarioGenerator.R: to obtain necessary scenario parameters.
* simulationResultsCleaning.R: data wrangling for visualization.
* visualization.R: reproduce result figures as presented in the manuscript and the supplementary document.
* model.stan: stan model.

We will consistently provide the necessary maintenance for these major R files.
