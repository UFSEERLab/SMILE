# SMILE

This repository is currently being developed by members of the SEER Lab. It stems from the published study [Decoupling environmental effects and host population dynamics for anthrax, a classic reservoir-driven disease](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208621) by Gomez et.al. 2018. You can check out the original code associated to that publication in the [early commits of this repository](https://github.com/UFSEERLab/SMILE/tree/e176302c90e16205b0211906e1ef93e641776407).

We are currently developing the model to incorporate the effect of vaccination. 

For contributing, please clone the repository and create a new branch. You can submit any changes with a pull request. For questions, please submit an issue.

The code in this repository contains the functions defined by the SMILE compartmental model. SMILE, model was created to model environmentally mediated diseases with a focus on anthrax dynamics based on outbreaks on the Northwestern plains of the United states. The model consists of five compartments, Suceptible, Immune, Infected, Local Infectious Zone, Environment. The model assumes that the transmission of the disease is indirect and the intensity of the transmission is driven by a pure birth stochastic process that accounts for heterogeneity in the dispersion effort of the pathogen. It also allows specifying a sinusoidal function that drives the seasonality in the outbreaks and decouples it from basic recruitment due to population dynamics. We also have derived a form of R0 for the disease based on the assumptions deriving the infection probability as a stochastic process. The las set of functions are used to estimate the paramters that define the distribution of the dispersion effort and the shape of the seasonality of the envrionmental driver.

In addition to the sorce function, we provide files with the data from a Montana Outbreak as well as the scripts written for the analysis ad figure reproduction provided in the text.
