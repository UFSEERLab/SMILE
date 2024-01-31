# SMILE
The code in this repository contains the functions defined by the SMILE compartmental model. SMILE, model was created to model environmentally mediated diseases with a focus on anthrax dynamics based on outbreaks on the Northwestern plains of the United states. The model consists of five compartments, Suceptible, Immune, Infected, Local Infectious Zone, Environment. The model assumes that the transmission of the disease is indirect and the intensity of the transmission is driven by a pure birth stochastic process that accounts for heterogeneity in the dispersion effort of the pathogen. It also allows specifying a sinusoidal function that drives the seasonality in the outbreaks and decouples it from basic recruitment due to population dynamics. We also have derived a form of R0 for the disease based on the assumptions deriving the infection probability as a stochastic process. The las set of functions are used to estimate the paramters that define the distribution of the dispersion effort and the shape of the seasonality of the envrionmental driver.

In addition to the sorce function, we provide files with the data from a Montana Outbreak as well as the scripts written for the analysis ad figure reproduction provided in the text.


Reference:
Gomez, J.P., D. Nekorchuk, S.J. Ryan, L. Mao, J.M. Ponciano, J.K. Blackburn. Compartmental model for environmentally mediated indirect disease transmission. To be Submitted.
