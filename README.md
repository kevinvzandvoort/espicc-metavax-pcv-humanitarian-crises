# **Effectiveness of pneumococcal vaccination campaigns in humanitarian settings: a modelling study. Data, model code, and analyses.**

This repository contains the data and code used for all analyses described in our manuscript:

Van Zandvoort K, Bobe MO, Hassan AI, Cummings R, Warsame A, Pell CP, Abdi MI, Ibrahim S, McGowan CR, Mulholland EK, Satzke C, Eggo RM, Hergeye MA, Checchi F, Flasche S. *Effectiveness of pneumococcal vaccination campaigns in humanitarian settings: a modelling study.* Available at <https://doi.org/10.1101/2025.05.16.25327803>.

This work is part of a larger study: Evaluating Strategies for Pneumococcal Immunization Campaigns in Crises [(ESPICC)](https://www.elrha.org/project/pneumococcal-vaccination-strategies-for-crisis-affected-populations/).

### **Model**

We developed *Metavax*, a flexible framework that can be used to run compartmental metapopulation models with an arbitrary number of populations, age-groups, vaccine strategies (and vaccinated strata), and compartmental model structures. *Metavax* is still in active development (<https://github.com/kevinvzandvoort/pcvm>), but the version used for this analysis has been copied to the `./model/metavax` folder.

### **Data**

The `./data` folder has two subdirectories:

-   `./data/demography` for demography related data (contact data and population data)

-   `./data/epidemiology` for epidemiological data (pneumococcal carriage prevalence, pneumococcal clearance rate

The `./data_load_all.R` script creates a list object named data, which is filled by sourcing several other R-scripts in the `./data/*` folders, that read and pre-process the data.

### **Running the model**

To run the model, we need to specify the location of the *Metavax* root folder in a `METAVAX_FOLDER` variable. We can then source the correct *Metavax* model implementation, which automatically compiles and loads the C++ code, and registers the number of compartments in R. The base model in our analysis uses the diamond-shaped streppneumo model:

``` r
METAVAX_FOLDER = "./model/metavax"
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))
```

Once the model is compiled and loaded, we need to source three files:

1.  `./functions.R`, which defines a number of project-specific functions used in the model and analyses (e.g. the log-likelihood function used when fitting the model)
2.  `./data_load_all.R`, as mentioned in the previous section
3.  `./model_setup.R`, which defines the model parameters list which is needed as the input for the model.
    -   .`/model_setup.R` defines a `modelSetup()` function, which returns model parameters set-up to run the base Digaale model with default values, and without vaccination.
    -   the model_parameters can be adapted using the `setParameter()` function, e.g. to add a vaccination campaign at a certain coverage value at a given timepoint (see the main *metavax* GitHub repository).

``` r
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

model_params = modelSetup()
```

The model can then be ran once using the `runModel()` function, which takes a number of required arguments:

-   `model_params`, a list with the model parameters used in the model

-   `initial_state`, a numeric vector with initial states for each compartment in the model (can be the output of a previous model run)

-   `steady_state`, a boolean to indicate whether the model should be ran to steady state (opposed to a number of time-steps)

    -   If `steady_state = FALSE`, the ODEs will be solved for the specified number of time-steps using the `{deSolve}` R-package

    -   If `steady_state = TRUE`, the model will be ran to steady state using the `{rootSolve}` R-package

        -   The default `{rootSolve}` package does not allow a list as input for compiled models. An adjusted package that does allow this, and which is essential for this model, can be installed from <https://github.com/kevinvzandvoort/rootSolve>

-   `times`, a numeric vector with the number of time-steps for which model results should be returned

-   `incidence`, a boolean to indicate whether model outputs should include incidence (in addition to prevalence)

``` r
result_prevacc = runModel(
  model_params = model_params_prevac_current,
  initial_state = model_params$global_settings$initial_states,
  steady_state = TRUE)

result_prevacc = result_prevacc$value
```

### **Replicating the analysis**

A number of scripts can be used to replicate the analysis:

To fit the models to the carriage prevalence data, `./model_fit_main.R`, `./model_fit_sensitivity.R`, and `model_fit_neutral.R` can be used to fit the base model for the four settings, the base model with different parameter values that affect transmission, and an alternative model structure that is structurally neutral. This uses MCMC (using the `{BayesianTools}` package) to fit different model parameters, as specified in `./model_setup.R`. The `./model_fit_assess.R` script can be used to assess this fit.

The fitted models are then used by the `./model_postvacc*` scripts to generate model results under different vaccination strategies:

-   `./model_postvacc_main.R`: runs the five main single-dose PCV campaigns (and sixth no-vaccination scenario) for each setting.

-   `./model_postvacc_with_routine.R`: runs the single-dose PCV campaigns in the presence of routine-vaccination prior to the PCV campaign

-   `./model_postvacc_sensitivity.R`: runs the PCV campaigns for the parametric sensitivity analyses.

-   `./model_postvacc_2dose_infants.R`: runs a sensitivity analysis where we assumed two PCV doses in infants for all PCV campaigns.

Finally, the modelled results are processed and tables and figures generated using the `./results*` scripts.
