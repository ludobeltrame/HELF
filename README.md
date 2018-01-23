# HYDRO-EPIDEMIOLOGICAL MODEL FOR LIVER FLUKE (HELF)

HELF is a mechanistic model which explicitly simulates the life-cycle of the liver fluke (Fasciola hepatica) in connection with temperature and soil moisture conditions, with a daily time step. The model integrates a hydrological and an epidemiological component. The former builds on TOPMODEL [Beven et al., 1995]. This is a widely-used catchment-scale rainfall-runoff model, employed here to simulate soil moisture dynamics. The latter was developed based on current understanding of the parasite life-cycle and its dependence upon temperature and soil moisture.

The model is coded in matlab. The function that runs the model is HELF.m. This in turn calls 4 other functions:

- calculateDailyPetHargreaves.m, which calculates daily potential evapotranspiration using Hargreaves equation.

- makeTIclasses.m, which discretizes the distribution of Topographic Index (TI) values over the catchment into a user-defined number of classes.

- HELF_hydro_model_component.m, which evaluates soil moisture for each TI class and calculates streamflow at the catchment outlet, based on concepts from TOPMODEL.

- HELF_fluke_model_component.m, which simulates the liver fluke life-cycle driven by temperature and soil moisture conditions. Running of this component relies on the class “FlukeEggToMeta.m”. This is employed to calculate the stage-specific development/mortality rates (using classes rateCalcLinear, rateLinearFuzzyMembership, rateCalcConstant and rateCombine), and, in turn, relies on the class “stage_cycle_cohort.m”, which describes the cohort-based dynamics of each life-cycle stage.
