# [Integrating count and detection-nondetection data to model population dynamics](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.1831)

### Elise F. Zipkin, Sam Rossman, Charles B. Yackulic, J. David Wiens, James T. Thorson, Raymond J. Davis, Evan H. Campbell Grant

### Ecology

### Please contact the first author for questions about the code or data: Elise F. Zipkin (ezipkin@msu.edu)
________________________________________________________________________________________________________________________________________

## Abstract
There is increasing need for methods that integrate multiple data types into a single analytical framework as the spatial and temporal scale of ecological research expands. Current work on this topic primarily focuses on combining capture–recapture data from marked individuals with other data types into integrated population models. Yet, studies of species distributions and trends often rely on data from unmarked individuals across broad scales where local abundance and environmental variables may vary. We present a modeling framework for integrating detection–nondetection and count data into a single analysis to estimate population dynamics, abundance, and individual detection probabilities during sampling. Our dynamic population model assumes that site-specific abundance can change over time according to survival of individuals and gains through reproduction and immigration. The observation process for each data type is modeled by assuming that every individual present at a site has an equal probability of being detected during sampling processes. We examine our modeling approach through a series of simulations illustrating the relative value of count vs. detection–nondetection data under a variety of parameter values and survey configurations. We also provide an empirical example of the model by combining long-term detection–nondetection data (1995–2014) with newly collected count data (2015–2016) from a growing population of Barred Owl (Strix varia) in the Pacific Northwest to examine the factors influencing population abundance over time. Our model provides a foundation for incorporating unmarked data within a single framework, even in cases where sampling processes yield different detection probabilities. This approach will be useful for survey design and to researchers interested in incorporating historical or citizen science data into analyses focused on understanding how demographic rates drive population abundance.

## **Data**
[BDOW_combo_data.csv](https://github.com/zipkinlab/Zipkin_etal_2017_Ecol/blob/master/BDOW_combo_data.csv) - Barred owl and covariate data. The rows contain the data for each of the 106 sites in the study area. The columns contain 1) detection/nondetection data (columns beginning with BO95-BO14) collected during the first 20 years of sampling (1995-2014); 2) count data (columns beginning with BO15-BO15) collected during the last two years of sampling (2015-2016); 3) detection covariate data for the count data indicating the proportion of the site that was surveyed (columns beginning with OV); and 4) covariate data on the amount of older growth forest in each site (columns starting with AHAB). More information on the data can be found here: https://www.sciencebase.gov/catalog/item/599db7e4e4b012c075b965d0


## **Code**
[Combo_simulations.R](https://github.com/zipkinlab/Zipkin_etal_2017_Ecol/blob/master/Combo_simulations.R) - R code to run the three versions of the combined count and detection-nondetection model described in the main text of the published paper. For each modeling scenario, we provide the code to: 1) simulate data, 2) generate a script file with the JAGS model (including likelihood and prior distributions for all parameters), and 3) create initial values and run the JAGS file.

[BDOW_combo_wrapper.R](https://github.com/zipkinlab/Zipkin_etal_2017_Ecol/blob/master/BDOW_combo_wrapper.R) - R code to run the barred owl application (uses the data in the csv file and the JAGS code in the model file). Contains code to import and reshape the data and run the model file in JAGS.

[BDOW_combo_model_cloglog.R](https://github.com/zipkinlab/Zipkin_etal_2017_Ecol/blob/master/BDOW_combo_model_cloglog.R) - R file containing the JAGS code.

[BDOW_figure4_graphs.R](https://github.com/zipkinlab/Zipkin_etal_2017_Ecol/blob/master/BDOW_figure4_graphs.R) - R file with code to create panels b-e in Figure 4 of the paper.
