# Appendix: Lisovski et al. 2016 Oikos

Original publication: http://onlinelibrary.wiley.com/doi/10.1111/oik.03796/full

## This repository contains R code used in the analysis of environmental and AIV surveillance data


## Avian Influenza Surveillance Data

Individual-level AIV infection data for dabbling ducks sampled across North America was obtained from the NIAID Influenza Research Database (www.fludb.org) on 2nd February 2016. Due to the low number of samples prior to 2005 and the potential for incomplete reports after 2013, we used only samples collected from 2005 to 2013. All records lacking spatial information, or from birds that were not live and free-living at the time of sampling, were removed, resulting in 65,358 records from individual birds representing 12 species sampled across 174 sites (sites within a radius of 50 km were treated as one site).
To quantify annual infection dynamics within each cluster, we fitted a generalized additive model (GAM with cubic spline and binomial response) using the R Package mgcv to weekly infection data and calculated the area under this infection curve as a quantitative measure for the total number of cases of infection over time.

##  Classification of Seasonality

We used the remotely-sensed Normalized Difference Vegetation Index (NDVI) to characterize seasonal variation in primary productivity across the North American continent. In this global dataset, NDVI is recorded weekly at a spatial resolution of 16x16 km, available from the National Oceanic and Atmospheric Administration (NOAA) [http://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php]. NDVI values range from -1 to +1. Negative values indicate water bodies, values between 0 and 0.2 indicate almost a complete lack of vegetation and values close to 1 indicate extremely dense green vegetation cover. We calculated seasonal amplitude (s_amp) and seasonal duration (s_dur) for each 16x16 km grid cell in North America based on annual maximum (max_NDIV) and minimum (min_NDVI) NDVI values from 2005 to 2012, after excluding those that overlapped with the coastline (to remove artefacts of ocean reflectance).  The remaining NDVI values were processed in the following sequence: 1) calculation of a smoothed seasonal curve across the 6 years of aggregated data using the loess function in R package base with span = 0.7; 2) deletion of values from albedo effects of ice/snow cover during winter (in instances where NDVI showed a trough in spring and/or autumn values before the minimum in spring and after the minimum in autumn were removed); 3) derivation of s_amp as max_NDIV – min_NDVI from the smoothed seasonal curve; 4) calculation of the duration of the annual productivity peak (s_dur), defined as the number of days that the smoothed seasonal NDVI curve is greater than or equal to 25% of s_amp.
