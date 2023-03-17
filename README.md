# Distributed lag non-linear models

This repository contain codes and files related to DLNM applied to COVID-19 data (counts of cases) in Belgium (per municipality, all provinces) in 2020-2022.

The raw data is located in the **data** folder, and the output figures and datasets are in the **output** folder.

Data for Covid-19 were downloaded [here](https://epistat.sciensano.be/covid/) and includes two datasets:

1) Confirmed cases by date and municipality
2) Administered vaccines by week, municipality, age and dose


Below is the description of each R file

| File name  | Description |
| ------------- | ------------- |
|  **covid19_exploration.Rmd** | Basic exploration of COVID19 cases in 2020-2023 in Belgium, such as number of month and weeks with data per municipality, missing data  |
| **pollution_exploration.Rmd** | Basic  exploration of the data on the administered vaccines in 2020-2023 in Belgium  |

? Pollution data on weighted daily mean pollutant concentrations (BC, NO2, PM10, PM2.5, O3)  are not deposited here due to their large size, they can be made available on request or downloaded [here](http://ftp.irceline.be/rio4x4/gemeente/).
