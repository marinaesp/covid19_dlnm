## ----setup, include=FALSE----------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE)


## ----------------------------------------------------------------------------------------------------
library(tidyverse)
library(here)
library(INLA)
library(tsModel)
library(sf)
library(spdep)
library(sp)
library(rgdal)
library(dlnm)


## ----------------------------------------------------------------------------------------------------
data <- read.csv(here("output", "master_dataset.csv"))


## ----------------------------------------------------------------------------------------------------
data_sort <- data |> 
  arrange(year, week_num)


## ----------------------------------------------------------------------------------------------------
nlag = 6


## ----------------------------------------------------------------------------------------------------
lag_ozone <- tsModel::Lag(data_sort$pollution_ozone_wk, group = data_sort$mcp_code, k = 0:nlag)
lag_no2 <- tsModel::Lag(data_sort$pollution_no2_wk, group = data_sort$mcp_code, k = 0:nlag)
lag_pm25 <- tsModel::Lag(data_sort$pollution_pm25_wk, group = data_sort$mcp_code, k = 0:nlag)
lag_pm10 <- tsModel::Lag(data_sort$pollution_pm10_wk, group = data_sort$mcp_code, k = 0:nlag)
lag_bc <- tsModel::Lag(data_sort$pollution_bc_wk, group = data_sort$mcp_code, k = 0:nlag)
lag_vaccine <-tsModel::Lag(data_sort$cumul_week, group = data_sort$mcp_code, k = 0:nlag)


## ----------------------------------------------------------------------------------------------------
colSums(is.na(lag_ozone))
#lag_ozone <- lag_ozone[!(data_sort$week_num < 47 & data_sort$year == 2020),]
lag_ozone <- lag_ozone[data_sort$year > 2020,]
lag_no2 <- lag_no2[data_sort$year > 2020, ]
lag_bc <- lag_bc[data_sort$year > 2020, ]

lag_pm25 <- lag_pm25[data_sort$year > 2020, ]
lag_pm10 <- lag_pm10[data_sort$year > 2020, ]


lag_vaccine <- lag_vaccine[data_sort$year > 2020, ]

# remove weeks year 2020 from the main dataset
data_prep <- data_sort[data_sort$year > 2020,]


## ----------------------------------------------------------------------------------------------------
basis_ozone <- crossbasis(lag_ozone,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$pollution_ozone_wk, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))

basis_no2 <-crossbasis(lag_no2,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$pollution_no2, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))
  
  
basis_pm25 <- crossbasis(lag_pm25,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$pollution_pm25, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))
  
basis_pm10 <- crossbasis(lag_pm10,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$pollution_pm10, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))
  
basis_bc <-crossbasis(lag_bc,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$pollution_bc, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))

basis_vaccine <- crossbasis(lag_vaccine,
                    argvar = list(fun = "ns", knots = equalknots(data_prep$cumul_week, 2)),
                    arglag = list(fun = "ns", knots = nlag/2))


## ----------------------------------------------------------------------------------------------------
colnames(basis_ozone) = paste0("basis_ozone.", colnames(basis_ozone))
colnames(basis_no2) = paste0("basis_no2.", colnames(basis_no2))
colnames(basis_pm25) = paste0("basis_pm25.", colnames(basis_pm25))
colnames(basis_pm10) = paste0("basis_pm10.", colnames(basis_pm10))
colnames(basis_bc) = paste0("basis_bc.", colnames(basis_bc))
colnames(basis_vaccine) = paste0("basis_vaccine.", colnames(basis_vaccine))


## ----------------------------------------------------------------------------------------------------
precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))


## ---- message=F--------------------------------------------------------------------------------------
# load shape file for Belgium
map <- readOGR("./data/belgium_shape/Apn_AdMu.shp")

neig.map <- poly2nb(map,row.names = map$AdMuKey)


## ---- echo=FALSE-------------------------------------------------------------------------------------
tiff(file = "./figs/adjacency_mat.tiff",width = 7, height = 6, units = 'in', res = 250)
plot(map, border = grey(0.5))   #Don't close the graph
Coords <- coordinates(map)
# Connect neighbours with a line
plot(neig.map, 
     coords = Coords, 
     add = TRUE,
     pch = 20,
     lwd = 1,
     col = 2)
dev.off()


## ----------------------------------------------------------------------------------------------------
adj.file <- "output/adjacency.mat"
if (!file.exists(adj.file)) nb2INLA(adj.file, neig.map)



## ----------------------------------------------------------------------------------------------------
mymodel <- function(formula, data = df, family = "nbinomial", config = FALSE)

  {
  model <- inla(formula = formula, data = data, family = family, offset = log(E),
       control.inla = list(strategy = 'adaptive'), 
       control.compute = list(dic = TRUE, config = config, 
                              cpo = TRUE, return.marginals = FALSE),
       control.fixed = list(correlation.matrix = TRUE, 
                            prec.intercept = 1, prec = 1),
       control.predictor = list(link = 1, compute = TRUE), 
       verbose = FALSE)
  model <- inla.rerun(model)
  return(model)
}


## ----------------------------------------------------------------------------------------------------
# total number of weeks
nweeks <- length(unique(data_prep$week))
#total number of muncipalities
nmcp <- length(unique(data_prep$mcp_code))
# total number of provinces
nprov <- length(unique(data_prep$province))


## ----------------------------------------------------------------------------------------------------
#index for municipality (run chunk 15 first!)
mcp_index <- data.frame(
  mcp_code = map$AdMuKey,
  mcp_id = seq(1, length(map$AdMuKey))
)
mcp_index$mcp_code <- as.numeric(mcp_index$mcp_code)

data_prep <- left_join(data_prep, mcp_index, by = "mcp_code")


## ----------------------------------------------------------------------------------------------------
Y <- data_prep$cases_per_week  #observed cases of COVID-19 per week
N <- length(Y)  #total length of the dataset
E <- data_prep$pop_tot/10^3  #for the incidence per 1000 population
T1 <- data_prep$week_num  #week indicator = week number
S1 <- data_prep$mcp_id  #municipality indicator (1 to 581)


## ----------------------------------------------------------------------------------------------------
df <- data.frame(Y, E, T1,S1)

