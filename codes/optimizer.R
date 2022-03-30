# Loading libraries
library(nloptr) # For COBYLA, Nelder-Mead and Subplex

# Original temperature values
original_temperature <- c(130,
                          150,
                          150.1,
                          150.15,
                          150.2,
                          150.25,
                          150.3,
                          150.4,
                          150.45,
                          150.5,
                          150.55)
# no. of parameters
n <- 15

# Objective function
objective <- function(k)
{
  
  # Running simulation
  system(paste0("isom.hysys_comparision.exe -override=E[1]=",k[1],
                ",E[2]=",k[2],",E[5]=",k[3],",E[6]=",k[4],",E[16]=",k[5],
                ",E[27]=",k[6],",E[29]=",k[7],",E[30]=",k[8],",E[35]=",k[9],
                ",E[37]=",k[10],",E[41]=",k[11],",E[44]=",k[12],",E[45]=",k[13],
                ",E[48]=",k[14],",E[49]=",k[15]))
  
  # Reading data
  data <- read.csv("isom.hysys_comparision_res.csv")
  temperature <- data$T_in_c[-nrow(data)]
  cat("temperature is : ",temperature,"\n")
  
  # Mean Square Error (MSE)
  error <- mean(c(original_temperature-temperature)^2)
  cat("Error value is: ",error,"\n")
  return(error)
}


initial_par = c(148.93,154.28, 150.98, 155.92, 105.4, 135.45, 154.54,
                98.63, 177.32, 59.8, 284.97, 112.05, 265, 265, 295.62)

opt=cobyla(x0 = initial_par,
       fn = objective,
       lower = c(140, 148, 140, 140, 100, 125, 135, 93, 160, 54, 265,
                 106, 250, 250, 280),
       upper = c(160, 170, 160, 170, 113, 150, 170, 110, 195, 65, 310, 
                 125, 270, 270, 330),
       nl.info=TRUE,
       control = list(xtol_rel= 1e-16, maxeval=10000))