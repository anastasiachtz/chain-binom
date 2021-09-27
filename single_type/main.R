library(rstan)
library("R2OpenBUGS")
library(readr)
rstan_options(auto_write = TRUE)           
options(mc.cores = parallel::detectCores())

population_GR <- 10423054    #Population numbers: https://www.worldometers.info/world-population/population-by-country/
time_series_covid19_confirmed_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv")
time_series_covid19_deaths_global <- read_csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv")
full_cum_confGR <- subset(time_series_covid19_confirmed_global, `Country/Region`=="Greece")[c(5:length(time_series_covid19_confirmed_global))]
full_cum_deadGR <- subset(time_series_covid19_deaths_global, `Country/Region`=="Greece")[c(5:length(time_series_covid19_deaths_global))]
timeGR <- as.vector(names(full_cum_confGR))
cum_confGR <- as.numeric(full_cum_confGR)[(match('3/7/20', timeGR)):(match('3/19/21', timeGR))]
cum_deadGR <- as.numeric(full_cum_deadGR)[(match('3/7/20', timeGR)):(match('3/19/21', timeGR))]
fit_timeGR <- timeGR[(match('3/7/20', timeGR)):(match('3/19/21', timeGR))]
sample_daysGR <- length(fit_timeGR)
sample_weeksGR <- sample_daysGR/7
sample_time=1:sample_daysGR
new_casesGR <- rep(0,sample_daysGR)
new_casesGR[1] <- cum_confGR[1]
for (t in 2:sample_daysGR){ 
  new_casesGR[t] <- cum_confGR[t]-cum_confGR[t-1]
}
deadGR <- rep(0,sample_daysGR)
deadGR[1] <- cum_deadGR[1]
for (t in 2:sample_daysGR){
  deadGR[t] <- cum_deadGR[t]-cum_deadGR[t-1]
  if (deadGR[t]<0) {
    deadGR[t]=0
  } 
}
# Discretize infection to death distribution
infection_death = rep(0,sample_daysGR)
infection_death[1] = pgamma(1.5,shape=6.29,rate =0.26) - pgamma(0,shape=6.29,rate =0.26)
for(i in 2:length(infection_death)) {
  infection_death[i] = pgamma(i+.5,shape=6.29,rate =0.26) - pgamma(i-.5,shape=6.29,rate =0.26)
}
mean_ifr <- c(0.01142,0.007994,0.0115)

# Modify data into a form suitable for Stan/BUGS
cov_data = list(n_obs = sample_daysGR,
                n_weeks = sample_weeksGR,
                n_pop = population_GR,
                yC = new_casesGR,
                yD = deadGR,
                S0 = population_GR - 100,
                I0 = 100,
                I_D = infection_death,
                IFRmu = mean_ifr,
                sigma_cp1=22,
                sigma_cp2=45,
                ifr_cp1 = match('7/28/20', fit_timeGR),
                ifr_cp2 = match('11/1/20', fit_timeGR),
                ifr_cp3 = match('2/9/21', fit_timeGR))

#############################################################################
                
# Specify parameters to monitor
parameters_nuts = c("beta", "rho", "phiD", "ifr", "sigma", "R_t",
               "eta", "eta_w", "c_tot", "muD",  "c_unrep",  
               "log_lik", "dev")


set.seed(1234)

#init_cases <- c(821, 888, 879, 809, 721, 652, 734, 900, 800, 810, 670, 610, 560, 425, 400, 280, 250, 250, 180, 160, 144, 138, 136, 132, 129, 127, 122, 119, 117, 114, 110, 107, 104, 102, 99, 95, 94, 93, 90, 90, 88, 89, 90, 95, 92, 95, 97, 98, 99, 99, 98, 98, 97, 95, 92, 89, 95, 100, 100, 91, 75, 75, 80, 65, 70, 70, 55, 54, 50, 41, 45, 38, 50, 45, 33, 40, 31, 40, 40, 35, 26, 24, 30, 30, 30, 30, 30, 15, 25, 13, 12, 12, 11, 10, 10, 9, 9, 8, 8, 8, 8, 9, 9, 9, 10, 11, 11, 12, 13, 15, 17, 19, 22, 26, 31, 36, 42, 49, 58, 66, 76, 85, 94, 101, 107, 111, 113, 116, 118, 119, 125, 128, 132, 137, 147, 160, 171, 186, 203, 222, 243, 266, 291, 308, 318, 326, 325, 325, 334, 341, 348, 356, 363, 371, 377, 384, 392, 397, 405, 410, 420, 426, 434, 441, 450, 460, 469, 479, 489, 503, 516, 527, 543, 555, 570, 585, 597, 607, 618, 630, 639, 645, 649, 655, 661, 666, 667, 671, 682, 690, 701, 715, 733, 751, 779, 806, 834, 865, 901, 938, 972, 1007, 1036, 1064, 1089, 1107, 1117, 1120, 1112, 1102, 1087, 1076, 1071, 1058, 1058, 1056, 1062, 1078, 1103, 1140, 1193, 1256, 1350, 1472, 1616, 1810, 2043, 2343, 2701, 3157, 3712, 4354, 5090, 5862, 6717, 7592, 8405, 9154, 9811, 10343, 10745, 11043, 11279, 11397, 11421, 11313, 11161, 10902, 10587, 10237, 9866, 9442, 9030, 9250, 8500, 7950, 7673, 7399, 7119, 6927, 6800, 6700, 6389, 6215, 6034, 5868, 5710, 5538, 5369, 5221, 5052, 4882, 4703, 4509, 4500, 4166, 4006, 3836, 3685, 3544, 3412, 3500, 3500, 3019, 2905, 2801, 2687, 2585, 2600, 2500, 2294, 2211, 2137, 2069, 2300, 2100, 1900, 1863, 1833, 1806, 1787, 1771, 1764, 1763, 1772, 1785, 1803, 1823, 1849, 1868, 1873, 1879, 1878, 1882, 1882, 1890, 1899, 1909, 1914, 1925, 1938, 1945, 1967, 1984, 2006, 2029, 2055, 2082, 2123, 2166, 2226, 2284, 2356, 2434, 2509, 2598, 2698, 2806, 2921, 3046, 3175, 3314, 3462, 3613, 3761, 3920, 4077, 4237, 4396, 4528, 4655, 4800, 4915, 5022, 5121, 5201, 5280, 5346, 5392, 5459, 5519, 5560, 5604, 5659, 5716, 5781, 5872, 5967, 6042, 6130, 6221, 6330, 6452, 6567, 6691, 6849, 7001, 7137)

inits_nuts = function(){
  list(eta_w = rep(runif(1,-2,-1.4),sample_weeksGR), 
    rho=runif(1,-0.99,0.99),
    reciprocal_phiD=runif(1,0.5,1.5),
    ifr=c(runif(1,0.0114,0.01145),runif(1,0.00797,0.00801), runif(1,0.01145,0.01155)),
    sigma=c(runif(1,0.1,4), runif(1,0.1,4), runif(1,0.1,4))
    #,c_tot = init_cases + runif(1,-5, 100)
    )
}

nuts_fit_CB = stan("single_type_CB_stan.stan",
                  data = cov_data, 
                  pars = parameters_nuts, 
                  init = inits, 
                  chains = 4, 
                  warmup = 2000, 
                  iter = 12000, 
                  thin=10, 
                  control=list(adapt_delta=0.99, max_treedepth=14), seed=4321)
####################################################################################################################
parameters_bugs = c("beta", "rho", "NBr", "ifr", "tau", "R_t",
               "eta", "eta_w", "c_tot", "muD",  "c_unrep",  
               "log_lik", "dev")

inits_bugs = function(){
  list(eta_w = rep(runif(1,-2,-1.4),sample_weeksGR),
       rho=runif(1,-0.99,0.99),
       NBr =1,
       ifr = c(runif(1,0.0114,0.01145),runif(1,0.00797,0.00801), runif(1,0.01145,0.01155)),
       tau =c (runif(1,0.2,2.5), runif(1,10,30), runif(1,15,90))
       #,c_tot = init_cases+ runif(1,-5, 100)
       )
}

mcmc_fit_CB = bugs(model.file = "single_type_CB_bugs.txt",
                   cov_data, 
                   inits = inits_bugs, 
                   parameters = parameters_bugs,
                   n.chains = 5, 
                   n.iter = 20000, 
                   n.thin=10, 
                   n.burnin = 5000,
                   debug=TRUE
                   )
