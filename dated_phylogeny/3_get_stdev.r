# Rscript to get stdev (in real space) from 95% HPD:

#1 Get sdlog (stdev in LOG space) by TRYING DIFFERENT VALUES OF STDEV until you get the desired HPD

library (TeachingDemos)

my_mean = 2.084
stdev = try_values_here

real_mean_HPD <- hpd(qlnorm, meanlog=log(my_mean)-stdev^2/2, sdlog=stdev)

print(real_mean_HPD)

# 2 Once you have the sdlog that gives you the desired HPD, transform into stdev

# To do so, first you have to CALCULATE THE meanlog (or TRY DIFFERENT VALUES until it gives you the real mean)

# Here you define function to calculate stdev i mean from sdlog and meanlog

logno_moments <- function(meanlog, sdlog) {m <- exp(meanlog + (1/2)*sdlog^2)
	s <- exp(meanlog + (1/2)*sdlog^2)*sqrt(exp(sdlog^2) - 1)
	return(list(mean = m, sd = s))}

# Try here values of meanlog until you get the desired mean, then keep the resulting stdev

meanlog <- try_values
sdlog <- write_previously_obtained_value
logno_moments(meanlog, sdlog)

#--------------------------------------------------------------------------------------------------------------------

# Potential way to use a normal 

library (TeachingDemos)

my_mean = 2.084
stdev = 0.36

real_mean_HPD <- hpd(qnorm, mean=my_mean, sd=stdev)

print(real_mean_HPD)
