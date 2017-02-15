#this is for exercise 1, 'Quantifying uncertainty' section, 'basic', (B)

# Load the library
# you might have to install this the first time
library(mlbench)

# Load the data
ozone = data(Ozone, package='mlbench')

# Scrub the missing values
# Extract the relevant columns 
ozone = na.omit(Ozone)[,4:13]

y = ozone[,1]
x = as.matrix(ozone[,2:10])

# add an intercept
x = cbind(1,x)

# compute the estimator
betahat = solve(t(x) %*% x) %*% t(x) %*% y

# Fill in the blank - ChenfengHe
# use variance calculated from the sample data to represent the actual unknown variance
# rho^2 ~ RSS/(n-p-1)
rss = sum((y-x %*% betahat)^2)
rho_samp = rss/(length(y)-dim(x)[2]-1)
# only compare the diagnal of covariance matrix, that is the variance
betacov = sqrt(rho_samp * diag(solve(t(x) %*% x)))



# Now compare to lm
# the 'minus 1' notation says not to fit an intercept (we've already hard-coded it as an extra column)
lm1 = lm(y~x-1)

summary(lm1)
betacovlm = vcov(lm1)
sqrt(diag(betacovlm))
