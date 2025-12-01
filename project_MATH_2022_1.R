# MATH-2022-1
# Mvomo Eto Wilfried
# s22625


##################################################################################
# 2.  Plot the distribution Be(α, β) for different values of its shape parameters#
##################################################################################

set.seed(1)                       # for reproducibility

# We proceed by the different values.

# Set the alpha and beta parameter values and colors

parameters <- list(
  list(alpha = 2, beta = 5, color = "blue"),
  list(alpha = 3, beta = 7, color = "red"), 
  list(alpha = 1, beta = 1, color = "saddlebrown"), 
  list(alpha = 7, beta = 2, color = "black"),
  list(alpha = 10, beta = 10, color = "yellow"),
  list(alpha = 1, beta = 8, color = "green"),
  list(alpha = 50, beta = 100, color = "orange"),
  list(alpha = 0.5, beta = 0.6, color = "darkred")
  
)

# Generate a sequence of values from 0 to 1

x <- seq(0, 1, length.out = 100)


# Create a blank plot

plot(NULL, xlim = c(0, 1), ylim = c(0, 10), xlab = "x", ylab = "pdf", main = "Beta Distribution")


# Loop through each set of parameters and plot the beta distribution

for (param in parameters) {
  
  alpha <- param$alpha
  beta <- param$beta
  color <- param$color
  
  # Calculate the PDF of the beta distribution for each value of x
  
  y <- dbeta(x, alpha, beta)
  
  # Plot the beta distribution with the specified color
  
  lines(x, y, lwd = 2, col = color)
}

# Add a legend

legend("topright", legend = c("alpha = 2, beta = 5", "alpha = 3, beta = 7", "alpha = 1, beta = 1",
                              "alpha = 7, beta = 2", "alpha = 10, beta = 1", "alpha = 1, beta = 8",
                              "alpha = 50, beta = 100","alpha = 0.5, beta = 0.6"), lwd = 2, col = c("blue", "red","saddlebrown","black","yellow","green","orange","darkred"))



########################################################
# 3. Simulate 10, 000 values of a Be(4, 6) distribution#
########################################################


set.seed(1)                       # for reproducibility

n <- 10^4   # Number of random numbers to generate
alpha <- 4  # Alpha parameter
beta <- 6   # Beta parameter
x <- rbeta(n, alpha, beta)
xseq<- seq(0,1,0.01)
hist(x, freq = FALSE, xlim = c(0, 1), ylim = c(0, 3), col = "skyblue",
     xlab = "x", ylab = "Density", main = "Beta Distribution")
lines(xseq, dbeta(xseq,alpha, beta), col="red")
legend("topright", legend = "Be(2,6)", lwd = 2, col = "red")



###############################################################################
# 4. Implement this transformation method                                     #
# algorithm and apply it to simulate 10, 000 values of a Be(2, 6) distribution#
###############################################################################

set.seed(1)                       # for reproducibility

alpha <- 2  
beta <- 6 
N <- 10^4

start_time1 <- Sys.time()

Y <- NULL  # create an empty vector to store the values of Y

for (i in 1:N) {
  
  log_U <- NULL  
  
  for (j in 1:(beta + alpha)) {
    
    U <- runif(1) 
    log_U[j] <- log(U)  
    
  }
  
  Y[i] <- sum(log_U[1:alpha]) / sum(log_U[1:(beta + alpha)])  # compute Y
}

end_time1 <- Sys.time()

xseq<- seq(0,1,0.01)
hist(Y, freq = FALSE, xlim = c(0, 1), ylim = c(0,3), col = "skyblue",
     xlab = "x", ylab = "Density", main = "Beta Distribution (using transform method)")
lines(xseq, dbeta(xseq,alpha, beta), col="red")
legend("topright", legend = "Be(2,6)", lwd = 2, col = "red")


######################################################################
# 5. Use(5.a) to generate a beta random variable and simulate 10, 000#
# values from a Be(2, 6).                                            #
######################################################################

set.seed(1)                       # for reproducibility

alpha <- 2
beta <- 6
N <- 10^4

# Generate random values from Gamma distributions

start_time2 <- Sys.time()

Y_1 <- rgamma(N, shape = alpha, rate = 1)
Y_2 <- rgamma(N, shape = beta, rate = 1)

# Compute X = Y1 / (Y1 + Y2)

X <- Y_1/ (Y_1 + Y_2)

end_time2 <- Sys.time()

xseq<- seq(0,1,0.01)
hist(X, freq = FALSE, xlim = c(0, 1), ylim = c(0,3), col = "skyblue",
     xlab = "x", ylab = "Density", main = "Beta Distribution")
lines(xseq, dbeta(xseq,alpha, beta), col="red")
legend("topright", legend = "Be(2,6)", lwd = 2, col = "red")


###################################################################################
# 6. Comparison of both approaches according to their run time and their precision#                                                             #
###################################################################################


# Print runtime and precision first approach (answer 4)

runtime1 <- end_time1 - start_time1
print(runtime1)

precision1 <- sum(Y >= 0 & Y <= 1) / N
print(precision1)


# Print runtime and precision second approach (answer 5)

runtime2 <- end_time2 - start_time2
print(runtime2)

precision2 <- sum(X >= 0 & X <= 1) / N
print(precision2)


#######################################################################################
# 7.b Determine (using R) this limiting constant C used in the Accept-Reject algorithm#
# if the target distribution is Be(1.6, 5.8)                                          #                                                             #
#######################################################################################


library(gtools)

set.seed(1)                       
target_density <- function(x) {
  dbeta(x, 1.6,  5.8)
}

instrumental_density <- function(x) {
  dunif(x, 0, 1)
}

x <- seq(0, 1, length = 1000)  # Range of x values

ratios <- target_density(x) / instrumental_density(x)

C <- max(ratios)     # set the limiting constant C

print(C)


######################################################################################################
# 7.c Implement the algorithm and perform 10, 000 simulations to obtain values from a Be(1.6, 5.8)   #
# distribution (using pdf of the uniform distribution) and 7.d Create a plot that shows all simulated#
# values, distinguishing between the accepted (in green)                                             # 
# and rejected values (in red).                                                                      #
#                                                                                                    #
######################################################################################################

library(AR)

set.seed(1)  
start_time3 <- Sys.time()
simulation_beta_unif <-  AR.Sim(n = 10^4, 
                                f_X = function(y){dbeta(y,1.6,5.8)},
                                Y.dist = "unif", Y.dist.par = c(0,1),
                                Rej.Num = TRUE,
                                Rej.Rate = TRUE, 
                                Acc.Rate = TRUE)
end_time3<- Sys.time()

simulation_beta_unif
runtime3 <- end_time3 - start_time3


##################
# 7.e Comparison #                                    
##################


C <- 0.32786
C_general <- 1/ 3.06
difference <- C - C_general
print(difference)           # the both results are sensitively equal



#####################################
# 8.b Redo 7.d using pdf of Be(1,5) #                                                             
#####################################

set.seed(1) 
start_time4 <- Sys.time()
simulation_beta_beta <- AR.Sim(n = 10^4,
                               f_X = function(y){dbeta(y, 1.6, 5.8)},
                               Y.dist = "beta", Y.dist.par = c(1, 5),
                               Rej.Num = TRUE,
                               Rej.Rate = TRUE,
                               Acc.Rate = TRUE)
end_time4 <- Sys.time()

simulation_beta_beta
runtime4 <- end_time4 - start_time4

print(runtime3)
print(runtime4)


######################################################################################
# 10.d ) Use the Monte Carlo approach to estimate the shape parameters of a Be(25,7) #
# distribution (set n = 10, 000) (Method of Moments (MOM) estimation)                #                     #
######################################################################################   

set.seed(1)

# Set the number of samples for Monte Carlo simulation
n <- 10^4

# Generate n random samples from the Beta distribution
alpha <- 25
beta <- 7
samples <- rbeta(n, alpha, beta)

# Calculate the empirical mean

mean_X <- mean(samples)

# Calculate the empirical variance

var_X <- var(samples)

# Estimate alpha and beta

hat_alpha_mom <- mean_X * ((1 - mean_X) / var_X - 1)
hat_beta_mom <- (1 - mean_X) * ((1 - mean_X) / var_X - 1)

# Print the estimated parameters
cat("Estimated alpha:", hat_alpha_mom, "\n")
cat("Estimated beta:", hat_beta_mom, "\n")


#######################################################################################################
# 11.acUse the previous implementation for the MOM and use the function ebeta contained in the package#
# EnvStats for the MLE estimator for the beta distribution to create histograms illustrating the      #
# distribution of each of the shape parameters                                                        #
#######################################################################################################

library(EnvStats)

set.seed(1)
N<- 10^4 
nSim <- 100  
alpha <- 25
beta <- 7


hat_alpha_mom <- NULL;hat_beta_mom <- NULL
hat_alpha_mle <- NULL;hat_beta_mle <- NULL


for (i in 1:N){    
  
  samples <- rbeta(nSim, alpha, beta)
  
  # Method of Moments (MOM) estimation
  
  mean_X <- mean(samples)
  var_X <- var(samples)
  hat_alpha_mom[i] <- mean_X * ((1 - mean_X) / var_X - 1)
  hat_beta_mom[i] <- (1 - mean_X) * ((1 - mean_X) / var_X - 1)
  
  
  # Method of Moments (MLE) estimation
  
  mle <- ebeta(samples)
  hat_alpha_mle[i] <- mle$parameters[1]    # estimated parameters
  hat_beta_mle[i] <-mle$parameters[2]
  
}

# create histograms illustrating the distribution of each of the shape parameters

par(mfrow = c(2, 2))
hist(hat_alpha_mom, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MOM)")
hist(hat_alpha_mle, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MLE)")
hist(hat_beta_mom, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated beta (MOM)")
hist(hat_beta_mle, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MLE)")



#################################################################################################
# 11.b Evaluate the efficiency of each method in terms of bias and variance. Compare the results#                                                      #
#################################################################################################


# Calculate bias and variance for alpha estimates

alpha_mom <- mean(hat_alpha_mom)
alpha_mle <- mean(hat_alpha_mle)

bias_alpha_mom <- alpha_mom - alpha
bias_alpha_mle <- alpha_mle - alpha

var_alpha_mom <- var(hat_alpha_mom)
var_alpha_mle <- var(hat_alpha_mle)

# Print bias and variance for alpha estimates

cat("Bias of alpha (MOM):", bias_alpha_mom, "\n")
cat("Bias of alpha (MLE):", bias_alpha_mle, "\n")
cat("Variance of alpha (MOM):", var_alpha_mom, "\n")
cat("Variance of alpha (MLE):", var_alpha_mle, "\n")


# Calculate bias and variance for beta estimates

beta_mom <- mean(hat_beta_mom)
beta_mle <- mean(hat_beta_mle)

bias_beta_mom <- beta_mom - beta
bias_beta_mle <- beta_mle - beta

var_beta_mom <- var(hat_beta_mom)
var_beta_mle <- var(hat_beta_mle)

# Print bias and variance for beta estimates

cat("Bias of beta (MOM):", bias_beta_mom, "\n")
cat("Bias of beta (MLE):", bias_beta_mle, "\n")
cat("Variance of beta (MOM):", var_beta_mom, "\n")
cat("Variance of beta (MLE):", var_beta_mle, "\n")


###########################################################################
# 11.c Illustrate the asymptotic normality of the obtained MLE estimators.#                                                      #
###########################################################################


xseq <- seq(0,60,0.01)

par(mfrow = c(2, 2))
hist(hat_alpha_mom, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MOM)")
lines(xseq, dnorm(xseq,mean(hat_alpha_mom),sd(hat_alpha_mom)), col="red")
legend("topright", legend = "Normal's pdf", lwd = 2, col = "red")
hist(hat_alpha_mle, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MLE)")
lines(xseq, dnorm(xseq,mean(hat_alpha_mle),sd(hat_alpha_mle)), col="red")
legend("topright", legend = "Normal's pdf", lwd = 2, col = "red")
hist(hat_beta_mom, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated beta (MOM)")
lines(xseq, dnorm(xseq,mean(hat_beta_mom),sd(hat_beta_mom)), col="red")
legend("topright", legend = "Normal's pdf", lwd = 2, col = "red")
hist(hat_beta_mle, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "Estimated aplha (MLE)")
lines(xseq, dnorm(xseq,mean(hat_beta_mle),sd(hat_beta_mle)), col="red")
legend("topright", legend = "Normal's pdf", lwd = 2, col = "red")













