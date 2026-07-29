#Moment estimation for difference gamma distribution
sigma <- 5

#p asymptotic confidence interval
p<- 0.1
#iterations
M <- 50000
#sample size
n <- 50000
#alpha
shape.true <- 0.5
#1/beta	 		
scale.true <- 5

A <- function(alpha){
  2 / sqrt(pi) *exp(lgamma(alpha + 0.5) - lgamma(alpha))
}

alpha_hat <- numeric(M)
beta_hat  <- numeric(M)
conf_alpha <- array(0,c(M,2))

alpha_hat_2 <- numeric(M)
beta_hat_2  <- numeric(M)

for (m in 1:M){
  epsilon <- rgamma(n, shape = shape.true, scale = scale.true) + c(0,cumsum(rnorm(n-1, mean = 0, sd = sqrt(sigma^2/n))))
  
  Ed<-abs(diff(epsilon))
  m1 <- mean(Ed)
  m2 <- mean(Ed^2)
  m4 <- mean(Ed^4)
  
  f <- function(alpha) {
    (m2 * A(alpha)^2) / (2 * alpha) - m1^2
  }
  grid <- exp(seq(log(1e-3), log(1e3), length.out = 200))
  vals <- sapply(grid, f)
  
  idx <- which(diff(sign(vals)) != 0)
  
  est <- tryCatch({
    if (length(idx) > 0) {
      lower <- grid[idx[1]]
      upper <- grid[idx[1] + 1]
      alpha <- uniroot(f, lower = lower, upper = upper,tol = 1e-11)$root
    } else {
      alpha <- optimize(function(a) f(a)^2,
                        interval = c(1e-3, 1e3))$minimum
    }
    beta <- A(alpha) / m1
    c(alpha, beta)
  }, error = function(e) {
    c(NA, NA)})
  
  alpha_hat[m] <- est[1]
  beta_hat[m]  <- est[2]
  
  alpha_hat_2[m] <- (3*m2^2)/(m4 - 3*m2^2)
  beta_hat_2[m] <- sqrt(2*alpha_hat_2[m]/m2)
  conf_alpha[m,1] <- est[1] + qnorm(p/2)/sqrt(n)* sqrt((8/3*pi*est[1]*exp(2*lgamma(est[1]) - 2*lgamma(est[1] + 0.5))) -
                                                   15 +
                                                   4*sqrt(3)*exp(lgamma(est[1] + 1/3) + lgamma(est[1] + 2/3) - 2*lgamma(est[1] + 0.5)) +
                                                   4*sqrt(pi)*exp(lgamma(est[1]) - lgamma(est[1] + 0.5)) * (2*pbeta(1/3, est[1], 2*est[1]) - 1) -
                                                   2/est[1])*(2*digamma(est[1] + 0.5) - 2*digamma(est[1]) - 1/est[1])^(-1)
  conf_alpha[m,2] <- est[1] + qnorm(1-p/2)/sqrt(n)* sqrt((8/3*pi*est[1]*exp(2*lgamma(est[1]) - 2*lgamma(est[1] + 0.5))) -
                                                   15 +
                                                   4*sqrt(3)*exp(lgamma(est[1] + 1/3) + lgamma(est[1] + 2/3) - 2*lgamma(est[1] + 0.5)) +
                                                   4*sqrt(pi)*exp(lgamma(est[1]) - lgamma(est[1] + 0.5)) * (2*pbeta(1/3, est[1], 2*est[1]) - 1) -
                                                   2/est[1])*(2*digamma(est[1] + 0.5) - 2*digamma(est[1]) - 1/est[1])^(-1)
}


#histograms of alpha_hat and alpha_hat_2
breaks <- seq(
  min(c(alpha_hat, alpha_hat_2)),
  max(c(alpha_hat, alpha_hat_2)),
  length.out = 100
)

hist(alpha_hat_2,
     breaks = breaks,
     freq = FALSE,
     col = rgb(1, 0, 0, 0.4),
     border = "white",
     xlab = "",
     ylab = "",
     xlim = c(0.36,0.64),
     ylim = c(0,40),
     main = "",
     cex.axis = 1.7
     )

hist(alpha_hat,
     breaks = breaks,
     freq = FALSE,
     col = rgb(0, 0, 1, 0.4),
     border = "white",
     add = TRUE)
curve(dnorm(x, mean = shape.true, sd = sqrt(((8/3*pi*shape.true*exp(2*lgamma(shape.true) - 2*lgamma(shape.true + 0.5))) -
                                                    15 +
                                                    4*sqrt(3)*exp(lgamma(shape.true + 1/3) + lgamma(shape.true + 2/3) - 2*lgamma(shape.true + 0.5)) +
                                                    4*sqrt(pi)*exp(lgamma(shape.true) - lgamma(shape.true + 0.5)) * (2*pbeta(1/3, shape.true, 2*shape.true) - 1) -
                                                    2/shape.true)*(2*digamma(shape.true + 0.5) - 2*digamma(shape.true) - 1/shape.true)^(-2)/n)), 
      add = TRUE, 
      col = "black", 
      lwd = 2,
      lty = 1)



#histogram of beta_hat
hist(beta_hat, main = bquote(
  "Histogram of " ~ hat(beta) ~
    " with EW =" ~ .(round(mean(beta_hat),3)) ~
    " and SD =" ~ .(round(sd(beta_hat), 3))),
  freq = FALSE,
  breaks = 100,
  xlab = bquote(hat(beta)))


