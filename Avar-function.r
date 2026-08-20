#Moment estimation for difference gamma distribution
sigma <- 5

#iterations
M <- 5000
#sample size
n_s <- c(10000,50000)
#alpha
steps <- 10
shape.true <- (1:steps) * 1.5 / (steps + 1)
#1/beta
scale.true <- 1                
A <- function(alpha){
  2 / sqrt(pi) *exp(lgamma(alpha + 0.5) - lgamma(alpha))
}
alpha_hat <- numeric(M)
AVar <- array(0,c(length(shape.true),length(n_s)))

i<- 0

for (n in n_s) {
  i<-i+1
  for(k in 1:length(shape.true)){
    for (m in 1:M){
    epsilon <- rgamma(n, shape = shape.true[k], scale = scale.true) + c(0,cumsum(rnorm(n-1, mean = 0, sd = sqrt(sigma^2/n))))
    #independent random variables version
    #Ed <- abs(epsilon[seq(1, n - 1, by = 2)] - epsilon[seq(2, n, by = 2)])
    #using all data 
    Ed<-abs(diff(epsilon))
    m1 <- mean(Ed)
    m2 <- mean(Ed^2)
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
        print("hi")
        alpha <- optimize(function(a) f(a)^2,
                          interval = c(1e-3, 1e3))$minimum
      }
      beta <- A(alpha) / m1
      c(alpha, beta)
    }, error = function(e) {
      c(NA, NA)})
    
    alpha_hat[m] <- est[1]
    }
    
  AVar[k,i] <- var(alpha_hat)*n
  print(k)
  }
}
  
varianz <- function (shape.true){
  ((8/3*pi*shape.true*exp(2*lgamma(shape.true) - 2*lgamma(shape.true + 0.5))) -
     15 +
     4*sqrt(3)*exp(lgamma(shape.true + 1/3) + lgamma(shape.true + 2/3) - 2*lgamma(shape.true + 0.5)) +
     4*sqrt(pi)*exp(lgamma(shape.true) - lgamma(shape.true + 0.5)) * (2*pbeta(1/3, shape.true, 2*shape.true) - 1) -
     2/shape.true)*(2*digamma(shape.true + 0.5) - 2*digamma(shape.true) - 1/shape.true)^(-2)
}

par(mgp = c(2.4, 1, 0))
plot(seq(0.00001,1.5,length.out = 200),varianz(seq(0.00001,1.5,length.out = 200)),
     type = "l",
     lwd = 2,
     xlab = "",
     ylab = "",
     cex.lab = 1.4,
     cex.axis = 1.2)

grid()
lines(shape.true,AVar[,2], col = "green3", type = "p", cex = 1.5, pch = 19)
lines(shape.true,AVar[,1], col = "red3", type = "p", pch = 19)
