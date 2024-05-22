nadaray_watson <- function(x,X,Y,h,K = dnorm){
  Kx <- rbind(sapply(X, function(Xi) K((x - Xi) / h) / h))
  W <- Kx / rowSums(Kx)
  drop(W %*% Y)
}

cvNW <- function(X, Y, h, K = dnorm) {
  sum(((Y - nadaray_watson(x = X, X = X, Y = Y, h = h, K = K)) /
         (1 - K(0) / colSums(K(outer(X, X, "-") / h))))^2)
}

cross_val <- function(X, Y, K = dnorm, plot.cv = FALSE) {
  h <- seq(0,5,0.01)
  CV_error <- sapply(h, function(h) cvNW(X = X, Y = Y, h = h, K = K))
  h_min <- h[which.min(CV_error)]
  if (plot.cv) {
    plot(h, CV_error, type = "o")
    abline(v = h_min, col = 2, lwd = 2)
    title("CV error trend in relation to bandwidth")
  }
  return(h_min)
}

# Generate some data to test the implementation GENERATION 1
n <- 500
eps <- rnorm(n, sd = 1)
m <- function(x) x^3 * sin(x)
X <- rnorm(n, sd = 1)
Y <- m(X) + eps
x_grid <- seq(-10, 10, l = 500)

hCV <- cross_val(X = X, Y = Y, plot.cv = TRUE)
h_gen1 <- hCV
plot(X, Y)
title("x^3 * sin(x)")
rug(X, side = 1); rug(Y, side = 2)
lines(x_grid, m(x_grid), col = 1)
lines(x_grid, nadaray_watson(x = x_grid, X = X, Y = Y, h = hCV), col = 2)
legend("top", legend = c("True regression", "Nadaraya-Watson"),
       lwd = 2, col = 1:2)

# Generate some data to test the implementation GENERATION 2
eps <- rnorm(n, sd = 1)
m <- function(x) x^2 * cos(x)
X <- rnorm(n, sd = 1)
Y <- m(X) + eps
x_grid <- seq(-10, 10, l = 500)

hCV <- cross_val(X = X, Y = Y, plot.cv = TRUE)
h_gen2 <- hCV
plot(X, Y)
title("x^2 * cos(x)")
rug(X, side = 1); rug(Y, side = 2)
lines(x_grid, m(x_grid), col = 1)
lines(x_grid, nadaray_watson(x = x_grid, X = X, Y = Y, h = hCV), col = 2)
legend("top", legend = c("True regression", "Nadaraya-Watson"),
       lwd = 2, col = 1:2)

#Loading and setting boston data
data(Boston,package ="MASS")

X <- sort(Boston$lstat)
Y <- Boston$medv[order(Boston$lstat)] 
hCV <- cross_val(X = X, Y = Y, plot.cv = TRUE)
h_Boston = hCV
plot(X,Y)
title("Boston data:")
lines(X,nadaray_watson(x = X,X=X,Y=Y,h=hCV),col = 2)

#TODO Cercare un altro dataset

data(airquality)
airquality <- na.omit(airquality)
pairs(airquality)
X <- sort(airquality$Temp)
Y <- airquality$Ozone[order(airquality$Temp)]
hCV <- cross_val(X = X, Y = Y, plot.cv = TRUE)
h_airquality <- hCV
plot(X,Y)
title("NY Air quality data")
lines(X,nadaray_watson(x = X,X=X,Y=Y,h=hCV),col = 2)
