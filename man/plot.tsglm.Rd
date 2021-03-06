% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.tsglm.R
\name{plot.tsglm}
\alias{plot.tsglm}
\title{Plot Effective Image Pixels for A \kbd{"tsglm"} Object}
\usage{
\method{plot}{tsglm}(x, method = p.adjust.methods, background = NULL, ...)
}
\arguments{
\item{x}{an object of class \kbd{"tsglm"}.}

\item{method}{p-value correction method. See \code{\link[stats]{p.adjust}}.}

\item{background}{an image data that used as the background of the effectiveness markers. 
If \code{background = NULL}, the background color shows the effect size of the each pixel.}

\item{...}{further arguments passed to the \code{\link[graphics]{image}} function.}
}
\description{
\kbd{plot} method for self-defined class \kbd{"tsglm"}.
}
\examples{
# Predefined function: sum of hadamard product in each array
`\%i\%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))

# Simulation data
n <- 1000 # number of observations
n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
beta_True <- rep(1, n_d)
B_True <- c(1,1,1)\%*\%t(rnorm(n_G)) + c(0, .5, .5)\%*\%t(rnorm(n_G))
B_True <- B_True / 10
W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))

# Binomial Responses
p_B <- exp(W\%*\%beta_True + X\%i\%B_True); p_B <- p_B/(1+p_B)
y_B <- rbinom(n, 1, p_B)
DATA_B <- list(y = y_B, W = W, X = X)

# Binomial Model
result_B <- TRtest.omics(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
opt = 1, max_ite = 100, tol = 10^(-7) )

# Plot the effect size of the effective pixels
plot(result_B, "fdr")

# Plot the effective pixels with data image as the background
x0 <- DATA_B$X[,,which(DATA_B$y == 0)]
m0 <- matrix(0, dim(DATA_B$X)[1], dim(DATA_B$X)[2])
for (i in 1:dim(x0)[3]) m0 <- m0 + x0[,,i]/dim(x0)[3]
plot(result_B, "fdr", m0, col = gray(seq(0, 1, 0.05)))

}
\seealso{
\code{\link{TRtest.omics}, \link{drawpixelmarks}}
}
\author{
Ping-Yang Chen
}
