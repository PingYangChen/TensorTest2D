#' Predict by Second-order Tensor Generalized Regression 
#'
#' \kbd{predict} method for self-defined class \kbd{"tsglm"}.
#'
#' @details \kbd{predict.tsglm} is 
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param object an object of class \kbd{"tsglm"}.
#' @param newdata a data array in which to look for variables with which to predict. 
#' @param type the type of prediction required.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \kbd{predict.tsglm} returns xxx.
#'
#' @seealso \code{\link{TRtest.omics}, \link{summary.tsglm}}
#'
#' @examples
#' # Predefined function: sum of hadamard product in each array
#' `%i%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
#'
#' # Simulation data
#' n <- 1000 # number of observations
#' n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
#' n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
#' beta_True <- rep(1, n_d)
#' B_True <- c(1,1,1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
#' B_True <- B_True / 10
#' W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
#' X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))
#' ## Regression
#' y_R<- as.vector(W%*%beta_True + X%i%B_True + rnorm(n))
#' DATA_R <- list(y = y_R, X = X, W = W)
#' ## Binomial
#' p_B <- exp(W%*%beta_True + X%i%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#' ## Poisson
#' p_P <- exp(W%*%beta_True + X%i%B_True)
#' y_P <- rpois(n, p_P)
#' y_P[which(y_P > 170)] <- 170 # If y_P > 170, factorial(y_P) == inf.
#' DATA_P <- list(y = y_P, W = W, X = X)
#'
#' # Execution
#' ## Regression
#' result_R <- TRtest.omics(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' predict(result_R, DATA_R$X)
#'
#' ## Binomial
#' result_B <- TRtest.omics(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_B)
#'
#' ## Poisson
#' result_P <- TRtest.omics(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_P)
#'
#'
#' @author Ping-Yang Chen
#'
#' @export
predict.tsglm <- function(object, newx = NULL, type = c("link", "response"), ...){
  
  `%i%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
  eta <- rep(1, dim(newx)[3]) %*% object$b_EST + newx %i% object$B_EST
  
  if(object$family == "gaussian"){
    return(eta)
  }else if(object$family == "binomial"){
    if(type == "link") {
      return(eta)
    }else if(type == "response") {
      mu <- exp(eta)/(1 + exp(eta))
      return(mu)
    }
  }else if(object$family == "poisson"){
    if(type == "link") {
      return(eta)
    }else if(type == "response") {
      mu <- exp(eta)
      return(mu)
    }
  }
}



#' Plots
#'
#' \kbd{plot} method for self-defined class \kbd{"tsglm"}.
#'
#' @details \kbd{plot.tsglm} is 
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param object an object of class \kbd{"tsglm"}.
#' @param X an image data. 
#' @param ... further arguments passed to or from other methods.
#'
#' @return \kbd{plot.tsglm} returns xxx.
#'
#' @seealso \code{\link{TRtest.omics}, \link{summary.tsglm}}
#'
#' @examples
#' # Predefined function: sum of hadamard product in each array
#' `%i%` <- function(X, B) sapply(1:dim(X)[3], function(i) sum(X[,,i]*B))
#'
#' # Simulation data
#' n <- 1000 # number of observations
#' n_P <- 3; n_G <- 64 # dimension of 3-D tensor variables.
#' n_d <- 1 # number of numerical variable, if n_d == 1,  numerical variable equals to intercept.
#' beta_True <- rep(1, n_d)
#' B_True <- c(1,1,1)%*%t(rnorm(n_G)) + c(0, .5, .5)%*%t(rnorm(n_G))
#' B_True <- B_True / 10
#' W <- matrix(rnorm(n*n_d), n, n_d); W[,1] <- 1
#' X <- array(rnorm(n*n_P*n_G), dim=c(n_P, n_G, n))
#' ## Regression
#' y_R<- as.vector(W%*%beta_True + X%i%B_True + rnorm(n))
#' DATA_R <- list(y = y_R, X = X, W = W)
#' ## Binomial
#' p_B <- exp(W%*%beta_True + X%i%B_True); p_B <- p_B/(1+p_B)
#' y_B <- rbinom(n, 1, p_B)
#' DATA_B <- list(y = y_B, W = W, X = X)
#' ## Poisson
#' p_P <- exp(W%*%beta_True + X%i%B_True)
#' y_P <- rpois(n, p_P)
#' y_P[which(y_P > 170)] <- 170 # If y_P > 170, factorial(y_P) == inf.
#' DATA_P <- list(y = y_P, W = W, X = X)
#'
#' # Execution
#' ## Regression
#' result_R <- TRtest.omics(y = DATA_R$y, X = DATA_R$X, W=NULL, n_R = 1, family = "gaussian",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' predict(result_R, DATA_R$X)
#'
#' ## Binomial
#' result_B <- TRtest.omics(y = DATA_B$y, X = DATA_B$X, W=NULL, n_R = 1, family = "binomial",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_B)
#'
#' ## Poisson
#' result_P <- TRtest.omics(y = DATA_P$y, X = DATA_P$X, W=NULL, n_R = 1, family = "poisson",
#' opt = 1, max_ite = 100, tol = 10^(-7) )
#' summary(result_P)
#'
#'
#' @author Ping-Yang Chen
#'
#' @export
plot.tsglm <- function(object, X, p.adjust.methods){
  
  adjp <- matrix(p.adjust(as.vector(object$B_PV), method = p.adjust.methods),
                 nrow(object$B_PV), ncol(object$B_PV))
  
  marks <- object$B_EST*(adjp < 0.05)
  effectivepixelmarks(X, marks)
  
}


#' drawpixelmarks
#'
#' \kbd{drawpixelmarks} .
#'
#' @details \kbd{drawpixelmarks} is 
#'
#' @importFrom stats coefficients glm pnorm pt rnorm symnum
#'
#' @param img an image data.
#' @param marks an 0, -1, 1 marking on the image data. 
#'
#' @return \kbd{plot.tsglm} returns xxx.
#'
#' @seealso \code{\link{TRtest.omics}, \link{plot.tsglm}}
#'
#' @examples
#' # 
#'
#' @author Ping-Yang Chen
#'
#' @export
drawpixelmarks <- function(img, marks){
  
  image(1:ncol(img),1:nrow(img),img,
        col=gray(seq(0,1,0.05)),xlab="",ylab="",axes=FALSE)
  
  for (i in 1:nrow(img)) {
    for (j in 1:ncol(img)) {
      bcol <- ifelse(marks[i,j] == 0, NA, ifelse(marks[i,j] < 0, "blue", "red"))
      rect(i-.5, j-.5, i+.5, j+.5, col = NA, border = bcol, lwd = 2)
    }
  }
  
}


