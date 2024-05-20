# Function to simulate ZINB data
#' Function to simulate ZINB data
#'
#' @param seed
#' @param K numbers of taxa
#' @param n1 sample size of group1
#' @param n2 sample size of group2
#' @param p the propotion of the differential taxa
#' @param bias bias size, "small" or "large"
#' @param zi the propotion of zero inflation, default is 0.3
#' @param overdisp default is 3
#' @param confounder "TRUE" or "FALSE"
#'
#' @return
#' @export
#'
#' @examples
#' ZINB_sim(seed=1,K=50,n1=20,n2=20,p=0.1,bias ="small",zi=0.3,confounder =F)



ZINB_sim <- function(seed, K, n1, n2, p, bias = c("large", "small"), zi = 0.3, overdisp = 3, confounder = c("TRUE", "FALSE")) {
  # Total number of samples
  library(MASS)
  n <- n1 + n2
  x <- c(rep(0, n1), rep(1, n2))
  w <- c(0.6, 0.3, 0.1)

  # Generate baseline beta0, ebeta0 means exp(beta0)
  set.seed(seed)
  index <- sample(c(1, 2, 3), K, replace = TRUE, prob = w)
  ebeta0 <- rep(NA, K)
  set.seed(seed)
  ebeta0[which(index == 1)] <- rgamma(length(which(index == 1)), shape = 50, rate = 1)
  ebeta0[which(index == 2)] <- rgamma(length(which(index == 2)), shape = 200, rate = 1)
  ebeta0[which(index == 3)] <- rgamma(length(which(index == 3)), shape = 10000, rate = 1)
  beta0 <- log(ebeta0)

  # Generate H
  set.seed(seed)
  H <- rep(0, K)
  C1 <- 1:floor(K * p)
  wt <- runif(1, 0, 1)
  C2 <- sample(C1, wt * length(C1), replace = FALSE)
  H[C1] <- 1

  # Generate delta, edelta means exp(delta)
  set.seed(seed)
  edelta <- rep(1, K)
  edelta[C1] <- runif(length(C1), 1.5, 10)
  edelta[C2] <- 1 / runif(length(C2), 1.5, 10)
  delta <- log(edelta)

  # Generate d
  if (bias == "small") {
    d_bias <- mean(abs(delta[C1])) / 2
  } else {
    d_bias <- mean(abs(delta[C1])) * 5
  }
  set.seed(seed)
  d <- c(rnorm(n1, 0, 1), rnorm(n2, d_bias, 1))

  # Generate Y, lY means log(Y)
  lY <- matrix(NA, n, K)
  for (i in 1:n) {
    lY[i, ] <- d[i] + beta0 + x[i] * delta
  }

  # Confounders
  if (confounder) {
    set.seed(seed)
    Gamma <- cbind(runif(n, 0, 1), rbinom(n, 1, 0.5))
    beta <- rbind(rep(-1, K), rep(1, K))
    lY <- lY + Gamma %*% beta
    rownames(Gamma) <- paste0("sub", seq(n))
  } else {
    Gamma <- NULL
  }
  Y <- exp(lY)

  # Generate Z
  set.seed(seed)
  eta <- zi
  Z <- matrix(0, n, K)
  set.seed(seed)
  for (j in 1:n) {
    Z[j, ] <- rbinom(K, size = 1, prob = eta)
  }

  # Generate C
  C <- matrix(NA, ncol = K, nrow = n)
  set.seed(seed)
  for (i in 1:n) {
    for (k in 1:K) {
      C[i, k] <- rnegbin(1, Y[i, k], theta = overdisp)
    }
  }
  C[Z == 1] <- 0
  C <- as.data.frame(C)
  Z <- as.data.frame(Z)

  # Prepare outputs
  rownames(Z) <- paste0("sub", seq(n))
  colnames(Z) <- paste0("taxon", seq(K))
  rownames(C) <- paste0("sub", seq(n))
  colnames(C) <- paste0("taxon", seq(K))
  names(H) <- paste0("taxon", seq(K))
  names(x) <- paste0("sub", seq(n))
  names(delta) <- paste0("taxon", seq(K))
  names(d) <- paste0("sub", seq(n))

  test.data <- list(C, delta, x, Gamma, H, d, Z, eta)
  names(test.data) <- c("count", "delta", "grp", "confounder", "diff.taxa", "d", "zero.pos", "zero.prop")
  return(test.data)
}
