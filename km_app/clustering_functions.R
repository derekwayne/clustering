mv_data <- function(mu_1, mu_2, Sigma_1, Sigma_2) {
  sample_1 <- mvrnorm(n=200, mu=mu_1, Sigma = Sigma_1)
  sample_2 <- mvrnorm(n=200, mu=mu_2, Sigma = Sigma_2)
  
  sample <- rbind(sample_1, sample_2)
  colnames(sample) <- c("X1", "X2")
  classes <- c(rep(1, 200), rep(2, 200))
  
  sample
}

# K means cost function
cost <- function(x, centroid, R) {
  sum((x - t(centroid) %*% t(R)) ** 2)
}

km_e_step <- function(x, centroid, A) {
  norm_vec <- function(y) sqrt(sum(y^2)) # euclidean norm
  for (i in 1:nrow(A)) {
    dist_1 <- norm_vec(x[,i] - centroid[1,])
    dist_2 <- norm_vec(x[,i] - centroid[2,])
    k <- which.min(c(dist_1, dist_2))
    # assign r[i,k]
    if (k == 1) {
      A[i,] <- c(1,0)
    } else {
      A[i,] <- c(0,1)
    }
  }
  A
}

km_m_step <- function(x, centroid, A) {
  C1 <- x %*% A[,1] / sum(A[,1])
  C2 <- x %*% A[,2] / sum(A[,2])
  centroid <- matrix(c(C1,C2),ncol=2,byrow = TRUE)
  centroid
}

km_ <- function(x, centroid, R, max.iter) {
  cost_old <- 1L
  cost_new <- 0
  i <- 0
  cost_vec <- c()
  while ((cost_new != cost_old) & (i < max.iter)) {
    R <- R <- km_e_step(x, centroid, R)
    cost_old <- cost(x, centroid, R)
    cost_vec <- c(cost_vec, cost_old)
    centroid <- km_m_step(x, centroid, R)
    cost_new <- cost(x, centroid, R)
    cost_vec <- c(cost_vec, cost_new)
    i <- i + 1
  }
  cluster <- apply(t(R), FUN = which.max, MARGIN = 2)
  r <- list(centroid=centroid, R=R, cluster=cluster, cost=cost_vec)
  r
}

# Numerically stable computation of mvnormal desnity
normal_density <- function(x, mean, sigma) {
  cf <- chol(sigma)
  tmp <- backsolve(cf, t(x) - mean, transpose = TRUE)
  rss <- colSums(tmp^2)
  lv <- -sum(log(diag(cf))) - 0.5 * ncol(x) * log(2 * 
                                                    pi) - 0.5 * rss
  exp(lv)
}

# EM

e_step <- function(x, theta) {
  mix.p <- c(theta[[1]], theta[[2]])
  mu_1 <- theta[[3]]; mu_2 <- theta[[4]]
  Sigma_1 <- theta[[5]]
  Sigma_2 <- theta[[6]]
  
  a <- mix.p[1] * normal_density(x, mean=mu_1, sigma=Sigma_1)
  b <- mix.p[2] * normal_density(x, mean=mu_2, sigma=Sigma_2)
  resp_1 <- a / (a+b)
  resp_2 <- b / (a+b)
  r <- list("loglik"=sum(log(a+b)), "responsibilities"=cbind(resp_1, resp_2))
  r
}

m_step <- function(x, resps) {
  mu_1.dagger <- as.vector(crossprod(resps[,1],x) / sum(resps[,1]))
  mu_2.dagger <- as.vector(crossprod(resps[,2],x) / sum(resps[,2]))
  
  Sigma_1.dagger <- t(resps[,1] * t(apply(x, 1, function(x) x - mu_1.dagger))) %*% 
    (resps[,1] * t(apply(x, 1, function(x) x - mu_1.dagger))) * 1/sum(resps[,1])
  Sigma_2.dagger <- t(resps[,2] * t(apply(x, 1, function(x) x - mu_2.dagger))) %*% 
    (resps[,2] * t(apply(x, 1, function(x) x - mu_2.dagger))) * 1/sum(resps[,2])
  
  pi_1.dagger <- sum(resps[,1]) / 400
  pi_2.dagger <- sum(resps[,2]) / 400
  
  r <- list(p1=pi_1.dagger, p2=pi_2.dagger, mu1=mu_1.dagger, mu2=mu_2.dagger,
            sig1=Sigma_1.dagger, sig2=Sigma_2.dagger)
  r
}

EM <- function(x, init, max.iter, tol=10e-4) {
  for (i in 1:max.iter) {
    if (i == 1) {
      E_ <- e_step(x, init)
      M_ <- m_step(x, E_$responsibilities)
      old.loglik <- E_$loglik
      loglik.vec <- E_$loglik
    } else {
      E_ <- e_step(x, M_)
      M_ <- m_step(x, E_$responsibilities)
      loglik.vec <- c(loglik.vec, E_$loglik)
      condition <- abs((E_$loglik - old.loglik))
      if (condition < tol) {
        break
      } else {
        old.loglik <- E_$loglik
      }
    }
  }
  list(p1=M_$p1, p2=M_$p2, mu1=M_$mu1, mu2=M_$mu2,
       sig1=M_$sig1, sig2=M_$sig2,
       responsibilities = E_$responsibilities, loglik = loglik.vec)
}