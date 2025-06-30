

#' Poisson mixture data generator
#'
#' This generates a poisson mixture data
#' @export
#' @param M number of rows to generate
#' @param x covariate
#' @param a vector of 'a' components
#' @param b vector of 'b' components
#' @param pi_a prior probability of 'a'
#' @param pi_b prior probability of 'b'
#' @examples
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'            pi_a = c(1), pi_b = c(0.8, 0.2))
#' head(data)
generate.data <- function(M, x, a, b, pi_a, pi_b){
  pi.v = c(t(outer(pi_a, pi_b, "*")))
  nsamp = length(x)
  I = length(a)
  J = length(b)
  K = I*J
  b0 = b[1]
  p_mat <- matrix(NA, nrow = K, ncol = 2)
  k = 1
  for (i in 1:I) {
    for (j in 1:J) {
      param <- c(a[i], b[j])
      p_mat[k, ] <- param
      k <- k + 1
    }
  }
  bees <- p_mat[, 2]  # means of the two normal distributions
  sds <- c(0.9, 0.5)    # standard deviations of the two normal distributions

  ## Generate the component labels for each gene
  comps <- sample(1:K, M, replace = TRUE, prob = pi.v)
  sub0 = seq(which.min(abs(b)), (I*J), J)
  comp0 <- which(comps %in% sub0)
  complnt_0 <- setdiff(seq_along(comps), comp0)

  # 3. Generate RNA-seq counts using the log-linear model
  seq_data <- matrix(nrow = M, ncol = nsamp)
  mu_mat <- matrix(nrow = M, ncol = nsamp)

  for (i in comp0) {
    # Select the intercept and coefficients for the current gene based on the component
    chosen_int <- p_mat[comps[i], 1]
    chosen_coef <- p_mat[comps[i], 2]


    # Compute the log-linear model for each sample
    log_mu <- chosen_int + chosen_coef*x

    # Convert log mean to the actual mean (Poisson parameter)
    mu <- exp(log_mu)
    mu_mat[i, ] <- mu
    # Generate Poisson-distributed counts for each sample
    seq_data[i, ] <- rpois(nsamp, mu)
  }
  for (i in complnt_0 ) {
    # Select the intercept and coefficients for the current gene based on the component
    chosen_int <- p_mat[comps[i], 1]
    chosen_coef <- p_mat[comps[i], 2]

    # Compute the log-linear model for each sample
    log_mu <- chosen_int + chosen_coef*x

    # Convert log mean to the actual mean (Poisson parameter)
    mu <- exp(log_mu)
    mu_mat[i, ] <- mu
    # Generate Poisson-distributed counts for each sample
    seq_data[i, ] <- rpois(nsamp, mu)
  }
  return(seq_data)
}

#' Starting values for EM algorithm
#'
#' Provides starting values to implement the EM algorithm by running a k-means clustering on the GLM estimates of each feature
#' @export
#' @param Y count data
#' @param I number of 'a' components
#' @param J number of 'b' components
#' @param x covariate
#' @examples
#' # Generate data
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'             pi_a = c(1), pi_b = c(0.8, 0.2))
#'
#' # Choose starting values
#' sv <- starting.values(data, x = c(0,0,0,0,1,1,1,1))
#' sv
#'
starting.values = function(Y, x, I = 1, J = 2){
  cnnt <- Y
  s = colSums(cnnt)/sum(cnnt)
  off.mat <- matrix(nrow = nrow(cnnt), ncol = length(x))
  for (i in 1:nrow(off.mat)) {
    off.mat[i, ] = s
  }
  A <- numeric(0)
  B <- numeric(0)
  for (i in 1:nrow(cnnt)) {
    model <- glm(cnnt[i, ] ~ x, family = "poisson", offset = log(off.mat[i, ]))
    A[i] <- model$coefficients[1]
    B[i] <- model$coefficients[2]
  }

  ##load data
  a.df <- data.frame(a = A)

  #perform k-means clustering with k clusters
  kma <- kmeans(a.df, centers = I, nstart = 25)

  #results of final k-means model
  a.df1 = as.data.frame(cbind(a.df, a.df))
  a.df1 = data.frame(a.df1, cluster = as.factor(kma$cluster))
  a.sv = aggregate(a ~ cluster, data = a.df1, mean)
  a.sv$prop = as.numeric(prop.table(table(a.df1$cluster)))
  a = a.sv[order(a.sv$a), 2]
  pi_a = a.sv[order(a.sv$a), 3]

  #load data
  b.df <- data.frame(b = B)
  kmb <- kmeans(b.df, centers = J, nstart = 25)
  b.df1 = as.data.frame(cbind(b.df,b.df))
  b.df1 = data.frame(b.df1, cluster = as.factor(kmb$cluster))
  b.sv = aggregate(b ~ cluster, data = b.df1, mean)
  b.sv$prop = as.numeric(prop.table(table(b.df1$cluster)))
  closest_index <- which.min(abs(b.sv[, 2]))
  remaining_indices <- setdiff(seq_len(nrow(b.sv)), closest_index)  # Exclude closest index
  ordered_indices <- c(closest_index, remaining_indices[order(b.sv[remaining_indices, 2])])
  b.sv <- b.sv[ordered_indices, ]
  b = b.sv[, 2]
  b[1] = 0
  pi_b = b.sv[, 3]

  pi.v = c(t(outer(pi_a, pi_b, "*")))
  s.vals = list(a = a, b = b, pi.v = pi.v, pi_b = pi_b, pi_a = pi_a, A = A, B = B)
  return(s.vals)
}


loglikelihood = function(Y, x, pi.v, a, b){
  s = colSums(Y)/sum(Y)
  M = nrow(Y)
  N = length(x)
  I = length(a)
  J = length(b)
  K = I*J
  ij_grid <- expand.grid(j = 1:J, i = 1:I)  # inner loop is j, outer is i
  mu_block <- apply(ij_grid, 1, function(idx) {
    i <- idx["i"]
    j <- idx["j"]
    exp(log(s) + a[i] + b[j]*x)
  })
  muij <- do.call(rbind, replicate(M, t(mu_block), simplify = FALSE))
  gpi.v <- rep(pi.v, times = M)
  rep.times <- rep(K, M)
  # Replicate each row according to replication_times
  rep.rows <- unlist(
    lapply(seq_len(nrow(Y)), function(i) {
      rep(list(Y[i, , drop = FALSE]), rep.times[i])
    }),
    recursive = FALSE
  )
  # Combine the replicated rows into a matrix
  bigY <- do.call(rbind, rep.rows)
  pvc = c()
  for(i in 1:nrow(muij)){
    pvc[i] = log(gpi.v[i]) + sum(dpois(bigY[i,], lambda = muij[i,], log = T))
  }
  num.pmat = matrix(pvc, nrow = M, byrow = T)
  a = apply(num.pmat, 1, max)
  a <- apply(num.pmat, 1, max)  # max of each row
  lprob <- a + log(rowSums(exp(num.pmat - a)))
  return(sum(lprob))
}


E <- function(Y, x, pi.v, a, b){
  s = colSums(Y)/sum(Y)
  M = nrow(Y)
  N = length(x)
  I = length(a)
  J = length(b)
  K = I*J
  ij_grid <- expand.grid(j = 1:J, i = 1:I)  # inner loop is j, outer is i
  mu_block <- apply(ij_grid, 1, function(idx) {
    i <- idx["i"]
    j <- idx["j"]
    exp(log(s) + a[i] + b[j]*x)
  })
  muij <- do.call(rbind, replicate(M, t(mu_block), simplify = FALSE))
  gpi.v <- rep(pi.v, times = M)
  rep.times <- rep(K, M)
  # Replicate each row according to replication_times
  rep.rows <- unlist(
    lapply(seq_len(nrow(Y)), function(i) {
      rep(list(Y[i, , drop = FALSE]), rep.times[i])
    }),
    recursive = FALSE
  )
  # Combine the replicated rows into a matrix
  bigY <- do.call(rbind, rep.rows)
  pvc = c()
  for(i in 1:nrow(muij)){
    pvc[i] = log(gpi.v[i]) + sum(dpois(bigY[i,], lambda = muij[i,], log = T))
  }
  num.pmat = matrix(pvc, nrow = M, byrow = T)
  a <- apply(num.pmat, 1, max)  # max of each row
  lprob <- a + log(rowSums(exp(num.pmat - a)))
  num = num.pmat - lprob
  e = exp(num)
  return(e)
}



split_tau <- function(mat, small_cols) {
  # The code below produces a list where the number of elements in the list is I = length(a).
  # This is to help when we want to pick out tau_{mij} for each i = 1,...,I
  # Initialize a list to store the smaller matrices
  small_matrices <- list()

  # Loop through the columns
  for (i in seq(1, ncol(mat), small_cols)) {
    # Calculate the ending column index for the subset
    end_col <- min(i + small_cols - 1, ncol(mat))

    # Subset the original matrix
    small_matrix <- mat[, i:end_col, drop = FALSE]  # Ensure that matrix remains a matrix even if only one column is selected

    # Store the smaller matrix in the list
    small_matrices[[length(small_matrices) + 1]] <- small_matrix
  }
  return(small_matrices)
}


tbj <- function(atmat) {
  # The code below produces a list where the number of elements in the list is J = length(b).
  # This is to help when we want to pick out tau_{mij} for each j = 1,...,J
  num_cols <- ncol(atmat[[1]])  # Get the number of columns in the matrices
  transposed_list <- vector("list", length = num_cols)  # Initialize an empty list to store transposed matrices

  # Loop through each column index
  for (i in 1:num_cols) {
    # Extract the i-th column from each matrix and bind them together into one matrix
    transposed_matrix <- do.call(cbind, lapply(atmat, function(mat) mat[, i]))
    # Add the transposed matrix to the transposed_list
    transposed_list[[i]] <- transposed_matrix
  }

  return(transposed_list)
}


M.b <- function(Y, Z, x, a = 2.1, p.b, limits = c(-8, 45)){
  s <- colSums(Y)/sum(Y)
  K <- length(Z[1,])
  new.b <- p.b
  M <- dim(Y)[1]
  N <- length(x)
  I <- length(a)
  J <- length(p.b)
  btmat <- tbj(split_tau(Z, J))
  lims <- c(limits[1], p.b[-1], limits[2])
  for(k in 2:J){
    obj <- function(b) {
      A <- outer(a, b * x, `+`)
      S <- matrix(log(s), nrow = I, ncol = N, byrow = TRUE)
      ExpPart <- -exp(A + S)
      Ypart <- A + S
      result <- sum(
        sapply(1:M, function(m) {
          W <- btmat[[k]][m, ]
          W_mat <- matrix(W, nrow = I, ncol = N)
          Y_vec <- Y[m, ]
          Y_mat <- matrix(Y_vec, nrow = I, ncol = N, byrow = TRUE)
          sum(W_mat * (ExpPart + Y_mat * Ypart))
        })
      )
      return(result)
    }
    new.b[k] <- optimize(obj, c(lims[k-1], lims[k+1]), maximum = T)[[1]]
    lims[k] <- new.b[k]
  }
  return(new.b)
}


M.a <- function(Y, Z, x, p.a = c(2), b, limits = c(-10,10)){
  s <- colSums(Y)/sum(Y)
  K <- length(Z[1,])
  log_s <- log(s)
  new.a <- p.a
  M = dim(Y)[1]
  N = length(x)
  I = length(p.a)
  J = length(b)
  atmat = split_tau(Z, J)
  for(k in 1:I){
    obj <- function(a) {
      linpart <- outer(b, x, function(b, x) a + b*x)  # J × N
      exp_part <- exp(t(t(linpart) + log_s))  # J × N
      term1 <- -exp_part  # J × N
      term2 <- t(t(linpart) + log_s)  # J × N
      # Multiply each (j, n) element by every (m, j) in atmat[[k]] and Y
      total <- 0
      for (m in 1:M) {
        # replicate across n
        y_m <- Y[m, ]  # length N
        for (j in 1:J) {
          a_kmj <- atmat[[k]][m, j]
          total <- total + sum(a_kmj * (term1[j, ] + y_m * term2[j, ]))
        }
      }
      return(total)
    }
    new.a[k] <- optimize(obj, limits, maximum = T)[[1]]
  }
  return(new.a)
}



OneIteration <- function(Y, x, pi.old, a.old, b.old, limits = c(-10, 45)){
  Z <- E(Y, x, pi.old, a.old, b.old)
  pi.new <- colMeans(Z)
  b.new <- M.b(Y, Z, x, a.old, b.old, limits)
  a.new <- M.a(Y, Z, x, a.old, b.old, limits)
  ###gets b.new using optim ###
  return(c(pi.new, a.new, b.new))
}


#' EM algorithm
#'
#' Implements the EM algorithm to get MLE estimates
#' @export
#' @param Y count data
#' @param x covariate
#' @param start.pi starting values for mixing proportions
#' @param start.a starting values for 'a'
#' @param start.b starting values for 'b'
#' @param max.it maximum number of iterations
#' @param tol tolerance level
#' @param limits interval to locate estimates
#' @examples
#' # Generate data
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'            pi_a = c(1), pi_b = c(0.8, 0.2))
#'
#' # Choose starting values
#' sv <- starting.values(Y = data, x = c(0,0,0,0,1,1,1,1))
#' a <- sv$a; b <- sv$b; pi.v <- sv$pi.v
#'
#' # Get MLE estimates
#' EMestimates <- EM.algorithm(Y = data, x = c(0,0,0,0,1,1,1,1), start.pi = pi.v,
#'                    start.a = a, start.b = b)
#' EMestimates
EM.algorithm <- function(Y, x, start.pi, start.a, start.b, max.it = 300, tol = 0.00001, limits = c(-10, 45))
{STOP=F
K <- length(start.pi)
I <- length(start.a)
J <- length(start.b)
new.pi = start.pi
new.a = start.a
new.b = start.b
parameters <- matrix(c(new.pi, new.a, new.b), nrow=1)
loglike = loglikelihood(Y, x, new.pi, new.a, new.b)
iteration = 1
print(c("interation ", " parameter estimates ", " log likelihood"))
while((STOP == F)&iteration< (max.it+1)){
  parameters<-rbind(parameters, OneIteration(Y, x, parameters[iteration, 1:K], parameters[iteration, (K+1):(K+I)], parameters[iteration, (K+I+1):(K+I+J)], limits))
  iteration<-iteration + 1
  loglike<-c(loglike, loglikelihood(Y, x, parameters[iteration, 1:K], parameters[iteration, (K+1):(K+I)], parameters[iteration, (K+I+1):(K+I+J)]))
  STOP<-max(loglike[iteration] - loglike[(iteration - 1)], abs(parameters[iteration] - parameters[iteration-1]))<tol
  print(c(round(iteration,0), round(parameters[iteration,],3),loglike[iteration]))
}
pi <- parameters[iteration, 1:K]
as <- parameters[iteration, (K+1):(K+I)]
bs <- parameters[iteration, (K+I+1):(K+I+J)]
ll <- loglike[iteration]
#fp <- (length(start.pi) - 1)*2
fp <- (I*J) + I + J - 1
AIC = -2*ll + 2*fp
BIC = -2*ll + log(length(Y[,1]))*fp
summaries = c(ll = ll, AIC = AIC, BIC = BIC)
return(list(pi = pi, as = as, bs = bs, summaries = summaries))
}



pi_split <- function(mat, small_cols) {
  ## Initialize a list to store the smaller matrices
  small_matrices <- list()

  ## Loop through the columns
  for (i in seq(1, ncol(mat), small_cols)) {
    ## Calculate the ending column index for the subset
    end_col <- min(i + small_cols - 1, ncol(mat))

    ## Subset the original matrix
    small_matrix <- mat[, i:end_col, drop = FALSE]  # Ensure that matrix remains a matrix even if only one column is selected

    ## Store the smaller matrix in the list
    small_matrices[[length(small_matrices) + 1]] <- small_matrix
  }
  return(small_matrices)
}


Lfdr <- function(Y, x, pi.v, a, b, eps = 0.36){
  s <- colSums(Y)/sum(Y)
  K = length(pi.v)
  M = nrow(Y)
  N = length(x)
  I = length(a)
  J = length(b)
  K = I*J
  b0 = b[b >= -eps & b <= eps]
  j0 = which(b >= -eps & b <= eps)
  K0 = I*length(b0)

  ### pi.est splits the pi vector so we can select single components e.g. pi_{12}
  pi.est = pi_split(matrix(pi.v, nrow = 1), J)
  pi.v0 = unlist(lapply(pi.est, `[`, j0))
  gpi.v0 <- rep(pi.v0, times = M)

  ### creating a matrix of means (b_0 component) for the poisson distn

  ij_grid0 <- expand.grid(j = 1:length(b0), i = 1:I)  # inner loop is j, outer is i
  mu_block0 <- apply(ij_grid0, 1, function(idx0) {
    i <- idx0["i"]
    j <- idx0["j"]
    exp(log(s) + a[i] + b0[j]*x)
  })
  mui0 <- do.call(rbind, replicate(M, t(mu_block0), simplify = FALSE))
  bigY0 <- Y[rep(1:nrow(Y), each = K0), ]

  ### creating a matrix of means [exp(a[i] + b[j]*x)] for the poisson distn
  ij_grid <- expand.grid(j = 1:J, i = 1:I)  # inner loop is j, outer is i
  mu_block <- apply(ij_grid, 1, function(idx) {
    i <- idx["i"]
    j <- idx["j"]
    exp(log(s) + a[i] + b[j]*x)
  })
  muij <- do.call(rbind, replicate(M, t(mu_block), simplify = FALSE))
  gpi.v <- rep(pi.v, times = M)
  rep.times <- rep(K, M)
  # Replicate each row according to replication_times
  rep.rows <- unlist(
    lapply(seq_len(nrow(Y)), function(i) {
      rep(list(Y[i, , drop = FALSE]), rep.times[i])
    }),
    recursive = FALSE
  )
  # Combine the replicated rows into a matrix
  bigY <- do.call(rbind, rep.rows)

  ### code for the numerator of the Lfdr_m including adjustments made for prob = 0 values
  pvc0 = c()
  for(i in 1:nrow(mui0)){
    pvc0[i] = log(gpi.v0[i]) + sum(dpois(bigY0[i,], lambda = mui0[i,], log = T))
  }
  num.pmat0 = matrix(pvc0, nrow = M, byrow = T)
  a0 <- apply(num.pmat0, 1, max)  # max of each row
  lprob0 <- a0 + log(rowSums(exp(num.pmat0 - a0)))

  ## code for the denominator of the Lfdr_m including adjustments made for prob = 0 values
  pvc = c()
  for(i in 1:nrow(muij)){
    pvc[i] = log(gpi.v[i]) + sum(dpois(bigY[i,], lambda = muij[i,], log = T))
  }
  num.pmat = matrix(pvc, nrow = M, byrow = T)
  a <- apply(num.pmat, 1, max)  # max of each row
  lprob <- a + log(rowSums(exp(num.pmat - a)))
  num = lprob0 - lprob
  e = exp(num)
  return(e)
}

#' Adaptive Lfdr Procedure
#'
#' Implements the adaptive Lfdr procedure to test for significant features
#' @export
#' @param Y count data
#' @param x covariate
#' @param pi.v mixing proportions/prior probabilities
#' @param a vector of 'a' components
#' @param b vector of 'b' components
#' @param alpha alpha level
#' @param epsilon epsilon value for the empirical null hypothesis
#' @examples
#' # Generate data
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'              pi_a = c(1), pi_b = c(0.8, 0.2))
#'
#' # Choose starting values
#' sv <- starting.values(Y = data, x = c(0,0,0,0,1,1,1,1))
#' a <- sv$a; b <- sv$b; pi.v <- sv$pi.v
#'
#' # Get MLE estimates
#' EMestimates <- EM.algorithm(Y = data, x = c(0,0,0,0,1,1,1,1), start.pi = pi.v,
#'                   start.a = a, start.b = b)
#' est.a <- EMestimates$as; est.b <- EMestimates$bs; est.pi <- EMestimates$pi
#'
#' # Implement procedure
#' procdata <- adaptive.procedure(Y = data, x = c(0,0,0,0,1,1,1,1), pi.v = est.pi,
#'                    a = est.a, b = est.b, epsilon = 0)
#' head(procdata)
adaptive.procedure <- function(Y, x, pi.v, a, b, alpha = 0.05, epsilon = 0){
  sv <- starting.values(Y, x = x)
  A <- sv$A
  B <- sv$B
  Y <- as.matrix(cbind(Y, A, B))
  Y1 = Y[, -c((ncol(Y)-1), ncol(Y))]
  lFDRs = Lfdr(Y1, x, pi.v, a, b, epsilon)
  B = Y[,ncol(Y)]
  B1 = B[order(lFDRs)]
  index = order(lFDRs)
  rlFDRs = lFDRs[order(lFDRs)]
  k = length(rlFDRs)
  crits = alpha*1:k
  summ = cumsum(rlFDRs)
  reject = rep(1,k)
  while(summ[k]>crits[k]&&k>0){
    reject[k] = 0
    k = k - 1
  }
  Ylfdr = Y[index,]
  result = data.frame(Ylfdr, crits, summ, rlFDRs, reject)
  num = rlFDRs*reject
  s.num = sum(num)
  s.den = sum(reject)
  est.FDR = s.num/s.den
  return(list(result.data = result, est.FDR = est.FDR))
}


#' MA plot for data with binary covariate
#'
#' Creates an MA plot for data with binary covariate
#' @export
#' @param lfdr.data result.data table from adaptive.procedure() function
#' @param row_name label of feature names
#' @importFrom dplyr %>%
#' @importFrom ggplot2 aes
#' @examples
#' # Generate data
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'                        pi_a = c(1), pi_b = c(0.8, 0.2))
#'
#' # Choose starting values
#' sv <- starting.values(Y = data, x = c(0,0,0,0,1,1,1,1))
#' a <- sv$a; b <- sv$b; pi.v <- sv$pi.v
#'
#' # Get MLE estimates
#' EMestimates <- EM.algorithm(Y = data, x = c(0,0,0,0,1,1,1,1), start.pi = pi.v,
#'                   start.a = a, start.b = b)
#' est.a <- EMestimates$as; est.b <- EMestimates$bs; est.pi <- EMestimates$pi
#'
#' # Implement procedure
#' procdata <- adaptive.procedure(Y = data, x = c(0,0,0,0,1,1,1,1), pi.v = est.pi,
#'                a = est.a, b = est.b)
#'
#' # Generate plot
#' MA.Plot(lfdr.data = procdata$result.data)
MA.Plot <- function(lfdr.data, row_name = "Features") {
  N <- ncol(lfdr.data[, -((ncol(lfdr.data)-5):ncol(lfdr.data))])
  mean2 <- rowMeans(lfdr.data[,1:(N/2)])
  mean1 <- rowMeans(lfdr.data[,((N/2)+1):N])
  A <- (mean1 + mean2)/2
  M <- log2((mean2)/(mean1))
  ma_df <- data.frame(row_name = rownames(lfdr.data), A = A, M = M, reject = lfdr.data$reject)

  ma_df <- ma_df %>%
    dplyr::mutate(
      is_labeled = reject == 1,
      label = ifelse(is_labeled, row_name, NA),
      color_group = ifelse(is_labeled, "sig", "non-sig")
    )

  ggplot2::ggplot(ma_df, aes(x = A, y = M, color = color_group)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(aes(label = label), max.overlaps = 15, size = 3) +
    ggplot2::scale_color_manual(values = c("sig" = "red", "non-sig" = "blue")) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      title = paste0("MA Plot of ",row_name),
      x = "Average Expression",
      y = "log2 Fold Change",
      color = "Decision"
    )
}

#' Plot for data with continuous covariate
#'
#' Creates a visual plot of MLE estimates of each feature in a data with continuous covariate
#' @export
#' @param lfdr.data output from adaptive.procedure() function
#' @importFrom ggplot2 aes
#' @examples
#' # Generate data
#' data <- generate.data(M = 100, x = c(0,0,0,0,1,1,1,1), a = 2, b = c(0, 0.7),
#'                        pi_a = c(1), pi_b = c(0.8, 0.2))
#'
#' # Choose starting values
#' sv <- starting.values(Y = data, x = c(0,0,0,0,1,1,1,1))
#' a <- sv$a; b <- sv$b; pi.v <- sv$pi.v
#'
#' # Get MLE estimates
#' EMestimates <- EM.algorithm(Y = data, x = c(0,0,0,0,1,1,1,1), start.pi = pi.v,
#'                   start.a = a, start.b = b)
#' est.a <- EMestimates$as; est.b <- EMestimates$bs; est.pi <- EMestimates$pi
#'
#' # Implement procedure
#' procdata <- adaptive.procedure(Y = data, x = c(0,0,0,0,1,1,1,1), pi.v = est.pi,
#'                a = est.a, b = est.b)
#'
#' # Generate plot
#' Continuous.plot(lfdr.data = procdata$result.data)
Continuous.plot = function(lfdr.data){
  ggplot2::ggplot(lfdr.data, aes(x = A, y = B, colour = as.factor(reject))) +
    ggplot2::geom_point(size = 1.3) +
    #geom_hline(yintercept = c(-epsilon, epsilon), color = "blue", linetype = "dashed", linewidth = 0.5)+
    ggplot2::labs(
      title = "Scatterplot of regression estimates",
      x = "a estimates",
      y = "b estimates"
    ) +
    ggplot2::theme_minimal()
}











