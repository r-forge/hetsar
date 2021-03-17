#-- Hessian of the function: phi(vtheta):=-(logL(vtheta) / T) --#
fn_hessian <- function(v_theta, l_data) {
     # retrieve data
     m_y  <- l_data[["m_y" ]]
     m_ys <- l_data[["m_ys"]]
     l_x  <- l_data[["l_x" ]]
     m_W  <- l_data[["m_W" ]]

     # sample size
     N  <- nrow(m_y)
     TT <- ncol(m_y)
     K  <- length(l_x)

     v_psi   <- v_theta[1:N] # (N,)
     v_beta  <- v_theta[(N + 1):(N + (K * N))]; # (KN,) x 1
     v_sgmsq <- v_theta[(N + (K * N) + 1):(N + (N * K) + N)] # (N,)

     # reshape beta
     m_beta_tr <- matrix(v_beta, K, N) # (K,N)
     m_beta <- t(m_beta_tr) # (N,K)

     v_sgm4h <- v_sgmsq^2
     v_sgm6h <- v_sgmsq^3

     # compute residuals
     m_beta_times_x <- matrix(0, N, TT) # 0 needed for recursive sum; (N x T) with generic element \vbeta_{i}'\vx_{it}
     for (k in 1:K) {
          m_x_k <- l_x[[k]] # (N,T)
          v_beta_k <- m_beta[, k] # (N,)
          m_beta_times_x_k <- m_x_k * v_beta_k #!! (N,T) "times" (N,1) = (N,T)
          m_beta_times_x <- m_beta_times_x + m_beta_times_x_k
     }
     m_psi_times_ys <- m_ys * v_psi #!! (N,T) "times" (N,1) = (N,T)
     m_eps <- m_y - m_psi_times_ys - m_beta_times_x # N x T
     v_ssr <- rowSums(m_eps^2) # (N,) sum of squared residuals: sum_t (eps_it)^2 !!crossprod()??

     m_Psi <- diag(v_psi, nrow = N) # (N,N)
     m_A <- diag(N) - (m_Psi %*% m_W)

     det_mA <- det(m_A)
     if (det_mA <= 0) {
          #stop('  !!TODO: replace with a meaningful error message')
          det_mA <- 1 ##################################################
     }

     # first derivative
     m_Q <- m_W %*% solve(m_A) #!! IS THERE A BETTER WAY TO WRITE THIS LINE?

     # second derivative
     m_H11 <- (m_Q * t(m_Q)) + diag(rowSums(m_ys^2) / v_sgmsq / TT, nrow = N) # N x N
     m_H13 <- diag(rowSums(m_ys * m_eps) / v_sgm4h / TT) # N x N
     m_H33 <- diag(-(1 / 2 / v_sgm4h) + (v_ssr / v_sgm6h / TT)) # N x N

     m_H12 <- matrix(0, N, N * K)
     m_H22 <- matrix(0, N * K, N * K)
     m_H23 <- matrix(0, N * K, N)

     # Note 1: the loop below can be probably eliminated in the same way v_dphi_dvbeta is computed above
     # Note 2: the code can be made probably faster by using sparse matrices
     for (i in 1:N) {
          ind <- ((i - 1) * K + 1):(i * K) #!! should I wrap this expression in c()?
          v_ysi <- m_ys[i, ] # (T,)
          m_Xi <- matrix(NA_real_, nrow = TT, ncol = K)
          for (k in 1:K) {
               m_x_k <- l_x[[k]] # (N,T)
               m_Xi[, k] <- m_x_k[i, ] # (T,)
          }
          v_epsi <- m_eps[i, ] # (T,)

          sgmsqi <- v_sgmsq[[i]]
          sgm4hi <- v_sgm4h[[i]]
          stopifnot(K > 1) #!! I think this bit of the code may not work when K=1 
          m_H12[i, ind] <- (rbind(v_ysi) %*% m_Xi) / sgmsqi / TT # (1,K)
          m_H22[ind, ind] <- (t(m_Xi) %*% m_Xi) / sgmsqi / TT # (K,K)
          m_H23[ind, i] <- (t(m_Xi) %*% cbind(v_epsi)) / sgm4hi / TT # (K,1)
     }

     m_H <- rbind(cbind(m_H11, m_H12, m_H13),
                  cbind(t(m_H12), m_H22, m_H23),
                  cbind(t(m_H13), t(m_H23), m_H33))

     #browser()
     #m_H <- Matrix::Matrix(m_H)
     #m_H <- Matrix::Matrix(m_H, sparse = TRUE)
     # impose sparsity as requested by trustOptim::trust.optim()
     m_H <- as(m_H,"CsparseMatrix") #!! ASK GIOVANNI: isn't it a waste to make a dense matrix sparse. Can we construct it sparse already?
     #print(class(m_H))
     #browser()

     # return
     m_H
}
